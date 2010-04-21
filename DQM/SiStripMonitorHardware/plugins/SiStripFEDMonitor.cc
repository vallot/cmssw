// -*- C++ -*-
//
// Package:    DQM/SiStripMonitorHardware
// Class:      SiStripFEDMonitorPlugin
// 
/**\class SiStripFEDMonitorPlugin SiStripFEDMonitor.cc DQM/SiStripMonitorHardware/plugins/SiStripFEDMonitor.cc

 Description: DQM source application to produce data integrety histograms for SiStrip data
*/
//
// Original Author:  Nicholas Cripps
//         Created:  2008/09/16
// $Id: SiStripFEDMonitor.cc,v 1.37 2010/04/02 15:41:32 amagnan Exp $
//
//Modified        :  Anne-Marie Magnan
//   ---- 2009/04/21 : histogram management put in separate class
//                     struct helper to simplify arguments of functions
//   ---- 2009/04/22 : add TkHistoMap with % of bad channels per module
//   ---- 2009/04/27 : create FEDErrors class 

#include <sstream>
#include <memory>
#include <list>
#include <algorithm>
#include <cassert>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"

#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"

#include "DQMServices/Core/interface/DQMStore.h"

#include "DQM/SiStripMonitorHardware/interface/FEDHistograms.hh"
#include "DQM/SiStripMonitorHardware/interface/FEDErrors.hh"


//
// Class declaration
//

class SiStripFEDMonitorPlugin : public edm::EDAnalyzer
{
 public:
  explicit SiStripFEDMonitorPlugin(const edm::ParameterSet&);
  ~SiStripFEDMonitorPlugin();
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
				    const edm::EventSetup& context);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
				  const edm::EventSetup& context);

  //update the cabling if necessary
  void updateCabling(const edm::EventSetup& eventSetup);


  //tag of FEDRawData collection
  edm::InputTag rawDataTag_;
  //histogram helper class
  FEDHistograms fedHists_;
  //folder name for histograms in DQMStore
  std::string folderName_;
  //book detailed histograms even if they will be empty (for merging)
  bool fillAllDetailedHistograms_;
  //do histos vs time with time=event number. Default time = orbit number (s)
  bool fillWithEvtNum_;
  //print debug messages when problems are found: 1=error debug, 2=light debug, 3=full debug
  unsigned int printDebug_;
  //bool printDebug_;
  //write the DQMStore to a root file at the end of the job
  bool writeDQMStore_;
  std::string dqmStoreFileName_;
  //the DQMStore
  DQMStore* dqm_;
  //FED cabling
  uint32_t cablingCacheId_;
  const SiStripFedCabling* cabling_;

  //add parameter to save computing time if TkHistoMap are not filled
  bool doTkHistoMap_;
  bool doMedHists_;

  unsigned int nEvt_;

  //FED errors
  //need class member for lumi histograms
  FEDErrors fedErrors_;

};


//
// Constructors and destructor
//

SiStripFEDMonitorPlugin::SiStripFEDMonitorPlugin(const edm::ParameterSet& iConfig)
  : rawDataTag_(iConfig.getUntrackedParameter<edm::InputTag>("RawDataTag",edm::InputTag("source",""))),
    folderName_(iConfig.getUntrackedParameter<std::string>("HistogramFolderName","SiStrip/ReadoutView/FedMonitoringSummary")),
    fillAllDetailedHistograms_(iConfig.getUntrackedParameter<bool>("FillAllDetailedHistograms",false)),
    fillWithEvtNum_(iConfig.getUntrackedParameter<bool>("FillWithEventNumber",false)),
    printDebug_(iConfig.getUntrackedParameter<unsigned int>("PrintDebugMessages",1)),
    //printDebug_(iConfig.getUntrackedParameter<bool>("PrintDebugMessages",false)),
    writeDQMStore_(iConfig.getUntrackedParameter<bool>("WriteDQMStore",false)),
    dqmStoreFileName_(iConfig.getUntrackedParameter<std::string>("DQMStoreFileName","DQMStore.root")),
    dqm_(0),
    cablingCacheId_(0)
{
  //print config to debug log
  std::ostringstream debugStream;
  if (printDebug_>1) {
    debugStream << "[SiStripFEDMonitorPlugin]Configuration for SiStripFEDMonitorPlugin: " << std::endl
                << "[SiStripFEDMonitorPlugin]\tRawDataTag: " << rawDataTag_ << std::endl
                << "[SiStripFEDMonitorPlugin]\tHistogramFolderName: " << folderName_ << std::endl
                << "[SiStripFEDMonitorPlugin]\tFillAllDetailedHistograms? " << (fillAllDetailedHistograms_ ? "yes" : "no") << std::endl
		<< "[SiStripFEDMonitorPlugin]\tFillWithEventNumber?" << (fillWithEvtNum_ ? "yes" : "no") << std::endl
                << "[SiStripFEDMonitorPlugin]\tPrintDebugMessages? " << (printDebug_ ? "yes" : "no") << std::endl
                << "[SiStripFEDMonitorPlugin]\tWriteDQMStore? " << (writeDQMStore_ ? "yes" : "no") << std::endl;
    if (writeDQMStore_) debugStream << "[SiStripFEDMonitorPlugin]\tDQMStoreFileName: " << dqmStoreFileName_ << std::endl;
  }
  
  //don;t generate debug mesages if debug is disabled
  std::ostringstream* pDebugStream = (printDebug_>1 ? &debugStream : NULL);
  
  fedHists_.initialise(iConfig,pDebugStream);

  doTkHistoMap_ = fedHists_.tkHistoMapEnabled();

  doMedHists_ = fedHists_.cmHistosEnabled();

  if (printDebug_) {
    LogTrace("SiStripMonitorHardware") << debugStream.str();

    //debugStream.str("");

    //debugStream << " -- Quelle est la difference entre un canard ? " << std::endl 
    //	<< " -- Reponse: c'est qu'il a les deux pattes de la meme longueur, surtout la gauche." << std::endl;

    //edm::LogError("SiStripMonitorHardware") << debugStream.str();
  }

  nEvt_ = 0;

}

SiStripFEDMonitorPlugin::~SiStripFEDMonitorPlugin()
{
}


//
// Member functions
//

// ------------ method called to for each event  ------------
void
SiStripFEDMonitorPlugin::analyze(const edm::Event& iEvent, 
				 const edm::EventSetup& iSetup)
{
  //update cabling
  updateCabling(iSetup);
  
  //get raw data
  edm::Handle<FEDRawDataCollection> rawDataCollectionHandle;
  iEvent.getByLabel(rawDataTag_,rawDataCollectionHandle);
  const FEDRawDataCollection& rawDataCollection = *rawDataCollectionHandle;
  
  fedErrors_.initialiseEvent();

  //initialise map of fedId/bad channel number
  std::map<unsigned int,std::pair<unsigned short,unsigned short> > badChannelFraction;

  unsigned int lNFEDMonitoring = 0;
  unsigned int lNFEDUnpacker = 0;
  unsigned int lNChannelMonitoring = 0;
  unsigned int lNChannelUnpacker = 0;

  unsigned int lNTotBadFeds = 0;
  unsigned int lNTotBadChannels = 0;
  unsigned int lNTotBadActiveChannels = 0;

  //loop over siStrip FED IDs
  for (unsigned int fedId = FEDNumbering::MINSiStripFEDID; 
       fedId <= FEDNumbering::MAXSiStripFEDID; 
       fedId++) {//loop over FED IDs
    const FEDRawData& fedData = rawDataCollection.FEDData(fedId);

    //create an object to fill all errors
    fedErrors_.initialiseFED(fedId,cabling_);
    bool lFullDebug = false;
 
    //Do detailed check
    //first check if data exists
    bool lDataExist = fedErrors_.checkDataPresent(fedData);
    if (!lDataExist) {
      fedHists_.fillFEDHistograms(fedErrors_,lFullDebug);
      continue;
    }


 
    //check for problems and fill detailed histograms
    fedErrors_.fillFEDErrors(fedData,
			     lFullDebug,
			     printDebug_,
			     lNChannelMonitoring,
			     lNChannelUnpacker,
			     doMedHists_,
			     fedHists_.cmHistPointer(false),
			     fedHists_.cmHistPointer(true)
			     );

    //check filled in previous method.
    bool lFailUnpackerFEDcheck = fedErrors_.failUnpackerFEDCheck();

    fedErrors_.incrementFEDCounters();
    fedHists_.fillFEDHistograms(fedErrors_,lFullDebug);


    bool lFailMonitoringFEDcheck = fedErrors_.failMonitoringFEDCheck();
    if (lFailMonitoringFEDcheck) lNTotBadFeds++;


    //sanity check: if something changed in the unpacking code 
    //but wasn't propagated here
    //print only the summary, and more info if printDebug>1
    if (lFailMonitoringFEDcheck != lFailUnpackerFEDcheck) {
      if (printDebug_>1) {
	std::ostringstream debugStream;
	debugStream << " --- WARNING: FED " << fedId << std::endl 
		    << " ------ Monitoring FED check " ;
	if (lFailMonitoringFEDcheck) debugStream << "failed." << std::endl;
	else debugStream << "passed." << std::endl ;
	debugStream << " ------ Unpacker FED check " ;
	if (lFailUnpackerFEDcheck) debugStream << "failed." << std::endl;
	else debugStream << "passed." << std::endl ;
	edm::LogError("SiStripMonitorHardware") << debugStream.str();
      }

      if (lFailMonitoringFEDcheck) lNFEDMonitoring++;
      else if (lFailUnpackerFEDcheck) lNFEDUnpacker++;
    }


    //Fill TkHistoMap:
    //add an entry for all channels (good = 0), 
    //so that tkHistoMap knows which channels should be there.
    if (doTkHistoMap_ && !fedHists_.tkHistoMapPointer()) {
      edm::LogWarning("SiStripMonitorHardware") << " -- Fedid " << fedId
						<< ", TkHistoMap enabled but pointer is null." << std::endl;
    }

    fedErrors_.fillBadChannelList(doTkHistoMap_,
				  fedHists_.tkHistoMapPointer(),
				  lNTotBadChannels,
				  lNTotBadActiveChannels);
  }//loop over FED IDs

  if ((lNTotBadFeds> 0 || lNTotBadChannels>0) && printDebug_>1) {
    std::ostringstream debugStream;
    debugStream << "[SiStripFEDMonitorPlugin] --- Total number of bad feds = " 
		<< lNTotBadFeds << std::endl
		<< "[SiStripFEDMonitorPlugin] --- Total number of bad channels = " 
		<< lNTotBadChannels << std::endl
		<< "[SiStripFEDMonitorPlugin] --- Total number of bad active channels = " 
		<< lNTotBadActiveChannels << std::endl;
    edm::LogInfo("SiStripMonitorHardware") << debugStream.str();
  }

  if ((lNFEDMonitoring > 0 || lNFEDUnpacker > 0 || lNChannelMonitoring > 0 || lNChannelUnpacker > 0) && printDebug_) {
    std::ostringstream debugStream;
    debugStream
      << "[SiStripFEDMonitorPlugin]-------------------------------------------------------------------------" << std::endl 
      << "[SiStripFEDMonitorPlugin]-------------------------------------------------------------------------" << std::endl 
      << "[SiStripFEDMonitorPlugin]-- Summary of differences between unpacker and monitoring at FED level : " << std::endl 
      << "[SiStripFEDMonitorPlugin] ---- Number of times monitoring fails but not unpacking = " << lNFEDMonitoring << std::endl 
      << "[SiStripFEDMonitorPlugin] ---- Number of times unpacking fails but not monitoring = " << lNFEDUnpacker << std::endl
      << "[SiStripFEDMonitorPlugin]-------------------------------------------------------------------------" << std::endl 
      << "[SiStripFEDMonitorPlugin]-- Summary of differences between unpacker and monitoring at Channel level : " << std::endl 
      << "[SiStripFEDMonitorPlugin] ---- Number of times monitoring fails but not unpacking = " << lNChannelMonitoring << std::endl 
      << "[SiStripFEDMonitorPlugin] ---- Number of times unpacking fails but not monitoring = " << lNChannelUnpacker << std::endl
      << "[SiStripFEDMonitorPlugin]-------------------------------------------------------------------------" << std::endl 
      << "[SiStripFEDMonitorPlugin]-------------------------------------------------------------------------" << std::endl ;
    edm::LogError("SiStripMonitorHardware") << debugStream.str();

  }

  FEDErrors::getFEDErrorsCounters().nTotalBadChannels = lNTotBadChannels;
  FEDErrors::getFEDErrorsCounters().nTotalBadActiveChannels = lNTotBadActiveChannels;

  //fedHists_.fillCountersHistograms(FEDErrors::getFEDErrorsCounters(), nEvt_);
  //time in seconds since beginning of the run or event number
  if (fillWithEvtNum_) fedHists_.fillCountersHistograms(FEDErrors::getFEDErrorsCounters(),FEDErrors::getChannelErrorsCounters(),iEvent.id().event());
  else fedHists_.fillCountersHistograms(FEDErrors::getFEDErrorsCounters(),FEDErrors::getChannelErrorsCounters(),iEvent.orbitNumber()/11223.);


  nEvt_++;

}//analyze method

// ------------ method called once each job just before starting event loop  ------------
void 
SiStripFEDMonitorPlugin::beginJob()
{
  //get DQM store
  dqm_ = &(*edm::Service<DQMStore>());
  dqm_->setCurrentFolder(folderName_);
  
  //this propagates dqm_ to the histoclass, must be called !
  fedHists_.bookTopLevelHistograms(dqm_);
  
  if (fillAllDetailedHistograms_) fedHists_.bookAllFEDHistograms();

  nEvt_ = 0;

  //const unsigned int siStripFedIdMin = FEDNumbering::MINSiStripFEDID;
  //const unsigned int siStripFedIdMax = FEDNumbering::MAXSiStripFEDID;

  //mark all channels as inactive until they have been 'locked' at least once
  //   activeChannels_.resize(siStripFedIdMax+1);
  //   for (unsigned int fedId = siStripFedIdMin; 
  //        fedId <= siStripFedIdMax; 
  //        fedId++) {
  //     activeChannels_[fedId].resize(sistrip::FEDCH_PER_FED,false);
  //   }
  



}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripFEDMonitorPlugin::endJob()
{
  if (writeDQMStore_) dqm_->save(dqmStoreFileName_);
}



void 
SiStripFEDMonitorPlugin::beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
					      const edm::EventSetup& context)
{

  fedErrors_.initialiseLumiBlock();

}


void 
SiStripFEDMonitorPlugin::endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
					    const edm::EventSetup& context)
{
  fedHists_.fillLumiHistograms(fedErrors_.getLumiErrors());
}




void SiStripFEDMonitorPlugin::updateCabling(const edm::EventSetup& eventSetup)
{
  uint32_t currentCacheId = eventSetup.get<SiStripFedCablingRcd>().cacheIdentifier();
  if (cablingCacheId_ != currentCacheId) {
    edm::ESHandle<SiStripFedCabling> cablingHandle;
    eventSetup.get<SiStripFedCablingRcd>().get(cablingHandle);
    cabling_ = cablingHandle.product();
    cablingCacheId_ = currentCacheId;
  }
}


//
// Define as a plug-in
//

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SiStripFEDMonitorPlugin);
