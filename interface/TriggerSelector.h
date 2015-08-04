// 
// Authors: Andres Florez:      Universidad de los Andes, Colombia. 
//          kaur amandeepkalsi: Panjab University, India. 
// 
#ifndef __TRIGGER_H_ 
#define __TRIGGER_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/Muon/interface/HLTMuonIsoFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "baseTree.h"
#include <TBranch.h>
using namespace std;
using namespace edm;
/////
//   Class declaration
/////
class TriggerSelector : public baseTree{
 public:
  TriggerSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~TriggerSelector();
  virtual void startTrigger (edm::EventSetup const& , edm::Run const &);
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
 private:
  TriggerSelector(){};
  /////
  //   Config variables
  /////
  HLTConfigProvider hltConfig_;
  edm::InputTag triggerResultsTag_;
  /////
  //   BSM methods/variables
  /////
  vector<int> Trigger_decision;
  vector<string> Trigger_names;
};
#endif
