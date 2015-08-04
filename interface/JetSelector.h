#ifndef __JET_MU_H_
#define __JET_MU_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>                                                                    
#include <TClonesArray.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "baseTree.h"
using namespace std;
using namespace pat;
using namespace edm;
/////
//   Class declaration
/////
class JetSelector : public  baseTree{
 public:
  JetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~JetSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
 private:
  JetSelector(){};
  /////
  //   Config variables
  /////
  edm::InputTag jetToken_;
  edm::InputTag _vertexInputTag;
  double _Jet_pt_min;
  /////
  //   BSM variables
  /////
  vector <double> Jet_pt,Jet_eta,Jet_phi,Jet_energy,Jet_bDiscriminator,Jet_mass,JetParton,JetjetId;
  vector <double> Jet_pileupId,JetIDPU,Jetpass_pileupJetId,Jet_neutralHadEnergyFraction,Jet_neutralEmEmEnergyFraction; 
  vector <double> Jet_chargedHadronEnergyFraction,Jet_chargedEmEnergyFraction,Jet_muonEnergyFraction; 
  vector <double> Jet_electronEnergy,Jet_photonEnergy,UncorrJet_pt; 
  vector <int> Jet_numberOfConstituents;
  vector <int> Jet_chargedMultiplicity;
  bool _super_TNT;
};
#endif
