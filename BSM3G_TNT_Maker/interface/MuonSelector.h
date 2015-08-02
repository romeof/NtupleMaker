// 
//Authors: Andres Florez:      Universidad de los Andes, Colombia. 
//         kaur amandeepkalsi: Panjab University, India. 
// 
#ifndef __MUON_MU_H_                                                                                                            
#define __MUON_MU_H_
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
#include <TClonesArray.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
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
class MuonSelector : public  baseTree{
 public:
  MuonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~MuonSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtx);
 private:
  MuonSelector(){};
  /////
  //   Config variables
  /////
  edm::InputTag _muonToken;
  edm::InputTag _vertexInputTag; 
  double _Muon_pt_min;
  double _Muon_eta_max;
  int    _Muon_vtx_ndof_min;
  int    _Muon_vtx_rho_max;
  double _Muon_vtx_position_z_max;
  /////
  //   BSM methods/variables
  /////
  bool _super_TNT; //super tiny ntuple?
  vector <double> Muon_isoSum,Muon_isoCharParPt;
  vector <double> Muon_pt,Muon_eta,Muon_phi,Muon_dz,Muon_energy,Muon_iso,Muon_relIso;
  vector <double> Muon_isoCharged,Muon_isoNeutralHadron,Muon_isoPhoton,Muon_isoPU;
  vector <double> Muon_charge,Muon_chi2,Muon_p,Muon_matchedStat,Muon_dxy,Muon_validHits,Muon_validHitsInner,Muon_TLayers; 
  vector <int>    Muon_loose,Muon_tight,Muon_soft,Muon_isHighPt,Muon_pf,Muon_isglobal;   
  vector <int>    Muon_pdgId;
};
#endif
