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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
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
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  //   BSM methods 
  bool isGoodVertex(const reco::Vertex& vtx);
  //   HN methods
  //   TTH methods
  double get_iso_rho(const pat::Muon& mu, double& rho); 
  void   get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& mujet_mindr, double& mujet_pt, double& muptSUmujetpt, double& mujet_btagdisc);
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
  bool   _super_TNT; //super tiny ntuple?
  /////
  //   BSM variables
  /////
  vector<double> Muon_pt, Muon_ptError, Muon_eta, Muon_phi, Muon_energy, Muon_p;
  vector<double> Muon_loose, Muon_soft, Muon_tight, Muon_isHighPt, Muon_pf, Muon_isglobal, Muon_pdgId; 
  vector<double> Muon_isoSum, Muon_isoCharParPt, Muon_isoCharged, Muon_relIso;
  vector<double> Muon_chi2, Muon_validHits, Muon_matchedStat, Muon_dxy, Muon_validHitsInner, Muon_TLayers, Muon_dz, Muon_isoNeutralHadron, Muon_isoPhoton, Muon_isoPU;
  vector<double> Muon_charge;
  /////
  //   HN variables
  /////
  vector<double> Muon_dB, Muon_bestTrack_pT, Muon_pTerrorOVERbestTrackpT, Muon_tunePBestTrack_pt;
  vector<double> Muon_isTrackerMuon, Muon_isMediumMuon, Muon_POGisGood;
  vector<double> Muon_chi2LocalPosition, Muon_trkKink, Muon_segmentCompatibility, Muon_validFraction, Muon_pixelLayersWithMeasurement, Muon_qualityhighPurity;
  vector<double> Muon_tunePBestTrackType;
  /////
  //   TTH variables
  /////
  edm::InputTag jetToken_;
  vector<double> Muon_pTErrorSUpT;
  vector<double> Muon_iso_rho, Muon_iso_nue_rel, Muon_iso_ch_rel;
  vector<double> Muon_dxy_it, Muon_dz_it, Muon_ip3d_val, Muon_ip3d_err, Muon_ip3d_sig, Muon_nchi2_gt;
  vector<double> Muon_mujet_mindr, Muon_mujet_pt, Muon_muptSUmujetpt, Muon_mujet_btagdisc;
};
#endif
