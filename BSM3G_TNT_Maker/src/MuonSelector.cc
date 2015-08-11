#include "NtupleMaker/BSM3G_TNT_Maker/interface/MuonSelector.h"
MuonSelector::MuonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
 _muonToken               = iConfig.getParameter<edm::InputTag>("muons");
 _vertexInputTag          = iConfig.getParameter<edm::InputTag>("vertices");
 _Muon_pt_min             = iConfig.getParameter<double>("Muon_pt_min");
 _Muon_eta_max            = iConfig.getParameter<double>("Muon_eta_max");
 _Muon_vtx_ndof_min       = iConfig.getParameter<int>("Muon_vtx_ndof_min");
 _Muon_vtx_rho_max        = iConfig.getParameter<int>("Muon_vtx_rho_max");
 _Muon_vtx_position_z_max = iConfig.getParameter<double>("Muon_vtx_position_z_max");
 _super_TNT               = iConfig.getParameter<bool>("super_TNT");
 jetToken_                = iConfig.getParameter<edm::InputTag>("jets");
 SetBranches();
}
MuonSelector::~MuonSelector(){
 delete tree_;
}
void MuonSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 Clear();
 /////
 //   Recall collections
 ///// 
 edm::Handle<edm::View<pat::Muon> > muon_h;
 iEvent.getByLabel(_muonToken, muon_h);
 edm::Handle<reco::VertexCollection> vtx_h;
 iEvent.getByLabel(_vertexInputTag, vtx_h);
 edm::Handle<double> rhoHandle;
 iEvent.getByLabel("fixedGridRhoFastjetAll",rhoHandle);
 double rho = *rhoHandle;
 /////
 //   Require a good vertex 
 ///// 
 reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
 for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++){
  isGoodVertex(*it);
  firstGoodVertex = it;
  break;
 }
 if(firstGoodVertex == vtx_h->end()) return;
 /////
 //   Get muon information
 /////
 int mucoun = 0; //Only muons in Acceptance 
 for(edm::View<pat::Muon>::const_iterator mu = muon_h->begin(); mu != muon_h->end(); mu++){
  /////
  //   BSM variables
  /////
  //Acceptance 
  if(mu->pt()<_Muon_pt_min)         continue;
  if(fabs(mu->eta())>_Muon_eta_max) continue;  
  //Kinematics
  Muon_pt.push_back(mu->pt());
  Muon_ptError.push_back(mu->pt());
  Muon_eta.push_back(mu->eta());
  Muon_phi.push_back(mu->phi());
  Muon_energy.push_back(mu->energy());
  Muon_p.push_back(mu->p());
  //ID
  Muon_loose.push_back(mu->isLooseMuon());
  Muon_soft.push_back(mu->isSoftMuon(*firstGoodVertex));
  Muon_tight.push_back(mu->isTightMuon(*firstGoodVertex));
  Muon_isHighPt.push_back(mu->isHighPtMuon(*firstGoodVertex));
  Muon_pf.push_back(mu->isPFMuon());   
  Muon_isglobal.push_back(mu->isGlobalMuon());   
  Muon_pdgId.push_back(mu->pdgId());
  //Isolation
  Muon_isoSum.push_back((mu->trackIso() + mu->ecalIso() + mu->hcalIso()));
  Muon_isoCharParPt.push_back((mu->pfIsolationR04().sumChargedParticlePt));
  Muon_isoCharged.push_back((mu->pfIsolationR04().sumChargedHadronPt));
  double SumChargedHadronPt = mu->pfIsolationR04().sumChargedHadronPt;
  double SumNeutralEt       = mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt;
  double SumPU              = 0.5*mu->pfIsolationR04().sumPUPt;
  double SumNeutralCorrEt   = std::max( 0.0, SumNeutralEt - SumPU );
  double relIso = (SumChargedHadronPt + SumNeutralCorrEt)/mu->pt();
  Muon_relIso.push_back(relIso);
  //Track related variables and neutral part isolation
  reco::TrackRef gtk = mu->globalTrack();
  if(!_super_TNT){
   if(gtk.isNonnull()){
    Muon_chi2.push_back(gtk->normalizedChi2());
    Muon_validHits.push_back(gtk->hitPattern().numberOfValidMuonHits()); 
   }else{                
    Muon_chi2.push_back(-9999);
    Muon_validHits.push_back(-9999);  
   }
   Muon_matchedStat.push_back(mu->numberOfMatchedStations());
   Muon_dxy.push_back(mu->muonBestTrack()->dxy(firstGoodVertex->position())); 
   if(mu->innerTrack().isNonnull()){
    Muon_validHitsInner.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
    Muon_TLayers.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
   }else{
    Muon_validHitsInner.push_back(-999);
    Muon_TLayers.push_back(-999);
   }
   Muon_dz.push_back(mu->muonBestTrack()->dz(firstGoodVertex->position()));
   Muon_isoNeutralHadron.push_back((mu->pfIsolationR04().sumNeutralHadronEt));
   Muon_isoPhoton.push_back((mu->pfIsolationR04().sumPhotonEt));
   Muon_isoPU.push_back((mu->pfIsolationR04().sumPUPt));
  }
  //Other prop
  Muon_charge.push_back(mu->charge());
  /////
  //   HN variables
  /////
  //Kinematics
  Muon_dB.push_back(mu->dB());
  Muon_bestTrack_pT.push_back(mu->muonBestTrack()->pt());
  Muon_pTerrorOVERbestTrackpT.push_back(mu->muonBestTrack()->ptError()/mu->muonBestTrack()->pt());
  reco::TrackRef tunePBestTrack = mu->tunePMuonBestTrack();
  Muon_tunePBestTrack_pt.push_back(tunePBestTrack->pt());
  //ID
  Muon_isTrackerMuon.push_back(mu->isTrackerMuon());
  Muon_isMediumMuon.push_back(mu->isMediumMuon());
  Muon_POGisGood.push_back(muon::isGoodMuon(*mu, muon::TMOneStationTight));
  //Track related variables
  Muon_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
  Muon_trkKink.push_back(mu->combinedQuality().trkKink);
  Muon_segmentCompatibility.push_back(mu->segmentCompatibility());
  if(mu->innerTrack().isNonnull()){
   Muon_validFraction.push_back(mu->innerTrack()->validFraction());
   Muon_pixelLayersWithMeasurement.push_back(mu->innerTrack()->hitPattern().pixelLayersWithMeasurement());
   Muon_qualityhighPurity.push_back(mu->innerTrack()->quality(reco::TrackBase::highPurity));
  }else{
   Muon_validFraction.push_back(-9999);
   Muon_pixelLayersWithMeasurement.push_back(-9999);
   Muon_qualityhighPurity.push_back(-9999);
  }
  //Other prop
  reco::Muon::MuonTrackType tunePBestTrackType = mu->tunePMuonBestTrackType();
  Muon_tunePBestTrackType.push_back(tunePBestTrackType); 
  //Counter
  mucoun++;
  /////
  //   TTH variables
  ///// 
  if(mu->innerTrack().isNonnull() && mu->globalTrack().isNonnull()){
   //Kinematics
   Muon_pTErrorSUpT.push_back(mu->innerTrack()->ptError()/mu->innerTrack()->pt());   
   //ID
   //Isolation
   double iso_rho = get_iso_rho(*mu,rho);
   Muon_iso_rho.push_back(iso_rho);
   Muon_iso_nue_rel.push_back(iso_rho-mu->pfIsolationR03().sumChargedHadronPt/mu->pt());
   Muon_iso_ch_rel.push_back(mu->pfIsolationR03().sumChargedHadronPt/mu->pt());
   //Track related variables
   Muon_dxy_it.push_back(fabs(mu->innerTrack()->dxy(firstGoodVertex->position())));
   Muon_dz_it.push_back(fabs(mu->innerTrack()->dz(firstGoodVertex->position())));
   Muon_ip3d_val.push_back(fabs(mu->dB(pat::Muon::PV3D)));
   Muon_ip3d_err.push_back(mu->edB(pat::Muon::PV3D));
   Muon_ip3d_sig.push_back(fabs(mu->dB(pat::Muon::PV3D)/mu->edB(pat::Muon::PV3D)));
   Muon_nchi2_gt.push_back(mu->globalTrack()->normalizedChi2());
   //Other prop
   double mujet_mindr    = 999;
   double mujet_pt       = -1;  
   double muptSUmujetpt  = -1;
   double mujet_btagdisc = -1;
   get_mujet_info(*mu,iEvent,iSetup,mujet_mindr,mujet_pt,muptSUmujetpt,mujet_btagdisc);
   Muon_mujet_mindr.push_back(mujet_mindr);
   Muon_mujet_pt.push_back(mujet_pt);
   Muon_muptSUmujetpt.push_back(muptSUmujetpt);
   Muon_mujet_btagdisc.push_back(mujet_btagdisc);
  }else{
   Muon_pTErrorSUpT.push_back(-9999);
   Muon_iso_rho.push_back(-9999);
   Muon_iso_nue_rel.push_back(-9999);
   Muon_iso_ch_rel.push_back(-9999);
   Muon_dxy_it.push_back(-9999);
   Muon_dz_it.push_back(-9999);
   Muon_ip3d_val.push_back(-9999);
   Muon_ip3d_err.push_back(-9999);
   Muon_ip3d_sig.push_back(-9999);
   Muon_nchi2_gt.push_back(-9999);
   Muon_mujet_mindr.push_back(-9999);
   Muon_mujet_pt.push_back(-9999);
   Muon_muptSUmujetpt.push_back(-9999);
   Muon_mujet_btagdisc.push_back(-9999);
  }
 }
}
void MuonSelector::SetBranches(){
 if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
 /////
 //   BSM variables
 /////
 //Kinematics
 AddBranch(&Muon_pt               ,"Muon_pt");
 AddBranch(&Muon_ptError          ,"Muon_ptError");
 AddBranch(&Muon_eta              ,"Muon_eta");
 AddBranch(&Muon_phi              ,"Muon_phi");
 AddBranch(&Muon_energy           ,"Muon_energy");
 AddBranch(&Muon_p                ,"Muon_p");
 //ID
 AddBranch(&Muon_loose            ,"Muon_loose");
 AddBranch(&Muon_soft             ,"Muon_soft");
 AddBranch(&Muon_tight            ,"Muon_tight");
 AddBranch(&Muon_isHighPt         ,"Muon_isHighPt");
 AddBranch(&Muon_pf               ,"Muon_pf");
 AddBranch(&Muon_isglobal         ,"Muon_isglobal");
 AddBranch(&Muon_pdgId            ,"Muon_pdgId");
 //Isolation
 AddBranch(&Muon_isoSum           ,"Muon_isoSum");
 AddBranch(&Muon_isoCharParPt     ,"Muon_isoCharParPt");
 AddBranch(&Muon_isoCharged       ,"Muon_isoCharged");
 AddBranch(&Muon_relIso           ,"Muon_relIso");
 //Track related variables and neutral part isolation
 if(!_super_TNT){  
  AddBranch(&Muon_chi2            ,"Muon_chi2");
  AddBranch(&Muon_validHits       ,"Muon_validHits");
  AddBranch(&Muon_matchedStat     ,"Muon_matchedStat");
  AddBranch(&Muon_dxy             ,"Muon_dxy");
  AddBranch(&Muon_validHitsInner  ,"Muon_validHitsInner");
  AddBranch(&Muon_TLayers         ,"Muon_TLayers");
  AddBranch(&Muon_dz              ,"Muon_dz");
  AddBranch(&Muon_isoNeutralHadron,"Muon_isoNeutralHadron");
  AddBranch(&Muon_isoPhoton       ,"Muon_isoPhoton");
  AddBranch(&Muon_isoPU           ,"Muon_isoPU");
 }
 //Other prop
 AddBranch(&Muon_charge           ,"Muon_charge");
 /////
 //   HN variables
 /////
 //Kinematics
 AddBranch(&Muon_dB               ,"Muon_dB");
 AddBranch(&Muon_bestTrack_pT     ,"Muon_bestTrack_pT");
 AddBranch(&Muon_pTerrorOVERbestTrackpT ,"Muon_pTerrorOVERbestTrackpT");
 AddBranch(&Muon_tunePBestTrack_pt      ,"Muon_tunePBestTrack_pt");
 //ID
 AddBranch(&Muon_isTrackerMuon    ,"Muon_isTrackerMuon");
 AddBranch(&Muon_isMediumMuon     ,"Muon_isMediumMuon");
 AddBranch(&Muon_POGisGood        ,"Muon_POGisGood");
 //Track related variables
 AddBranch(&Muon_chi2LocalPosition          ,"Muon_chi2LocalPosition");
 AddBranch(&Muon_trkKink                    ,"Muon_trkKink");
 AddBranch(&Muon_segmentCompatibility       ,"Muon_segmentCompatibility");
 AddBranch(&Muon_validFraction              ,"Muon_validFraction");
 AddBranch(&Muon_pixelLayersWithMeasurement ,"Muon_pixelLayersWithMeasurement");
 AddBranch(&Muon_qualityhighPurity          ,"Muon_qualityhighPurity"); 
 //Other prop
 AddBranch(&Muon_tunePBestTrackType         ,"Muon_tunePBestTrackType");
 /////
 //   TTH variables
 /////
 //Kinematics
 AddBranch(&Muon_pTErrorSUpT                ,"Muon_pTErrorSUpT");
 //ID
 //Isolation
 AddBranch(&Muon_iso_rho                    ,"Muon_iso_rho");
 AddBranch(&Muon_iso_nue_rel                ,"Muon_iso_nue_rel");
 AddBranch(&Muon_iso_ch_rel                 ,"Muon_iso_ch_rel");
 //Track related variables
 AddBranch(&Muon_dxy_it                     ,"Muon_dxy_it");
 AddBranch(&Muon_dz_it                      ,"Muon_dz_it");
 AddBranch(&Muon_ip3d_val                   ,"Muon_ip3d_val");
 AddBranch(&Muon_ip3d_err                   ,"Muon_ip3d_err");
 AddBranch(&Muon_ip3d_sig                   ,"Muon_ip3d_sig");
 AddBranch(&Muon_nchi2_gt                   ,"Muon_nchi2_gt");
 //Other prop
 AddBranch(&Muon_mujet_mindr                ,"Muon_mujet_mindr");
 AddBranch(&Muon_mujet_pt                   ,"Muon_mujet_pt");
 AddBranch(&Muon_muptSUmujetpt              ,"Muon_muptSUmujetpt");
 AddBranch(&Muon_mujet_btagdisc             ,"Muon_mujet_btagdisc");
 if(debug_) std::cout<<"set branches"<<std::endl;
}
void MuonSelector::Clear(){
 /////
 //   BSM variables
 /////
 //Kinematics  
 Muon_pt.clear();
 Muon_ptError.clear();
 Muon_eta.clear();
 Muon_phi.clear();
 Muon_energy.clear();
 Muon_p.clear(); 
 //ID
 Muon_loose.clear();
 Muon_soft.clear();
 Muon_tight.clear();
 Muon_isHighPt.clear();
 Muon_pf.clear();   
 Muon_isglobal.clear();   
 Muon_pdgId.clear();
 //Isolation
 Muon_isoSum.clear();
 Muon_isoCharParPt.clear();
 Muon_isoCharged.clear();
 Muon_relIso.clear();
 //Track related variables and neutral part isolation
 Muon_chi2.clear(); 
 Muon_validHits.clear();
 Muon_matchedStat.clear(); 
 Muon_dxy.clear(); 
 Muon_validHitsInner.clear(); 
 Muon_TLayers.clear(); 
 Muon_dz.clear();
 Muon_isoNeutralHadron.clear();
 Muon_isoPhoton.clear();
 Muon_isoPU.clear();
 //Other prop
 Muon_charge.clear(); 
 /////
 //   HN variables
 /////
 //Kinematics
 Muon_dB.clear();
 Muon_bestTrack_pT.clear();
 Muon_pTerrorOVERbestTrackpT.clear();
 Muon_tunePBestTrack_pt.clear();
 //ID
 Muon_isTrackerMuon.clear();
 Muon_isMediumMuon.clear();
 Muon_POGisGood.clear();
 //Track related variables
 Muon_chi2LocalPosition.clear();
 Muon_trkKink.clear();
 Muon_segmentCompatibility.clear();
 Muon_validFraction.clear();
 Muon_pixelLayersWithMeasurement.clear();
 Muon_qualityhighPurity.clear();
 //Other prop
 Muon_tunePBestTrackType.clear();
 /////
 //   TTH variables
 /////
 //Kinematics
 Muon_pTErrorSUpT.clear();
 //ID
 //Isolation
 Muon_iso_rho.clear();
 Muon_iso_nue_rel.clear();
 Muon_iso_ch_rel.clear();
 //Track related variables
 Muon_dxy_it.clear();
 Muon_dz_it.clear();
 Muon_ip3d_val.clear();
 Muon_ip3d_err.clear();
 Muon_ip3d_sig.clear();
 Muon_nchi2_gt.clear();
 //Other prop
 Muon_mujet_mindr.clear();
 Muon_mujet_pt.clear();
 Muon_muptSUmujetpt.clear();
 Muon_mujet_btagdisc.clear();
}
bool MuonSelector::isGoodVertex(const reco::Vertex& vtx){
 if(vtx.isFake())                                        return false;
 if(vtx.ndof()<_Muon_vtx_ndof_min)                       return false;
 if(vtx.position().Rho()>_Muon_vtx_rho_max)              return false;
 if(fabs(vtx.position().Z()) > _Muon_vtx_position_z_max) return false;
 return true;
}
double MuonSelector::get_iso_rho(const pat::Muon& mu, double& rho){
 double Eta = mu.eta();
 double pfIsoCharged = mu.pfIsolationR03().sumChargedHadronPt;
 double pfIsoNeutral = mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt;
 double EffArea = 9999.;
 if(abs(Eta) < 0.8)      EffArea = 0.0913;
 else if(abs(Eta) < 1.3) EffArea = 0.0765;
 else if(abs(Eta) < 2.0) EffArea = 0.0546;
 else if(abs(Eta) < 2.2) EffArea = 0.0728;
 else                    EffArea = 0.1177;
 double correction = rho*EffArea;
 double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction);
 double iso = (pfIsoCharged + pfIsoPUSubtracted)/mu.pt();
 return iso;
}
void MuonSelector::get_mujet_info(const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup, double& mujet_mindr, double& mujet_pt, double& muptSUmujetpt, double& mujet_btagdisc){
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByLabel(jetToken_, jets);
 pat::Jet mujet;
 const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );
 for(const pat::Jet &j : *jets){
  pat::Jet jet = j.correctedJet(0);
  double scale = corrector->correction(jet, iEvent, iSetup);
  jet.scaleEnergy(scale);
  double dr = deltaR(mu.p4(),jet.p4());
  if(dr<mujet_mindr){
   mujet_mindr = dr;
   mujet       = jet;
  }
 }
 mujet_pt        = mujet.pt(); 
 muptSUmujetpt   = min(mu.pt()/mujet.pt(), 1.5);
 mujet_btagdisc  = max(mujet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"), float(0.0)); 
}
