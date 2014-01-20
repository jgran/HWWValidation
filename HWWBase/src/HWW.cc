#include "HWWValidation/HWWBase/interface/HWW.h"

static HWW hww;
HWW& GetHWW() {return hww;}

namespace HWWVal {

//to access data members
  //vertex 
  std::vector<LorentzVector>  &vtxs_position(){
    return(GetHWW().vtxs_position());
  }
  std::vector<float>          &vtxs_ndof(){
    return(GetHWW().vtxs_ndof());
  }
  std::vector<float>          &vtxs_sumpt(){
    return(GetHWW().vtxs_sumpt());
  }
  std::vector<int>            &vtxs_isFake(){
    return(GetHWW().vtxs_isFake());
  }
  std::vector<float>          &vtxs_xError(){
    return(GetHWW().vtxs_xError());
  }
  std::vector<float>          &vtxs_yError(){
    return(GetHWW().vtxs_yError());
  }
  std::vector<float>          &vtxs_zError(){
    return(GetHWW().vtxs_zError());
  }
  std::vector<vector<float>>  &vtxs_covMatrix(){
    return(GetHWW().vtxs_covMatrix());
  }

  //tracks
  std::vector<LorentzVector>  &trks_trk_p4(){
    return(GetHWW().trks_trk_p4());
  }
  std::vector<float>          &trks_chi2(){
    return(GetHWW().trks_chi2());
  }
  std::vector<float>          &trks_ndof(){
    return(GetHWW().trks_ndof());
  }
  std::vector<float>          &trks_d0(){
    return(GetHWW().trks_d0());
  }
  std::vector<int>            &trks_nlayers(){
    return(GetHWW().trks_nlayers());
  }
  std::vector<int>            &trks_valid_pixelhits(){
    return(GetHWW().trks_valid_pixelhits());
  }
  std::vector<float>          &trks_z0(){
    return(GetHWW().trks_z0());
  }
  std::vector<float>          &trks_z0Err(){
    return(GetHWW().trks_z0Err());
  }
  std::vector<float>          &trks_etaErr(){
    return(GetHWW().trks_etaErr());
  }
  std::vector<float>          &trks_d0Err(){
    return(GetHWW().trks_d0Err());
  }
  std::vector<float>          &trks_phiErr(){
    return(GetHWW().trks_phiErr());
  }
  std::vector<float>          &trks_d0phiCov(){
    return(GetHWW().trks_d0phiCov());
  }
  std::vector<int>            &trks_qualityMask(){
    return(GetHWW().trks_qualityMask());
  }
  std::vector<int>            &trks_charge(){
    return(GetHWW().trks_charge());
  }

  //electrons
  std::vector<LorentzVector>  &els_p4(){
    return(GetHWW().els_p4());
  }
  std::vector<LorentzVector>  &els_trk_p4(){
    return(GetHWW().els_trk_p4());
  }
  std::vector<LorentzVector>  &els_vertex_p4(){
    return(GetHWW().els_vertex_p4());
  }
  std::vector<float>          &els_lh(){
    return(GetHWW().els_lh());
  }
  std::vector<float>          &els_etaSC(){
    return(GetHWW().els_etaSC());
  }
  std::vector<float>          &els_sigmaIEtaIEta(){
    return(GetHWW().els_sigmaIEtaIEta());
  }
  std::vector<float>          &els_dEtaIn(){
    return(GetHWW().els_dEtaIn());
  }
  std::vector<float>          &els_dPhiIn(){
    return(GetHWW().els_dPhiIn());
  }
  std::vector<float>          &els_hOverE(){
    return(GetHWW().els_hOverE());
  }
  std::vector<float>          &els_tkIso(){
    return(GetHWW().els_tkIso());
  }
  std::vector<float>          &els_d0corr(){
    return(GetHWW().els_d0corr());
  }
  std::vector<float>          &els_d0(){
    return(GetHWW().els_d0());
  }
  std::vector<float>          &els_z0corr(){
    return(GetHWW().els_z0corr());
  }
  std::vector<float>          &els_fbrem(){
    return(GetHWW().els_fbrem());
  }
  std::vector<float>          &els_eOverPIn(){
    return(GetHWW().els_eOverPIn());
  }
  std::vector<float>          &els_eSeedOverPOut(){
    return(GetHWW().els_eSeedOverPOut());
  }
  std::vector<float>          &els_eSeedOverPIn(){
    return(GetHWW().els_eSeedOverPIn());
  }
  std::vector<float>          &els_sigmaIPhiIPhi(){
    return(GetHWW().els_sigmaIPhiIPhi());
  }
  std::vector<float>          &els_eSC(){
    return(GetHWW().els_eSC());
  }
  std::vector<float>          &els_ip3d(){
    return(GetHWW().els_ip3d());
  }
  std::vector<float>          &els_ip3derr(){
    return(GetHWW().els_ip3derr());
  }
  std::vector<float>          &els_chi2(){
    return(GetHWW().els_chi2());
  }
  std::vector<float>          &els_ndof(){
    return(GetHWW().els_ndof());
  }
  std::vector<float>          &els_dEtaOut(){
    return(GetHWW().els_dEtaOut());
  }
  std::vector<float>          &els_dPhiOut(){
    return(GetHWW().els_dPhiOut());
  }
  std::vector<float>          &els_eSCRaw(){
    return(GetHWW().els_eSCRaw());
  }
  std::vector<float>          &els_etaSCwidth(){
    return(GetHWW().els_etaSCwidth());
  }
  std::vector<float>          &els_phiSCwidth(){
    return(GetHWW().els_phiSCwidth());
  }
  std::vector<float>          &els_eSCPresh(){
    return(GetHWW().els_eSCPresh());
  }
  std::vector<float>          &els_iso03_pf_ch(){
    return(GetHWW().els_iso03_pf_ch());
  }
  std::vector<float>          &els_iso03_pf_nhad05(){
    return(GetHWW().els_iso03_pf_nhad05());
  }
  std::vector<float>          &els_iso03_pf_gamma05(){
    return(GetHWW().els_iso03_pf_gamma05());
  }
  std::vector<float>          &els_iso04_pf_ch(){
    return(GetHWW().els_iso04_pf_ch());
  }
  std::vector<float>          &els_iso04_pf_nhad05(){
    return(GetHWW().els_iso04_pf_nhad05());
  }
  std::vector<float>          &els_iso04_pf_gamma05(){
    return(GetHWW().els_iso04_pf_gamma05());
  }
  std::vector<float>          &els_e5x5(){
    return(GetHWW().els_e5x5());
  }
  std::vector<float>          &els_e1x5(){
    return(GetHWW().els_e1x5());
  }
  std::vector<float>          &els_e3x3(){
    return(GetHWW().els_e3x3());
  }
  std::vector<float>          &els_ecalEnergy(){
    return(GetHWW().els_ecalEnergy());
  }
  std::vector<float>          &els_eOverPOut(){
    return(GetHWW().els_eOverPOut());
  }
  std::vector<float>          &els_ecalIso(){
    return(GetHWW().els_ecalIso());
  }
  std::vector<float>          &els_hcalIso(){
    return(GetHWW().els_hcalIso());
  }
  std::vector<float>          &els_trkshFrac(){
    return(GetHWW().els_trkshFrac());
  }
  std::vector<float>          &els_conv_dist(){
    return(GetHWW().els_conv_dist());
  }
  std::vector<float>          &els_conv_dcot(){
    return(GetHWW().els_conv_dcot());
  }
  std::vector<float>          &els_conv_old_dist(){
    return(GetHWW().els_conv_old_dist());
  }
  std::vector<float>          &els_conv_old_dcot(){
    return(GetHWW().els_conv_old_dcot());
  }
  std::vector<float>          &els_iso04_pf2012_ch(){
    return(GetHWW().els_iso04_pf2012_ch());
  }
  std::vector<float>          &els_iso04_pf2012_em(){
    return(GetHWW().els_iso04_pf2012_em());
  }
  std::vector<float>          &els_iso04_pf2012_nh(){
    return(GetHWW().els_iso04_pf2012_nh());
  }
  std::vector<float>          &els_iso03_pf2012_ch(){
    return(GetHWW().els_iso03_pf2012_ch());
  }
  std::vector<float>          &els_iso03_pf2012_em(){
    return(GetHWW().els_iso03_pf2012_em());
  }
  std::vector<float>          &els_iso03_pf2012_nh(){
    return(GetHWW().els_iso03_pf2012_nh());
  }
  std::vector<float>          &els_ecalIso04(){
    return(GetHWW().els_ecalIso04());
  }
  std::vector<float>          &els_hcalIso04(){
    return(GetHWW().els_hcalIso04());
  }
  std::vector<int>            &els_nSeed(){
    return(GetHWW().els_nSeed());
  }
  std::vector<int>            &els_scindex(){
    return(GetHWW().els_scindex());
  }
  std::vector<int>            &els_charge(){
    return(GetHWW().els_charge());
  }
  std::vector<int>            &els_gsftrkidx(){
    return(GetHWW().els_gsftrkidx());
  }
  std::vector<int>            &els_exp_innerlayers(){
    return(GetHWW().els_exp_innerlayers());
  }
  std::vector<int>            &els_trkidx(){
    return(GetHWW().els_trkidx());
  }
  std::vector<int>            &els_type(){
    return(GetHWW().els_type());
  }
  std::vector<int>            &els_fiduciality(){
    return(GetHWW().els_fiduciality());
  }
  std::vector<int>            &els_sccharge(){
    return(GetHWW().els_sccharge());
  }
  std::vector<int>            &els_trk_charge(){
    return(GetHWW().els_trk_charge());
  }
  std::vector<int>            &els_closestMuon(){
    return(GetHWW().els_closestMuon());
  }

  //muons
  std::vector<LorentzVector>  &mus_p4(){
    return(GetHWW().mus_p4());
  }
  std::vector<LorentzVector>  &mus_trk_p4(){
    return(GetHWW().mus_trk_p4());
  }
  std::vector<LorentzVector>  &mus_vertex_p4(){
    return(GetHWW().mus_vertex_p4());
  }
  std::vector<LorentzVector>  &mus_sta_p4(){
    return(GetHWW().mus_sta_p4());
  }
  std::vector<float>          &mus_gfit_chi2(){
    return(GetHWW().mus_gfit_chi2());
  }
  std::vector<float>          &mus_gfit_ndof(){
    return(GetHWW().mus_gfit_ndof());
  }
  std::vector<float>          &mus_ptErr(){
    return(GetHWW().mus_ptErr());
  }
  std::vector<float>          &mus_trkKink(){
    return(GetHWW().mus_trkKink());
  }
  std::vector<float>          &mus_d0corr(){
    return(GetHWW().mus_d0corr());
  }
  std::vector<float>          &mus_d0(){
    return(GetHWW().mus_d0());
  }
  std::vector<float>          &mus_z0corr(){
    return(GetHWW().mus_z0corr());
  }
  std::vector<float>          &mus_chi2(){
    return(GetHWW().mus_chi2());
  }
  std::vector<float>          &mus_ndof(){
    return(GetHWW().mus_ndof());
  }
  std::vector<float>          &mus_ip3d(){
    return(GetHWW().mus_ip3d());
  }
  std::vector<float>          &mus_ip3derr(){
    return(GetHWW().mus_ip3derr());
  }
  std::vector<float>          &mus_segmCompatibility(){
    return(GetHWW().mus_segmCompatibility());
  }
  std::vector<float>          &mus_caloCompatibility(){
    return(GetHWW().mus_caloCompatibility());
  }
  std::vector<float>          &mus_e_had(){
    return(GetHWW().mus_e_had());
  }
  std::vector<float>          &mus_e_ho(){
    return(GetHWW().mus_e_ho());
  }
  std::vector<float>          &mus_e_em(){
    return(GetHWW().mus_e_em());
  }
  std::vector<float>          &mus_e_hadS9(){
    return(GetHWW().mus_e_hadS9());
  }
  std::vector<float>          &mus_e_hoS9(){
    return(GetHWW().mus_e_hoS9());
  }
  std::vector<float>          &mus_e_emS9(){
    return(GetHWW().mus_e_emS9());
  }
  std::vector<float>          &mus_iso03_sumPt(){
    return(GetHWW().mus_iso03_sumPt());
  }
  std::vector<float>          &mus_iso03_emEt(){
    return(GetHWW().mus_iso03_emEt());
  }
  std::vector<float>          &mus_iso03_hadEt(){
    return(GetHWW().mus_iso03_hadEt());
  }
  std::vector<float>          &mus_iso05_sumPt(){
    return(GetHWW().mus_iso05_sumPt());
  }
  std::vector<float>          &mus_iso05_emEt(){
    return(GetHWW().mus_iso05_emEt());
  }
  std::vector<float>          &mus_iso05_hadEt(){
    return(GetHWW().mus_iso05_hadEt());
  }
  std::vector<float>          &mus_sta_d0(){
    return(GetHWW().mus_sta_d0());
  }
  std::vector<float>          &mus_sta_z0corr(){
    return(GetHWW().mus_sta_z0corr());
  }
  std::vector<float>          &mus_isoR03_pf_ChargedHadronPt(){
    return(GetHWW().mus_isoR03_pf_ChargedHadronPt());
  }
  std::vector<float>          &mus_isoR03_pf_NeutralHadronEt(){
    return(GetHWW().mus_isoR03_pf_NeutralHadronEt());
  }
  std::vector<float>          &mus_isoR03_pf_PhotonEt(){
    return(GetHWW().mus_isoR03_pf_PhotonEt());
  }
  std::vector<float>          &mus_isoR03_pf_PUPt(){
    return(GetHWW().mus_isoR03_pf_PUPt());
  }
  std::vector<float>          &mus_iso_ecalvetoDep(){
    return(GetHWW().mus_iso_ecalvetoDep());
  }
  std::vector<float>          &mus_iso_hcalvetoDep(){
    return(GetHWW().mus_iso_hcalvetoDep());
  }
  std::vector<int>            &mus_gfit_validSTAHits(){
    return(GetHWW().mus_gfit_validSTAHits());
  }
  std::vector<int>            &mus_numberOfMatchedStations(){
    return(GetHWW().mus_numberOfMatchedStations());
  }
  std::vector<int>            &mus_pfmusidx(){
    return(GetHWW().mus_pfmusidx());
  }
  std::vector<int>            &mus_charge(){
    return(GetHWW().mus_charge());
  }
  std::vector<int>            &mus_validHits(){
    return(GetHWW().mus_validHits());
  }
  std::vector<int>            &mus_trkidx(){
    return(GetHWW().mus_trkidx());
  }
  std::vector<int>            &mus_pid_PFMuon(){
    return(GetHWW().mus_pid_PFMuon());
  }
  std::vector<int>            &mus_pid_TMLastStationTight(){
    return(GetHWW().mus_pid_TMLastStationTight());
  }
  std::vector<int>            &mus_nmatches(){
    return(GetHWW().mus_nmatches());
  }
  std::vector<int>            &mus_goodmask(){
    return(GetHWW().mus_goodmask());
  }
  std::vector<int>            &mus_type(){
    return(GetHWW().mus_type());
  }

  //dilepton hypothesis
  std::vector<vector<LorentzVector> > &hyp_jets_p4(){
    return(GetHWW().hyp_jets_p4());
  }
  std::vector<LorentzVector>  &hyp_p4(){
    return(GetHWW().hyp_p4());
  }
  std::vector<LorentzVector>  &hyp_ll_p4(){
    return(GetHWW().hyp_ll_p4());
  }
  std::vector<LorentzVector>  &hyp_lt_p4(){
    return(GetHWW().hyp_lt_p4());
  }
  std::vector<int>            &hyp_ll_index(){
    return(GetHWW().hyp_ll_index());
  }
  std::vector<int>            &hyp_lt_index(){
    return(GetHWW().hyp_lt_index());
  }
  std::vector<int>            &hyp_ll_id(){
    return(GetHWW().hyp_ll_id());
  }
  std::vector<int>            &hyp_lt_id(){
    return(GetHWW().hyp_lt_id());
  }
  std::vector<int>            &hyp_ll_charge(){
    return(GetHWW().hyp_ll_charge());
  }
  std::vector<int>            &hyp_lt_charge(){
    return(GetHWW().hyp_lt_charge());
  }
  std::vector<int>            &hyp_type(){
    return(GetHWW().hyp_type());
  }

  //event variables
  unsigned int                &evt_run(){
    return(GetHWW().evt_run());
  }
  unsigned int                &evt_lumiBlock(){
    return(GetHWW().evt_lumiBlock());
  }
  unsigned int                &evt_event(){
    return(GetHWW().evt_event());
  }
  int                         &evt_isRealData(){
    return(GetHWW().evt_isRealData());
  }
  float                       &evt_ww_rho_vor(){
    return(GetHWW().evt_ww_rho_vor());
  }
  float                       &evt_ww_rho(){
    return(GetHWW().evt_ww_rho());
  }
  float                       &evt_rho(){
    return(GetHWW().evt_rho());
  }
  float                       &evt_kt6pf_foregiso_rho(){
    return(GetHWW().evt_kt6pf_foregiso_rho());
  }
  float                       &evt_pfmet(){
    return(GetHWW().evt_pfmet());
  }
  float                       &evt_pfmetPhi(){
    return(GetHWW().evt_pfmetPhi());
  }


  std::vector<float>          &convs_ndof(){
    return(GetHWW().convs_ndof());
  }
  std::vector<float>          &convs_chi2(){
    return(GetHWW().convs_chi2());
  }
  std::vector<float>          &convs_dl(){
    return(GetHWW().convs_dl());
  }
  std::vector<int>            &convs_isConverted(){
    return(GetHWW().convs_isConverted());
  }
  std::vector<vector<int>>    &convs_tkalgo(){
    return(GetHWW().convs_tkalgo());
  }
  std::vector<vector<int>>    &convs_tkidx(){
    return(GetHWW().convs_tkidx());
  }
  std::vector<vector<int>>    &convs_nHitsBeforeVtx(){
    return(GetHWW().convs_nHitsBeforeVtx());
  }
  std::vector<int>            &convs_quality(){
    return(GetHWW().convs_quality());
  }
  std::vector<float>          &scs_sigmaIEtaIPhi(){
    return(GetHWW().scs_sigmaIEtaIPhi());
  }
  std::vector<LorentzVector>  &scs_pos_p4(){
    return(GetHWW().scs_pos_p4());
  }
  std::vector<LorentzVector>  &gsftrks_p4(){
    return(GetHWW().gsftrks_p4());
  }
  std::vector<LorentzVector>  &gsftrks_vertex_p4(){
    return(GetHWW().gsftrks_vertex_p4());
  }
  std::vector<float>          &gsftrks_d0(){
    return(GetHWW().gsftrks_d0());
  }
  std::vector<float>          &gsftrks_d0Err(){
    return(GetHWW().gsftrks_d0Err());
  }
  std::vector<float>          &gsftrks_phiErr(){
    return(GetHWW().gsftrks_phiErr());
  }
  std::vector<float>          &gsftrks_d0phiCov(){
    return(GetHWW().gsftrks_d0phiCov());
  }
  std::vector<float>          &gsftrks_z0Err(){
    return(GetHWW().gsftrks_z0Err());
  }
  std::vector<float>          &gsftrks_z0(){
    return(GetHWW().gsftrks_z0());
  }
  std::vector<float>          &gsftrks_etaErr(){
    return(GetHWW().gsftrks_etaErr());
  }
  std::vector<LorentzVector>  &pfcands_p4(){
    return(GetHWW().pfcands_p4());
  }
  std::vector<int>            &pfcands_trkidx(){
    return(GetHWW().pfcands_trkidx());
  }
  std::vector<int>            &pfcands_particleId(){
    return(GetHWW().pfcands_particleId());
  }
  std::vector<int>            &pfcands_pfelsidx(){
    return(GetHWW().pfcands_pfelsidx());
  }
  std::vector<int>            &pfcands_vtxidx(){
    return(GetHWW().pfcands_vtxidx());
  }
  std::vector<int>            &pfcands_charge(){
    return(GetHWW().pfcands_charge());
  }
  std::vector<int>            &pfels_elsidx(){
    return(GetHWW().pfels_elsidx());
  }
  std::vector<LorentzVector>  &pfels_p4(){
    return(GetHWW().pfels_p4());
  }
  std::vector<LorentzVector>  &pfmus_p4(){
    return(GetHWW().pfmus_p4());
  }
  vector<float>               &trk_met(){
    return(GetHWW().trk_met());
  }
  vector<float>               &trk_metPhi(){
    return(GetHWW().trk_metPhi());
  }
  vector<LorentzVector>       &pfjets_p4(){
    return(GetHWW().pfjets_p4());
  }
  vector<LorentzVector>       &pfjets_corr_p4(){
    return(GetHWW().pfjets_corr_p4());
  }
  vector<float>               &pfjets_area(){
    return(GetHWW().pfjets_area());
  }
  vector<float>               &pfjets_JEC(){
    return(GetHWW().pfjets_JEC());
  }
  vector<float>               &pfjets_mvavalue(){
    return(GetHWW().pfjets_mvavalue());
  }
  vector<float>               &pfjets_trackCountingHighEffBJetTag(){
    return(GetHWW().pfjets_trackCountingHighEffBJetTag());
  }





//isLoaded
  void Load_vtxs_position(){
    GetHWW().vtxs_position_isLoaded = true;
  }
  void Load_vtxs_ndof(){
    GetHWW().vtxs_ndof_isLoaded = true;
  }
  void Load_vtxs_sumpt(){
    GetHWW().vtxs_sumpt_isLoaded = true;
  }
  void Load_vtxs_isFake(){
    GetHWW().vtxs_isFake_isLoaded = true;
  }
  void Load_vtxs_xError(){
    GetHWW().vtxs_xError_isLoaded = true;
  }
  void Load_vtxs_yError(){
    GetHWW().vtxs_yError_isLoaded = true;
  }
  void Load_vtxs_zError(){
    GetHWW().vtxs_zError_isLoaded = true;
  }
  void Load_vtxs_covMatrix(){
    GetHWW().vtxs_covMatrix_isLoaded = true;
  }
  void Load_trks_trk_p4(){
    GetHWW().trks_trk_p4_isLoaded = true;
  }
  void Load_trks_chi2(){
    GetHWW().trks_chi2_isLoaded = true;
  }
  void Load_trks_ndof(){
    GetHWW().trks_ndof_isLoaded = true;
  }
  void Load_trks_d0(){
    GetHWW().trks_d0_isLoaded = true;
  }
  void Load_trks_nlayers(){
    GetHWW().trks_nlayers_isLoaded = true;
  }
  void Load_trks_valid_pixelhits(){
    GetHWW().trks_valid_pixelhits_isLoaded = true;
  }
  void Load_trks_z0(){
    GetHWW().trks_z0_isLoaded = true;
  }
  void Load_trks_z0Err(){
    GetHWW().trks_z0Err_isLoaded = true;
  }
  void Load_trks_etaErr(){
    GetHWW().trks_etaErr_isLoaded = true;
  }
  void Load_trks_d0Err(){
    GetHWW().trks_d0Err_isLoaded = true;
  }
  void Load_trks_phiErr(){
    GetHWW().trks_phiErr_isLoaded = true;
  }
  void Load_trks_d0phiCov(){
    GetHWW().trks_d0phiCov_isLoaded = true;
  }
  void Load_trks_qualityMask(){
    GetHWW().trks_qualityMask_isLoaded = true;
  }
  void Load_trks_charge(){
    GetHWW().trks_charge_isLoaded = true;
  }
  void Load_els_p4(){
    GetHWW().els_p4_isLoaded = true;
  }
  void Load_els_trk_p4(){
    GetHWW().els_trk_p4_isLoaded = true;
  }
  void Load_els_vertex_p4(){
    GetHWW().els_vertex_p4_isLoaded = true;
  }
  void Load_els_lh(){
    GetHWW().els_lh_isLoaded = true;
  }
  void Load_els_etaSC(){
    GetHWW().els_etaSC_isLoaded = true;
  }
  void Load_els_sigmaIEtaIEta(){
    GetHWW().els_sigmaIEtaIEta_isLoaded = true;
  }
  void Load_els_dEtaIn(){
    GetHWW().els_dEtaIn_isLoaded = true;
  }
  void Load_els_dPhiIn(){
    GetHWW().els_dPhiIn_isLoaded = true;
  }
  void Load_els_hOverE(){
    GetHWW().els_hOverE_isLoaded = true;
  }
  void Load_els_tkIso(){
    GetHWW().els_tkIso_isLoaded = true;
  }
  void Load_els_d0corr(){
    GetHWW().els_d0corr_isLoaded = true;
  }
  void Load_els_d0(){
    GetHWW().els_d0_isLoaded = true;
  }
  void Load_els_z0corr(){
    GetHWW().els_z0corr_isLoaded = true;
  }
  void Load_els_fbrem(){
    GetHWW().els_fbrem_isLoaded = true;
  }
  void Load_els_eOverPIn(){
    GetHWW().els_eOverPIn_isLoaded = true;
  }
  void Load_els_eSeedOverPOut(){
    GetHWW().els_eSeedOverPOut_isLoaded = true;
  }
  void Load_els_eSeedOverPIn(){
    GetHWW().els_eSeedOverPIn_isLoaded = true;
  }
  void Load_els_sigmaIPhiIPhi(){
    GetHWW().els_sigmaIPhiIPhi_isLoaded = true;
  }
  void Load_els_eSC(){
    GetHWW().els_eSC_isLoaded = true;
  }
  void Load_els_ip3d(){
    GetHWW().els_ip3d_isLoaded = true;
  }
  void Load_els_ip3derr(){
    GetHWW().els_ip3derr_isLoaded = true;
  }
  void Load_els_chi2(){
    GetHWW().els_chi2_isLoaded = true;
  }
  void Load_els_ndof(){
    GetHWW().els_ndof_isLoaded = true;
  }
  void Load_els_dEtaOut(){
    GetHWW().els_dEtaOut_isLoaded = true;
  }
  void Load_els_dPhiOut(){
    GetHWW().els_dPhiOut_isLoaded = true;
  }
  void Load_els_eSCRaw(){
    GetHWW().els_eSCRaw_isLoaded = true;
  }
  void Load_els_etaSCwidth(){
    GetHWW().els_etaSCwidth_isLoaded = true;
  }
  void Load_els_phiSCwidth(){
    GetHWW().els_phiSCwidth_isLoaded = true;
  }
  void Load_els_eSCPresh(){
    GetHWW().els_eSCPresh_isLoaded = true;
  }
  void Load_els_iso03_pf_ch(){
    GetHWW().els_iso03_pf_ch_isLoaded = true;
  }
  void Load_els_iso03_pf_nhad05(){
    GetHWW().els_iso03_pf_nhad05_isLoaded = true;
  }
  void Load_els_iso03_pf_gamma05(){
    GetHWW().els_iso03_pf_gamma05_isLoaded = true;
  }
  void Load_els_iso04_pf_ch(){
    GetHWW().els_iso04_pf_ch_isLoaded = true;
  }
  void Load_els_iso04_pf_nhad05(){
    GetHWW().els_iso04_pf_nhad05_isLoaded = true;
  }
  void Load_els_iso04_pf_gamma05(){
    GetHWW().els_iso04_pf_gamma05_isLoaded = true;
  }
  void Load_els_e5x5(){
    GetHWW().els_e5x5_isLoaded = true;
  }
  void Load_els_e1x5(){
    GetHWW().els_e1x5_isLoaded = true;
  }
  void Load_els_e3x3(){
    GetHWW().els_e3x3_isLoaded = true;
  }
  void Load_els_ecalEnergy(){
    GetHWW().els_ecalEnergy_isLoaded = true;
  }
  void Load_els_eOverPOut(){
    GetHWW().els_eOverPOut_isLoaded = true;
  }
  void Load_els_ecalIso(){
    GetHWW().els_ecalIso_isLoaded = true;
  }
  void Load_els_hcalIso(){
    GetHWW().els_hcalIso_isLoaded = true;
  }
  void Load_els_trkshFrac(){
    GetHWW().els_trkshFrac_isLoaded = true;
  }
  void Load_els_conv_dist(){
    GetHWW().els_conv_dist_isLoaded = true;
  }
  void Load_els_conv_dcot(){
    GetHWW().els_conv_dcot_isLoaded = true;
  }
  void Load_els_conv_old_dist(){
    GetHWW().els_conv_old_dist_isLoaded = true;
  }
  void Load_els_conv_old_dcot(){
    GetHWW().els_conv_old_dcot_isLoaded = true;
  }
  void Load_els_iso04_pf2012_ch(){
    GetHWW().els_iso04_pf2012_ch_isLoaded = true;
  }
  void Load_els_iso04_pf2012_em(){
    GetHWW().els_iso04_pf2012_em_isLoaded = true;
  }
  void Load_els_iso04_pf2012_nh(){
    GetHWW().els_iso04_pf2012_nh_isLoaded = true;
  }
  void Load_els_iso03_pf2012_ch(){
    GetHWW().els_iso03_pf2012_ch_isLoaded = true;
  }
  void Load_els_iso03_pf2012_em(){
    GetHWW().els_iso03_pf2012_em_isLoaded = true;
  }
  void Load_els_iso03_pf2012_nh(){
    GetHWW().els_iso03_pf2012_nh_isLoaded = true;
  }
  void Load_els_ecalIso04(){
    GetHWW().els_ecalIso04_isLoaded = true;
  }
  void Load_els_hcalIso04(){
    GetHWW().els_hcalIso04_isLoaded = true;
  }
  void Load_els_nSeed(){
    GetHWW().els_nSeed_isLoaded = true;
  }
  void Load_els_scindex(){
    GetHWW().els_scindex_isLoaded = true;
  }
  void Load_els_charge(){
    GetHWW().els_charge_isLoaded = true;
  }
  void Load_els_gsftrkidx(){
    GetHWW().els_gsftrkidx_isLoaded = true;
  }
  void Load_els_exp_innerlayers(){
    GetHWW().els_exp_innerlayers_isLoaded = true;
  }
  void Load_els_trkidx(){
    GetHWW().els_trkidx_isLoaded = true;
  }
  void Load_els_type(){
    GetHWW().els_type_isLoaded = true;
  }
  void Load_els_fiduciality(){
    GetHWW().els_fiduciality_isLoaded = true;
  }
  void Load_els_sccharge(){
    GetHWW().els_sccharge_isLoaded = true;
  }
  void Load_els_trk_charge(){
    GetHWW().els_trk_charge_isLoaded = true;
  }
  void Load_els_closestMuon(){
    GetHWW().els_closestMuon_isLoaded = true;
  }
  void Load_mus_p4(){
    GetHWW().mus_p4_isLoaded = true;
  }
  void Load_mus_trk_p4(){
    GetHWW().mus_trk_p4_isLoaded = true;
  }
  void Load_mus_vertex_p4(){
    GetHWW().mus_vertex_p4_isLoaded = true;
  }
  void Load_mus_sta_p4(){
    GetHWW().mus_sta_p4_isLoaded = true;
  }
  void Load_mus_gfit_chi2(){
    GetHWW().mus_gfit_chi2_isLoaded = true;
  }
  void Load_mus_gfit_ndof(){
    GetHWW().mus_gfit_ndof_isLoaded = true;
  }
  void Load_mus_ptErr(){
    GetHWW().mus_ptErr_isLoaded = true;
  }
  void Load_mus_trkKink(){
    GetHWW().mus_trkKink_isLoaded = true;
  }
  void Load_mus_d0corr(){
    GetHWW().mus_d0corr_isLoaded = true;
  }
  void Load_mus_d0(){
    GetHWW().mus_d0_isLoaded = true;
  }
  void Load_mus_z0corr(){
    GetHWW().mus_z0corr_isLoaded = true;
  }
  void Load_mus_chi2(){
    GetHWW().mus_chi2_isLoaded = true;
  }
  void Load_mus_ndof(){
    GetHWW().mus_ndof_isLoaded = true;
  }
  void Load_mus_ip3d(){
    GetHWW().mus_ip3d_isLoaded = true;
  }
  void Load_mus_ip3derr(){
    GetHWW().mus_ip3derr_isLoaded = true;
  }
  void Load_mus_segmCompatibility(){
    GetHWW().mus_segmCompatibility_isLoaded = true;
  }
  void Load_mus_caloCompatibility(){
    GetHWW().mus_caloCompatibility_isLoaded = true;
  }
  void Load_mus_e_had(){
    GetHWW().mus_e_had_isLoaded = true;
  }
  void Load_mus_e_ho(){
    GetHWW().mus_e_ho_isLoaded = true;
  }
  void Load_mus_e_em(){
    GetHWW().mus_e_em_isLoaded = true;
  }
  void Load_mus_e_hadS9(){
    GetHWW().mus_e_hadS9_isLoaded = true;
  }
  void Load_mus_e_hoS9(){
    GetHWW().mus_e_hoS9_isLoaded = true;
  }
  void Load_mus_e_emS9(){
    GetHWW().mus_e_emS9_isLoaded = true;
  }
  void Load_mus_iso03_sumPt(){
    GetHWW().mus_iso03_sumPt_isLoaded = true;
  }
  void Load_mus_iso03_emEt(){
    GetHWW().mus_iso03_emEt_isLoaded = true;
  }
  void Load_mus_iso03_hadEt(){
    GetHWW().mus_iso03_hadEt_isLoaded = true;
  }
  void Load_mus_iso05_sumPt(){
    GetHWW().mus_iso05_sumPt_isLoaded = true;
  }
  void Load_mus_iso05_emEt(){
    GetHWW().mus_iso05_emEt_isLoaded = true;
  }
  void Load_mus_iso05_hadEt(){
    GetHWW().mus_iso05_hadEt_isLoaded = true;
  }
  void Load_mus_sta_d0(){
    GetHWW().mus_sta_d0_isLoaded = true;
  }
  void Load_mus_sta_z0corr(){
    GetHWW().mus_sta_z0corr_isLoaded = true;
  }
  void Load_mus_isoR03_pf_ChargedHadronPt(){
    GetHWW().mus_isoR03_pf_ChargedHadronPt_isLoaded = true;
  }
  void Load_mus_isoR03_pf_NeutralHadronEt(){
    GetHWW().mus_isoR03_pf_NeutralHadronEt_isLoaded = true;
  }
  void Load_mus_isoR03_pf_PhotonEt(){
    GetHWW().mus_isoR03_pf_PhotonEt_isLoaded = true;
  }
  void Load_mus_isoR03_pf_PUPt(){
    GetHWW().mus_isoR03_pf_PUPt_isLoaded = true;
  }
  void Load_mus_iso_ecalvetoDep(){
    GetHWW().mus_iso_ecalvetoDep_isLoaded = true;
  }
  void Load_mus_iso_hcalvetoDep(){
    GetHWW().mus_iso_hcalvetoDep_isLoaded = true;
  }
  void Load_mus_gfit_validSTAHits(){
    GetHWW().mus_gfit_validSTAHits_isLoaded = true;
  }
  void Load_mus_numberOfMatchedStations(){
    GetHWW().mus_numberOfMatchedStations_isLoaded = true;
  }
  void Load_mus_pfmusidx(){
    GetHWW().mus_pfmusidx_isLoaded = true;
  }
  void Load_mus_charge(){
    GetHWW().mus_charge_isLoaded = true;
  }
  void Load_mus_validHits(){
    GetHWW().mus_validHits_isLoaded = true;
  }
  void Load_mus_trkidx(){
    GetHWW().mus_trkidx_isLoaded = true;
  }
  void Load_mus_pid_PFMuon(){
    GetHWW().mus_pid_PFMuon_isLoaded = true;
  }
  void Load_mus_pid_TMLastStationTight(){
    GetHWW().mus_pid_TMLastStationTight_isLoaded = true;
  }
  void Load_mus_nmatches(){
    GetHWW().mus_nmatches_isLoaded = true;
  }
  void Load_mus_goodmask(){
    GetHWW().mus_goodmask_isLoaded = true;
  }
  void Load_mus_type(){
    GetHWW().mus_type_isLoaded = true;
  }
  void Load_hyp_jets_p4(){
    GetHWW().hyp_jets_p4_isLoaded = true;
  }
  void Load_hyp_p4(){
    GetHWW().hyp_p4_isLoaded = true;
  }
  void Load_hyp_ll_p4(){
    GetHWW().hyp_ll_p4_isLoaded = true;
  }
  void Load_hyp_lt_p4(){
    GetHWW().hyp_lt_p4_isLoaded = true;
  }
  void Load_hyp_ll_index(){
    GetHWW().hyp_ll_index_isLoaded = true;
  }
  void Load_hyp_lt_index(){
    GetHWW().hyp_lt_index_isLoaded = true;
  }
  void Load_hyp_ll_id(){
    GetHWW().hyp_ll_id_isLoaded = true;
  }
  void Load_hyp_lt_id(){
    GetHWW().hyp_lt_id_isLoaded = true;
  }
  void Load_hyp_ll_charge(){
    GetHWW().hyp_ll_charge_isLoaded = true;
  }
  void Load_hyp_lt_charge(){
    GetHWW().hyp_lt_charge_isLoaded = true;
  }
  void Load_hyp_type(){
    GetHWW().hyp_type_isLoaded = true;
  }
  void Load_evt_run(){
    GetHWW().evt_run_isLoaded = true;
  }
  void Load_evt_lumiBlock(){
    GetHWW().evt_lumiBlock_isLoaded = true;
  }
  void Load_evt_event(){
    GetHWW().evt_event_isLoaded = true;
  }
  void Load_evt_isRealData(){
    GetHWW().evt_isRealData_isLoaded = true;
  }
  void Load_evt_ww_rho_vor(){
    GetHWW().evt_ww_rho_vor_isLoaded = true;
  }
  void Load_evt_ww_rho(){
    GetHWW().evt_ww_rho_isLoaded = true;
  }
  void Load_evt_rho(){
    GetHWW().evt_rho_isLoaded = true;
  }
  void Load_evt_kt6pf_foregiso_rho(){
    GetHWW().evt_kt6pf_foregiso_rho_isLoaded = true;
  }
  void Load_evt_pfmet(){
    GetHWW().evt_pfmet_isLoaded = true;
  }
  void Load_evt_pfmetPhi(){
    GetHWW().evt_pfmetPhi_isLoaded = true;
  }
  void Load_convs_ndof(){
    GetHWW().convs_ndof_isLoaded = true;
  }
  void Load_convs_chi2(){
    GetHWW().convs_chi2_isLoaded = true;
  }
  void Load_convs_dl(){
    GetHWW().convs_dl_isLoaded = true;
  }
  void Load_convs_isConverted(){
    GetHWW().convs_isConverted_isLoaded = true;
  }
  void Load_convs_tkalgo(){
    GetHWW().convs_tkalgo_isLoaded = true;
  }
  void Load_convs_tkidx(){
    GetHWW().convs_tkidx_isLoaded = true;
  }
  void Load_convs_nHitsBeforeVtx(){
    GetHWW().convs_nHitsBeforeVtx_isLoaded = true;
  }
  void Load_convs_quality(){
    GetHWW().convs_quality_isLoaded = true;
  }
  void Load_scs_sigmaIEtaIPhi(){
    GetHWW().scs_sigmaIEtaIPhi_isLoaded = true;
  }
  void Load_scs_pos_p4(){
    GetHWW().scs_pos_p4_isLoaded = true;
  }
  void Load_gsftrks_p4(){
    GetHWW().gsftrks_p4_isLoaded = true;
  }
  void Load_gsftrks_vertex_p4(){
    GetHWW().gsftrks_vertex_p4_isLoaded = true;
  }
  void Load_gsftrks_d0(){
    GetHWW().gsftrks_d0_isLoaded = true;
  }
  void Load_gsftrks_d0Err(){
    GetHWW().gsftrks_d0Err_isLoaded = true;
  }
  void Load_gsftrks_phiErr(){
    GetHWW().gsftrks_phiErr_isLoaded = true;
  }
  void Load_gsftrks_d0phiCov(){
    GetHWW().gsftrks_d0phiCov_isLoaded = true;
  }
  void Load_gsftrks_z0Err(){
    GetHWW().gsftrks_z0Err_isLoaded = true;
  }
  void Load_gsftrks_z0(){
    GetHWW().gsftrks_z0_isLoaded = true;
  }
  void Load_gsftrks_etaErr(){
    GetHWW().gsftrks_etaErr_isLoaded = true;
  }
  void Load_pfcands_p4(){
    GetHWW().pfcands_p4_isLoaded = true;
  }
  void Load_pfcands_trkidx(){
    GetHWW().pfcands_trkidx_isLoaded = true;
  }
  void Load_pfcands_particleId(){
    GetHWW().pfcands_particleId_isLoaded = true;
  }
  void Load_pfcands_pfelsidx(){
    GetHWW().pfcands_pfelsidx_isLoaded = true;
  }
  void Load_pfcands_vtxidx(){
    GetHWW().pfcands_vtxidx_isLoaded = true;
  }
  void Load_pfcands_charge(){
    GetHWW().pfcands_charge_isLoaded = true;
  }
  void Load_pfels_elsidx(){
    GetHWW().pfels_elsidx_isLoaded = true;
  }
  void Load_pfels_p4(){
    GetHWW().pfels_p4_isLoaded = true;
  }
  void Load_pfmus_p4(){
    GetHWW().pfmus_p4_isLoaded = true;
  }
  void Load_trk_met(){
    GetHWW().trk_met_isLoaded = true;
  }
  void Load_trk_metPhi(){
    GetHWW().trk_metPhi_isLoaded = true;
  }
  void Load_pfjets_p4(){
    GetHWW().pfjets_p4_isLoaded = true;
  }
  void Load_pfjets_corr_p4(){
    GetHWW().pfjets_corr_p4_isLoaded = true;
  }
  void Load_pfjets_area(){
    GetHWW().pfjets_area_isLoaded = true;
  }
  void Load_pfjets_JEC(){
    GetHWW().pfjets_JEC_isLoaded = true;
  }
  void Load_pfjets_mvavalue(){
    GetHWW().pfjets_mvavalue_isLoaded = true;
  }
  void Load_pfjets_trackCountingHighEffBJetTag(){
    GetHWW().pfjets_trackCountingHighEffBJetTag_isLoaded = true;
  }


  void Reset() {GetHWW().Reset();}

}//end namespace


void HWW::Reset(){

  vtxs_position_.clear();
  vtxs_ndof_.clear();
  vtxs_sumpt_.clear();
  vtxs_isFake_.clear();
  vtxs_xError_.clear();
  vtxs_yError_.clear();
  vtxs_zError_.clear();
  vtxs_covMatrix_.clear();

  
  trks_trk_p4_.clear();
  trks_chi2_.clear();
  trks_ndof_.clear();
  trks_d0_.clear();
  trks_nlayers_.clear();
  trks_valid_pixelhits_.clear();
  trks_z0_.clear();
  trks_z0Err_.clear();
  trks_etaErr_.clear();
  trks_d0Err_.clear();
  trks_phiErr_.clear();
  trks_d0phiCov_.clear();
  trks_qualityMask_.clear();
  trks_charge_.clear();

  
  els_p4_.clear();
  els_trk_p4_.clear();
  els_vertex_p4_.clear();
  els_lh_.clear();
  els_etaSC_.clear();
  els_sigmaIEtaIEta_.clear();
  els_dEtaIn_.clear();
  els_dPhiIn_.clear();
  els_hOverE_.clear();
  els_tkIso_.clear();
  els_d0corr_.clear();
  els_d0_.clear();
  els_z0corr_.clear();
  els_fbrem_.clear();
  els_eOverPIn_.clear();
  els_eSeedOverPOut_.clear();
  els_eSeedOverPIn_.clear();
  els_sigmaIPhiIPhi_.clear();
  els_eSC_.clear();
  els_ip3d_.clear();
  els_ip3derr_.clear();
  els_chi2_.clear();
  els_ndof_.clear();
  els_dEtaOut_.clear();
  els_dPhiOut_.clear();
  els_eSCRaw_.clear();
  els_etaSCwidth_.clear();
  els_phiSCwidth_.clear();
  els_eSCPresh_.clear();
  els_iso03_pf_ch_.clear();
  els_iso03_pf_nhad05_.clear();
  els_iso03_pf_gamma05_.clear();
  els_iso04_pf_ch_.clear();
  els_iso04_pf_nhad05_.clear();
  els_iso04_pf_gamma05_.clear();
  els_e5x5_.clear();
  els_e1x5_.clear();
  els_e3x3_.clear();
  els_ecalEnergy_.clear();
  els_eOverPOut_.clear();
  els_ecalIso_.clear();
  els_hcalIso_.clear();
  els_trkshFrac_.clear();
  els_conv_dist_.clear();
  els_conv_dcot_.clear();
  els_conv_old_dist_.clear();
  els_conv_old_dcot_.clear();
  els_iso04_pf2012_ch_.clear();
  els_iso04_pf2012_em_.clear();
  els_iso04_pf2012_nh_.clear();
  els_iso03_pf2012_ch_.clear();
  els_iso03_pf2012_em_.clear();
  els_iso03_pf2012_nh_.clear();
  els_ecalIso04_.clear();
  els_hcalIso04_.clear();
  els_nSeed_.clear();
  els_scindex_.clear();
  els_charge_.clear();
  els_gsftrkidx_.clear();
  els_exp_innerlayers_.clear();
  els_trkidx_.clear();
  els_type_.clear();
  els_fiduciality_.clear();
  els_sccharge_.clear();
  els_trk_charge_.clear();
  els_closestMuon_.clear();

  
  mus_p4_.clear();
  mus_trk_p4_.clear();
  mus_vertex_p4_.clear();
  mus_sta_p4_.clear();
  mus_gfit_chi2_.clear();
  mus_gfit_ndof_.clear();
  mus_ptErr_.clear();
  mus_trkKink_.clear();
  mus_d0corr_.clear();
  mus_d0_.clear();
  mus_z0corr_.clear();
  mus_chi2_.clear();
  mus_ndof_.clear();
  mus_ip3d_.clear();
  mus_ip3derr_.clear();
  mus_segmCompatibility_.clear();
  mus_caloCompatibility_.clear();
  mus_e_had_.clear();
  mus_e_ho_.clear();
  mus_e_em_.clear();
  mus_e_hadS9_.clear();
  mus_e_hoS9_.clear();
  mus_e_emS9_.clear();
  mus_iso03_sumPt_.clear();
  mus_iso03_emEt_.clear();
  mus_iso03_hadEt_.clear();
  mus_iso05_sumPt_.clear();
  mus_iso05_emEt_.clear();
  mus_iso05_hadEt_.clear();
  mus_sta_d0_.clear();
  mus_sta_z0corr_.clear();
  mus_isoR03_pf_ChargedHadronPt_.clear();
  mus_isoR03_pf_NeutralHadronEt_.clear();
  mus_isoR03_pf_PhotonEt_.clear();
  mus_isoR03_pf_PUPt_.clear();
  mus_iso_ecalvetoDep_.clear();
  mus_iso_hcalvetoDep_.clear();
  mus_gfit_validSTAHits_.clear();
  mus_numberOfMatchedStations_.clear();
  mus_pfmusidx_.clear();
  mus_charge_.clear();
  mus_validHits_.clear();
  mus_trkidx_.clear();
  mus_pid_PFMuon_.clear();
  mus_pid_TMLastStationTight_.clear();
  mus_nmatches_.clear();
  mus_goodmask_.clear();
  mus_type_.clear();

  
  hyp_jets_p4_.clear();
  hyp_p4_.clear();
  hyp_ll_p4_.clear();
  hyp_lt_p4_.clear();
  hyp_ll_index_.clear();
  hyp_lt_index_.clear();
  hyp_ll_id_.clear();
  hyp_lt_id_.clear();
  hyp_ll_charge_.clear();
  hyp_lt_charge_.clear();
  hyp_type_.clear();

  
  evt_run_ = 999;
  evt_lumiBlock_ = 999;
  evt_event_ = 999;
  evt_isRealData_ = -999;
  evt_ww_rho_vor_ = -999.0;
  evt_ww_rho_ = -999.0;
  evt_rho_ = -999.0;
  evt_kt6pf_foregiso_rho_ = -999.0;
  evt_pfmet_ = -999.0;
  evt_pfmetPhi_ = -999.0;

  convs_ndof_.clear();
  convs_chi2_.clear();
  convs_dl_.clear();
  convs_isConverted_.clear();
  convs_tkalgo_.clear();
  convs_tkidx_.clear();
  convs_nHitsBeforeVtx_.clear();
  convs_quality_.clear();
  scs_sigmaIEtaIPhi_.clear();
  scs_pos_p4_.clear();
  gsftrks_p4_.clear();
  gsftrks_vertex_p4_.clear();
  gsftrks_d0_.clear();
  gsftrks_d0Err_.clear();
  gsftrks_phiErr_.clear();
  gsftrks_d0phiCov_.clear();
  gsftrks_z0Err_.clear();
  gsftrks_z0_.clear();
  gsftrks_etaErr_.clear();
  pfcands_p4_.clear();
  pfcands_trkidx_.clear();
  pfcands_particleId_.clear();
  pfcands_pfelsidx_.clear();
  pfcands_vtxidx_.clear();
  pfcands_charge_.clear();
  pfels_elsidx_.clear();
  pfels_p4_.clear();
  pfmus_p4_.clear();
  trk_met_.clear();
  trk_metPhi_.clear();
  pfjets_p4_.clear();
  pfjets_corr_p4_.clear();
  pfjets_area_.clear();
  pfjets_JEC_.clear();
  pfjets_mvavalue_.clear();
  pfjets_trackCountingHighEffBJetTag_.clear();



  vtxs_position_isLoaded = false;
  vtxs_ndof_isLoaded = false;
  vtxs_sumpt_isLoaded = false;
  vtxs_isFake_isLoaded = false;
  vtxs_xError_isLoaded = false;
  vtxs_yError_isLoaded = false;
  vtxs_zError_isLoaded = false;
  vtxs_covMatrix_isLoaded = false;

  trks_trk_p4_isLoaded = false;
  trks_chi2_isLoaded = false;
  trks_ndof_isLoaded = false;
  trks_d0_isLoaded = false;
  trks_nlayers_isLoaded = false;
  trks_valid_pixelhits_isLoaded = false;
  trks_z0_isLoaded = false;
  trks_z0Err_isLoaded = false;
  trks_etaErr_isLoaded = false;
  trks_d0Err_isLoaded = false;
  trks_phiErr_isLoaded = false;
  trks_d0phiCov_isLoaded = false;
  trks_qualityMask_isLoaded = false;
  trks_charge_isLoaded = false;

  els_p4_isLoaded = false;
  els_trk_p4_isLoaded = false;
  els_vertex_p4_isLoaded = false;
  els_lh_isLoaded = false;
  els_etaSC_isLoaded = false;
  els_sigmaIEtaIEta_isLoaded = false;
  els_dEtaIn_isLoaded = false;
  els_dPhiIn_isLoaded = false;
  els_hOverE_isLoaded = false;
  els_tkIso_isLoaded = false;
  els_d0corr_isLoaded = false;
  els_d0_isLoaded = false;
  els_z0corr_isLoaded = false;
  els_fbrem_isLoaded = false;
  els_eOverPIn_isLoaded = false;
  els_eSeedOverPOut_isLoaded = false;
  els_eSeedOverPIn_isLoaded = false;
  els_sigmaIPhiIPhi_isLoaded = false;
  els_eSC_isLoaded = false;
  els_ip3d_isLoaded = false;
  els_ip3derr_isLoaded = false;
  els_chi2_isLoaded = false;
  els_ndof_isLoaded = false;
  els_dEtaOut_isLoaded = false;
  els_dPhiOut_isLoaded = false;
  els_eSCRaw_isLoaded = false;
  els_etaSCwidth_isLoaded = false;
  els_phiSCwidth_isLoaded = false;
  els_eSCPresh_isLoaded = false;
  els_iso03_pf_ch_isLoaded = false;
  els_iso03_pf_nhad05_isLoaded = false;
  els_iso03_pf_gamma05_isLoaded = false;
  els_iso04_pf_ch_isLoaded = false;
  els_iso04_pf_nhad05_isLoaded = false;
  els_iso04_pf_gamma05_isLoaded = false;
  els_e5x5_isLoaded = false;
  els_e1x5_isLoaded = false;
  els_e3x3_isLoaded = false;
  els_ecalEnergy_isLoaded = false;
  els_eOverPOut_isLoaded = false;
  els_ecalIso_isLoaded = false;
  els_hcalIso_isLoaded = false;
  els_trkshFrac_isLoaded = false;
  els_conv_dist_isLoaded = false;
  els_conv_dcot_isLoaded = false;
  els_conv_old_dist_isLoaded = false;
  els_conv_old_dcot_isLoaded = false;
  els_iso04_pf2012_ch_isLoaded = false;
  els_iso04_pf2012_em_isLoaded = false;
  els_iso04_pf2012_nh_isLoaded = false;
  els_iso03_pf2012_ch_isLoaded = false;
  els_iso03_pf2012_em_isLoaded = false;
  els_iso03_pf2012_nh_isLoaded = false;
  els_ecalIso04_isLoaded = false;
  els_hcalIso04_isLoaded = false;
  els_nSeed_isLoaded = false;
  els_scindex_isLoaded = false;
  els_charge_isLoaded = false;
  els_gsftrkidx_isLoaded = false;
  els_exp_innerlayers_isLoaded = false;
  els_trkidx_isLoaded = false;
  els_type_isLoaded = false;
  els_fiduciality_isLoaded = false;
  els_sccharge_isLoaded = false;
  els_trk_charge_isLoaded = false;
  els_closestMuon_isLoaded = false;

  mus_p4_isLoaded = false;
  mus_trk_p4_isLoaded = false;
  mus_vertex_p4_isLoaded = false;
  mus_sta_p4_isLoaded = false;
  mus_gfit_chi2_isLoaded = false;
  mus_gfit_ndof_isLoaded = false;
  mus_ptErr_isLoaded = false;
  mus_trkKink_isLoaded = false;
  mus_d0corr_isLoaded = false;
  mus_d0_isLoaded = false;
  mus_z0corr_isLoaded = false;
  mus_chi2_isLoaded = false;
  mus_ndof_isLoaded = false;
  mus_ip3d_isLoaded = false;
  mus_ip3derr_isLoaded = false;
  mus_segmCompatibility_isLoaded = false;
  mus_caloCompatibility_isLoaded = false;
  mus_e_had_isLoaded = false;
  mus_e_ho_isLoaded = false;
  mus_e_em_isLoaded = false;
  mus_e_hadS9_isLoaded = false;
  mus_e_hoS9_isLoaded = false;
  mus_e_emS9_isLoaded = false;
  mus_iso03_sumPt_isLoaded = false;
  mus_iso03_emEt_isLoaded = false;
  mus_iso03_hadEt_isLoaded = false;
  mus_iso05_sumPt_isLoaded = false;
  mus_iso05_emEt_isLoaded = false;
  mus_iso05_hadEt_isLoaded = false;
  mus_sta_d0_isLoaded = false;
  mus_sta_z0corr_isLoaded = false;
  mus_isoR03_pf_ChargedHadronPt_isLoaded = false;
  mus_isoR03_pf_NeutralHadronEt_isLoaded = false;
  mus_isoR03_pf_PhotonEt_isLoaded = false;
  mus_isoR03_pf_PUPt_isLoaded = false;
  mus_iso_ecalvetoDep_isLoaded = false;
  mus_iso_hcalvetoDep_isLoaded = false;
  mus_gfit_validSTAHits_isLoaded = false;
  mus_numberOfMatchedStations_isLoaded = false;
  mus_pfmusidx_isLoaded = false;
  mus_charge_isLoaded = false;
  mus_validHits_isLoaded = false;
  mus_trkidx_isLoaded = false;
  mus_pid_PFMuon_isLoaded = false;
  mus_pid_TMLastStationTight_isLoaded = false;
  mus_nmatches_isLoaded = false;
  mus_goodmask_isLoaded = false;
  mus_type_isLoaded = false;

  hyp_jets_p4_isLoaded = false;
  hyp_p4_isLoaded = false;
  hyp_ll_p4_isLoaded = false;
  hyp_lt_p4_isLoaded = false;
  hyp_ll_index_isLoaded = false;
  hyp_lt_index_isLoaded = false;
  hyp_ll_id_isLoaded = false;
  hyp_lt_id_isLoaded = false;
  hyp_ll_charge_isLoaded = false;
  hyp_lt_charge_isLoaded = false;
  hyp_type_isLoaded = false;

  evt_run_isLoaded = false;
  evt_lumiBlock_isLoaded = false;
  evt_event_isLoaded = false;
  evt_isRealData_isLoaded = false;
  evt_ww_rho_vor_isLoaded = false;
  evt_ww_rho_isLoaded = false;
  evt_rho_isLoaded = false;
  evt_kt6pf_foregiso_rho_isLoaded = false;
  evt_pfmet_isLoaded = false;
  evt_pfmetPhi_isLoaded = false;

  convs_ndof_isLoaded = false;
  convs_chi2_isLoaded = false;
  convs_dl_isLoaded = false;
  convs_isConverted_isLoaded = false;
  convs_tkalgo_isLoaded = false;
  convs_tkidx_isLoaded = false;
  convs_nHitsBeforeVtx_isLoaded = false;
  convs_quality_isLoaded = false;
  scs_sigmaIEtaIPhi_isLoaded = false;
  scs_pos_p4_isLoaded = false;
  gsftrks_p4_isLoaded = false;
  gsftrks_vertex_p4_isLoaded = false;
  gsftrks_d0_isLoaded = false;
  gsftrks_d0Err_isLoaded = false;
  gsftrks_phiErr_isLoaded = false;
  gsftrks_d0phiCov_isLoaded = false;
  gsftrks_z0Err_isLoaded = false;
  gsftrks_z0_isLoaded = false;
  gsftrks_etaErr_isLoaded = false;
  pfcands_p4_isLoaded = false;
  pfcands_trkidx_isLoaded = false;
  pfcands_particleId_isLoaded = false;
  pfcands_pfelsidx_isLoaded = false;
  pfcands_vtxidx_isLoaded = false;
  pfcands_charge_isLoaded = false;
  pfels_elsidx_isLoaded = false;
  pfels_p4_isLoaded = false;
  pfmus_p4_isLoaded = false;
  trk_met_isLoaded = false;
  trk_metPhi_isLoaded = false;
  pfjets_p4_isLoaded = false;
  pfjets_corr_p4_isLoaded = false;
  pfjets_area_isLoaded = false;
  pfjets_JEC_isLoaded = false;
  pfjets_mvavalue_isLoaded = false;
  pfjets_trackCountingHighEffBJetTag_isLoaded = false;

}
