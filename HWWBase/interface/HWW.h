#ifndef HWW_H
#define HWW_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
using namespace std;

class HWW{
public:


  //vertex 
  std::vector<LorentzVector>  vtxs_position_;
  std::vector<float>          vtxs_ndof_;
  std::vector<float>          vtxs_sumpt_;
  std::vector<int>            vtxs_isFake_;
  std::vector<float>          vtxs_xError_;
  std::vector<float>          vtxs_yError_;
  std::vector<float>          vtxs_zError_;
  std::vector<vector<float>>  vtxs_covMatrix_;

  //tracks
  std::vector<LorentzVector>  trks_trk_p4_;
  std::vector<float>          trks_chi2_;
  std::vector<float>          trks_ndof_;
  std::vector<float>          trks_d0_;
  std::vector<int>            trks_nlayers_;
  std::vector<int>            trks_valid_pixelhits_;
  std::vector<float>          trks_z0_;
  std::vector<float>          trks_z0Err_;
  std::vector<float>          trks_etaErr_;
  std::vector<float>          trks_d0Err_;
  std::vector<float>          trks_phiErr_;
  std::vector<float>          trks_d0phiCov_;
  std::vector<int>            trks_qualityMask_;
  std::vector<int>            trks_charge_;

  //electrons
  std::vector<LorentzVector>  els_p4_;
  std::vector<LorentzVector>  els_trk_p4_;
  std::vector<LorentzVector>  els_vertex_p4_;
  std::vector<float>          els_lh_;
  std::vector<float>          els_etaSC_;
  std::vector<float>          els_sigmaIEtaIEta_;
  std::vector<float>          els_dEtaIn_;
  std::vector<float>          els_dPhiIn_;
  std::vector<float>          els_hOverE_;
  std::vector<float>          els_tkIso_;
  std::vector<float>          els_d0corr_;
  std::vector<float>          els_d0_;
  std::vector<float>          els_z0corr_;
  std::vector<float>          els_fbrem_;
  std::vector<float>          els_eOverPIn_;
  std::vector<float>          els_eSeedOverPOut_;
  std::vector<float>          els_eSeedOverPIn_;
  std::vector<float>          els_sigmaIPhiIPhi_;
  std::vector<float>          els_eSC_;
  std::vector<float>          els_ip3d_;
  std::vector<float>          els_ip3derr_;
  std::vector<float>          els_chi2_;
  std::vector<float>          els_ndof_;
  std::vector<float>          els_dEtaOut_;
  std::vector<float>          els_dPhiOut_;
  std::vector<float>          els_eSCRaw_;
  std::vector<float>          els_etaSCwidth_;
  std::vector<float>          els_phiSCwidth_;
  std::vector<float>          els_eSCPresh_;
  std::vector<float>          els_iso03_pf_ch_;
  std::vector<float>          els_iso03_pf_nhad05_;
  std::vector<float>          els_iso03_pf_gamma05_;
  std::vector<float>          els_iso04_pf_ch_;
  std::vector<float>          els_iso04_pf_nhad05_;
  std::vector<float>          els_iso04_pf_gamma05_;
  std::vector<float>          els_e5x5_;
  std::vector<float>          els_e1x5_;
  std::vector<float>          els_e3x3_;
  std::vector<float>          els_ecalEnergy_;
  std::vector<float>          els_eOverPOut_;
  std::vector<float>          els_ecalIso_;
  std::vector<float>          els_hcalIso_;
  std::vector<float>          els_trkshFrac_;
  std::vector<float>          els_conv_dist_;
  std::vector<float>          els_conv_dcot_;
  std::vector<float>          els_conv_old_dist_;
  std::vector<float>          els_conv_old_dcot_;
  std::vector<float>          els_iso04_pf2012_ch_;
  std::vector<float>          els_iso04_pf2012_em_;
  std::vector<float>          els_iso04_pf2012_nh_;
  std::vector<float>          els_iso03_pf2012_ch_;
  std::vector<float>          els_iso03_pf2012_em_;
  std::vector<float>          els_iso03_pf2012_nh_;
  std::vector<float>          els_ecalIso04_;
  std::vector<float>          els_hcalIso04_;
  std::vector<int>            els_nSeed_;
  std::vector<int>            els_scindex_;
  std::vector<int>            els_charge_;
  std::vector<int>            els_gsftrkidx_;
  std::vector<int>            els_exp_innerlayers_;
  std::vector<int>            els_trkidx_;
  std::vector<int>            els_type_;
  std::vector<int>            els_fiduciality_;
  std::vector<int>            els_sccharge_;
  std::vector<int>            els_trk_charge_;
  std::vector<int>            els_closestMuon_;

  //muons
  std::vector<LorentzVector>  mus_p4_;
  std::vector<LorentzVector>  mus_trk_p4_;
  std::vector<LorentzVector>  mus_vertex_p4_;
  std::vector<LorentzVector>  mus_sta_p4_;
  std::vector<float>          mus_gfit_chi2_;
  std::vector<float>          mus_gfit_ndof_;
  std::vector<float>          mus_ptErr_;
  std::vector<float>          mus_trkKink_;
  std::vector<float>          mus_d0corr_;
  std::vector<float>          mus_d0_;
  std::vector<float>          mus_z0corr_;
  std::vector<float>          mus_chi2_;
  std::vector<float>          mus_ndof_;
  std::vector<float>          mus_ip3d_;
  std::vector<float>          mus_ip3derr_;
  std::vector<float>          mus_segmCompatibility_;
  std::vector<float>          mus_caloCompatibility_;
  std::vector<float>          mus_e_had_;
  std::vector<float>          mus_e_ho_;
  std::vector<float>          mus_e_em_;
  std::vector<float>          mus_e_hadS9_;
  std::vector<float>          mus_e_hoS9_;
  std::vector<float>          mus_e_emS9_;
  std::vector<float>          mus_iso03_sumPt_;
  std::vector<float>          mus_iso03_emEt_;
  std::vector<float>          mus_iso03_hadEt_;
  std::vector<float>          mus_iso05_sumPt_;
  std::vector<float>          mus_iso05_emEt_;
  std::vector<float>          mus_iso05_hadEt_;
  std::vector<float>          mus_sta_d0_;
  std::vector<float>          mus_sta_z0corr_;
  std::vector<float>          mus_isoR03_pf_ChargedHadronPt_;
  std::vector<float>          mus_isoR03_pf_NeutralHadronEt_;
  std::vector<float>          mus_isoR03_pf_PhotonEt_;
  std::vector<float>          mus_isoR03_pf_PUPt_;
  std::vector<float>          mus_iso_ecalvetoDep_;
  std::vector<float>          mus_iso_hcalvetoDep_;
  std::vector<int>            mus_gfit_validSTAHits_;
  std::vector<int>            mus_numberOfMatchedStations_;
  std::vector<int>            mus_pfmusidx_;
  std::vector<int>            mus_charge_;
  std::vector<int>            mus_validHits_;
  std::vector<int>            mus_trkidx_;
  std::vector<int>            mus_pid_PFMuon_;
  std::vector<int>            mus_pid_TMLastStationTight_;
  std::vector<int>            mus_nmatches_;
  std::vector<int>            mus_goodmask_;
  std::vector<int>            mus_type_;

  //dilepton hypothesis
  std::vector<vector<LorentzVector> > hyp_jets_p4_;
  std::vector<LorentzVector>  hyp_p4_;
  std::vector<LorentzVector>  hyp_ll_p4_;
  std::vector<LorentzVector>  hyp_lt_p4_;
  std::vector<int>            hyp_ll_index_;
  std::vector<int>            hyp_lt_index_;
  std::vector<int>            hyp_ll_id_;
  std::vector<int>            hyp_lt_id_;
  std::vector<int>            hyp_ll_charge_;
  std::vector<int>            hyp_lt_charge_;
  std::vector<int>            hyp_type_;

  //event variables
  unsigned int                evt_run_;
  unsigned int                evt_lumiBlock_;
  unsigned int                evt_event_;
  int                         evt_isRealData_;
  float                       evt_ww_rho_vor_;
  float                       evt_ww_rho_;
  float                       evt_rho_;
  float                       evt_kt6pf_foregiso_rho_;
  float                       evt_pfmet_;
  float                       evt_pfmetPhi_;

  std::vector<float>          convs_ndof_;
  std::vector<float>          convs_chi2_;
  std::vector<float>          convs_dl_;
  std::vector<int>            convs_isConverted_;
  std::vector<vector<int>>    convs_tkalgo_;
  std::vector<vector<int>>    convs_tkidx_;
  std::vector<vector<int>>    convs_nHitsBeforeVtx_;
  std::vector<int>            convs_quality_;
  std::vector<float>          scs_sigmaIEtaIPhi_;
  std::vector<LorentzVector>  scs_pos_p4_;
  std::vector<LorentzVector>  gsftrks_p4_;
  std::vector<LorentzVector>  gsftrks_vertex_p4_;
  std::vector<float>          gsftrks_d0_;
  std::vector<float>          gsftrks_d0Err_;
  std::vector<float>          gsftrks_phiErr_;
  std::vector<float>          gsftrks_d0phiCov_;
  std::vector<float>          gsftrks_z0Err_;
  std::vector<float>          gsftrks_z0_;
  std::vector<float>          gsftrks_etaErr_;
  std::vector<LorentzVector>  pfcands_p4_;
  std::vector<int>            pfcands_trkidx_;
  std::vector<int>            pfcands_particleId_;
  std::vector<int>            pfcands_pfelsidx_;
  std::vector<int>            pfcands_vtxidx_;
  std::vector<int>            pfcands_charge_;
  std::vector<int>            pfels_elsidx_;
  std::vector<LorentzVector>  pfels_p4_;
  std::vector<LorentzVector>  pfmus_p4_;
  vector<float>               trk_met_;
  vector<float>               trk_metPhi_;
  vector<LorentzVector>       pfjets_p4_;
  vector<LorentzVector>       pfjets_corr_p4_;
  vector<float>               pfjets_area_;
  vector<float>               pfjets_JEC_;
  vector<float>               pfjets_mvavalue_;
  vector<float>               pfjets_trackCountingHighEffBJetTag_;


  ///////////////////////////////////////////
  ////  bools to keep track of variables ////
  ///////////////////////////////////////////

  bool vtxs_position_isLoaded;
  bool vtxs_ndof_isLoaded;
  bool vtxs_sumpt_isLoaded;
  bool vtxs_isFake_isLoaded;
  bool vtxs_xError_isLoaded;
  bool vtxs_yError_isLoaded;
  bool vtxs_zError_isLoaded;
  bool vtxs_covMatrix_isLoaded;

  bool trks_trk_p4_isLoaded;
  bool trks_chi2_isLoaded;
  bool trks_ndof_isLoaded;
  bool trks_d0_isLoaded;
  bool trks_nlayers_isLoaded;
  bool trks_valid_pixelhits_isLoaded;
  bool trks_z0_isLoaded;
  bool trks_z0Err_isLoaded;
  bool trks_etaErr_isLoaded;
  bool trks_d0Err_isLoaded;
  bool trks_phiErr_isLoaded;
  bool trks_d0phiCov_isLoaded;
  bool trks_qualityMask_isLoaded;
  bool trks_charge_isLoaded;

  bool els_p4_isLoaded;
  bool els_trk_p4_isLoaded;
  bool els_vertex_p4_isLoaded;
  bool els_lh_isLoaded;
  bool els_etaSC_isLoaded;
  bool els_sigmaIEtaIEta_isLoaded;
  bool els_dEtaIn_isLoaded;
  bool els_dPhiIn_isLoaded;
  bool els_hOverE_isLoaded;
  bool els_tkIso_isLoaded;
  bool els_d0corr_isLoaded;
  bool els_d0_isLoaded;
  bool els_z0corr_isLoaded;
  bool els_fbrem_isLoaded;
  bool els_eOverPIn_isLoaded;
  bool els_eSeedOverPOut_isLoaded;
  bool els_eSeedOverPIn_isLoaded;
  bool els_sigmaIPhiIPhi_isLoaded;
  bool els_eSC_isLoaded;
  bool els_ip3d_isLoaded;
  bool els_ip3derr_isLoaded;
  bool els_chi2_isLoaded;
  bool els_ndof_isLoaded;
  bool els_dEtaOut_isLoaded;
  bool els_dPhiOut_isLoaded;
  bool els_eSCRaw_isLoaded;
  bool els_etaSCwidth_isLoaded;
  bool els_phiSCwidth_isLoaded;
  bool els_eSCPresh_isLoaded;
  bool els_iso03_pf_ch_isLoaded;
  bool els_iso03_pf_nhad05_isLoaded;
  bool els_iso03_pf_gamma05_isLoaded;
  bool els_iso04_pf_ch_isLoaded;
  bool els_iso04_pf_nhad05_isLoaded;
  bool els_iso04_pf_gamma05_isLoaded;
  bool els_e5x5_isLoaded;
  bool els_e1x5_isLoaded;
  bool els_e3x3_isLoaded;
  bool els_ecalEnergy_isLoaded;
  bool els_eOverPOut_isLoaded;
  bool els_ecalIso_isLoaded;
  bool els_hcalIso_isLoaded;
  bool els_trkshFrac_isLoaded;
  bool els_conv_dist_isLoaded;
  bool els_conv_dcot_isLoaded;
  bool els_conv_old_dist_isLoaded;
  bool els_conv_old_dcot_isLoaded;
  bool els_iso04_pf2012_ch_isLoaded;
  bool els_iso04_pf2012_em_isLoaded;
  bool els_iso04_pf2012_nh_isLoaded;
  bool els_iso03_pf2012_ch_isLoaded;
  bool els_iso03_pf2012_em_isLoaded;
  bool els_iso03_pf2012_nh_isLoaded;
  bool els_ecalIso04_isLoaded;
  bool els_hcalIso04_isLoaded;
  bool els_nSeed_isLoaded;
  bool els_scindex_isLoaded;
  bool els_charge_isLoaded;
  bool els_gsftrkidx_isLoaded;
  bool els_exp_innerlayers_isLoaded;
  bool els_trkidx_isLoaded;
  bool els_type_isLoaded;
  bool els_fiduciality_isLoaded;
  bool els_sccharge_isLoaded;
  bool els_trk_charge_isLoaded;
  bool els_closestMuon_isLoaded;

  bool mus_p4_isLoaded;
  bool mus_trk_p4_isLoaded;
  bool mus_vertex_p4_isLoaded;
  bool mus_sta_p4_isLoaded;
  bool mus_gfit_chi2_isLoaded;
  bool mus_gfit_ndof_isLoaded;
  bool mus_ptErr_isLoaded;
  bool mus_trkKink_isLoaded;
  bool mus_d0corr_isLoaded;
  bool mus_d0_isLoaded;
  bool mus_z0corr_isLoaded;
  bool mus_chi2_isLoaded;
  bool mus_ndof_isLoaded;
  bool mus_ip3d_isLoaded;
  bool mus_ip3derr_isLoaded;
  bool mus_segmCompatibility_isLoaded;
  bool mus_caloCompatibility_isLoaded;
  bool mus_e_had_isLoaded;
  bool mus_e_ho_isLoaded;
  bool mus_e_em_isLoaded;
  bool mus_e_hadS9_isLoaded;
  bool mus_e_hoS9_isLoaded;
  bool mus_e_emS9_isLoaded;
  bool mus_iso03_sumPt_isLoaded;
  bool mus_iso03_emEt_isLoaded;
  bool mus_iso03_hadEt_isLoaded;
  bool mus_iso05_sumPt_isLoaded;
  bool mus_iso05_emEt_isLoaded;
  bool mus_iso05_hadEt_isLoaded;
  bool mus_sta_d0_isLoaded;
  bool mus_sta_z0corr_isLoaded;
  bool mus_isoR03_pf_ChargedHadronPt_isLoaded;
  bool mus_isoR03_pf_NeutralHadronEt_isLoaded;
  bool mus_isoR03_pf_PhotonEt_isLoaded;
  bool mus_isoR03_pf_PUPt_isLoaded;
  bool mus_iso_ecalvetoDep_isLoaded;
  bool mus_iso_hcalvetoDep_isLoaded;
  bool mus_gfit_validSTAHits_isLoaded;
  bool mus_numberOfMatchedStations_isLoaded;
  bool mus_pfmusidx_isLoaded;
  bool mus_charge_isLoaded;
  bool mus_validHits_isLoaded;
  bool mus_trkidx_isLoaded;
  bool mus_pid_PFMuon_isLoaded;
  bool mus_pid_TMLastStationTight_isLoaded;
  bool mus_nmatches_isLoaded;
  bool mus_goodmask_isLoaded;
  bool mus_type_isLoaded;

  bool hyp_jets_p4_isLoaded;
  bool hyp_p4_isLoaded;
  bool hyp_ll_p4_isLoaded;
  bool hyp_lt_p4_isLoaded;
  bool hyp_ll_index_isLoaded;
  bool hyp_lt_index_isLoaded;
  bool hyp_ll_id_isLoaded;
  bool hyp_lt_id_isLoaded;
  bool hyp_ll_charge_isLoaded;
  bool hyp_lt_charge_isLoaded;
  bool hyp_type_isLoaded;

  bool evt_run_isLoaded;
  bool evt_lumiBlock_isLoaded;
  bool evt_event_isLoaded;
  bool evt_isRealData_isLoaded;
  bool evt_ww_rho_vor_isLoaded;
  bool evt_ww_rho_isLoaded;
  bool evt_rho_isLoaded;
  bool evt_kt6pf_foregiso_rho_isLoaded;
  bool evt_pfmet_isLoaded;
  bool evt_pfmetPhi_isLoaded;

  bool convs_ndof_isLoaded;
  bool convs_chi2_isLoaded;
  bool convs_dl_isLoaded;
  bool convs_isConverted_isLoaded;
  bool convs_tkalgo_isLoaded;
  bool convs_tkidx_isLoaded;
  bool convs_nHitsBeforeVtx_isLoaded;
  bool convs_quality_isLoaded;
  bool scs_sigmaIEtaIPhi_isLoaded;
  bool scs_pos_p4_isLoaded;
  bool gsftrks_p4_isLoaded;
  bool gsftrks_vertex_p4_isLoaded;
  bool gsftrks_d0_isLoaded;
  bool gsftrks_d0Err_isLoaded;
  bool gsftrks_phiErr_isLoaded;
  bool gsftrks_d0phiCov_isLoaded;
  bool gsftrks_z0Err_isLoaded;
  bool gsftrks_z0_isLoaded;
  bool gsftrks_etaErr_isLoaded;
  bool pfcands_p4_isLoaded;
  bool pfcands_trkidx_isLoaded;
  bool pfcands_particleId_isLoaded;
  bool pfcands_pfelsidx_isLoaded;
  bool pfcands_vtxidx_isLoaded;
  bool pfcands_charge_isLoaded;
  bool pfels_elsidx_isLoaded;
  bool pfels_p4_isLoaded;
  bool pfmus_p4_isLoaded;
  bool trk_met_isLoaded;
  bool trk_metPhi_isLoaded;
  bool pfjets_p4_isLoaded;
  bool pfjets_corr_p4_isLoaded;
  bool pfjets_area_isLoaded;
  bool pfjets_JEC_isLoaded;
  bool pfjets_mvavalue_isLoaded;
  bool pfjets_trackCountingHighEffBJetTag_isLoaded;


  ///////////////////////////////////////////
  //define functions to access data members//
  ///////////////////////////////////////////


  //vertex 
  std::vector<LorentzVector>  &vtxs_position(){
    if(!vtxs_position_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_position not loaded!";
    return(vtxs_position_);
  }
  std::vector<float>          &vtxs_ndof(){
    if(!vtxs_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_ndof not loaded!";
    return(vtxs_ndof_);
  }
  std::vector<float>          &vtxs_sumpt(){
    if(!vtxs_sumpt_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_sumpt not loaded!";
    return(vtxs_sumpt_);
  }
  std::vector<int>            &vtxs_isFake(){
    if(!vtxs_isFake_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_isFake not loaded!";
    return(vtxs_isFake_);
  }
  std::vector<float>          &vtxs_xError(){
    if(!vtxs_xError_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_xError not loaded!";
    return(vtxs_xError_);
  }
  std::vector<float>          &vtxs_yError(){
    if(!vtxs_yError_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_yError not loaded!";
    return(vtxs_yError_);
  }
  std::vector<float>          &vtxs_zError(){
    if(!vtxs_zError_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_zError not loaded!";
    return(vtxs_zError_);
  }
  std::vector<vector<float>>  &vtxs_covMatrix(){
    if(!vtxs_covMatrix_isLoaded) edm::LogWarning("VariableNotSet") << "vtxs_covMatrix not loaded!";
    return(vtxs_covMatrix_);
  }

  //tracks
  std::vector<LorentzVector>  &trks_trk_p4(){
    if(!trks_trk_p4_isLoaded) edm::LogWarning("VariableNotSet") << "trks_trk_p4 not loaded!";
    return(trks_trk_p4_);
  }
  std::vector<float>          &trks_chi2(){
    if(!trks_chi2_isLoaded) edm::LogWarning("VariableNotSet") << "trks_chi2 not loaded!";
    return(trks_chi2_);
  }
  std::vector<float>          &trks_ndof(){
    if(!trks_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "trks_ndof not loaded!";
    return(trks_ndof_);
  }
  std::vector<float>          &trks_d0(){
    if(!trks_d0_isLoaded) edm::LogWarning("VariableNotSet") << "trks_d0 not loaded!";
    return(trks_d0_);
  }
  std::vector<int>            &trks_nlayers(){
    if(!trks_nlayers_isLoaded) edm::LogWarning("VariableNotSet") << "trks_nlayers not loaded!";
    return(trks_nlayers_);
  }
  std::vector<int>            &trks_valid_pixelhits(){
    if(!trks_valid_pixelhits_isLoaded) edm::LogWarning("VariableNotSet") << "trks_valid_pixelhits not loaded!";
    return(trks_valid_pixelhits_);
  }
  std::vector<float>          &trks_z0(){
    if(!trks_z0_isLoaded) edm::LogWarning("VariableNotSet") << "trks_z0 not loaded!";
    return(trks_z0_);
  }
  std::vector<float>          &trks_z0Err(){
    if(!trks_z0Err_isLoaded) edm::LogWarning("VariableNotSet") << "trks_z0Err not loaded!";
    return(trks_z0Err_);
  }
  std::vector<float>          &trks_etaErr(){
    if(!trks_etaErr_isLoaded) edm::LogWarning("VariableNotSet") << "trks_etaErr not loaded!";
    return(trks_etaErr_);
  }
  std::vector<float>          &trks_d0Err(){
    if(!trks_d0Err_isLoaded) edm::LogWarning("VariableNotSet") << "trks_d0Err not loaded!";
    return(trks_d0Err_);
  }
  std::vector<float>          &trks_phiErr(){
    if(!trks_phiErr_isLoaded) edm::LogWarning("VariableNotSet") << "trks_phiErr not loaded!";
    return(trks_phiErr_);
  }
  std::vector<float>          &trks_d0phiCov(){
    if(!trks_d0phiCov_isLoaded) edm::LogWarning("VariableNotSet") << "trks_d0phiCov not loaded!";
    return(trks_d0phiCov_);
  }
  std::vector<int>            &trks_qualityMask(){
    if(!trks_qualityMask_isLoaded) edm::LogWarning("VariableNotSet") << "trks_qualityMask not loaded!";
    return(trks_qualityMask_);
  }
  std::vector<int>            &trks_charge(){
    if(!trks_charge_isLoaded) edm::LogWarning("VariableNotSet") << "trks_charge not loaded!";
    return(trks_charge_);
  }

  //electrons
  std::vector<LorentzVector>  &els_p4(){
    if(!els_p4_isLoaded) edm::LogWarning("VariableNotSet") << "els_p4 not loaded!";
    return(els_p4_);
  }
  std::vector<LorentzVector>  &els_trk_p4(){
    if(!els_trk_p4_isLoaded) edm::LogWarning("VariableNotSet") << "els_trk_p4 not loaded!";
    return(els_trk_p4_);
  }
  std::vector<LorentzVector>  &els_vertex_p4(){
    if(!els_vertex_p4_isLoaded) edm::LogWarning("VariableNotSet") << "els_vertex_p4 not loaded!";
    return(els_vertex_p4_);
  }
  std::vector<float>          &els_lh(){
    if(!els_lh_isLoaded) edm::LogWarning("VariableNotSet") << "els_lh not loaded!";
    return(els_lh_);
  }
  std::vector<float>          &els_etaSC(){
    if(!els_etaSC_isLoaded) edm::LogWarning("VariableNotSet") << "els_etaSC not loaded!";
    return(els_etaSC_);
  }
  std::vector<float>          &els_sigmaIEtaIEta(){
    if(!els_sigmaIEtaIEta_isLoaded) edm::LogWarning("VariableNotSet") << "els_sigmaIEtaIEta not loaded!";
    return(els_sigmaIEtaIEta_);
  }
  std::vector<float>          &els_dEtaIn(){
    if(!els_dEtaIn_isLoaded) edm::LogWarning("VariableNotSet") << "els_dEtaIn not loaded!";
    return(els_dEtaIn_);
  }
  std::vector<float>          &els_dPhiIn(){
    if(!els_dPhiIn_isLoaded) edm::LogWarning("VariableNotSet") << "els_dPhiIn not loaded!";
    return(els_dPhiIn_);
  }
  std::vector<float>          &els_hOverE(){
    if(!els_hOverE_isLoaded) edm::LogWarning("VariableNotSet") << "els_hOverE not loaded!";
    return(els_hOverE_);
  }
  std::vector<float>          &els_tkIso(){
    if(!els_tkIso_isLoaded) edm::LogWarning("VariableNotSet") << "els_tkIso not loaded!";
    return(els_tkIso_);
  }
  std::vector<float>          &els_d0corr(){
    if(!els_d0corr_isLoaded) edm::LogWarning("VariableNotSet") << "els_d0corr not loaded!";
    return(els_d0corr_);
  }
  std::vector<float>          &els_d0(){
    if(!els_d0_isLoaded) edm::LogWarning("VariableNotSet") << "els_d0 not loaded!";
    return(els_d0_);
  }
  std::vector<float>          &els_z0corr(){
    if(!els_z0corr_isLoaded) edm::LogWarning("VariableNotSet") << "els_z0corr not loaded!";
    return(els_z0corr_);
  }
  std::vector<float>          &els_fbrem(){
    if(!els_fbrem_isLoaded) edm::LogWarning("VariableNotSet") << "els_fbrem not loaded!";
    return(els_fbrem_);
  }
  std::vector<float>          &els_eOverPIn(){
    if(!els_eOverPIn_isLoaded) edm::LogWarning("VariableNotSet") << "els_eOverPIn not loaded!";
    return(els_eOverPIn_);
  }
  std::vector<float>          &els_eSeedOverPOut(){
    if(!els_eSeedOverPOut_isLoaded) edm::LogWarning("VariableNotSet") << "els_eSeedOverPOut not loaded!";
    return(els_eSeedOverPOut_);
  }
  std::vector<float>          &els_eSeedOverPIn(){
    if(!els_eSeedOverPIn_isLoaded) edm::LogWarning("VariableNotSet") << "els_eSeedOverPIn not loaded!";
    return(els_eSeedOverPIn_);
  }
  std::vector<float>          &els_sigmaIPhiIPhi(){
    if(!els_sigmaIPhiIPhi_isLoaded) edm::LogWarning("VariableNotSet") << "els_sigmaIPhiIPhi not loaded!";
    return(els_sigmaIPhiIPhi_);
  }
  std::vector<float>          &els_eSC(){
    if(!els_eSC_isLoaded) edm::LogWarning("VariableNotSet") << "els_eSC not loaded!";
    return(els_eSC_);
  }
  std::vector<float>          &els_ip3d(){
    if(!els_ip3d_isLoaded) edm::LogWarning("VariableNotSet") << "els_ip3d not loaded!";
    return(els_ip3d_);
  }
  std::vector<float>          &els_ip3derr(){
    if(!els_ip3derr_isLoaded) edm::LogWarning("VariableNotSet") << "els_ip3derr not loaded!";
    return(els_ip3derr_);
  }
  std::vector<float>          &els_chi2(){
    if(!els_chi2_isLoaded) edm::LogWarning("VariableNotSet") << "els_chi2 not loaded!";
    return(els_chi2_);
  }
  std::vector<float>          &els_ndof(){
    if(!els_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "els_ndof not loaded!";
    return(els_ndof_);
  }
  std::vector<float>          &els_dEtaOut(){
    if(!els_dEtaOut_isLoaded) edm::LogWarning("VariableNotSet") << "els_dEtaOut not loaded!";
    return(els_dEtaOut_);
  }
  std::vector<float>          &els_dPhiOut(){
    if(!els_dPhiOut_isLoaded) edm::LogWarning("VariableNotSet") << "els_dPhiOut not loaded!";
    return(els_dPhiOut_);
  }
  std::vector<float>          &els_eSCRaw(){
    if(!els_eSCRaw_isLoaded) edm::LogWarning("VariableNotSet") << "els_eSCRaw not loaded!";
    return(els_eSCRaw_);
  }
  std::vector<float>          &els_etaSCwidth(){
    if(!els_etaSCwidth_isLoaded) edm::LogWarning("VariableNotSet") << "els_etaSCwidth not loaded!";
    return(els_etaSCwidth_);
  }
  std::vector<float>          &els_phiSCwidth(){
    if(!els_phiSCwidth_isLoaded) edm::LogWarning("VariableNotSet") << "els_phiSCwidth not loaded!";
    return(els_phiSCwidth_);
  }
  std::vector<float>          &els_eSCPresh(){
    if(!els_eSCPresh_isLoaded) edm::LogWarning("VariableNotSet") << "els_eSCPresh not loaded!";
    return(els_eSCPresh_);
  }
  std::vector<float>          &els_iso03_pf_ch(){
    if(!els_iso03_pf_ch_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf_ch not loaded!";
    return(els_iso03_pf_ch_);
  }
  std::vector<float>          &els_iso03_pf_nhad05(){
    if(!els_iso03_pf_nhad05_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf_nhad05 not loaded!";
    return(els_iso03_pf_nhad05_);
  }
  std::vector<float>          &els_iso03_pf_gamma05(){
    if(!els_iso03_pf_gamma05_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf_gamma05 not loaded!";
    return(els_iso03_pf_gamma05_);
  }
  std::vector<float>          &els_iso04_pf_ch(){
    if(!els_iso04_pf_ch_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf_ch not loaded!";
    return(els_iso04_pf_ch_);
  }
  std::vector<float>          &els_iso04_pf_nhad05(){
    if(!els_iso04_pf_nhad05_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf_nhad05 not loaded!";
    return(els_iso04_pf_nhad05_);
  }
  std::vector<float>          &els_iso04_pf_gamma05(){
    if(!els_iso04_pf_gamma05_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf_gamma05 not loaded!";
    return(els_iso04_pf_gamma05_);
  }
  std::vector<float>          &els_e5x5(){
    if(!els_e5x5_isLoaded) edm::LogWarning("VariableNotSet") << "els_e5x5 not loaded!";
    return(els_e5x5_);
  }
  std::vector<float>          &els_e1x5(){
    if(!els_e1x5_isLoaded) edm::LogWarning("VariableNotSet") << "els_e1x5 not loaded!";
    return(els_e1x5_);
  }
  std::vector<float>          &els_e3x3(){
    if(!els_e3x3_isLoaded) edm::LogWarning("VariableNotSet") << "els_e3x3 not loaded!";
    return(els_e3x3_);
  }
  std::vector<float>          &els_ecalEnergy(){
    if(!els_ecalEnergy_isLoaded) edm::LogWarning("VariableNotSet") << "els_ecalEnergy not loaded!";
    return(els_ecalEnergy_);
  }
  std::vector<float>          &els_eOverPOut(){
    if(!els_eOverPOut_isLoaded) edm::LogWarning("VariableNotSet") << "els_eOverPOut not loaded!";
    return(els_eOverPOut_);
  }
  std::vector<float>          &els_ecalIso(){
    if(!els_ecalIso_isLoaded) edm::LogWarning("VariableNotSet") << "els_ecalIso not loaded!";
    return(els_ecalIso_);
  }
  std::vector<float>          &els_hcalIso(){
    if(!els_hcalIso_isLoaded) edm::LogWarning("VariableNotSet") << "els_hcalIso not loaded!";
    return(els_hcalIso_);
  }
  std::vector<float>          &els_trkshFrac(){
    if(!els_trkshFrac_isLoaded) edm::LogWarning("VariableNotSet") << "els_trkshFrac not loaded!";
    return(els_trkshFrac_);
  }
  std::vector<float>          &els_conv_dist(){
    if(!els_conv_dist_isLoaded) edm::LogWarning("VariableNotSet") << "els_conv_dist not loaded!";
    return(els_conv_dist_);
  }
  std::vector<float>          &els_conv_dcot(){
    if(!els_conv_dcot_isLoaded) edm::LogWarning("VariableNotSet") << "els_conv_dcot not loaded!";
    return(els_conv_dcot_);
  }
  std::vector<float>          &els_conv_old_dist(){
    if(!els_conv_old_dist_isLoaded) edm::LogWarning("VariableNotSet") << "els_conv_old_dist not loaded!";
    return(els_conv_old_dist_);
  }
  std::vector<float>          &els_conv_old_dcot(){
    if(!els_conv_old_dcot_isLoaded) edm::LogWarning("VariableNotSet") << "els_conv_old_dcot not loaded!";
    return(els_conv_old_dcot_);
  }
  std::vector<float>          &els_iso04_pf2012_ch(){
    if(!els_iso04_pf2012_ch_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf2012_ch not loaded!";
    return(els_iso04_pf2012_ch_);
  }
  std::vector<float>          &els_iso04_pf2012_em(){
    if(!els_iso04_pf2012_em_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf2012_em not loaded!";
    return(els_iso04_pf2012_em_);
  }
  std::vector<float>          &els_iso04_pf2012_nh(){
    if(!els_iso04_pf2012_nh_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso04_pf2012_nh not loaded!";
    return(els_iso04_pf2012_nh_);
  }
  std::vector<float>          &els_iso03_pf2012_ch(){
    if(!els_iso03_pf2012_ch_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf2012_ch not loaded!";
    return(els_iso03_pf2012_ch_);
  }
  std::vector<float>          &els_iso03_pf2012_em(){
    if(!els_iso03_pf2012_em_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf2012_em not loaded!";
    return(els_iso03_pf2012_em_);
  }
  std::vector<float>          &els_iso03_pf2012_nh(){
    if(!els_iso03_pf2012_nh_isLoaded) edm::LogWarning("VariableNotSet") << "els_iso03_pf2012_nh not loaded!";
    return(els_iso03_pf2012_nh_);
  }
  std::vector<float>          &els_ecalIso04(){
    if(!els_ecalIso04_isLoaded) edm::LogWarning("VariableNotSet") << "els_ecalIso04 not loaded!";
    return(els_ecalIso04_);
  }
  std::vector<float>          &els_hcalIso04(){
    if(!els_hcalIso04_isLoaded) edm::LogWarning("VariableNotSet") << "els_hcalIso04 not loaded!";
    return(els_hcalIso04_);
  }
  std::vector<int>            &els_nSeed(){
    if(!els_nSeed_isLoaded) edm::LogWarning("VariableNotSet") << "els_nSeed not loaded!";
    return(els_nSeed_);
  }
  std::vector<int>            &els_scindex(){
    if(!els_scindex_isLoaded) edm::LogWarning("VariableNotSet") << "els_scindex not loaded!";
    return(els_scindex_);
  }
  std::vector<int>            &els_charge(){
    if(!els_charge_isLoaded) edm::LogWarning("VariableNotSet") << "els_charge not loaded!";
    return(els_charge_);
  }
  std::vector<int>            &els_gsftrkidx(){
    if(!els_gsftrkidx_isLoaded) edm::LogWarning("VariableNotSet") << "els_gsftrkidx not loaded!";
    return(els_gsftrkidx_);
  }
  std::vector<int>            &els_exp_innerlayers(){
    if(!els_exp_innerlayers_isLoaded) edm::LogWarning("VariableNotSet") << "els_exp_innerlayers not loaded!";
    return(els_exp_innerlayers_);
  }
  std::vector<int>            &els_trkidx(){
    if(!els_trkidx_isLoaded) edm::LogWarning("VariableNotSet") << "els_trkidx not loaded!";
    return(els_trkidx_);
  }
  std::vector<int>            &els_type(){
    if(!els_type_isLoaded) edm::LogWarning("VariableNotSet") << "els_type not loaded!";
    return(els_type_);
  }
  std::vector<int>            &els_fiduciality(){
    if(!els_fiduciality_isLoaded) edm::LogWarning("VariableNotSet") << "els_fiduciality not loaded!";
    return(els_fiduciality_);
  }
  std::vector<int>            &els_sccharge(){
    if(!els_sccharge_isLoaded) edm::LogWarning("VariableNotSet") << "els_sccharge not loaded!";
    return(els_sccharge_);
  }
  std::vector<int>            &els_trk_charge(){
    if(!els_trk_charge_isLoaded) edm::LogWarning("VariableNotSet") << "els_trk_charge not loaded!";
    return(els_trk_charge_);
  }
  std::vector<int>            &els_closestMuon(){
    if(!els_closestMuon_isLoaded) edm::LogWarning("VariableNotSet") << "els_closestMuon not loaded!";
    return(els_closestMuon_);
  }

  //muons
  std::vector<LorentzVector>  &mus_p4(){
    if(!mus_p4_isLoaded) edm::LogWarning("VariableNotSet") << "mus_p4 not loaded!";
    return(mus_p4_);
  }
  std::vector<LorentzVector>  &mus_trk_p4(){
    if(!mus_trk_p4_isLoaded) edm::LogWarning("VariableNotSet") << "mus_trk_p4 not loaded!";
    return(mus_trk_p4_);
  }
  std::vector<LorentzVector>  &mus_vertex_p4(){
    if(!mus_vertex_p4_isLoaded) edm::LogWarning("VariableNotSet") << "mus_vertex_p4 not loaded!";
    return(mus_vertex_p4_);
  }
  std::vector<LorentzVector>  &mus_sta_p4(){
    if(!mus_sta_p4_isLoaded) edm::LogWarning("VariableNotSet") << "mus_sta_p4 not loaded!";
    return(mus_sta_p4_);
  }
  std::vector<float>          &mus_gfit_chi2(){
    if(!mus_gfit_chi2_isLoaded) edm::LogWarning("VariableNotSet") << "mus_gfit_chi2 not loaded!";
    return(mus_gfit_chi2_);
  }
  std::vector<float>          &mus_gfit_ndof(){
    if(!mus_gfit_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "mus_gfit_ndof not loaded!";
    return(mus_gfit_ndof_);
  }
  std::vector<float>          &mus_ptErr(){
    if(!mus_ptErr_isLoaded) edm::LogWarning("VariableNotSet") << "mus_ptErr not loaded!";
    return(mus_ptErr_);
  }
  std::vector<float>          &mus_trkKink(){
    if(!mus_trkKink_isLoaded) edm::LogWarning("VariableNotSet") << "mus_trkKink not loaded!";
    return(mus_trkKink_);
  }
  std::vector<float>          &mus_d0corr(){
    if(!mus_d0corr_isLoaded) edm::LogWarning("VariableNotSet") << "mus_d0corr not loaded!";
    return(mus_d0corr_);
  }
  std::vector<float>          &mus_d0(){
    if(!mus_d0_isLoaded) edm::LogWarning("VariableNotSet") << "mus_d0 not loaded!";
    return(mus_d0_);
  }
  std::vector<float>          &mus_z0corr(){
    if(!mus_z0corr_isLoaded) edm::LogWarning("VariableNotSet") << "mus_z0corr not loaded!";
    return(mus_z0corr_);
  }
  std::vector<float>          &mus_chi2(){
    if(!mus_chi2_isLoaded) edm::LogWarning("VariableNotSet") << "mus_chi2 not loaded!";
    return(mus_chi2_);
  }
  std::vector<float>          &mus_ndof(){
    if(!mus_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "mus_ndof not loaded!";
    return(mus_ndof_);
  }
  std::vector<float>          &mus_ip3d(){
    if(!mus_ip3d_isLoaded) edm::LogWarning("VariableNotSet") << "mus_ip3d not loaded!";
    return(mus_ip3d_);
  }
  std::vector<float>          &mus_ip3derr(){
    if(!mus_ip3derr_isLoaded) edm::LogWarning("VariableNotSet") << "mus_ip3derr not loaded!";
    return(mus_ip3derr_);
  }
  std::vector<float>          &mus_segmCompatibility(){
    if(!mus_segmCompatibility_isLoaded) edm::LogWarning("VariableNotSet") << "mus_segmCompatibility not loaded!";
    return(mus_segmCompatibility_);
  }
  std::vector<float>          &mus_caloCompatibility(){
    if(!mus_caloCompatibility_isLoaded) edm::LogWarning("VariableNotSet") << "mus_caloCompatibility not loaded!";
    return(mus_caloCompatibility_);
  }
  std::vector<float>          &mus_e_had(){
    if(!mus_e_had_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_had not loaded!";
    return(mus_e_had_);
  }
  std::vector<float>          &mus_e_ho(){
    if(!mus_e_ho_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_ho not loaded!";
    return(mus_e_ho_);
  }
  std::vector<float>          &mus_e_em(){
    if(!mus_e_em_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_em not loaded!";
    return(mus_e_em_);
  }
  std::vector<float>          &mus_e_hadS9(){
    if(!mus_e_hadS9_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_hadS9 not loaded!";
    return(mus_e_hadS9_);
  }
  std::vector<float>          &mus_e_hoS9(){
    if(!mus_e_hoS9_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_hoS9 not loaded!";
    return(mus_e_hoS9_);
  }
  std::vector<float>          &mus_e_emS9(){
    if(!mus_e_emS9_isLoaded) edm::LogWarning("VariableNotSet") << "mus_e_emS9 not loaded!";
    return(mus_e_emS9_);
  }
  std::vector<float>          &mus_iso03_sumPt(){
    if(!mus_iso03_sumPt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso03_sumPt not loaded!";
    return(mus_iso03_sumPt_);
  }
  std::vector<float>          &mus_iso03_emEt(){
    if(!mus_iso03_emEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso03_emEt not loaded!";
    return(mus_iso03_emEt_);
  }
  std::vector<float>          &mus_iso03_hadEt(){
    if(!mus_iso03_hadEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso03_hadEt not loaded!";
    return(mus_iso03_hadEt_);
  }
  std::vector<float>          &mus_iso05_sumPt(){
    if(!mus_iso05_sumPt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso05_sumPt not loaded!";
    return(mus_iso05_sumPt_);
  }
  std::vector<float>          &mus_iso05_emEt(){
    if(!mus_iso05_emEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso05_emEt not loaded!";
    return(mus_iso05_emEt_);
  }
  std::vector<float>          &mus_iso05_hadEt(){
    if(!mus_iso05_hadEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso05_hadEt not loaded!";
    return(mus_iso05_hadEt_);
  }
  std::vector<float>          &mus_sta_d0(){
    if(!mus_sta_d0_isLoaded) edm::LogWarning("VariableNotSet") << "mus_sta_d0 not loaded!";
    return(mus_sta_d0_);
  }
  std::vector<float>          &mus_sta_z0corr(){
    if(!mus_sta_z0corr_isLoaded) edm::LogWarning("VariableNotSet") << "mus_sta_z0corr not loaded!";
    return(mus_sta_z0corr_);
  }
  std::vector<float>          &mus_isoR03_pf_ChargedHadronPt(){
    if(!mus_isoR03_pf_ChargedHadronPt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_isoR03_pf_ChargedHadronPt not loaded!";
    return(mus_isoR03_pf_ChargedHadronPt_);
  }
  std::vector<float>          &mus_isoR03_pf_NeutralHadronEt(){
    if(!mus_isoR03_pf_NeutralHadronEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_isoR03_pf_NeutralHadronEt not loaded!";
    return(mus_isoR03_pf_NeutralHadronEt_);
  }
  std::vector<float>          &mus_isoR03_pf_PhotonEt(){
    if(!mus_isoR03_pf_PhotonEt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_isoR03_pf_PhotonEt not loaded!";
    return(mus_isoR03_pf_PhotonEt_);
  }
  std::vector<float>          &mus_isoR03_pf_PUPt(){
    if(!mus_isoR03_pf_PUPt_isLoaded) edm::LogWarning("VariableNotSet") << "mus_isoR03_pf_PUPt not loaded!";
    return(mus_isoR03_pf_PUPt_);
  }
  std::vector<float>          &mus_iso_ecalvetoDep(){
    if(!mus_iso_ecalvetoDep_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso_ecalvetoDep not loaded!";
    return(mus_iso_ecalvetoDep_);
  }
  std::vector<float>          &mus_iso_hcalvetoDep(){
    if(!mus_iso_hcalvetoDep_isLoaded) edm::LogWarning("VariableNotSet") << "mus_iso_hcalvetoDep not loaded!";
    return(mus_iso_hcalvetoDep_);
  }
  std::vector<int>            &mus_gfit_validSTAHits(){
    if(!mus_gfit_validSTAHits_isLoaded) edm::LogWarning("VariableNotSet") << "mus_gfit_validSTAHits not loaded!";
    return(mus_gfit_validSTAHits_);
  }
  std::vector<int>            &mus_numberOfMatchedStations(){
    if(!mus_numberOfMatchedStations_isLoaded) edm::LogWarning("VariableNotSet") << "mus_numberOfMatchedStations not loaded!";
    return(mus_numberOfMatchedStations_);
  }
  std::vector<int>            &mus_pfmusidx(){
    if(!mus_pfmusidx_isLoaded) edm::LogWarning("VariableNotSet") << "mus_pfmusidx not loaded!";
    return(mus_pfmusidx_);
  }
  std::vector<int>            &mus_charge(){
    if(!mus_charge_isLoaded) edm::LogWarning("VariableNotSet") << "mus_charge not loaded!";
    return(mus_charge_);
  }
  std::vector<int>            &mus_validHits(){
    if(!mus_validHits_isLoaded) edm::LogWarning("VariableNotSet") << "mus_validHits not loaded!";
    return(mus_validHits_);
  }
  std::vector<int>            &mus_trkidx(){
    if(!mus_trkidx_isLoaded) edm::LogWarning("VariableNotSet") << "mus_trkidx not loaded!";
    return(mus_trkidx_);
  }
  std::vector<int>            &mus_pid_PFMuon(){
    if(!mus_pid_PFMuon_isLoaded) edm::LogWarning("VariableNotSet") << "mus_pid_PFMuon not loaded!";
    return(mus_pid_PFMuon_);
  }
  std::vector<int>            &mus_pid_TMLastStationTight(){
    if(!mus_pid_TMLastStationTight_isLoaded) edm::LogWarning("VariableNotSet") << "mus_pid_TMLastStationTight not loaded!";
    return(mus_pid_TMLastStationTight_);
  }
  std::vector<int>            &mus_nmatches(){
    if(!mus_nmatches_isLoaded) edm::LogWarning("VariableNotSet") << "mus_nmatches not loaded!";
    return(mus_nmatches_);
  }
  std::vector<int>            &mus_goodmask(){
    if(!mus_goodmask_isLoaded) edm::LogWarning("VariableNotSet") << "mus_goodmask not loaded!";
    return(mus_goodmask_);
  }
  std::vector<int>            &mus_type(){
    if(!mus_type_isLoaded) edm::LogWarning("VariableNotSet") << "mus_type not loaded!";
    return(mus_type_);
  }

  //dilepton hypothesis
  std::vector<vector<LorentzVector> > &hyp_jets_p4(){
    if(!hyp_jets_p4_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_jets_p4 not loaded!";
    return(hyp_jets_p4_);
  }
  std::vector<LorentzVector>  &hyp_p4(){
    if(!hyp_p4_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_p4 not loaded!";
    return(hyp_p4_);
  }
  std::vector<LorentzVector>  &hyp_ll_p4(){
    if(!hyp_ll_p4_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_ll_p4 not loaded!";
    return(hyp_ll_p4_);
  }
  std::vector<LorentzVector>  &hyp_lt_p4(){
    if(!hyp_lt_p4_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_lt_p4 not loaded!";
    return(hyp_lt_p4_);
  }
  std::vector<int>            &hyp_ll_index(){
    if(!hyp_ll_index_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_ll_index not loaded!";
    return(hyp_ll_index_);
  }
  std::vector<int>            &hyp_lt_index(){
    if(!hyp_lt_index_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_lt_index not loaded!";
    return(hyp_lt_index_);
  }
  std::vector<int>            &hyp_ll_id(){
    if(!hyp_ll_id_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_ll_id not loaded!";
    return(hyp_ll_id_);
  }
  std::vector<int>            &hyp_lt_id(){
    if(!hyp_lt_id_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_lt_id not loaded!";
    return(hyp_lt_id_);
  }
  std::vector<int>            &hyp_ll_charge(){
    if(!hyp_ll_charge_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_ll_charge not loaded!";
    return(hyp_ll_charge_);
  }
  std::vector<int>            &hyp_lt_charge(){
    if(!hyp_lt_charge_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_lt_charge not loaded!";
    return(hyp_lt_charge_);
  }
  std::vector<int>            &hyp_type(){
    if(!hyp_type_isLoaded) edm::LogWarning("VariableNotSet") << "hyp_type not loaded!";
    return(hyp_type_);
  }

  //event variables
  unsigned int                &evt_run(){
    if(!evt_run_isLoaded) edm::LogWarning("VariableNotSet") << "evt_run not loaded!";
    return(evt_run_);
  }
  unsigned int                &evt_lumiBlock(){
    if(!evt_lumiBlock_isLoaded) edm::LogWarning("VariableNotSet") << "evt_lumiBlock not loaded!";
    return(evt_lumiBlock_);
  }
  unsigned int                &evt_event(){
    if(!evt_event_isLoaded) edm::LogWarning("VariableNotSet") << "evt_event not loaded!";
    return(evt_event_);
  }
  int                         &evt_isRealData(){
    if(!evt_isRealData_isLoaded) edm::LogWarning("VariableNotSet") << "evt_isRealData not loaded!";
    return(evt_isRealData_);
  }
  float                       &evt_ww_rho_vor(){
    if(!evt_ww_rho_vor_isLoaded) edm::LogWarning("VariableNotSet") << "evt_ww_rho_vor not loaded!";
    return(evt_ww_rho_vor_);
  }
  float                       &evt_ww_rho(){
    if(!evt_ww_rho_isLoaded) edm::LogWarning("VariableNotSet") << "evt_ww_rho not loaded!";
    return(evt_ww_rho_);
  }
  float                       &evt_rho(){
    if(!evt_rho_isLoaded) edm::LogWarning("VariableNotSet") << "evt_rho not loaded!";
    return(evt_rho_);
  }
  float                       &evt_kt6pf_foregiso_rho(){
    if(!evt_kt6pf_foregiso_rho_isLoaded) edm::LogWarning("VariableNotSet") << "evt_kt6pf_foregiso_rho not loaded!";
    return(evt_kt6pf_foregiso_rho_);
  }
  float                       &evt_pfmet(){
    if(!evt_pfmet_isLoaded) edm::LogWarning("VariableNotSet") << "evt_pfmet not loaded!";
    return(evt_pfmet_);
  }
  float                       &evt_pfmetPhi(){
    if(!evt_pfmetPhi_isLoaded) edm::LogWarning("VariableNotSet") << "evt_pfmetPhi not loaded!";
    return(evt_pfmetPhi_);
  }


  std::vector<float>          &convs_ndof(){
    if(!convs_ndof_isLoaded) edm::LogWarning("VariableNotSet") << "convs_ndof not loaded!";
    return(convs_ndof_);
  }
  std::vector<float>          &convs_chi2(){
    if(!convs_chi2_isLoaded) edm::LogWarning("VariableNotSet") << "convs_chi2 not loaded!";
    return(convs_chi2_);
  }
  std::vector<float>          &convs_dl(){
    if(!convs_dl_isLoaded) edm::LogWarning("VariableNotSet") << "convs_dl not loaded!";
    return(convs_dl_);
  }
  std::vector<int>            &convs_isConverted(){
    if(!convs_isConverted_isLoaded) edm::LogWarning("VariableNotSet") << "convs_isConverted not loaded!";
    return(convs_isConverted_);
  }
  std::vector<vector<int>>    &convs_tkalgo(){
    if(!convs_tkalgo_isLoaded) edm::LogWarning("VariableNotSet") << "convs_tkalgo not loaded!";
    return(convs_tkalgo_);
  }
  std::vector<vector<int>>    &convs_tkidx(){
    if(!convs_tkidx_isLoaded) edm::LogWarning("VariableNotSet") << "convs_tkidx not loaded!";
    return(convs_tkidx_);
  }
  std::vector<vector<int>>    &convs_nHitsBeforeVtx(){
    if(!convs_nHitsBeforeVtx_isLoaded) edm::LogWarning("VariableNotSet") << "convs_nHitsBeforeVtx not loaded!";
    return(convs_nHitsBeforeVtx_);
  }
  std::vector<int>            &convs_quality(){
    if(!convs_quality_isLoaded) edm::LogWarning("VariableNotSet") << "convs_quality not loaded!";
    return(convs_quality_);
  }
  std::vector<float>          &scs_sigmaIEtaIPhi(){
    if(!scs_sigmaIEtaIPhi_isLoaded) edm::LogWarning("VariableNotSet") << "scs_sigmaIEtaIPhi not loaded!";
    return(scs_sigmaIEtaIPhi_);
  }
  std::vector<LorentzVector>  &scs_pos_p4(){
    if(!scs_pos_p4_isLoaded) edm::LogWarning("VariableNotSet") << "scs_pos_p4 not loaded!";
    return(scs_pos_p4_);
  }
  std::vector<LorentzVector>  &gsftrks_p4(){
    if(!gsftrks_p4_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_p4 not loaded!";
    return(gsftrks_p4_);
  }
  std::vector<LorentzVector>  &gsftrks_vertex_p4(){
    if(!gsftrks_vertex_p4_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_vertex_p4 not loaded!";
    return(gsftrks_vertex_p4_);
  }
  std::vector<float>          &gsftrks_d0(){
    if(!gsftrks_d0_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_d0 not loaded!";
    return(gsftrks_d0_);
  }
  std::vector<float>          &gsftrks_d0Err(){
    if(!gsftrks_d0Err_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_d0Err not loaded!";
    return(gsftrks_d0Err_);
  }
  std::vector<float>          &gsftrks_phiErr(){
    if(!gsftrks_phiErr_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_phiErr not loaded!";
    return(gsftrks_phiErr_);
  }
  std::vector<float>          &gsftrks_d0phiCov(){
    if(!gsftrks_d0phiCov_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_d0phiCov not loaded!";
    return(gsftrks_d0phiCov_);
  }
  std::vector<float>          &gsftrks_z0Err(){
    if(!gsftrks_z0Err_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_z0Err not loaded!";
    return(gsftrks_z0Err_);
  }
  std::vector<float>          &gsftrks_z0(){
    if(!gsftrks_z0_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_z0 not loaded!";
    return(gsftrks_z0_);
  }
  std::vector<float>          &gsftrks_etaErr(){
    if(!gsftrks_etaErr_isLoaded) edm::LogWarning("VariableNotSet") << "gsftrks_etaErr not loaded!";
    return(gsftrks_etaErr_);
  }
  std::vector<LorentzVector>  &pfcands_p4(){
    if(!pfcands_p4_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_p4 not loaded!";
    return(pfcands_p4_);
  }
  std::vector<int>            &pfcands_trkidx(){
    if(!pfcands_trkidx_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_trkidx not loaded!";
    return(pfcands_trkidx_);
  }
  std::vector<int>            &pfcands_particleId(){
    if(!pfcands_particleId_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_particleId not loaded!";
    return(pfcands_particleId_);
  }
  std::vector<int>            &pfcands_pfelsidx(){
    if(!pfcands_pfelsidx_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_pfelsidx not loaded!";
    return(pfcands_pfelsidx_);
  }
  std::vector<int>            &pfcands_vtxidx(){
    if(!pfcands_vtxidx_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_vtxidx not loaded!";
    return(pfcands_vtxidx_);
  }
  std::vector<int>            &pfcands_charge(){
    if(!pfcands_charge_isLoaded) edm::LogWarning("VariableNotSet") << "pfcands_charge not loaded!";
    return(pfcands_charge_);
  }
  std::vector<int>            &pfels_elsidx(){
    if(!pfels_elsidx_isLoaded) edm::LogWarning("VariableNotSet") << "pfels_elsidx not loaded!";
    return(pfels_elsidx_);
  }
  std::vector<LorentzVector>  &pfels_p4(){
    if(!pfels_p4_isLoaded) edm::LogWarning("VariableNotSet") << "pfels_p4 not loaded!";
    return(pfels_p4_);
  }
  std::vector<LorentzVector>  &pfmus_p4(){
    if(!pfmus_p4_isLoaded) edm::LogWarning("VariableNotSet") << "pfmus_p4 not loaded!";
    return(pfmus_p4_);
  }
  vector<float>               &trk_met(){
    if(!trk_met_isLoaded) edm::LogWarning("VariableNotSet") << "trk_met not loaded!";
    return(trk_met_);
  }
  vector<float>               &trk_metPhi(){
    if(!trk_metPhi_isLoaded) edm::LogWarning("VariableNotSet") << "trk_metPhi not loaded!";
    return(trk_metPhi_);
  }
  vector<LorentzVector>       &pfjets_p4(){
    if(!pfjets_p4_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_p4 not loaded!";
    return(pfjets_p4_);
  }
  vector<LorentzVector>       &pfjets_corr_p4(){
    if(!pfjets_corr_p4_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_corr_p4 not loaded!";
    return(pfjets_corr_p4_);
  }
  vector<float>               &pfjets_area(){
    if(!pfjets_area_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_area not loaded!";
    return(pfjets_area_);
  }
  vector<float>               &pfjets_JEC(){
    if(!pfjets_JEC_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_JEC not loaded!";
    return(pfjets_JEC_);
  }
  vector<float>               &pfjets_mvavalue(){
    if(!pfjets_mvavalue_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_mvavalue not loaded!";
    return(pfjets_mvavalue_);
  }
  vector<float>               &pfjets_trackCountingHighEffBJetTag(){
    if(!pfjets_trackCountingHighEffBJetTag_isLoaded) edm::LogWarning("VariableNotSet") << "pfjets_trackCountingHighEffBJetTag not loaded!";
    return(pfjets_trackCountingHighEffBJetTag_);
  }

  void Reset();

};

HWW& GetHWW();

namespace HWWVal {

//to access data members
  //vertex 
  std::vector<LorentzVector>  &vtxs_position();
  std::vector<float>          &vtxs_ndof();
  std::vector<float>          &vtxs_sumpt();
  std::vector<int>            &vtxs_isFake();
  std::vector<float>          &vtxs_xError();
  std::vector<float>          &vtxs_yError();
  std::vector<float>          &vtxs_zError();
  std::vector<vector<float>>  &vtxs_covMatrix();

  //tracks
  std::vector<LorentzVector>  &trks_trk_p4();
  std::vector<float>          &trks_chi2();
  std::vector<float>          &trks_ndof();
  std::vector<float>          &trks_d0();
  std::vector<int>            &trks_nlayers();
  std::vector<int>            &trks_valid_pixelhits();
  std::vector<float>          &trks_z0();
  std::vector<float>          &trks_z0Err();
  std::vector<float>          &trks_etaErr();
  std::vector<float>          &trks_d0Err();
  std::vector<float>          &trks_phiErr();
  std::vector<float>          &trks_d0phiCov();
  std::vector<int>            &trks_qualityMask();
  std::vector<int>            &trks_charge();

  //electrons
  std::vector<LorentzVector>  &els_p4();
  std::vector<LorentzVector>  &els_trk_p4();
  std::vector<LorentzVector>  &els_vertex_p4();
  std::vector<float>          &els_lh();
  std::vector<float>          &els_etaSC();
  std::vector<float>          &els_sigmaIEtaIEta();
  std::vector<float>          &els_dEtaIn();
  std::vector<float>          &els_dPhiIn();
  std::vector<float>          &els_hOverE();
  std::vector<float>          &els_tkIso();
  std::vector<float>          &els_d0corr();
  std::vector<float>          &els_d0();
  std::vector<float>          &els_z0corr();
  std::vector<float>          &els_fbrem();
  std::vector<float>          &els_eOverPIn();
  std::vector<float>          &els_eSeedOverPOut();
  std::vector<float>          &els_eSeedOverPIn();
  std::vector<float>          &els_sigmaIPhiIPhi();
  std::vector<float>          &els_eSC();
  std::vector<float>          &els_ip3d();
  std::vector<float>          &els_ip3derr();
  std::vector<float>          &els_chi2();
  std::vector<float>          &els_ndof();
  std::vector<float>          &els_dEtaOut();
  std::vector<float>          &els_dPhiOut();
  std::vector<float>          &els_eSCRaw();
  std::vector<float>          &els_etaSCwidth();
  std::vector<float>          &els_phiSCwidth();
  std::vector<float>          &els_eSCPresh();
  std::vector<float>          &els_iso03_pf_ch();
  std::vector<float>          &els_iso03_pf_nhad05();
  std::vector<float>          &els_iso03_pf_gamma05();
  std::vector<float>          &els_iso04_pf_ch();
  std::vector<float>          &els_iso04_pf_nhad05();
  std::vector<float>          &els_iso04_pf_gamma05();
  std::vector<float>          &els_e5x5();
  std::vector<float>          &els_e1x5();
  std::vector<float>          &els_e3x3();
  std::vector<float>          &els_ecalEnergy();
  std::vector<float>          &els_eOverPOut();
  std::vector<float>          &els_ecalIso();
  std::vector<float>          &els_hcalIso();
  std::vector<float>          &els_trkshFrac();
  std::vector<float>          &els_conv_dist();
  std::vector<float>          &els_conv_dcot();
  std::vector<float>          &els_conv_old_dist();
  std::vector<float>          &els_conv_old_dcot();
  std::vector<float>          &els_iso04_pf2012_ch();
  std::vector<float>          &els_iso04_pf2012_em();
  std::vector<float>          &els_iso04_pf2012_nh();
  std::vector<float>          &els_iso03_pf2012_ch();
  std::vector<float>          &els_iso03_pf2012_em();
  std::vector<float>          &els_iso03_pf2012_nh();
  std::vector<float>          &els_ecalIso04();
  std::vector<float>          &els_hcalIso04();
  std::vector<int>            &els_nSeed();
  std::vector<int>            &els_scindex();
  std::vector<int>            &els_charge();
  std::vector<int>            &els_gsftrkidx();
  std::vector<int>            &els_exp_innerlayers();
  std::vector<int>            &els_trkidx();
  std::vector<int>            &els_type();
  std::vector<int>            &els_fiduciality();
  std::vector<int>            &els_sccharge();
  std::vector<int>            &els_trk_charge();
  std::vector<int>            &els_closestMuon();

  //muons
  std::vector<LorentzVector>  &mus_p4();
  std::vector<LorentzVector>  &mus_trk_p4();
  std::vector<LorentzVector>  &mus_vertex_p4();
  std::vector<LorentzVector>  &mus_sta_p4();
  std::vector<float>          &mus_gfit_chi2();
  std::vector<float>          &mus_gfit_ndof();
  std::vector<float>          &mus_ptErr();
  std::vector<float>          &mus_trkKink();
  std::vector<float>          &mus_d0corr();
  std::vector<float>          &mus_d0();
  std::vector<float>          &mus_z0corr();
  std::vector<float>          &mus_chi2();
  std::vector<float>          &mus_ndof();
  std::vector<float>          &mus_ip3d();
  std::vector<float>          &mus_ip3derr();
  std::vector<float>          &mus_segmCompatibility();
  std::vector<float>          &mus_caloCompatibility();
  std::vector<float>          &mus_e_had();
  std::vector<float>          &mus_e_ho();
  std::vector<float>          &mus_e_em();
  std::vector<float>          &mus_e_hadS9();
  std::vector<float>          &mus_e_hoS9();
  std::vector<float>          &mus_e_emS9();
  std::vector<float>          &mus_iso03_sumPt();
  std::vector<float>          &mus_iso03_emEt();
  std::vector<float>          &mus_iso03_hadEt();
  std::vector<float>          &mus_iso05_sumPt();
  std::vector<float>          &mus_iso05_emEt();
  std::vector<float>          &mus_iso05_hadEt();
  std::vector<float>          &mus_sta_d0();
  std::vector<float>          &mus_sta_z0corr();
  std::vector<float>          &mus_isoR03_pf_ChargedHadronPt();
  std::vector<float>          &mus_isoR03_pf_NeutralHadronEt();
  std::vector<float>          &mus_isoR03_pf_PhotonEt();
  std::vector<float>          &mus_isoR03_pf_PUPt();
  std::vector<float>          &mus_iso_ecalvetoDep();
  std::vector<float>          &mus_iso_hcalvetoDep();
  std::vector<int>            &mus_gfit_validSTAHits();
  std::vector<int>            &mus_numberOfMatchedStations();
  std::vector<int>            &mus_pfmusidx();
  std::vector<int>            &mus_charge();
  std::vector<int>            &mus_validHits();
  std::vector<int>            &mus_trkidx();
  std::vector<int>            &mus_pid_PFMuon();
  std::vector<int>            &mus_pid_TMLastStationTight();
  std::vector<int>            &mus_nmatches();
  std::vector<int>            &mus_goodmask();
  std::vector<int>            &mus_type();

  //dilepton hypothesis
  std::vector<vector<LorentzVector> > &hyp_jets_p4();
  std::vector<LorentzVector>  &hyp_p4();
  std::vector<LorentzVector>  &hyp_ll_p4();
  std::vector<LorentzVector>  &hyp_lt_p4();
  std::vector<int>            &hyp_ll_index();
  std::vector<int>            &hyp_lt_index();
  std::vector<int>            &hyp_ll_id();
  std::vector<int>            &hyp_lt_id();
  std::vector<int>            &hyp_ll_charge();
  std::vector<int>            &hyp_lt_charge();
  std::vector<int>            &hyp_type();

  //event variables
  unsigned int                &evt_run();
  unsigned int                &evt_lumiBlock();
  unsigned int                &evt_event();
  int                         &evt_isRealData();
  float                       &evt_ww_rho_vor();
  float                       &evt_ww_rho();
  float                       &evt_rho();
  float                       &evt_kt6pf_foregiso_rho();
  float                       &evt_pfmet();
  float                       &evt_pfmetPhi();


  std::vector<float>          &convs_ndof();
  std::vector<float>          &convs_chi2();
  std::vector<float>          &convs_dl();
  std::vector<int>            &convs_isConverted();
  std::vector<vector<int>>    &convs_tkalgo();
  std::vector<vector<int>>    &convs_tkidx();
  std::vector<vector<int>>    &convs_nHitsBeforeVtx();
  std::vector<int>            &convs_quality();
  std::vector<float>          &scs_sigmaIEtaIPhi();
  std::vector<LorentzVector>  &scs_pos_p4();
  std::vector<LorentzVector>  &gsftrks_p4();
  std::vector<LorentzVector>  &gsftrks_vertex_p4();
  std::vector<float>          &gsftrks_d0();
  std::vector<float>          &gsftrks_d0Err();
  std::vector<float>          &gsftrks_phiErr();
  std::vector<float>          &gsftrks_d0phiCov();
  std::vector<float>          &gsftrks_z0Err();
  std::vector<float>          &gsftrks_z0();
  std::vector<float>          &gsftrks_etaErr();
  std::vector<LorentzVector>  &pfcands_p4();
  std::vector<int>            &pfcands_trkidx();
  std::vector<int>            &pfcands_particleId();
  std::vector<int>            &pfcands_pfelsidx();
  std::vector<int>            &pfcands_vtxidx();
  std::vector<int>            &pfcands_charge();
  std::vector<int>            &pfels_elsidx();
  std::vector<LorentzVector>  &pfels_p4();
  std::vector<LorentzVector>  &pfmus_p4();
  vector<float>               &trk_met();
  vector<float>               &trk_metPhi();
  vector<LorentzVector>       &pfjets_p4();
  vector<LorentzVector>       &pfjets_corr_p4();
  vector<float>               &pfjets_area();
  vector<float>               &pfjets_JEC();
  vector<float>               &pfjets_mvavalue();
  vector<float>               &pfjets_trackCountingHighEffBJetTag();


//isLoaded
  void Load_vtxs_position();
  void Load_vtxs_ndof();
  void Load_vtxs_sumpt();
  void Load_vtxs_isFake();
  void Load_vtxs_xError();
  void Load_vtxs_yError();
  void Load_vtxs_zError();
  void Load_vtxs_covMatrix();
  void Load_trks_trk_p4();
  void Load_trks_chi2();
  void Load_trks_ndof();
  void Load_trks_d0();
  void Load_trks_nlayers();
  void Load_trks_valid_pixelhits();
  void Load_trks_z0();
  void Load_trks_z0Err();
  void Load_trks_etaErr();
  void Load_trks_d0Err();
  void Load_trks_phiErr();
  void Load_trks_d0phiCov();
  void Load_trks_qualityMask();
  void Load_trks_charge();
  void Load_els_p4();
  void Load_els_trk_p4();
  void Load_els_vertex_p4();
  void Load_els_lh();
  void Load_els_etaSC();
  void Load_els_sigmaIEtaIEta();
  void Load_els_dEtaIn();
  void Load_els_dPhiIn();
  void Load_els_hOverE();
  void Load_els_tkIso();
  void Load_els_d0corr();
  void Load_els_d0();
  void Load_els_z0corr();
  void Load_els_fbrem();
  void Load_els_eOverPIn();
  void Load_els_eSeedOverPOut();
  void Load_els_eSeedOverPIn();
  void Load_els_sigmaIPhiIPhi();
  void Load_els_eSC();
  void Load_els_ip3d();
  void Load_els_ip3derr();
  void Load_els_chi2();
  void Load_els_ndof();
  void Load_els_dEtaOut();
  void Load_els_dPhiOut();
  void Load_els_eSCRaw();
  void Load_els_etaSCwidth();
  void Load_els_phiSCwidth();
  void Load_els_eSCPresh();
  void Load_els_iso03_pf_ch();
  void Load_els_iso03_pf_nhad05();
  void Load_els_iso03_pf_gamma05();
  void Load_els_iso04_pf_ch();
  void Load_els_iso04_pf_nhad05();
  void Load_els_iso04_pf_gamma05();
  void Load_els_e5x5();
  void Load_els_e1x5();
  void Load_els_e3x3();
  void Load_els_ecalEnergy();
  void Load_els_eOverPOut();
  void Load_els_ecalIso();
  void Load_els_hcalIso();
  void Load_els_trkshFrac();
  void Load_els_conv_dist();
  void Load_els_conv_dcot();
  void Load_els_conv_old_dist();
  void Load_els_conv_old_dcot();
  void Load_els_iso04_pf2012_ch();
  void Load_els_iso04_pf2012_em();
  void Load_els_iso04_pf2012_nh();
  void Load_els_iso03_pf2012_ch();
  void Load_els_iso03_pf2012_em();
  void Load_els_iso03_pf2012_nh();
  void Load_els_ecalIso04();
  void Load_els_hcalIso04();
  void Load_els_nSeed();
  void Load_els_scindex();
  void Load_els_charge();
  void Load_els_gsftrkidx();
  void Load_els_exp_innerlayers();
  void Load_els_trkidx();
  void Load_els_type();
  void Load_els_fiduciality();
  void Load_els_sccharge();
  void Load_els_trk_charge();
  void Load_els_closestMuon();
  void Load_mus_p4();
  void Load_mus_trk_p4();
  void Load_mus_vertex_p4();
  void Load_mus_sta_p4();
  void Load_mus_gfit_chi2();
  void Load_mus_gfit_ndof();
  void Load_mus_ptErr();
  void Load_mus_trkKink();
  void Load_mus_d0corr();
  void Load_mus_d0();
  void Load_mus_z0corr();
  void Load_mus_chi2();
  void Load_mus_ndof();
  void Load_mus_ip3d();
  void Load_mus_ip3derr();
  void Load_mus_segmCompatibility();
  void Load_mus_caloCompatibility();
  void Load_mus_e_had();
  void Load_mus_e_ho();
  void Load_mus_e_em();
  void Load_mus_e_hadS9();
  void Load_mus_e_hoS9();
  void Load_mus_e_emS9();
  void Load_mus_iso03_sumPt();
  void Load_mus_iso03_emEt();
  void Load_mus_iso03_hadEt();
  void Load_mus_iso05_sumPt();
  void Load_mus_iso05_emEt();
  void Load_mus_iso05_hadEt();
  void Load_mus_sta_d0();
  void Load_mus_sta_z0corr();
  void Load_mus_isoR03_pf_ChargedHadronPt();
  void Load_mus_isoR03_pf_NeutralHadronEt();
  void Load_mus_isoR03_pf_PhotonEt();
  void Load_mus_isoR03_pf_PUPt();
  void Load_mus_iso_ecalvetoDep();
  void Load_mus_iso_hcalvetoDep();
  void Load_mus_gfit_validSTAHits();
  void Load_mus_numberOfMatchedStations();
  void Load_mus_pfmusidx();
  void Load_mus_charge();
  void Load_mus_validHits();
  void Load_mus_trkidx();
  void Load_mus_pid_PFMuon();
  void Load_mus_pid_TMLastStationTight();
  void Load_mus_nmatches();
  void Load_mus_goodmask();
  void Load_mus_type();
  void Load_hyp_jets_p4();
  void Load_hyp_p4();
  void Load_hyp_ll_p4();
  void Load_hyp_lt_p4();
  void Load_hyp_ll_index();
  void Load_hyp_lt_index();
  void Load_hyp_ll_id();
  void Load_hyp_lt_id();
  void Load_hyp_ll_charge();
  void Load_hyp_lt_charge();
  void Load_hyp_type();
  void Load_evt_run();
  void Load_evt_lumiBlock();
  void Load_evt_event();
  void Load_evt_isRealData();
  void Load_evt_ww_rho_vor();
  void Load_evt_ww_rho();
  void Load_evt_rho();
  void Load_evt_kt6pf_foregiso_rho();
  void Load_evt_pfmet();
  void Load_evt_pfmetPhi();
  void Load_convs_ndof();
  void Load_convs_chi2();
  void Load_convs_dl();
  void Load_convs_isConverted();
  void Load_convs_tkalgo();
  void Load_convs_tkidx();
  void Load_convs_nHitsBeforeVtx();
  void Load_convs_quality();
  void Load_scs_sigmaIEtaIPhi();
  void Load_scs_pos_p4();
  void Load_gsftrks_p4();
  void Load_gsftrks_vertex_p4();
  void Load_gsftrks_d0();
  void Load_gsftrks_d0Err();
  void Load_gsftrks_phiErr();
  void Load_gsftrks_d0phiCov();
  void Load_gsftrks_z0Err();
  void Load_gsftrks_z0();
  void Load_gsftrks_etaErr();
  void Load_pfcands_p4();
  void Load_pfcands_trkidx();
  void Load_pfcands_particleId();
  void Load_pfcands_pfelsidx();
  void Load_pfcands_vtxidx();
  void Load_pfcands_charge();
  void Load_pfels_elsidx();
  void Load_pfels_p4();
  void Load_pfmus_p4();
  void Load_trk_met();
  void Load_trk_metPhi();
  void Load_pfjets_p4();
  void Load_pfjets_corr_p4();
  void Load_pfjets_area();
  void Load_pfjets_JEC();
  void Load_pfjets_mvavalue();
  void Load_pfjets_trackCountingHighEffBJetTag();

  void Reset();

}//end namespace

#endif
