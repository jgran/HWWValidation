#ifndef HWW_H
#define HWW_H

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
  std::vector<float>          scs_e1x3_;
  std::vector<float>          scs_e3x1_;
  std::vector<float>          scs_eMax_;
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
  //define functions to access data members//
  ///////////////////////////////////////////


  //vertex 
  std::vector<LorentzVector>  &vtxs_position(){
    return(vtxs_position_);
  }
  std::vector<float>          &vtxs_ndof(){
    return(vtxs_ndof_);
  }
  std::vector<float>          &vtxs_sumpt(){
    return(vtxs_sumpt_);
  }
  std::vector<int>            &vtxs_isFake(){
    return(vtxs_isFake_);
  }
  std::vector<float>          &vtxs_xError(){
    return(vtxs_xError_);
  }
  std::vector<float>          &vtxs_yError(){
    return(vtxs_yError_);
  }
  std::vector<float>          &vtxs_zError(){
    return(vtxs_zError_);
  }
  std::vector<vector<float>>  &vtxs_covMatrix(){
    return(vtxs_covMatrix_);
  }

  //tracks
  std::vector<LorentzVector>  &trks_trk_p4(){
    return(trks_trk_p4_);
  }
  std::vector<float>          &trks_chi2(){
    return(trks_chi2_);
  }
  std::vector<float>          &trks_ndof(){
    return(trks_ndof_);
  }
  std::vector<float>          &trks_d0(){
    return(trks_d0_);
  }
  std::vector<int>            &trks_nlayers(){
    return(trks_nlayers_);
  }
  std::vector<int>            &trks_valid_pixelhits(){
    return(trks_valid_pixelhits_);
  }
  std::vector<float>          &trks_z0(){
    return(trks_z0_);
  }
  std::vector<float>          &trks_z0Err(){
    return(trks_z0Err_);
  }
  std::vector<float>          &trks_etaErr(){
    return(trks_etaErr_);
  }
  std::vector<float>          &trks_d0Err(){
    return(trks_d0Err_);
  }
  std::vector<float>          &trks_phiErr(){
    return(trks_phiErr_);
  }
  std::vector<float>          &trks_d0phiCov(){
    return(trks_d0phiCov_);
  }
  std::vector<int>            &trks_qualityMask(){
    return(trks_qualityMask_);
  }
  std::vector<int>            &trks_charge(){
    return(trks_charge_);
  }

  //electrons
  std::vector<LorentzVector>  &els_p4(){
    return(els_p4_);
  }
  std::vector<LorentzVector>  &els_trk_p4(){
    return(els_trk_p4_);
  }
  std::vector<LorentzVector>  &els_vertex_p4(){
    return(els_vertex_p4_);
  }
  std::vector<float>          &els_lh(){
    return(els_lh_);
  }
  std::vector<float>          &els_etaSC(){
    return(els_etaSC_);
  }
  std::vector<float>          &els_sigmaIEtaIEta(){
    return(els_sigmaIEtaIEta_);
  }
  std::vector<float>          &els_dEtaIn(){
    return(els_dEtaIn_);
  }
  std::vector<float>          &els_dPhiIn(){
    return(els_dPhiIn_);
  }
  std::vector<float>          &els_hOverE(){
    return(els_hOverE_);
  }
  std::vector<float>          &els_tkIso(){
    return(els_tkIso_);
  }
  std::vector<float>          &els_d0corr(){
    return(els_d0corr_);
  }
  std::vector<float>          &els_d0(){
    return(els_d0_);
  }
  std::vector<float>          &els_z0corr(){
    return(els_z0corr_);
  }
  std::vector<float>          &els_fbrem(){
    return(els_fbrem_);
  }
  std::vector<float>          &els_eOverPIn(){
    return(els_eOverPIn_);
  }
  std::vector<float>          &els_eSeedOverPOut(){
    return(els_eSeedOverPOut_);
  }
  std::vector<float>          &els_eSeedOverPIn(){
    return(els_eSeedOverPIn_);
  }
  std::vector<float>          &els_sigmaIPhiIPhi(){
    return(els_sigmaIPhiIPhi_);
  }
  std::vector<float>          &els_eSC(){
    return(els_eSC_);
  }
  std::vector<float>          &els_ip3d(){
    return(els_ip3d_);
  }
  std::vector<float>          &els_ip3derr(){
    return(els_ip3derr_);
  }
  std::vector<float>          &els_chi2(){
    return(els_chi2_);
  }
  std::vector<float>          &els_ndof(){
    return(els_ndof_);
  }
  std::vector<float>          &els_dEtaOut(){
    return(els_dEtaOut_);
  }
  std::vector<float>          &els_dPhiOut(){
    return(els_dPhiOut_);
  }
  std::vector<float>          &els_eSCRaw(){
    return(els_eSCRaw_);
  }
  std::vector<float>          &els_etaSCwidth(){
    return(els_etaSCwidth_);
  }
  std::vector<float>          &els_phiSCwidth(){
    return(els_phiSCwidth_);
  }
  std::vector<float>          &els_eSCPresh(){
    return(els_eSCPresh_);
  }
  std::vector<float>          &els_iso03_pf_ch(){
    return(els_iso03_pf_ch_);
  }
  std::vector<float>          &els_iso03_pf_nhad05(){
    return(els_iso03_pf_nhad05_);
  }
  std::vector<float>          &els_iso03_pf_gamma05(){
    return(els_iso03_pf_gamma05_);
  }
  std::vector<float>          &els_iso04_pf_ch(){
    return(els_iso04_pf_ch_);
  }
  std::vector<float>          &els_iso04_pf_nhad05(){
    return(els_iso04_pf_nhad05_);
  }
  std::vector<float>          &els_iso04_pf_gamma05(){
    return(els_iso04_pf_gamma05_);
  }
  std::vector<float>          &els_e5x5(){
    return(els_e5x5_);
  }
  std::vector<float>          &els_e1x5(){
    return(els_e1x5_);
  }
  std::vector<float>          &els_e3x3(){
    return(els_e3x3_);
  }
  std::vector<float>          &els_ecalEnergy(){
    return(els_ecalEnergy_);
  }
  std::vector<float>          &els_eOverPOut(){
    return(els_eOverPOut_);
  }
  std::vector<float>          &els_ecalIso(){
    return(els_ecalIso_);
  }
  std::vector<float>          &els_hcalIso(){
    return(els_hcalIso_);
  }
  std::vector<float>          &els_trkshFrac(){
    return(els_trkshFrac_);
  }
  std::vector<float>          &els_conv_dist(){
    return(els_conv_dist_);
  }
  std::vector<float>          &els_conv_dcot(){
    return(els_conv_dcot_);
  }
  std::vector<float>          &els_conv_old_dist(){
    return(els_conv_old_dist_);
  }
  std::vector<float>          &els_conv_old_dcot(){
    return(els_conv_old_dcot_);
  }
  std::vector<float>          &els_iso04_pf2012_ch(){
    return(els_iso04_pf2012_ch_);
  }
  std::vector<float>          &els_iso04_pf2012_em(){
    return(els_iso04_pf2012_em_);
  }
  std::vector<float>          &els_iso04_pf2012_nh(){
    return(els_iso04_pf2012_nh_);
  }
  std::vector<float>          &els_iso03_pf2012_ch(){
    return(els_iso03_pf2012_ch_);
  }
  std::vector<float>          &els_iso03_pf2012_em(){
    return(els_iso03_pf2012_em_);
  }
  std::vector<float>          &els_iso03_pf2012_nh(){
    return(els_iso03_pf2012_nh_);
  }
  std::vector<float>          &els_ecalIso04(){
    return(els_ecalIso04_);
  }
  std::vector<float>          &els_hcalIso04(){
    return(els_hcalIso04_);
  }
  std::vector<int>            &els_nSeed(){
    return(els_nSeed_);
  }
  std::vector<int>            &els_scindex(){
    return(els_scindex_);
  }
  std::vector<int>            &els_charge(){
    return(els_charge_);
  }
  std::vector<int>            &els_gsftrkidx(){
    return(els_gsftrkidx_);
  }
  std::vector<int>            &els_exp_innerlayers(){
    return(els_exp_innerlayers_);
  }
  std::vector<int>            &els_trkidx(){
    return(els_trkidx_);
  }
  std::vector<int>            &els_type(){
    return(els_type_);
  }
  std::vector<int>            &els_fiduciality(){
    return(els_fiduciality_);
  }
  std::vector<int>            &els_sccharge(){
    return(els_sccharge_);
  }
  std::vector<int>            &els_trk_charge(){
    return(els_trk_charge_);
  }
  std::vector<int>            &els_closestMuon(){
    return(els_closestMuon_);
  }

  //muons
  std::vector<LorentzVector>  &mus_p4(){
    return(mus_p4_);
  }
  std::vector<LorentzVector>  &mus_trk_p4(){
    return(mus_trk_p4_);
  }
  std::vector<LorentzVector>  &mus_vertex_p4(){
    return(mus_vertex_p4_);
  }
  std::vector<LorentzVector>  &mus_sta_p4(){
    return(mus_sta_p4_);
  }
  std::vector<float>          &mus_gfit_chi2(){
    return(mus_gfit_chi2_);
  }
  std::vector<float>          &mus_gfit_ndof(){
    return(mus_gfit_ndof_);
  }
  std::vector<float>          &mus_ptErr(){
    return(mus_ptErr_);
  }
  std::vector<float>          &mus_trkKink(){
    return(mus_trkKink_);
  }
  std::vector<float>          &mus_d0corr(){
    return(mus_d0corr_);
  }
  std::vector<float>          &mus_d0(){
    return(mus_d0_);
  }
  std::vector<float>          &mus_z0corr(){
    return(mus_z0corr_);
  }
  std::vector<float>          &mus_chi2(){
    return(mus_chi2_);
  }
  std::vector<float>          &mus_ndof(){
    return(mus_ndof_);
  }
  std::vector<float>          &mus_ip3d(){
    return(mus_ip3d_);
  }
  std::vector<float>          &mus_ip3derr(){
    return(mus_ip3derr_);
  }
  std::vector<float>          &mus_segmCompatibility(){
    return(mus_segmCompatibility_);
  }
  std::vector<float>          &mus_caloCompatibility(){
    return(mus_caloCompatibility_);
  }
  std::vector<float>          &mus_e_had(){
    return(mus_e_had_);
  }
  std::vector<float>          &mus_e_ho(){
    return(mus_e_ho_);
  }
  std::vector<float>          &mus_e_em(){
    return(mus_e_em_);
  }
  std::vector<float>          &mus_e_hadS9(){
    return(mus_e_hadS9_);
  }
  std::vector<float>          &mus_e_hoS9(){
    return(mus_e_hoS9_);
  }
  std::vector<float>          &mus_e_emS9(){
    return(mus_e_emS9_);
  }
  std::vector<float>          &mus_iso03_sumPt(){
    return(mus_iso03_sumPt_);
  }
  std::vector<float>          &mus_iso03_emEt(){
    return(mus_iso03_emEt_);
  }
  std::vector<float>          &mus_iso03_hadEt(){
    return(mus_iso03_hadEt_);
  }
  std::vector<float>          &mus_iso05_sumPt(){
    return(mus_iso05_sumPt_);
  }
  std::vector<float>          &mus_iso05_emEt(){
    return(mus_iso05_emEt_);
  }
  std::vector<float>          &mus_iso05_hadEt(){
    return(mus_iso05_hadEt_);
  }
  std::vector<float>          &mus_sta_d0(){
    return(mus_sta_d0_);
  }
  std::vector<float>          &mus_sta_z0corr(){
    return(mus_sta_z0corr_);
  }
  std::vector<float>          &mus_isoR03_pf_ChargedHadronPt(){
    return(mus_isoR03_pf_ChargedHadronPt_);
  }
  std::vector<float>          &mus_isoR03_pf_NeutralHadronEt(){
    return(mus_isoR03_pf_NeutralHadronEt_);
  }
  std::vector<float>          &mus_isoR03_pf_PhotonEt(){
    return(mus_isoR03_pf_PhotonEt_);
  }
  std::vector<float>          &mus_isoR03_pf_PUPt(){
    return(mus_isoR03_pf_PUPt_);
  }
  std::vector<float>          &mus_iso_ecalvetoDep(){
    return(mus_iso_ecalvetoDep_);
  }
  std::vector<float>          &mus_iso_hcalvetoDep(){
    return(mus_iso_hcalvetoDep_);
  }
  std::vector<int>            &mus_gfit_validSTAHits(){
    return(mus_gfit_validSTAHits_);
  }
  std::vector<int>            &mus_numberOfMatchedStations(){
    return(mus_numberOfMatchedStations_);
  }
  std::vector<int>            &mus_pfmusidx(){
    return(mus_pfmusidx_);
  }
  std::vector<int>            &mus_charge(){
    return(mus_charge_);
  }
  std::vector<int>            &mus_validHits(){
    return(mus_validHits_);
  }
  std::vector<int>            &mus_trkidx(){
    return(mus_trkidx_);
  }
  std::vector<int>            &mus_pid_PFMuon(){
    return(mus_pid_PFMuon_);
  }
  std::vector<int>            &mus_pid_TMLastStationTight(){
    return(mus_pid_TMLastStationTight_);
  }
  std::vector<int>            &mus_nmatches(){
    return(mus_nmatches_);
  }
  std::vector<int>            &mus_goodmask(){
    return(mus_goodmask_);
  }
  std::vector<int>            &mus_type(){
    return(mus_type_);
  }

  //dilepton hypothesis
  std::vector<LorentzVector>  &hyp_p4(){
    return(hyp_p4_);
  }
  std::vector<LorentzVector>  &hyp_ll_p4(){
    return(hyp_ll_p4_);
  }
  std::vector<LorentzVector>  &hyp_lt_p4(){
    return(hyp_lt_p4_);
  }
  std::vector<int>            &hyp_ll_index(){
    return(hyp_ll_index_);
  }
  std::vector<int>            &hyp_lt_index(){
    return(hyp_lt_index_);
  }
  std::vector<int>            &hyp_ll_id(){
    return(hyp_ll_id_);
  }
  std::vector<int>            &hyp_lt_id(){
    return(hyp_lt_id_);
  }
  std::vector<int>            &hyp_ll_charge(){
    return(hyp_ll_charge_);
  }
  std::vector<int>            &hyp_lt_charge(){
    return(hyp_lt_charge_);
  }
  std::vector<int>            &hyp_type(){
    return(hyp_type_);
  }

  //event variables
  unsigned int                &evt_run(){
    return(evt_run_);
  }
  unsigned int                &evt_lumiBlock(){
    return(evt_lumiBlock_);
  }
  unsigned int                &evt_event(){
    return(evt_event_);
  }
  int                         &evt_isRealData(){
    return(evt_isRealData_);
  }
  float                       &evt_ww_rho_vor(){
    return(evt_ww_rho_vor_);
  }
  float                       &evt_ww_rho(){
    return(evt_ww_rho_);
  }
  float                       &evt_rho(){
    return(evt_rho_);
  }
  float                       &evt_kt6pf_foregiso_rho(){
    return(evt_kt6pf_foregiso_rho_);
  }
  float                       &evt_pfmet(){
    return(evt_pfmet_);
  }
  float                       &evt_pfmetPhi(){
    return(evt_pfmetPhi_);
  }

  std::vector<float>          &convs_ndof(){
    return(convs_ndof_);
  }
  std::vector<float>          &convs_chi2(){
    return(convs_chi2_);
  }
  std::vector<float>          &convs_dl(){
    return(convs_dl_);
  }
  std::vector<int>            &convs_isConverted(){
    return(convs_isConverted_);
  }
  std::vector<vector<int>>    &convs_tkalgo(){
    return(convs_tkalgo_);
  }
  std::vector<vector<int>>    &convs_tkidx(){
    return(convs_tkidx_);
  }
  std::vector<vector<int>>    &convs_nHitsBeforeVtx(){
    return(convs_nHitsBeforeVtx_);
  }
  std::vector<int>            &convs_quality(){
    return(convs_quality_);
  }
  std::vector<float>          &scs_sigmaIEtaIPhi(){
    return(scs_sigmaIEtaIPhi_);
  }
  std::vector<float>          &scs_e1x3(){
    return(scs_e1x3_);
  }
  std::vector<float>          &scs_e3x1(){
    return(scs_e3x1_);
  }
  std::vector<float>          &scs_eMax(){
    return(scs_eMax_);
  }
  std::vector<LorentzVector>  &scs_pos_p4(){
    return(scs_pos_p4_);
  }
  std::vector<LorentzVector>  &gsftrks_p4(){
    return(gsftrks_p4_);
  }
  std::vector<LorentzVector>  &gsftrks_vertex_p4(){
    return(gsftrks_vertex_p4_);
  }
  std::vector<float>          &gsftrks_d0(){
    return(gsftrks_d0_);
  }
  std::vector<float>          &gsftrks_d0Err(){
    return(gsftrks_d0Err_);
  }
  std::vector<float>          &gsftrks_phiErr(){
    return(gsftrks_phiErr_);
  }
  std::vector<float>          &gsftrks_d0phiCov(){
    return(gsftrks_d0phiCov_);
  }
  std::vector<float>          &gsftrks_z0Err(){
    return(gsftrks_z0Err_);
  }
  std::vector<float>          &gsftrks_z0(){
    return(gsftrks_z0_);
  }
  std::vector<float>          &gsftrks_etaErr(){
    return(gsftrks_etaErr_);
  }
  std::vector<LorentzVector>  &pfcands_p4(){
    return(pfcands_p4_);
  }
  std::vector<int>            &pfcands_trkidx(){
    return(pfcands_trkidx_);
  }
  std::vector<int>            &pfcands_particleId(){
    return(pfcands_particleId_);
  }
  std::vector<int>            &pfcands_pfelsidx(){
    return(pfcands_pfelsidx_);
  }
  std::vector<int>            &pfcands_vtxidx(){
    return(pfcands_vtxidx_);
  }
  std::vector<int>            &pfcands_charge(){
    return(pfcands_charge_);
  }
  std::vector<int>            &pfels_elsidx(){
    return(pfels_elsidx_);
  }
  std::vector<LorentzVector>  &pfmus_p4(){
    return(pfmus_p4_);
  }
  vector<float>               &trk_met(){
    return(trk_met_);
  }
  vector<float>               &trk_metPhi(){
    return(trk_metPhi_);
  }
  vector<LorentzVector>       &pfjets_p4(){
    return(pfjets_p4_);
  }
  vector<LorentzVector>       &pfjets_corr_p4(){
    return(pfjets_corr_p4_);
  }
  vector<float>               &pfjets_area(){
    return(pfjets_area_);
  }
  vector<float>               &pfjets_JEC(){
    return(pfjets_JEC_);
  }
  vector<float>               &pfjets_mvavalue(){
    return(pfjets_mvavalue_);
  }
  vector<float>               &pfjets_trackCountingHighEffBJetTag(){
    return(pfjets_trackCountingHighEffBJetTag_);
  }


};

extern HWW hww;

#endif
