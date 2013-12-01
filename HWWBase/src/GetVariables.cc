#include "HWWValidation/HWWBase/interface/HWWAnalyzer.h"
#include "HWWValidation/HWWBase/interface/HWW.h"



void HWWAnalyzer::GetFirstVariables(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::InputTag hyp_p4_tag(hypInputTag.label(),"hypp4");
  edm::Handle<std::vector<LorentzVector> > hyp_p4_h;
  iEvent.getByLabel(hyp_p4_tag, hyp_p4_h);
  hww.hyp_p4_ = *hyp_p4_h.product();

  edm::InputTag hyp_ll_p4_tag(hypInputTag.label(),"hypllp4");
  edm::Handle<std::vector<LorentzVector> > hyp_ll_p4_h;
  iEvent.getByLabel(hyp_ll_p4_tag, hyp_ll_p4_h);
  hww.hyp_ll_p4_ = *hyp_ll_p4_h.product();

  edm::InputTag hyp_lt_p4_tag(hypInputTag.label(),"hypltp4");
  edm::Handle<std::vector<LorentzVector> > hyp_lt_p4_h;
  iEvent.getByLabel(hyp_lt_p4_tag, hyp_lt_p4_h);
  hww.hyp_lt_p4_ = *hyp_lt_p4_h.product();

  edm::InputTag hyp_ll_index_tag(hypInputTag.label(),"hypllindex");
  edm::Handle<std::vector<int> > hyp_ll_index_h;
  iEvent.getByLabel(hyp_ll_index_tag, hyp_ll_index_h);
  hww.hyp_ll_index_ = *hyp_ll_index_h.product();

  edm::InputTag hyp_lt_index_tag(hypInputTag.label(),"hypltindex");
  edm::Handle<std::vector<int> > hyp_lt_index_h;
  iEvent.getByLabel(hyp_lt_index_tag, hyp_lt_index_h);
  hww.hyp_lt_index_ = *hyp_lt_index_h.product();

  edm::InputTag hyp_ll_id_tag(hypInputTag.label(),"hypllid");
  edm::Handle<std::vector<int> > hyp_ll_id_h;
  iEvent.getByLabel(hyp_ll_id_tag, hyp_ll_id_h);
  hww.hyp_ll_id_ = *hyp_ll_id_h.product();

  edm::InputTag hyp_lt_id_tag(hypInputTag.label(),"hypltid");
  edm::Handle<std::vector<int> > hyp_lt_id_h;
  iEvent.getByLabel(hyp_lt_id_tag, hyp_lt_id_h);
  hww.hyp_lt_id_ = *hyp_lt_id_h.product();

  edm::InputTag vtxs_position_tag(vertexInputTag.label(),"vtxsposition");
  edm::Handle<std::vector<LorentzVector> > vtxs_position_h;
  iEvent.getByLabel(vtxs_position_tag, vtxs_position_h);
  hww.vtxs_position_ = *vtxs_position_h.product();

  edm::InputTag vtxs_ndof_tag(vertexInputTag.label(),"vtxsndof");
  edm::Handle<std::vector<float> > vtxs_ndof_h;
  iEvent.getByLabel(vtxs_ndof_tag, vtxs_ndof_h);
  hww.vtxs_ndof_ = *vtxs_ndof_h.product();

  edm::InputTag vtxs_isFake_tag(vertexInputTag.label(),"vtxsisFake");
  edm::Handle<std::vector<int> > vtxs_isFake_h;
  iEvent.getByLabel(vtxs_isFake_tag, vtxs_isFake_h);
  hww.vtxs_isFake_ = *vtxs_isFake_h.product();

  edm::InputTag vtxs_sumpt_tag(vertexInputTag.label(),"vtxssumpt");
  edm::Handle<std::vector<float> > vtxs_sumpt_h;
  iEvent.getByLabel(vtxs_sumpt_tag, vtxs_sumpt_h);
  hww.vtxs_sumpt_ = *vtxs_sumpt_h.product();

  edm::InputTag trks_d0_tag(trackInputTag.label(),"trksd0");
  edm::Handle<std::vector<float> > trks_d0_h;
  iEvent.getByLabel(trks_d0_tag, trks_d0_h);
  hww.trks_d0_ = *trks_d0_h.product();

}





void HWWAnalyzer::GetVariables(const edm::Event& iEvent, const edm::EventSetup& iSetup){

//get event variables
  edm::InputTag evt_run_tag(eventInputTag.label(),"evtrun");
  edm::Handle<unsigned int> evt_run_h;
  iEvent.getByLabel(evt_run_tag, evt_run_h);
  hww.evt_run_ = *evt_run_h.product();

  edm::InputTag evt_lumiBlock_tag(eventInputTag.label(),"evtlumiBlock");
  edm::Handle<unsigned int> evt_lumiBlock_h;
  iEvent.getByLabel(evt_lumiBlock_tag, evt_lumiBlock_h);
  hww.evt_lumiBlock_ = *evt_lumiBlock_h.product();

  edm::InputTag evt_event_tag(eventInputTag.label(),"evtevent");
  edm::Handle<unsigned int> evt_event_h;
  iEvent.getByLabel(evt_event_tag, evt_event_h);
  hww.evt_event_ = *evt_event_h.product();

  edm::InputTag evt_isRealData_tag(eventInputTag.label(),"evtisRealData");
  edm::Handle<int> evt_isRealData_h;
  iEvent.getByLabel(evt_isRealData_tag, evt_isRealData_h);
  hww.evt_isRealData_ = *evt_isRealData_h.product();

//get vertex variables
  edm::InputTag vtxs_xError_tag(vertexInputTag.label(),"vtxsxError");
  edm::Handle<std::vector<float> > vtxs_xError_h;
  iEvent.getByLabel(vtxs_xError_tag, vtxs_xError_h);
  hww.vtxs_xError_ = *vtxs_xError_h.product();

  edm::InputTag vtxs_yError_tag(vertexInputTag.label(),"vtxsyError");
  edm::Handle<std::vector<float> > vtxs_yError_h;
  iEvent.getByLabel(vtxs_yError_tag, vtxs_yError_h);
  hww.vtxs_yError_ = *vtxs_yError_h.product();

  edm::InputTag vtxs_zError_tag(vertexInputTag.label(),"vtxszError");
  edm::Handle<std::vector<float> > vtxs_zError_h;
  iEvent.getByLabel(vtxs_zError_tag, vtxs_zError_h);
  hww.vtxs_zError_ = *vtxs_zError_h.product();

  edm::InputTag vtxs_covMatrix_tag(vertexInputTag.label(),"vtxscovMatrix");
  edm::Handle<std::vector<vector<float>> > vtxs_covMatrix_h;
  iEvent.getByLabel(vtxs_covMatrix_tag, vtxs_covMatrix_h);
  hww.vtxs_covMatrix_ = *vtxs_covMatrix_h.product();

//get track variables
  edm::InputTag trks_trk_p4_tag(trackInputTag.label(),"trkstrkp4");
  edm::Handle<std::vector<LorentzVector> > trks_trk_p4_h;
  iEvent.getByLabel(trks_trk_p4_tag, trks_trk_p4_h);
  hww.trks_trk_p4_ = *trks_trk_p4_h.product();

  edm::InputTag trks_chi2_tag(trackInputTag.label(),"trkschi2");
  edm::Handle<std::vector<float> > trks_chi2_h;
  iEvent.getByLabel(trks_chi2_tag, trks_chi2_h);
  hww.trks_chi2_ = *trks_chi2_h.product();

  edm::InputTag trks_ndof_tag(trackInputTag.label(),"trksndof");
  edm::Handle<std::vector<float> > trks_ndof_h;
  iEvent.getByLabel(trks_ndof_tag, trks_ndof_h);
  hww.trks_ndof_ = *trks_ndof_h.product();

  edm::InputTag trks_nlayers_tag(trackInputTag.label(),"trksnlayers");
  edm::Handle<std::vector<int> > trks_nlayers_h;
  iEvent.getByLabel(trks_nlayers_tag, trks_nlayers_h);
  hww.trks_nlayers_ = *trks_nlayers_h.product();

  edm::InputTag trks_valid_pixelhits_tag(trackInputTag.label(),"trksvalidpixelhits");
  edm::Handle<std::vector<int> > trks_valid_pixelhits_h;
  iEvent.getByLabel(trks_valid_pixelhits_tag, trks_valid_pixelhits_h);
  hww.trks_valid_pixelhits_ = *trks_valid_pixelhits_h.product();

  edm::InputTag trks_z0_tag(trackInputTag.label(),"trksz0");
  edm::Handle<std::vector<float> > trks_z0_h;
  iEvent.getByLabel(trks_z0_tag, trks_z0_h);
  hww.trks_z0_ = *trks_z0_h.product();

  edm::InputTag trks_z0Err_tag(trackInputTag.label(),"trksz0Err");
  edm::Handle<std::vector<float> > trks_z0Err_h;
  iEvent.getByLabel(trks_z0Err_tag, trks_z0Err_h);
  hww.trks_z0Err_ = *trks_z0Err_h.product();

  edm::InputTag trks_etaErr_tag(trackInputTag.label(),"trksetaErr");
  edm::Handle<std::vector<float> > trks_etaErr_h;
  iEvent.getByLabel(trks_etaErr_tag, trks_etaErr_h);
  hww.trks_etaErr_ = *trks_etaErr_h.product();

  edm::InputTag trks_d0Err_tag(trackInputTag.label(),"trksd0Err");
  edm::Handle<std::vector<float> > trks_d0Err_h;
  iEvent.getByLabel(trks_d0Err_tag, trks_d0Err_h);
  hww.trks_d0Err_ = *trks_d0Err_h.product();

  edm::InputTag trks_phiErr_tag(trackInputTag.label(),"trksphiErr");
  edm::Handle<std::vector<float> > trks_phiErr_h;
  iEvent.getByLabel(trks_phiErr_tag, trks_phiErr_h);
  hww.trks_phiErr_ = *trks_phiErr_h.product();

  edm::InputTag trks_d0phiCov_tag(trackInputTag.label(),"trksd0phiCov");
  edm::Handle<std::vector<float> > trks_d0phiCov_h;
  iEvent.getByLabel(trks_d0phiCov_tag, trks_d0phiCov_h);
  hww.trks_d0phiCov_ = *trks_d0phiCov_h.product();

  edm::InputTag trks_qualityMask_tag(trackInputTag.label(),"trksqualityMask");
  edm::Handle<std::vector<int> > trks_qualityMask_h;
  iEvent.getByLabel(trks_qualityMask_tag, trks_qualityMask_h);
  hww.trks_qualityMask_ = *trks_qualityMask_h.product();

  edm::InputTag trks_charge_tag(trackInputTag.label(),"trkscharge");
  edm::Handle<std::vector<int> > trks_charge_h;
  iEvent.getByLabel(trks_charge_tag, trks_charge_h);
  hww.trks_charge_ = *trks_charge_h.product();


//get electron variables
  edm::InputTag els_p4_tag(electronInputTag.label(),"elsp4");
  edm::Handle<std::vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel(els_p4_tag, els_p4_h);
  hww.els_p4_ = *els_p4_h.product();

  edm::InputTag els_trk_p4_tag(electronInputTag.label(),"elstrkp4");
  edm::Handle<std::vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel(els_trk_p4_tag, els_trk_p4_h);
  hww.els_trk_p4_ = *els_trk_p4_h.product();

  edm::InputTag els_trkidx_tag(electronInputTag.label(),"elstrkidx");
  edm::Handle<std::vector<int> > els_trkidx_h;
  iEvent.getByLabel(els_trkidx_tag, els_trkidx_h);
  hww.els_trkidx_ = *els_trkidx_h.product();

  edm::InputTag els_vertex_p4_tag(electronInputTag.label(),"elsvertexp4");
  edm::Handle<std::vector<LorentzVector> > els_vertex_p4_h;
  iEvent.getByLabel(els_vertex_p4_tag, els_vertex_p4_h);
  hww.els_vertex_p4_ = *els_vertex_p4_h.product();

  edm::InputTag els_lh_tag(electronInputTag.label(),"elslh");
  edm::Handle<std::vector<float> > els_lh_h;
  iEvent.getByLabel(els_lh_tag, els_lh_h);
  hww.els_lh_ = *els_lh_h.product();

  edm::InputTag els_etaSC_tag(electronInputTag.label(),"elsetaSC");
  edm::Handle<std::vector<float> > els_etaSC_h;
  iEvent.getByLabel(els_etaSC_tag, els_etaSC_h);
  hww.els_etaSC_ = *els_etaSC_h.product();

  edm::InputTag els_sigmaIEtaIEta_tag(electronInputTag.label(),"elssigmaIEtaIEta");
  edm::Handle<std::vector<float> > els_sigmaIEtaIEta_h;
  iEvent.getByLabel(els_sigmaIEtaIEta_tag, els_sigmaIEtaIEta_h);
  hww.els_sigmaIEtaIEta_ = *els_sigmaIEtaIEta_h.product();

  edm::InputTag els_dEtaIn_tag(electronInputTag.label(),"elsdEtaIn");
  edm::Handle<std::vector<float> > els_dEtaIn_h;
  iEvent.getByLabel(els_dEtaIn_tag, els_dEtaIn_h);
  hww.els_dEtaIn_ = *els_dEtaIn_h.product();

  edm::InputTag els_dPhiIn_tag(electronInputTag.label(),"elsdPhiIn");
  edm::Handle<std::vector<float> > els_dPhiIn_h;
  iEvent.getByLabel(els_dPhiIn_tag, els_dPhiIn_h);
  hww.els_dPhiIn_ = *els_dPhiIn_h.product();

  edm::InputTag els_hOverE_tag(electronInputTag.label(),"elshOverE");
  edm::Handle<std::vector<float> > els_hOverE_h;
  iEvent.getByLabel(els_hOverE_tag, els_hOverE_h);
  hww.els_hOverE_ = *els_hOverE_h.product();

  edm::InputTag els_tkIso_tag(electronInputTag.label(),"elstkIso");
  edm::Handle<std::vector<float> > els_tkIso_h;
  iEvent.getByLabel(els_tkIso_tag, els_tkIso_h);
  hww.els_tkIso_ = *els_tkIso_h.product();

  edm::InputTag els_d0corr_tag(electronInputTag.label(),"elsd0corr");
  edm::Handle<std::vector<float> > els_d0corr_h;
  iEvent.getByLabel(els_d0corr_tag, els_d0corr_h);
  hww.els_d0corr_ = *els_d0corr_h.product();

  edm::InputTag els_d0_tag(electronInputTag.label(),"elsd0");
  edm::Handle<std::vector<float> > els_d0_h;
  iEvent.getByLabel(els_d0_tag, els_d0_h);
  hww.els_d0_ = *els_d0_h.product();

  edm::InputTag els_z0corr_tag(electronInputTag.label(),"elsz0corr");
  edm::Handle<std::vector<float> > els_z0corr_h;
  iEvent.getByLabel(els_z0corr_tag, els_z0corr_h);
  hww.els_z0corr_ = *els_z0corr_h.product();

  edm::InputTag els_fbrem_tag(electronInputTag.label(),"elsfbrem");
  edm::Handle<std::vector<float> > els_fbrem_h;
  iEvent.getByLabel(els_fbrem_tag, els_fbrem_h);
  hww.els_fbrem_ = *els_fbrem_h.product();

  edm::InputTag els_eOverPIn_tag(electronInputTag.label(),"elseOverPIn");
  edm::Handle<std::vector<float> > els_eOverPIn_h;
  iEvent.getByLabel(els_eOverPIn_tag, els_eOverPIn_h);
  hww.els_eOverPIn_ = *els_eOverPIn_h.product();

  edm::InputTag els_eSeedOverPOut_tag(electronInputTag.label(),"elseSeedOverPOut");
  edm::Handle<std::vector<float> > els_eSeedOverPOut_h;
  iEvent.getByLabel(els_eSeedOverPOut_tag, els_eSeedOverPOut_h);
  hww.els_eSeedOverPOut_ = *els_eSeedOverPOut_h.product();

  edm::InputTag els_eSeedOverPIn_tag(electronInputTag.label(),"elseSeedOverPIn");
  edm::Handle<std::vector<float> > els_eSeedOverPIn_h;
  iEvent.getByLabel(els_eSeedOverPIn_tag, els_eSeedOverPIn_h);
  hww.els_eSeedOverPIn_ = *els_eSeedOverPIn_h.product();

  edm::InputTag els_sigmaIPhiIPhi_tag(electronInputTag.label(),"elssigmaIPhiIPhi");
  edm::Handle<std::vector<float> > els_sigmaIPhiIPhi_h;
  iEvent.getByLabel(els_sigmaIPhiIPhi_tag, els_sigmaIPhiIPhi_h);
  hww.els_sigmaIPhiIPhi_ = *els_sigmaIPhiIPhi_h.product();

  edm::InputTag els_eSC_tag(electronInputTag.label(),"elseSC");
  edm::Handle<std::vector<float> > els_eSC_h;
  iEvent.getByLabel(els_eSC_tag, els_eSC_h);
  hww.els_eSC_ = *els_eSC_h.product();

  edm::InputTag els_ip3d_tag(electronInputTag.label(),"elsip3d");
  edm::Handle<std::vector<float> > els_ip3d_h;
  iEvent.getByLabel(els_ip3d_tag, els_ip3d_h);
  hww.els_ip3d_ = *els_ip3d_h.product();

  edm::InputTag els_ip3derr_tag(electronInputTag.label(),"elsip3derr");
  edm::Handle<std::vector<float> > els_ip3derr_h;
  iEvent.getByLabel(els_ip3derr_tag, els_ip3derr_h);
  hww.els_ip3derr_ = *els_ip3derr_h.product();

  edm::InputTag els_chi2_tag(electronInputTag.label(),"elschi2");
  edm::Handle<std::vector<float> > els_chi2_h;
  iEvent.getByLabel(els_chi2_tag, els_chi2_h);
  hww.els_chi2_ = *els_chi2_h.product();

  edm::InputTag els_ndof_tag(electronInputTag.label(),"elsndof");
  edm::Handle<std::vector<float> > els_ndof_h;
  iEvent.getByLabel(els_ndof_tag, els_ndof_h);
  hww.els_ndof_ = *els_ndof_h.product();

  edm::InputTag els_dEtaOut_tag(electronInputTag.label(),"elsdEtaOut");
  edm::Handle<std::vector<float> > els_dEtaOut_h;
  iEvent.getByLabel(els_dEtaOut_tag, els_dEtaOut_h);
  hww.els_dEtaOut_ = *els_dEtaOut_h.product();

  edm::InputTag els_dPhiOut_tag(electronInputTag.label(),"elsdPhiOut");
  edm::Handle<std::vector<float> > els_dPhiOut_h;
  iEvent.getByLabel(els_dPhiOut_tag, els_dPhiOut_h);
  hww.els_dPhiOut_ = *els_dPhiOut_h.product();

  edm::InputTag els_eSCRaw_tag(electronInputTag.label(),"elseSCRaw");
  edm::Handle<std::vector<float> > els_eSCRaw_h;
  iEvent.getByLabel(els_eSCRaw_tag, els_eSCRaw_h);
  hww.els_eSCRaw_ = *els_eSCRaw_h.product();

  edm::InputTag els_etaSCwidth_tag(electronInputTag.label(),"elsetaSCwidth");
  edm::Handle<std::vector<float> > els_etaSCwidth_h;
  iEvent.getByLabel(els_etaSCwidth_tag, els_etaSCwidth_h);
  hww.els_etaSCwidth_ = *els_etaSCwidth_h.product();

  edm::InputTag els_phiSCwidth_tag(electronInputTag.label(),"elsphiSCwidth");
  edm::Handle<std::vector<float> > els_phiSCwidth_h;
  iEvent.getByLabel(els_phiSCwidth_tag, els_phiSCwidth_h);
  hww.els_phiSCwidth_ = *els_phiSCwidth_h.product();

  edm::InputTag els_eSCPresh_tag(electronInputTag.label(),"elseSCPresh");
  edm::Handle<std::vector<float> > els_eSCPresh_h;
  iEvent.getByLabel(els_eSCPresh_tag, els_eSCPresh_h);
  hww.els_eSCPresh_ = *els_eSCPresh_h.product();

  edm::InputTag els_iso03_pf_ch_tag(electronInputTag.label(),"elsiso03pfch");
  edm::Handle<std::vector<float> > els_iso03_pf_ch_h;
  iEvent.getByLabel(els_iso03_pf_ch_tag, els_iso03_pf_ch_h);
  hww.els_iso03_pf_ch_ = *els_iso03_pf_ch_h.product();

  edm::InputTag els_iso03_pf_nhad05_tag(electronInputTag.label(),"elsiso03pfnhad05");
  edm::Handle<std::vector<float> > els_iso03_pf_nhad05_h;
  iEvent.getByLabel(els_iso03_pf_nhad05_tag, els_iso03_pf_nhad05_h);
  hww.els_iso03_pf_nhad05_ = *els_iso03_pf_nhad05_h.product();

  edm::InputTag els_iso03_pf_gamma05_tag(electronInputTag.label(),"elsiso03pfgamma05");
  edm::Handle<std::vector<float> > els_iso03_pf_gamma05_h;
  iEvent.getByLabel(els_iso03_pf_gamma05_tag, els_iso03_pf_gamma05_h);
  hww.els_iso03_pf_gamma05_ = *els_iso03_pf_gamma05_h.product();

  edm::InputTag els_iso04_pf_ch_tag(electronInputTag.label(),"elsiso04pfch");
  edm::Handle<std::vector<float> > els_iso04_pf_ch_h;
  iEvent.getByLabel(els_iso04_pf_ch_tag, els_iso04_pf_ch_h);
  hww.els_iso04_pf_ch_ = *els_iso04_pf_ch_h.product();

  edm::InputTag els_iso04_pf_nhad05_tag(electronInputTag.label(),"elsiso04pfnhad05");
  edm::Handle<std::vector<float> > els_iso04_pf_nhad05_h;
  iEvent.getByLabel(els_iso04_pf_nhad05_tag, els_iso04_pf_nhad05_h);
  hww.els_iso04_pf_nhad05_ = *els_iso04_pf_nhad05_h.product();

  edm::InputTag els_iso04_pf_gamma05_tag(electronInputTag.label(),"elsiso04pfgamma05");
  edm::Handle<std::vector<float> > els_iso04_pf_gamma05_h;
  iEvent.getByLabel(els_iso04_pf_gamma05_tag, els_iso04_pf_gamma05_h);
  hww.els_iso04_pf_gamma05_ = *els_iso04_pf_gamma05_h.product();

  edm::InputTag els_e5x5_tag(electronInputTag.label(),"else5x5");
  edm::Handle<std::vector<float> > els_e5x5_h;
  iEvent.getByLabel(els_e5x5_tag, els_e5x5_h);
  hww.els_e5x5_ = *els_e5x5_h.product();

  edm::InputTag els_e1x5_tag(electronInputTag.label(),"else1x5");
  edm::Handle<std::vector<float> > els_e1x5_h;
  iEvent.getByLabel(els_e1x5_tag, els_e1x5_h);
  hww.els_e1x5_ = *els_e1x5_h.product();

  edm::InputTag els_e3x3_tag(electronInputTag.label(),"else3x3");
  edm::Handle<std::vector<float> > els_e3x3_h;
  iEvent.getByLabel(els_e3x3_tag, els_e3x3_h);
  hww.els_e3x3_ = *els_e3x3_h.product();

  edm::InputTag els_ecalEnergy_tag(electronInputTag.label(),"elsecalEnergy");
  edm::Handle<std::vector<float> > els_ecalEnergy_h;
  iEvent.getByLabel(els_ecalEnergy_tag, els_ecalEnergy_h);
  hww.els_ecalEnergy_ = *els_ecalEnergy_h.product();

  edm::InputTag els_eOverPOut_tag(electronInputTag.label(),"elseOverPOut");
  edm::Handle<std::vector<float> > els_eOverPOut_h;
  iEvent.getByLabel(els_eOverPOut_tag, els_eOverPOut_h);
  hww.els_eOverPOut_ = *els_eOverPOut_h.product();

  edm::InputTag els_ecalIso_tag(electronInputTag.label(),"elsecalIso");
  edm::Handle<std::vector<float> > els_ecalIso_h;
  iEvent.getByLabel(els_ecalIso_tag, els_ecalIso_h);
  hww.els_ecalIso_ = *els_ecalIso_h.product();

  edm::InputTag els_hcalIso_tag(electronInputTag.label(),"elshcalIso");
  edm::Handle<std::vector<float> > els_hcalIso_h;
  iEvent.getByLabel(els_hcalIso_tag, els_hcalIso_h);
  hww.els_hcalIso_ = *els_hcalIso_h.product();

  edm::InputTag els_nSeed_tag(electronInputTag.label(),"elsnSeed");
  edm::Handle<std::vector<int> > els_nSeed_h;
  iEvent.getByLabel(els_nSeed_tag, els_nSeed_h);
  hww.els_nSeed_ = *els_nSeed_h.product();

  edm::InputTag els_scindex_tag(electronInputTag.label(),"elsscindex");
  edm::Handle<std::vector<int> > els_scindex_h;
  iEvent.getByLabel(els_scindex_tag, els_scindex_h);
  hww.els_scindex_ = *els_scindex_h.product();

  edm::InputTag els_charge_tag(electronInputTag.label(),"elscharge");
  edm::Handle<std::vector<int> > els_charge_h;
  iEvent.getByLabel(els_charge_tag, els_charge_h);
  hww.els_charge_ = *els_charge_h.product();

  edm::InputTag els_gsftrkidx_tag(electronInputTag.label(),"elsgsftrkidx");
  edm::Handle<std::vector<int> > els_gsftrkidx_h;
  iEvent.getByLabel(els_gsftrkidx_tag, els_gsftrkidx_h);
  hww.els_gsftrkidx_ = *els_gsftrkidx_h.product();

  edm::InputTag els_exp_innerlayers_tag(electronInputTag.label(),"elsexpinnerlayers");
  edm::Handle<std::vector<int> > els_exp_innerlayers_h;
  iEvent.getByLabel(els_exp_innerlayers_tag, els_exp_innerlayers_h);
  hww.els_exp_innerlayers_ = *els_exp_innerlayers_h.product();

  edm::InputTag els_type_tag(electronInputTag.label(),"elstype");
  edm::Handle<std::vector<int> > els_type_h;
  iEvent.getByLabel(els_type_tag, els_type_h);
  hww.els_type_ = *els_type_h.product();

  edm::InputTag els_fiduciality_tag(electronInputTag.label(),"elsfiduciality");
  edm::Handle<std::vector<int> > els_fiduciality_h;
  iEvent.getByLabel(els_fiduciality_tag, els_fiduciality_h);
  hww.els_fiduciality_ = *els_fiduciality_h.product();

  edm::InputTag els_conv_dist_tag(electronInputTag.label(),"elsconvdist");
  edm::Handle<std::vector<float> > els_conv_dist_h;
  iEvent.getByLabel(els_conv_dist_tag, els_conv_dist_h);
  hww.els_conv_dist_ = *els_conv_dist_h.product();

  edm::InputTag els_conv_dcot_tag(electronInputTag.label(),"elsconvdcot");
  edm::Handle<std::vector<float> > els_conv_dcot_h;
  iEvent.getByLabel(els_conv_dcot_tag, els_conv_dcot_h);
  hww.els_conv_dcot_ = *els_conv_dcot_h.product();

  edm::InputTag els_conv_old_dist_tag(electronInputTag.label(),"elsconvolddist");
  edm::Handle<std::vector<float> > els_conv_old_dist_h;
  iEvent.getByLabel(els_conv_old_dist_tag, els_conv_old_dist_h);
  hww.els_conv_old_dist_ = *els_conv_old_dist_h.product();

  edm::InputTag els_conv_old_dcot_tag(electronInputTag.label(),"elsconvolddcot");
  edm::Handle<std::vector<float> > els_conv_old_dcot_h;
  iEvent.getByLabel(els_conv_old_dcot_tag, els_conv_old_dcot_h);
  hww.els_conv_old_dcot_ = *els_conv_old_dcot_h.product();

  edm::InputTag els_sccharge_tag(electronInputTag.label(),"elssccharge");
  edm::Handle<std::vector<int> > els_sccharge_h;
  iEvent.getByLabel(els_sccharge_tag, els_sccharge_h);
  hww.els_sccharge_ = *els_sccharge_h.product();

  edm::InputTag els_trk_charge_tag(electronInputTag.label(),"elstrkcharge");
  edm::Handle<std::vector<int> > els_trk_charge_h;
  iEvent.getByLabel(els_trk_charge_tag, els_trk_charge_h);
  hww.els_trk_charge_ = *els_trk_charge_h.product();

  edm::InputTag els_trkshFrac_tag(electronInputTag.label(),"elstrkshFrac");
  edm::Handle<std::vector<float> > els_trkshFrac_h;
  iEvent.getByLabel(els_trkshFrac_tag, els_trkshFrac_h);
  hww.els_trkshFrac_ = *els_trkshFrac_h.product();

  edm::InputTag els_iso04_pf2012_ch_tag(electronInputTag.label(),"elsiso04pf2012ch");
  edm::Handle<std::vector<float> > els_iso04_pf2012_ch_h;
  iEvent.getByLabel(els_iso04_pf2012_ch_tag, els_iso04_pf2012_ch_h);
  hww.els_iso04_pf2012_ch_ = *els_iso04_pf2012_ch_h.product();

  edm::InputTag els_iso04_pf2012_em_tag(electronInputTag.label(),"elsiso04pf2012em");
  edm::Handle<std::vector<float> > els_iso04_pf2012_em_h;
  iEvent.getByLabel(els_iso04_pf2012_em_tag, els_iso04_pf2012_em_h);
  hww.els_iso04_pf2012_em_ = *els_iso04_pf2012_em_h.product();

  edm::InputTag els_iso04_pf2012_nh_tag(electronInputTag.label(),"elsiso04pf2012nh");
  edm::Handle<std::vector<float> > els_iso04_pf2012_nh_h;
  iEvent.getByLabel(els_iso04_pf2012_nh_tag, els_iso04_pf2012_nh_h);
  hww.els_iso04_pf2012_nh_ = *els_iso04_pf2012_nh_h.product();

  edm::InputTag els_iso03_pf2012_ch_tag(electronInputTag.label(),"elsiso03pf2012ch");
  edm::Handle<std::vector<float> > els_iso03_pf2012_ch_h;
  iEvent.getByLabel(els_iso03_pf2012_ch_tag, els_iso03_pf2012_ch_h);
  hww.els_iso03_pf2012_ch_ = *els_iso03_pf2012_ch_h.product();

  edm::InputTag els_iso03_pf2012_em_tag(electronInputTag.label(),"elsiso03pf2012em");
  edm::Handle<std::vector<float> > els_iso03_pf2012_em_h;
  iEvent.getByLabel(els_iso03_pf2012_em_tag, els_iso03_pf2012_em_h);
  hww.els_iso03_pf2012_em_ = *els_iso03_pf2012_em_h.product();

  edm::InputTag els_iso03_pf2012_nh_tag(electronInputTag.label(),"elsiso03pf2012nh");
  edm::Handle<std::vector<float> > els_iso03_pf2012_nh_h;
  iEvent.getByLabel(els_iso03_pf2012_nh_tag, els_iso03_pf2012_nh_h);
  hww.els_iso03_pf2012_nh_ = *els_iso03_pf2012_nh_h.product();

  edm::InputTag els_ecalIso04_tag(electronInputTag.label(),"elsecalIso04");
  edm::Handle<std::vector<float> > els_ecalIso04_h;
  iEvent.getByLabel(els_ecalIso04_tag, els_ecalIso04_h);
  hww.els_ecalIso04_ = *els_ecalIso04_h.product();

  edm::InputTag els_hcalIso04_tag(electronInputTag.label(),"elshcalIso04");
  edm::Handle<std::vector<float> > els_hcalIso04_h;
  iEvent.getByLabel(els_hcalIso04_tag, els_hcalIso04_h);
  hww.els_hcalIso04_ = *els_hcalIso04_h.product();


//get muon variables
  edm::InputTag mus_p4_tag(muonInputTag.label(),"musp4");
  edm::Handle<std::vector<LorentzVector> > mus_p4_h;
  iEvent.getByLabel(mus_p4_tag, mus_p4_h);
  hww.mus_p4_ = *mus_p4_h.product();

  edm::InputTag mus_trk_p4_tag(muonInputTag.label(),"mustrkp4");
  edm::Handle<std::vector<LorentzVector> > mus_trk_p4_h;
  iEvent.getByLabel(mus_trk_p4_tag, mus_trk_p4_h);
  hww.mus_trk_p4_ = *mus_trk_p4_h.product();

  edm::InputTag mus_vertex_p4_tag(muonInputTag.label(),"musvertexp4");
  edm::Handle<std::vector<LorentzVector> > mus_vertex_p4_h;
  iEvent.getByLabel(mus_vertex_p4_tag, mus_vertex_p4_h);
  hww.mus_vertex_p4_ = *mus_vertex_p4_h.product();

  edm::InputTag mus_sta_p4_tag(muonInputTag.label(),"musstap4");
  edm::Handle<std::vector<LorentzVector> > mus_sta_p4_h;
  iEvent.getByLabel(mus_sta_p4_tag, mus_sta_p4_h);
  hww.mus_sta_p4_ = *mus_sta_p4_h.product();

  edm::InputTag mus_gfit_chi2_tag(muonInputTag.label(),"musgfitchi2");
  edm::Handle<std::vector<float> > mus_gfit_chi2_h;
  iEvent.getByLabel(mus_gfit_chi2_tag, mus_gfit_chi2_h);
  hww.mus_gfit_chi2_ = *mus_gfit_chi2_h.product();

  edm::InputTag mus_gfit_ndof_tag(muonInputTag.label(),"musgfitndof");
  edm::Handle<std::vector<float> > mus_gfit_ndof_h;
  iEvent.getByLabel(mus_gfit_ndof_tag, mus_gfit_ndof_h);
  hww.mus_gfit_ndof_ = *mus_gfit_ndof_h.product();

  edm::InputTag mus_ptErr_tag(muonInputTag.label(),"musptErr");
  edm::Handle<std::vector<float> > mus_ptErr_h;
  iEvent.getByLabel(mus_ptErr_tag, mus_ptErr_h);
  hww.mus_ptErr_ = *mus_ptErr_h.product();

  edm::InputTag mus_trkKink_tag(muonInputTag.label(),"mustrkKink");
  edm::Handle<std::vector<float> > mus_trkKink_h;
  iEvent.getByLabel(mus_trkKink_tag, mus_trkKink_h);
  hww.mus_trkKink_ = *mus_trkKink_h.product();

  edm::InputTag mus_d0corr_tag(muonInputTag.label(),"musd0corr");
  edm::Handle<std::vector<float> > mus_d0corr_h;
  iEvent.getByLabel(mus_d0corr_tag, mus_d0corr_h);
  hww.mus_d0corr_ = *mus_d0corr_h.product();

  edm::InputTag mus_d0_tag(muonInputTag.label(),"musd0");
  edm::Handle<std::vector<float> > mus_d0_h;
  iEvent.getByLabel(mus_d0_tag, mus_d0_h);
  hww.mus_d0_ = *mus_d0_h.product();

  edm::InputTag mus_z0corr_tag(muonInputTag.label(),"musz0corr");
  edm::Handle<std::vector<float> > mus_z0corr_h;
  iEvent.getByLabel(mus_z0corr_tag, mus_z0corr_h);
  hww.mus_z0corr_ = *mus_z0corr_h.product();

  edm::InputTag mus_chi2_tag(muonInputTag.label(),"muschi2");
  edm::Handle<std::vector<float> > mus_chi2_h;
  iEvent.getByLabel(mus_chi2_tag, mus_chi2_h);
  hww.mus_chi2_ = *mus_chi2_h.product();

  edm::InputTag mus_ndof_tag(muonInputTag.label(),"musndof");
  edm::Handle<std::vector<float> > mus_ndof_h;
  iEvent.getByLabel(mus_ndof_tag, mus_ndof_h);
  hww.mus_ndof_ = *mus_ndof_h.product();

  edm::InputTag mus_ip3d_tag(muonInputTag.label(),"musip3d");
  edm::Handle<std::vector<float> > mus_ip3d_h;
  iEvent.getByLabel(mus_ip3d_tag, mus_ip3d_h);
  hww.mus_ip3d_ = *mus_ip3d_h.product();

  edm::InputTag mus_ip3derr_tag(muonInputTag.label(),"musip3derr");
  edm::Handle<std::vector<float> > mus_ip3derr_h;
  iEvent.getByLabel(mus_ip3derr_tag, mus_ip3derr_h);
  hww.mus_ip3derr_ = *mus_ip3derr_h.product();

  edm::InputTag mus_segmCompatibility_tag(muonInputTag.label(),"mussegmCompatibility");
  edm::Handle<std::vector<float> > mus_segmCompatibility_h;
  iEvent.getByLabel(mus_segmCompatibility_tag, mus_segmCompatibility_h);
  hww.mus_segmCompatibility_ = *mus_segmCompatibility_h.product();

  edm::InputTag mus_caloCompatibility_tag(muonInputTag.label(),"muscaloCompatibility");
  edm::Handle<std::vector<float> > mus_caloCompatibility_h;
  iEvent.getByLabel(mus_caloCompatibility_tag, mus_caloCompatibility_h);
  hww.mus_caloCompatibility_ = *mus_caloCompatibility_h.product();

  edm::InputTag mus_e_had_tag(muonInputTag.label(),"musehad");
  edm::Handle<std::vector<float> > mus_e_had_h;
  iEvent.getByLabel(mus_e_had_tag, mus_e_had_h);
  hww.mus_e_had_ = *mus_e_had_h.product();

  edm::InputTag mus_e_ho_tag(muonInputTag.label(),"museho");
  edm::Handle<std::vector<float> > mus_e_ho_h;
  iEvent.getByLabel(mus_e_ho_tag, mus_e_ho_h);
  hww.mus_e_ho_ = *mus_e_ho_h.product();

  edm::InputTag mus_e_em_tag(muonInputTag.label(),"museem");
  edm::Handle<std::vector<float> > mus_e_em_h;
  iEvent.getByLabel(mus_e_em_tag, mus_e_em_h);
  hww.mus_e_em_ = *mus_e_em_h.product();

  edm::InputTag mus_e_hadS9_tag(muonInputTag.label(),"musehadS9");
  edm::Handle<std::vector<float> > mus_e_hadS9_h;
  iEvent.getByLabel(mus_e_hadS9_tag, mus_e_hadS9_h);
  hww.mus_e_hadS9_ = *mus_e_hadS9_h.product();

  edm::InputTag mus_e_hoS9_tag(muonInputTag.label(),"musehoS9");
  edm::Handle<std::vector<float> > mus_e_hoS9_h;
  iEvent.getByLabel(mus_e_hoS9_tag, mus_e_hoS9_h);
  hww.mus_e_hoS9_ = *mus_e_hoS9_h.product();

  edm::InputTag mus_e_emS9_tag(muonInputTag.label(),"museemS9");
  edm::Handle<std::vector<float> > mus_e_emS9_h;
  iEvent.getByLabel(mus_e_emS9_tag, mus_e_emS9_h);
  hww.mus_e_emS9_ = *mus_e_emS9_h.product();

  edm::InputTag mus_iso03_sumPt_tag(muonInputTag.label(),"musiso03sumPt");
  edm::Handle<std::vector<float> > mus_iso03_sumPt_h;
  iEvent.getByLabel(mus_iso03_sumPt_tag, mus_iso03_sumPt_h);
  hww.mus_iso03_sumPt_ = *mus_iso03_sumPt_h.product();

  edm::InputTag mus_iso03_emEt_tag(muonInputTag.label(),"musiso03emEt");
  edm::Handle<std::vector<float> > mus_iso03_emEt_h;
  iEvent.getByLabel(mus_iso03_emEt_tag, mus_iso03_emEt_h);
  hww.mus_iso03_emEt_ = *mus_iso03_emEt_h.product();

  edm::InputTag mus_iso03_hadEt_tag(muonInputTag.label(),"musiso03hadEt");
  edm::Handle<std::vector<float> > mus_iso03_hadEt_h;
  iEvent.getByLabel(mus_iso03_hadEt_tag, mus_iso03_hadEt_h);
  hww.mus_iso03_hadEt_ = *mus_iso03_hadEt_h.product();

  edm::InputTag mus_iso05_sumPt_tag(muonInputTag.label(),"musiso05sumPt");
  edm::Handle<std::vector<float> > mus_iso05_sumPt_h;
  iEvent.getByLabel(mus_iso05_sumPt_tag, mus_iso05_sumPt_h);
  hww.mus_iso05_sumPt_ = *mus_iso05_sumPt_h.product();

  edm::InputTag mus_iso05_emEt_tag(muonInputTag.label(),"musiso05emEt");
  edm::Handle<std::vector<float> > mus_iso05_emEt_h;
  iEvent.getByLabel(mus_iso05_emEt_tag, mus_iso05_emEt_h);
  hww.mus_iso05_emEt_ = *mus_iso05_emEt_h.product();

  edm::InputTag mus_iso05_hadEt_tag(muonInputTag.label(),"musiso05hadEt");
  edm::Handle<std::vector<float> > mus_iso05_hadEt_h;
  iEvent.getByLabel(mus_iso05_hadEt_tag, mus_iso05_hadEt_h);
  hww.mus_iso05_hadEt_ = *mus_iso05_hadEt_h.product();

  edm::InputTag mus_gfit_validSTAHits_tag(muonInputTag.label(),"musgfitvalidSTAHits");
  edm::Handle<std::vector<int> > mus_gfit_validSTAHits_h;
  iEvent.getByLabel(mus_gfit_validSTAHits_tag, mus_gfit_validSTAHits_h);
  hww.mus_gfit_validSTAHits_ = *mus_gfit_validSTAHits_h.product();

  edm::InputTag mus_charge_tag(muonInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByLabel(mus_charge_tag, mus_charge_h);
  hww.mus_charge_ = *mus_charge_h.product();

  edm::InputTag mus_validHits_tag(muonInputTag.label(),"musvalidHits");
  edm::Handle<std::vector<int> > mus_validHits_h;
  iEvent.getByLabel(mus_validHits_tag, mus_validHits_h);
  hww.mus_validHits_ = *mus_validHits_h.product();

  edm::InputTag mus_trkidx_tag(muonInputTag.label(),"mustrkidx");
  edm::Handle<std::vector<int> > mus_trkidx_h;
  iEvent.getByLabel(mus_trkidx_tag, mus_trkidx_h);
  hww.mus_trkidx_ = *mus_trkidx_h.product();

  edm::InputTag mus_pid_PFMuon_tag(muonInputTag.label(),"muspidPFMuon");
  edm::Handle<std::vector<int> > mus_pid_PFMuon_h;
  iEvent.getByLabel(mus_pid_PFMuon_tag, mus_pid_PFMuon_h);
  hww.mus_pid_PFMuon_ = *mus_pid_PFMuon_h.product();

  edm::InputTag mus_pid_TMLastStationTight_tag(muonInputTag.label(),"muspidTMLastStationTight");
  edm::Handle<std::vector<int> > mus_pid_TMLastStationTight_h;
  iEvent.getByLabel(mus_pid_TMLastStationTight_tag, mus_pid_TMLastStationTight_h);
  hww.mus_pid_TMLastStationTight_ = *mus_pid_TMLastStationTight_h.product();

  edm::InputTag mus_nmatches_tag(muonInputTag.label(),"musnmatches");
  edm::Handle<std::vector<int> > mus_nmatches_h;
  iEvent.getByLabel(mus_nmatches_tag, mus_nmatches_h);
  hww.mus_nmatches_ = *mus_nmatches_h.product();

  edm::InputTag mus_goodmask_tag(muonInputTag.label(),"musgoodmask");
  edm::Handle<std::vector<int> > mus_goodmask_h;
  iEvent.getByLabel(mus_goodmask_tag, mus_goodmask_h);
  hww.mus_goodmask_ = *mus_goodmask_h.product();

  edm::InputTag mus_type_tag(muonInputTag.label(),"mustype");
  edm::Handle<std::vector<int> > mus_type_h;
  iEvent.getByLabel(mus_type_tag, mus_type_h);
  hww.mus_type_ = *mus_type_h.product();

  edm::InputTag mus_sta_d0_tag(muonInputTag.label(),"musstad0");
  edm::Handle<std::vector<float> > mus_sta_d0_h;
  iEvent.getByLabel(mus_sta_d0_tag, mus_sta_d0_h);
  hww.mus_sta_d0_ = *mus_sta_d0_h.product();

  edm::InputTag mus_sta_z0corr_tag(muonInputTag.label(),"musstaz0corr");
  edm::Handle<std::vector<float> > mus_sta_z0corr_h;
  iEvent.getByLabel(mus_sta_z0corr_tag, mus_sta_z0corr_h);
  hww.mus_sta_z0corr_ = *mus_sta_z0corr_h.product();

  edm::InputTag mus_isoR03_pf_ChargedHadronPt_tag(muonInputTag.label(),"musisoR03pfChargedHadronPt");
  edm::Handle<std::vector<float> > mus_isoR03_pf_ChargedHadronPt_h;
  iEvent.getByLabel(mus_isoR03_pf_ChargedHadronPt_tag, mus_isoR03_pf_ChargedHadronPt_h);
  hww.mus_isoR03_pf_ChargedHadronPt_ = *mus_isoR03_pf_ChargedHadronPt_h.product();

  edm::InputTag mus_isoR03_pf_NeutralHadronEt_tag(muonInputTag.label(),"musisoR03pfNeutralHadronEt");
  edm::Handle<std::vector<float> > mus_isoR03_pf_NeutralHadronEt_h;
  iEvent.getByLabel(mus_isoR03_pf_NeutralHadronEt_tag, mus_isoR03_pf_NeutralHadronEt_h);
  hww.mus_isoR03_pf_NeutralHadronEt_ = *mus_isoR03_pf_NeutralHadronEt_h.product();

  edm::InputTag mus_isoR03_pf_PhotonEt_tag(muonInputTag.label(),"musisoR03pfPhotonEt");
  edm::Handle<std::vector<float> > mus_isoR03_pf_PhotonEt_h;
  iEvent.getByLabel(mus_isoR03_pf_PhotonEt_tag, mus_isoR03_pf_PhotonEt_h);
  hww.mus_isoR03_pf_PhotonEt_ = *mus_isoR03_pf_PhotonEt_h.product();

  edm::InputTag mus_isoR03_pf_PUPt_tag(muonInputTag.label(),"musisoR03pfPUPt");
  edm::Handle<std::vector<float> > mus_isoR03_pf_PUPt_h;
  iEvent.getByLabel(mus_isoR03_pf_PUPt_tag, mus_isoR03_pf_PUPt_h);
  hww.mus_isoR03_pf_PUPt_ = *mus_isoR03_pf_PUPt_h.product();

  edm::InputTag mus_iso_ecalvetoDep_tag(muonInputTag.label(),"musisoecalvetoDep");
  edm::Handle<std::vector<float> > mus_iso_ecalvetoDep_h;
  iEvent.getByLabel(mus_iso_ecalvetoDep_tag, mus_iso_ecalvetoDep_h);
  hww.mus_iso_ecalvetoDep_ = *mus_iso_ecalvetoDep_h.product();

  edm::InputTag mus_iso_hcalvetoDep_tag(muonInputTag.label(),"musisohcalvetoDep");
  edm::Handle<std::vector<float> > mus_iso_hcalvetoDep_h;
  iEvent.getByLabel(mus_iso_hcalvetoDep_tag, mus_iso_hcalvetoDep_h);
  hww.mus_iso_hcalvetoDep_ = *mus_iso_hcalvetoDep_h.product();

  edm::InputTag mus_numberOfMatchedStations_tag(muonInputTag.label(),"musnumberOfMatchedStations");
  edm::Handle<std::vector<int> > mus_numberOfMatchedStations_h;
  iEvent.getByLabel(mus_numberOfMatchedStations_tag, mus_numberOfMatchedStations_h);
  hww.mus_numberOfMatchedStations_ = *mus_numberOfMatchedStations_h.product();


//get dilepton hypothesis variables
  edm::InputTag hyp_type_tag(hypInputTag.label(),"hyptype");
  edm::Handle<std::vector<int> > hyp_type_h;
  iEvent.getByLabel(hyp_type_tag, hyp_type_h);
  hww.hyp_type_ = *hyp_type_h.product();

  edm::InputTag hyp_ll_charge_tag(hypInputTag.label(),"hypllcharge");
  edm::Handle<std::vector<int> > hyp_ll_charge_h;
  iEvent.getByLabel(hyp_ll_charge_tag, hyp_ll_charge_h);
  hww.hyp_ll_charge_ = *hyp_ll_charge_h.product();

  edm::InputTag hyp_lt_charge_tag(hypInputTag.label(),"hypltcharge");
  edm::Handle<std::vector<int> > hyp_lt_charge_h;
  iEvent.getByLabel(hyp_lt_charge_tag, hyp_lt_charge_h);
  hww.hyp_lt_charge_ = *hyp_lt_charge_h.product();


//get conversion variables
  edm::InputTag convs_ndof_tag(conversionInputTag.label(),"convsndof");
  edm::Handle<std::vector<float> > convs_ndof_h;
  iEvent.getByLabel(convs_ndof_tag, convs_ndof_h);
  hww.convs_ndof_ = *convs_ndof_h.product();

  edm::InputTag convs_chi2_tag(conversionInputTag.label(),"convschi2");
  edm::Handle<std::vector<float> > convs_chi2_h;
  iEvent.getByLabel(convs_chi2_tag, convs_chi2_h);
  hww.convs_chi2_ = *convs_chi2_h.product();

  edm::InputTag convs_dl_tag(conversionInputTag.label(),"convsdl");
  edm::Handle<std::vector<float> > convs_dl_h;
  iEvent.getByLabel(convs_dl_tag, convs_dl_h);
  hww.convs_dl_ = *convs_dl_h.product();

  edm::InputTag convs_isConverted_tag(conversionInputTag.label(),"convsisConverted");
  edm::Handle<std::vector<int> > convs_isConverted_h;
  iEvent.getByLabel(convs_isConverted_tag, convs_isConverted_h);
  hww.convs_isConverted_ = *convs_isConverted_h.product();

  edm::InputTag convs_tkalgo_tag(conversionInputTag.label(),"convstkalgo");
  edm::Handle<std::vector<vector<int>> > convs_tkalgo_h;
  iEvent.getByLabel(convs_tkalgo_tag, convs_tkalgo_h);
  hww.convs_tkalgo_ = *convs_tkalgo_h.product();

  edm::InputTag convs_tkidx_tag(conversionInputTag.label(),"convstkidx");
  edm::Handle<std::vector<vector<int>> > convs_tkidx_h;
  iEvent.getByLabel(convs_tkidx_tag, convs_tkidx_h);
  hww.convs_tkidx_ = *convs_tkidx_h.product();

  edm::InputTag convs_nHitsBeforeVtx_tag(conversionInputTag.label(),"convsnHitsBeforeVtx");
  edm::Handle<std::vector<vector<int>> > convs_nHitsBeforeVtx_h;
  iEvent.getByLabel(convs_nHitsBeforeVtx_tag, convs_nHitsBeforeVtx_h);
  hww.convs_nHitsBeforeVtx_ = *convs_nHitsBeforeVtx_h.product();

  edm::InputTag convs_quality_tag(conversionInputTag.label(),"convsquality");
  edm::Handle<std::vector<int> > convs_quality_h;
  iEvent.getByLabel(convs_quality_tag, convs_quality_h);
  hww.convs_quality_ = *convs_quality_h.product();

//get SC variables
  edm::InputTag scs_sigmaIEtaIPhi_tag(scInputTag.label(),"scssigmaIEtaIPhi");
  edm::Handle<std::vector<float> > scs_sigmaIEtaIPhi_h;
  iEvent.getByLabel(scs_sigmaIEtaIPhi_tag, scs_sigmaIEtaIPhi_h);
  hww.scs_sigmaIEtaIPhi_ = *scs_sigmaIEtaIPhi_h.product();

  edm::InputTag scs_e1x3_tag(scInputTag.label(),"scse1x3");
  edm::Handle<std::vector<float> > scs_e1x3_h;
  iEvent.getByLabel(scs_e1x3_tag, scs_e1x3_h);
  hww.scs_e1x3_ = *scs_e1x3_h.product();

  edm::InputTag scs_e3x1_tag(scInputTag.label(),"scse3x1");
  edm::Handle<std::vector<float> > scs_e3x1_h;
  iEvent.getByLabel(scs_e3x1_tag, scs_e3x1_h);
  hww.scs_e3x1_ = *scs_e3x1_h.product();

  edm::InputTag scs_eMax_tag(scInputTag.label(),"scseMax");
  edm::Handle<std::vector<float> > scs_eMax_h;
  iEvent.getByLabel(scs_eMax_tag, scs_eMax_h);
  hww.scs_eMax_ = *scs_eMax_h.product();

  edm::InputTag scs_pos_p4_tag(scInputTag.label(),"scsposp4");
  edm::Handle<std::vector<LorentzVector> > scs_pos_p4_h;
  iEvent.getByLabel(scs_pos_p4_tag, scs_pos_p4_h);
  hww.scs_pos_p4_ = *scs_pos_p4_h.product();

//get pf candidate variables
  edm::InputTag pfcands_p4_tag(pfCandidateInputTag.label(),"pfcandsp4");
  edm::Handle<std::vector<LorentzVector> > pfcands_p4_h;
  iEvent.getByLabel(pfcands_p4_tag, pfcands_p4_h);
  hww.pfcands_p4_ = *pfcands_p4_h.product();

  edm::InputTag pfcands_trkidx_tag(pfCandidateInputTag.label(),"pfcandstrkidx");
  edm::Handle<std::vector<int> > pfcands_trkidx_h;
  iEvent.getByLabel(pfcands_trkidx_tag, pfcands_trkidx_h);
  hww.pfcands_trkidx_ = *pfcands_trkidx_h.product();

  edm::InputTag pfcands_particleId_tag(pfCandidateInputTag.label(),"pfcandsparticleId");
  edm::Handle<std::vector<int> > pfcands_particleId_h;
  iEvent.getByLabel(pfcands_particleId_tag, pfcands_particleId_h);
  hww.pfcands_particleId_ = *pfcands_particleId_h.product();

  edm::InputTag pfcands_pfelsidx_tag(pfCandidateInputTag.label(),"pfcandspfelsidx");
  edm::Handle<std::vector<int> > pfcands_pfelsidx_h;
  iEvent.getByLabel(pfcands_pfelsidx_tag, pfcands_pfelsidx_h);
  hww.pfcands_pfelsidx_ = *pfcands_pfelsidx_h.product();

  edm::InputTag pfcands_vtxidx_tag(pfCandidateInputTag.label(),"pfcandsvtxidx");
  edm::Handle<std::vector<int> > pfcands_vtxidx_h;
  iEvent.getByLabel(pfcands_vtxidx_tag, pfcands_vtxidx_h);
  hww.pfcands_vtxidx_ = *pfcands_vtxidx_h.product();

  edm::InputTag pfcands_charge_tag(pfCandidateInputTag.label(),"pfcandscharge");
  edm::Handle<std::vector<int> > pfcands_charge_h;
  iEvent.getByLabel(pfcands_charge_tag, pfcands_charge_h);
  hww.pfcands_charge_ = *pfcands_charge_h.product();

//get gsf track variables
  edm::InputTag gsftrks_p4_tag(gsfTrackInputTag.label(),"gsftrksp4");
  edm::Handle<std::vector<LorentzVector> > gsftrks_p4_h;
  iEvent.getByLabel(gsftrks_p4_tag, gsftrks_p4_h);
  hww.gsftrks_p4_ = *gsftrks_p4_h.product();

  edm::InputTag gsftrks_vertex_p4_tag(gsfTrackInputTag.label(),"gsftrksvertexp4");
  edm::Handle<std::vector<LorentzVector> > gsftrks_vertex_p4_h;
  iEvent.getByLabel(gsftrks_vertex_p4_tag, gsftrks_vertex_p4_h);
  hww.gsftrks_vertex_p4_ = *gsftrks_vertex_p4_h.product();

  edm::InputTag gsftrks_d0_tag(gsfTrackInputTag.label(),"gsftrksd0");
  edm::Handle<std::vector<float> > gsftrks_d0_h;
  iEvent.getByLabel(gsftrks_d0_tag, gsftrks_d0_h);
  hww.gsftrks_d0_ = *gsftrks_d0_h.product();

  edm::InputTag gsftrks_d0Err_tag(gsfTrackInputTag.label(),"gsftrksd0Err");
  edm::Handle<std::vector<float> > gsftrks_d0Err_h;
  iEvent.getByLabel(gsftrks_d0Err_tag, gsftrks_d0Err_h);
  hww.gsftrks_d0Err_ = *gsftrks_d0Err_h.product();

  edm::InputTag gsftrks_phiErr_tag(gsfTrackInputTag.label(),"gsftrksphiErr");
  edm::Handle<std::vector<float> > gsftrks_phiErr_h;
  iEvent.getByLabel(gsftrks_phiErr_tag, gsftrks_phiErr_h);
  hww.gsftrks_phiErr_ = *gsftrks_phiErr_h.product();

  edm::InputTag gsftrks_d0phiCov_tag(gsfTrackInputTag.label(),"gsftrksd0phiCov");
  edm::Handle<std::vector<float> > gsftrks_d0phiCov_h;
  iEvent.getByLabel(gsftrks_d0phiCov_tag, gsftrks_d0phiCov_h);
  hww.gsftrks_d0phiCov_ = *gsftrks_d0phiCov_h.product();

  edm::InputTag gsftrks_z0_tag(gsfTrackInputTag.label(),"gsftrksz0");
  edm::Handle<std::vector<float> > gsftrks_z0_h;
  iEvent.getByLabel(gsftrks_z0_tag, gsftrks_z0_h);
  hww.gsftrks_z0_ = *gsftrks_z0_h.product();

  edm::InputTag gsftrks_z0Err_tag(gsfTrackInputTag.label(),"gsftrksz0Err");
  edm::Handle<std::vector<float> > gsftrks_z0Err_h;
  iEvent.getByLabel(gsftrks_z0Err_tag, gsftrks_z0Err_h);
  hww.gsftrks_z0Err_ = *gsftrks_z0Err_h.product();

  edm::InputTag gsftrks_etaErr_tag(gsfTrackInputTag.label(),"gsftrksetaErr");
  edm::Handle<std::vector<float> > gsftrks_etaErr_h;
  iEvent.getByLabel(gsftrks_etaErr_tag, gsftrks_etaErr_h);
  hww.gsftrks_etaErr_ = *gsftrks_etaErr_h.product();

//get other variables
  edm::InputTag evt_kt6pf_foregiso_rho_tag(kt6PFInputTag.label(),"evtkt6pfforegisorho");
  edm::Handle<float> evt_kt6pf_foregiso_rho_h;
  iEvent.getByLabel(evt_kt6pf_foregiso_rho_tag, evt_kt6pf_foregiso_rho_h);
  hww.evt_kt6pf_foregiso_rho_ = *evt_kt6pf_foregiso_rho_h.product();

  edm::InputTag evt_rho_tag(fastJetInputTag.label(),"evtrho");
  edm::Handle<float> evt_rho_h;
  iEvent.getByLabel(evt_rho_tag, evt_rho_h);
  hww.evt_rho_ = *evt_rho_h.product();

  edm::InputTag evt_ww_rho_tag(wwRhoDefaultInputTag.label(),"evtwwrho");
  edm::Handle<float> evt_ww_rho_h;
  iEvent.getByLabel(evt_ww_rho_tag, evt_ww_rho_h);
  hww.evt_ww_rho_ = *evt_ww_rho_h.product();

  edm::InputTag evt_ww_rho_vor_tag(wwRhoVoronoiInputTag.label(),"evtwwrhovor");
  edm::Handle<float > evt_ww_rho_vor_h;
  iEvent.getByLabel(evt_ww_rho_vor_tag, evt_ww_rho_vor_h);
  hww.evt_ww_rho_vor_ = *evt_ww_rho_vor_h.product();

  edm::InputTag els_closestMuon_tag(elToMuAssInputTag.label(),"elsclosestMuon");
  edm::Handle<std::vector<int> > els_closestMuon_h;
  iEvent.getByLabel(els_closestMuon_tag, els_closestMuon_h);
  hww.els_closestMuon_ = *els_closestMuon_h.product();

  edm::InputTag pfels_elsidx_tag(pfElToElAssInputTag.label(),"pfelselsidx");
  edm::Handle<std::vector<int> > pfels_elsidx_h;
  iEvent.getByLabel(pfels_elsidx_tag, pfels_elsidx_h);
  hww.pfels_elsidx_ = *pfels_elsidx_h.product();

  edm::InputTag evt_pfmet_tag(pfmetInputTag.label(),"evtpfmet");
  edm::Handle<float> evt_pfmet_h;
  iEvent.getByLabel(evt_pfmet_tag, evt_pfmet_h);
  hww.evt_pfmet_ = *evt_pfmet_h.product();

  edm::InputTag evt_pfmetPhi_tag(pfmetInputTag.label(),"evtpfmetPhi");
  edm::Handle<float> evt_pfmetPhi_h;
  iEvent.getByLabel(evt_pfmetPhi_tag, evt_pfmetPhi_h);
  hww.evt_pfmetPhi_ = *evt_pfmetPhi_h.product();

  edm::InputTag trk_met_tag(trkMetInputTag.label(),"trkmet");
  edm::Handle<vector<float> > trk_met_h;
  iEvent.getByLabel(trk_met_tag, trk_met_h);
  hww.trk_met_ = *trk_met_h.product();

  edm::InputTag trk_metPhi_tag(trkMetInputTag.label(),"trkmetPhi");
  edm::Handle<vector<float> > trk_metPhi_h;
  iEvent.getByLabel(trk_metPhi_tag, trk_metPhi_h);
  hww.trk_metPhi_ = *trk_metPhi_h.product();

  edm::InputTag pfjets_p4_tag(pfJetInputTag.label(),"pfjetsp4");
  edm::Handle<vector<LorentzVector> > pfjets_p4_h;
  iEvent.getByLabel(pfjets_p4_tag, pfjets_p4_h);
  hww.pfjets_p4_ = *pfjets_p4_h.product();

  edm::InputTag pfjets_area_tag(pfJetInputTag.label(),"pfjetsarea");
  edm::Handle<vector<float> > pfjets_area_h;
  iEvent.getByLabel(pfjets_area_tag, pfjets_area_h);
  hww.pfjets_area_ = *pfjets_area_h.product();

  edm::InputTag pfjets_JEC_tag(mvaJetIdInputTag.label(),"pfjetsJEC");
  edm::Handle<vector<float> > pfjets_JEC_h;
  iEvent.getByLabel(pfjets_JEC_tag, pfjets_JEC_h);
  hww.pfjets_JEC_ = *pfjets_JEC_h.product();

  edm::InputTag pfjets_mvavalue_tag(mvaJetIdInputTag.label(),"pfjetsmvavalue");
  edm::Handle<vector<float> > pfjets_mvavalue_h;
  iEvent.getByLabel(pfjets_mvavalue_tag, pfjets_mvavalue_h);
  hww.pfjets_mvavalue_ = *pfjets_mvavalue_h.product();

  edm::InputTag pfjets_corr_p4_tag(mvaJetIdInputTag.label(),"pfjetscorrp4");
  edm::Handle<vector<LorentzVector> > pfjets_corr_p4_h;
  iEvent.getByLabel(pfjets_corr_p4_tag, pfjets_corr_p4_h);
  hww.pfjets_corr_p4_ = *pfjets_corr_p4_h.product();

  edm::InputTag pfjets_trackCountingHighEffBJetTag_tag(bTagPFJetInputTag.label(),"pfjetstrackCountingHighEffBJetTag");
  edm::Handle<vector<float> > pfjets_trackCountingHighEffBJetTag_h;
  iEvent.getByLabel(pfjets_trackCountingHighEffBJetTag_tag, pfjets_trackCountingHighEffBJetTag_h);
  hww.pfjets_trackCountingHighEffBJetTag_ = *pfjets_trackCountingHighEffBJetTag_h.product();

}
