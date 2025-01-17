#include "Math/VectorUtil.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HWWValidation/HWWBase/interface/electronSelections.h"
#include "HWWValidation/HWWBase/interface/MITConversionUtilities.h"
#include "HWWValidation/HWWBase/interface/muonSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"
#include "HWWValidation/HWWBase/interface/analysisSelections.h"

using namespace std;
using namespace wp2012;

namespace HWWFunctions {

  bool pass_electronSelectionCompareMask( const cuts_t cuts_passed, const cuts_t selectionType ) {
      if ((cuts_passed & selectionType) == selectionType) return true;
      return false;
  }

  bool pass_electronSelection( const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, bool useGsfTrack) {
    checkElectronSelections();
    cuts_t cuts_passed = electronSelection(index, applyAlignmentCorrection, removedEtaCutInEndcap, useGsfTrack);
    if ( (cuts_passed & selectionType) == selectionType ) return true;
    return false;
  }

  cuts_t electronSelection( const unsigned int index, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, bool useGsfTrack) {

      // keep track of which cuts passed
      cuts_t cuts_passed = 0;

      /////////////// 
      // Isolation //
      ///////////////


      //pf iso
      float pfiso = electronIsoValuePF(index,0);
      if (fabs(HWWVal::els_p4().at(index).eta()) < 1.479){
        if (pfiso<0.15) cuts_passed |= (1ll<<ELEISO_SMURFV4);
        if (pfiso<0.13) cuts_passed |= (1ll<<ELEISO_SMURFV5);
      } else if (pfiso<0.09) {
        cuts_passed |= (1ll<<ELEISO_SMURFV4);
        cuts_passed |= (1ll<<ELEISO_SMURFV5);
      }

      ////////
      // d0 //
      ////////
      if (fabs(electron_d0PV_smurfV3(index)) < 0.02 && fabs(electron_dzPV_smurfV3(index)) < 0.2 ) cuts_passed |= (1ll<<ELEIP_PV_SMURFV3);
      if (fabs(electron_dzPV_smurfV3(index)) < 0.1 ) cuts_passed |= (1ll<<ELEIP_PV_DZ_1MM);
      if (fabs(electron_d0PV_smurfV3(index)) < 0.04 && fabs(electron_dzPV_smurfV3(index)) < 1.0 ) cuts_passed |= (1ll<<ELEIP_PV_OSV2);
      if (fabs(electron_d0PV_smurfV3(index)) < 0.20 && fabs(electron_dzPV_smurfV3(index)) < 1.0 ) cuts_passed |= (1ll<<ELEIP_PV_OSV2_FO);

      ////////////////////
      // Identification //
      ////////////////////

      // SMURF ID
      if (electronId_smurf_v1(index)) cuts_passed |= (1ll<<ELEID_SMURFV1_EXTRA);
      if (electronId_smurf_v2(index)) cuts_passed |= (1ll<<ELEID_SMURFV2_EXTRA);
      if (electronId_smurf_v3(index)) cuts_passed |= (1ll<<ELEID_SMURFV3_EXTRA);


      electronIdComponent_t answer_med_2012 = electronId_WP2012(index, MEDIUM);
      if ((answer_med_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) cuts_passed |= (1ll<<ELEID_WP2012_MEDIUM_NOISO);
      if ((answer_med_2012 & PassWP2012CutsNoIsoNoIP) == PassWP2012CutsNoIsoNoIP) cuts_passed |= (1ll<<ELEID_WP2012_MEDIUM_NOISO_NOIP);

      electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, LOOSE);
      if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) cuts_passed |= (1ll<<ELEID_WP2012_LOOSE_NOISO);

      // VBTF ID
      electronIdComponent_t answer_vbtf = 0;
      // VBTF95 (optimised in 35X)
      answer_vbtf = electronId_VBTF(index, VBTF_35X_95, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_95);
      // VBTF90 (optimised in 35X)
      answer_vbtf = electronId_VBTF(index, VBTF_35X_90, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_90);
      // VBTF80 (optimised in 35X)
      answer_vbtf = electronId_VBTF(index, VBTF_35X_80, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_80);
      // VBTF85 no H/E in endcap
      answer_vbtf = electronId_VBTF(index, VBTF_85_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_85_NOHOEEND);
      // VBTF85
      answer_vbtf = electronId_VBTF(index, VBTF_85 , applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_85);
      // VBTF80 no H/E in endcap
      answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_80_NOHOEEND);
      // VBTF70 no H/E in endcap
      answer_vbtf = electronId_VBTF(index, VBTF_70_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_70_NOHOEEND);
      // VBTF90 with H/E and dPhiIn tuned to match HLT
      answer_vbtf = electronId_VBTF(index, VBTF_90_HLT, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_90_HLT);
      // VBTF90 with H/E and dPhiIn tuned to match HLT (CaloIdT+TrkIdVL)
      answer_vbtf = electronId_VBTF(index, VBTF_90_HLT_CALOIDT_TRKIDVL, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL);
      // VBTF95 no H/E in endcap
      answer_vbtf = electronId_VBTF(index, VBTF_95_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
      if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_95_NOHOEEND);


      //////////////////////////
      // Conversion Rejection //
      //////////////////////////
      if (!isFromConversionPartnerTrack(index)) cuts_passed |= (1ll<<ELENOTCONV_DISTDCOT002);
      if (HWWVal::els_exp_innerlayers().at(index) == 0) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN_0MHITS);
      if(!isFromConversionMIT(index)) cuts_passed |= (1ll<<ELENOTCONV_MIT);

      /////////////////
      // Fiduciality //
      /////////////////
      if ((HWWVal::els_type().at(index) & (1ll<<ISECALDRIVEN))) cuts_passed |= (1ll<<ELESEED_ECAL);
      if (fabs(HWWVal::els_p4().at(index).eta()) < 2.5) cuts_passed |= (1ll<<ELEETA_250);
      if (fabs(HWWVal::els_p4().at(index).eta()) < 2.4) cuts_passed |= (1ll<<ELEETA_240);

      ////////
      // Pt //
      ////////
      if( HWWVal::els_p4().at(index).pt() > 10.0 ) cuts_passed |= (1ll<<ELEPT_010);

      // Veto electron in transition region
      if( fabs(HWWVal::els_etaSC().at(index)) < 1.4442 || fabs(HWWVal::els_etaSC().at(index)) > 1.566 )  cuts_passed |= (1ll<<ELE_NOT_TRANSITION);

      /////////////////
      // Charge Flip //
      /////////////////
      if (!isChargeFlip3agree(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP3AGREE);

      // return which selections passed
      return cuts_passed;

  }

  ////////////////////
  // Identification //
  ////////////////////

  bool electronId_smurf_v1(const unsigned int index)
  {

    if (HWWVal::els_p4().at(index).pt() > 20.0) return true;

    if (HWWVal::els_fbrem().at(index) > 0.15) return true;
    
    if (fabs(HWWVal::els_etaSC().at(index)) < 1.) {
      if (HWWVal::els_eOverPIn().at(index) > 0.95 && fabs(HWWVal::els_dPhiIn().at(index)*HWWVal::els_charge().at(index)) < 0.006) return true; 
    }

    return false;
  }

  bool electronId_smurf_v2(const unsigned int index)	 
   {	 
     
     if (HWWVal::els_p4().at(index).pt() > 20.0) return true;	 
     
     if (HWWVal::els_fbrem().at(index) > 0.15) return true;	 
     
     if (fabs(HWWVal::els_etaSC().at(index)) < 1.) {	 
       if (HWWVal::els_eOverPIn().at(index) > 0.95) return true;	 
     }	 
     
     return false;	 
   }	 
     
   bool electronId_smurf_v3(const unsigned int index)	 
   {	 
     if (fabs(HWWVal::els_etaSC().at(index)) > 1.479 && HWWVal::els_hOverE().at(index) > 0.1) return false;
     
     if (HWWVal::els_p4().at(index).pt() > 20.0) return true;	 
     
     electronIdComponent_t answer_vbtf = 0;	 
     answer_vbtf = electronId_VBTF(index, VBTF_70_NOHOEEND, false, false);	 
     if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {	 
     
       if (HWWVal::els_fbrem().at(index) > 0.15) return true;	 
     
       if (fabs(HWWVal::els_etaSC().at(index)) < 1.) {	 
         if (HWWVal::els_eOverPIn().at(index) > 0.95) return true;	 
       }	 
     
     }	 
     
     return false;	 
   }


  // WP2012
  electronIdComponent_t electronId_WP2012(const unsigned int index, const wp2012_tightness tightness)
  {

      // set return value
      unsigned int mask = 0;

      // cut values
      std::vector<double> dEtaInThresholds;
      std::vector<double> dPhiInThresholds;
      std::vector<double> sigmaIEtaIEtaThresholds;
      std::vector<double> hoeThresholds;
      std::vector<double> ooemoopThresholds;
      std::vector<double> d0VtxThresholds;
      std::vector<double> dzVtxThresholds;
      std::vector<bool> vtxFitThresholds;
      std::vector<int> mHitsThresholds;
      std::vector<double> isoHiThresholds;
      std::vector<double> isoLoThresholds;

      // set cut values
      eidGetWP2012(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, sigmaIEtaIEtaThresholds, 
                      ooemoopThresholds, d0VtxThresholds, dzVtxThresholds, vtxFitThresholds, mHitsThresholds, 
                      isoHiThresholds, isoLoThresholds);

      // useful kinematic variables
      unsigned int det = ((HWWVal::els_fiduciality().at(index) & (1<<ISEB)) == (1<<ISEB)) ? 0 : 1;
      float etaAbs = fabs(HWWVal::els_etaSC().at(index));
      float pt     = HWWVal::els_p4().at(index).pt();

      // get effective area
      float AEff = 0.;
      if (etaAbs <= 1.0) AEff = 0.10;
      else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
      else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
      else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
      else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
      else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
      else if (etaAbs > 2.4) AEff = 0.13;

      // pf iso
      // calculate from the ntuple for now...
      float pfiso_ch = HWWVal::els_iso03_pf2012_ch().at(index);
      float pfiso_em = HWWVal::els_iso03_pf2012_em().at(index);
      float pfiso_nh = HWWVal::els_iso03_pf2012_nh().at(index);

      // rho
      float rhoPrime = std::max(HWWVal::evt_kt6pf_foregiso_rho(), float(0.0));
      float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
      float pfiso = (pfiso_ch + pfiso_n) / pt;

      // |1/E - 1/p|
      float ooemoop = fabs( (1.0/HWWVal::els_ecalEnergy().at(index)) - (HWWVal::els_eOverPIn().at(index)/HWWVal::els_ecalEnergy().at(index)) );

      // MIT conversion vtx fit
      bool vtxFitConversion = isMITConversion(index, 0,   1e-6,   2.0,   true,  false);

      // d0
      float d0vtx = electron_d0PV_smurfV3(index);
      float dzvtx = electron_dzPV_smurfV3(index);
   
      // test cuts
      if (fabs(HWWVal::els_dEtaIn().at(index)) < dEtaInThresholds[det])             mask |= wp2012::DETAIN;
      if (fabs(HWWVal::els_dPhiIn().at(index)) < dPhiInThresholds[det])             mask |= wp2012::DPHIIN;
      if (HWWVal::els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[det])     mask |= wp2012::SIGMAIETAIETA;
      if (HWWVal::els_hOverE().at(index) < hoeThresholds[det])                      mask |= wp2012::HOE;
      if (ooemoop < ooemoopThresholds[det])                                   mask |= wp2012::OOEMOOP;
      if (fabs(d0vtx) < d0VtxThresholds[det])                                 mask |= wp2012::D0VTX;
      if (fabs(dzvtx) < dzVtxThresholds[det])                                 mask |= wp2012::DZVTX;
      if (!vtxFitThresholds[det] || !vtxFitConversion)                        mask |= wp2012::VTXFIT;
      if (HWWVal::els_exp_innerlayers().at(index) <= mHitsThresholds[det])          mask |= wp2012::MHITS;
      if (pt >= 20.0 && pfiso < isoHiThresholds[det])                         mask |= wp2012::ISO;
      if (pt < 20.0 && pfiso < isoLoThresholds[det])                          mask |= wp2012::ISO;

      // return the mask
      return mask;

  }

  float fastJetEffArea04_v1(const float eta)
  {
    // use absolute eta
      const float etaAbs = fabs(eta);

      // get effective area
      if      (etaAbs <= 1.0                  ) {return 0.19;}
      else if (etaAbs > 1.0 && etaAbs <= 1.479) {return 0.25;}
      else if (etaAbs > 1.479 && etaAbs <= 2.0) {return 0.12;}
      else if (etaAbs > 2.0 && etaAbs <= 2.2  ) {return 0.21;}
      else if (etaAbs > 2.2 && etaAbs <= 2.3  ) {return 0.27;}
      else if (etaAbs > 2.3 && etaAbs <= 2.4  ) {return 0.44;}
      else if (etaAbs > 2.4                   ) {return 0.52;}
      return -9999.0f;
  }

  // VBTF stuff
  electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness, bool applyAlignementCorrection, bool removedEtaCutInEndcap)
  {

      unsigned int answer = 0;

      std::vector<double> relisoThresholds;
      std::vector<double> dEtaInThresholds;
      std::vector<double> dPhiInThresholds;
      std::vector<double> hoeThresholds;
      std::vector<double> sigmaIEtaIEtaThresholds;

      eidGetVBTF(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, 
              sigmaIEtaIEtaThresholds, relisoThresholds);

      //
      // get corrected dEtaIn and dPhiIn
      //

      float dEtaIn = HWWVal::els_dEtaIn().at(index);
      float dPhiIn = HWWVal::els_dPhiIn().at(index);
      if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

      // barrel
      if (fabs(HWWVal::els_etaSC().at(index)) < 1.479) {

          if (electronIsolation_rel(index, true) < relisoThresholds[0])
              answer |= (1<<ELEID_ISO);

          if (fabs(dEtaIn) < dEtaInThresholds[0] &&
                  fabs(dPhiIn) < dPhiInThresholds[0] &&
                  HWWVal::els_hOverE().at(index) < hoeThresholds[0] &&
                  HWWVal::els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[0])
              answer |= (1<<ELEID_ID);
      }

      // endcap
      if (fabs(HWWVal::els_etaSC().at(index)) > 1.479) {
          if (electronIsolation_rel(index, true) < relisoThresholds[1])
              answer |= (1<<ELEID_ISO);
          bool passdEtaCut = fabs(dEtaIn) < dEtaInThresholds[1];
          if(removedEtaCutInEndcap) passdEtaCut = true;
          if ( passdEtaCut &&
                  fabs(dPhiIn) < dPhiInThresholds[1] &&
                  HWWVal::els_hOverE().at(index) < hoeThresholds[1] &&
                  HWWVal::els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[1])
              answer |= (1<<ELEID_ID);
      }

      return answer;

  }

  bool passLikelihoodId_v2(unsigned int index, float lhValue, int workingPoint)
  {

      float etaSC = HWWVal::els_etaSC().at(index);
      float pt = HWWVal::els_p4().at(index).Pt();
      unsigned int nbrem = HWWVal::els_nSeed().at(index);

      if ( workingPoint == 0 ) {

          if (pt > 20.0) {
              if ( ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > 3.5 ) ||
                   ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > 4.0 ) ||
                   ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > 4.0 ) ||
                   ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > 4.0 )) return true;
          } else if (pt > 10.0 && pt <= 20.0) {
              if ( ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > 4.0 ) ||
                   ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > 4.5 ) ||
                   ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > 4.0 ) ||
                   ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > 4.0 )) return true;
              }

          } else {
          edm::LogError("InvalidInput") << "Error! Likelihood WP not supported: " 
                                        << workingPoint << ". Please choose 0 for Emanuele 8th September";
      }

      return false;
  }

  ///////////////
  // Isolation //
  ///////////////

  // relative truncated
  float electronIsolation_rel( const unsigned int index, bool use_calo_iso ) {
      float sum = HWWVal::els_tkIso().at(index);
      if (use_calo_iso) {
          if (fabs(HWWVal::els_etaSC().at(index)) > 1.479) sum += HWWVal::els_ecalIso().at(index);
          if (fabs(HWWVal::els_etaSC().at(index)) <= 1.479) sum += max(0., (HWWVal::els_ecalIso().at(index) -1.));
          sum += HWWVal::els_hcalIso().at(index);
      }
      double pt = HWWVal::els_p4().at(index).pt();
      return sum/max(pt, 20.);
  }

  float electronIsoValuePF( const unsigned int iel, unsigned int ivtx, float coner, float minptn, float dzcut, float footprintdr, float gammastripveto, float elestripveto, int filterId ) {

    int elgsftkid = HWWVal::els_gsftrkidx().at(iel);
    int eltkid = HWWVal::els_trkidx().at(iel);
    float eldz = elgsftkid>=0 ? gsftrks_dz_pv( elgsftkid,ivtx ).first : trks_dz_pv(eltkid,ivtx).first;
    float eleta = HWWVal::els_p4().at(iel).eta();

    float pfciso = 0.;
    float pfniso = 0.;
    float pffootprint = 0.;
    float pfjurveto = 0.;
    float pfjurvetoq = 0.;
    for (unsigned int ipf=0; ipf<HWWVal::pfcands_p4().size(); ++ipf){

      float dR = ROOT::Math::VectorUtil::DeltaR( HWWVal::pfcands_p4().at(ipf), HWWVal::els_p4().at(iel) );

      if (dR>coner) continue;

      float pfpt = HWWVal::pfcands_p4().at(ipf).pt();    
      float pfeta = HWWVal::pfcands_p4().at(ipf).eta();    
      float deta = fabs(pfeta - eleta);
      int pfid = abs(HWWVal::pfcands_particleId().at(ipf));

      if (filterId!=0 && filterId!=pfid) continue;

      if (HWWVal::pfcands_charge().at(ipf)==0) {
        //neutrals
        if (pfpt>minptn) {
    pfniso+=pfpt;
    if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
    if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
        }
      } else {
        //charged  
        //avoid double counting of electron itself
        //if either the gsf or the ctf track are shared with the candidate, skip it
        int pftkid = HWWVal::pfcands_trkidx().at(ipf);
        if (eltkid>=0 && pftkid>=0 && eltkid==pftkid) continue;
        if (pfid==11 && HWWVal::pfcands_pfelsidx().at(ipf)>=0 && HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))>=0) {
    int pfgsfid = HWWVal::els_gsftrkidx().at(HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))); 
    if (elgsftkid>=0 && pfgsfid>=0 && elgsftkid==pfgsfid) continue;
        }
        //check electrons with gsf track
        if (pfid==11 && HWWVal::pfcands_pfelsidx().at(ipf)>=0 && HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))>=0) {
    int gsfid = HWWVal::els_gsftrkidx().at(HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))); 
    if (gsfid>=0) { 
      if(fabs(gsftrks_dz_pv( gsfid,ivtx ).first - eldz )<dzcut) {//dz cut
        pfciso+=pfpt;
        if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
      }   
      continue;//and avoid double counting
    }
        }      
        //then check anything that has a ctf track
        if (pftkid>=0) {//charged (with a ctf track)
    if(fabs( trks_dz_pv(HWWVal::pfcands_trkidx().at(ipf),ivtx).first - eldz )<dzcut) {//dz cut
      pfciso+=pfpt;
      if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
    }
        }
      } 
    }


    return (pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq)/HWWVal::els_p4().at(iel).pt();

  }

  //////////////////////////
  // Conversion Rejection //
  //////////////////////////

  bool isFromConversionPartnerTrack(const unsigned int index) {
      if( fabs(HWWVal::els_conv_dist().at(index)) < 0.02 && fabs(HWWVal::els_conv_dcot().at(index)) < 0.02 ) return true;
      return false;
  }

  bool isFromConversionMIT(const unsigned int index){
    return isMITConversion(index, 0,   1e-6,   2.0,   true,  false);
  }


  bool isChargeFlip3agree(int elIndex){
    if (HWWVal::els_trkidx().at(elIndex) >= 0) {
    // false if 3 charge measurements agree
      if(
          (HWWVal::els_trk_charge().at(elIndex)                        // gsf
          == HWWVal::trks_charge().at(HWWVal::els_trkidx().at(elIndex)))  // ctf 
          &&
          (HWWVal::trks_charge().at(HWWVal::els_trkidx().at(elIndex))     // ctf 
          == HWWVal::els_sccharge().at(elIndex)) )                     // sc
      return false;  
    }
    return true;
  }

  /////////////////////////
  // position correction //
  /////////////////////////
  void electronCorrection_pos( const unsigned int index, float &dEtaIn, float &dPhiIn ) {

      //
      // uncorrected dEtaIn and dPhiIn
      //

      dEtaIn = HWWVal::els_dEtaIn().at(index);
      dPhiIn = HWWVal::els_dPhiIn().at(index);

      //
      // if configered not to apply correction
      // or in barrel or no valid super cluster
      // return uncorrected values
      //

      if (!(HWWVal::els_fiduciality().at(index) & 1<<ISEE)) return;
      if (HWWVal::els_scindex().at(index) == -1) return;

      //
      // set up correction parameters for EE+ and EE-
      // RecoEgamma/EgammaTools/python/correctedElectronsProducer_cfi.py?revision=1.2
      //

      //                                      X',     Y',     Z'
      float scPositionCorrectionEEP[3] = {   0.52,   -0.81,  0.81};
      float scPositionCorrectionEEM[3] = {    -0.02,  -0.81,  -0.94};

      LorentzVector initial_pos = HWWVal::scs_pos_p4().at(HWWVal::els_scindex().at(index));
      LorentzVector corrected_pos;

      //
      // work out corrected position
      //

      if (HWWVal::els_etaSC().at(index) < 0) {
          corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEM[0],
                  initial_pos.y() + scPositionCorrectionEEM[1],
                  initial_pos.z() + scPositionCorrectionEEM[2], 0.0);
      }
      if (HWWVal::els_etaSC().at(index) > 0) {
          corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEP[0],
                  initial_pos.y() + scPositionCorrectionEEP[1],
                  initial_pos.z() + scPositionCorrectionEEP[2], 0.0);
      }

      //
      // work out correction to dEtaIn and dPhiIn
      //

      float deta_sc = corrected_pos.Eta() - initial_pos.Eta();
      float dphi_sc = acos(cos(corrected_pos.Phi() - initial_pos.Phi()));
      dEtaIn = deta_sc + HWWVal::els_dEtaIn().at(index);
      dPhiIn = acos(cos(dphi_sc + HWWVal::els_dPhiIn().at(index)));

  }

  ////////
  // d0 //
  ////////

  double electron_d0PV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dxyPV = HWWVal::els_d0().at(index)-
      HWWVal::vtxs_position().at(vtxIndex).x()*sin(HWWVal::els_trk_p4().at(index).phi())+
      HWWVal::vtxs_position().at(vtxIndex).y()*cos(HWWVal::els_trk_p4().at(index).phi());
    return dxyPV;
  }

  double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
  }
  double electron_dzPV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dzpv = dzPV(HWWVal::els_vertex_p4().at(index), HWWVal::els_trk_p4().at(index), HWWVal::vtxs_position().at(vtxIndex));
    return dzpv;
  }

  double electron_dzPV_wwV1(unsigned int index){ 
      if ( HWWVal::vtxs_sumpt().empty() ) return 9999.;
      double sumPtMax = -1;
      int iMax = -1;
      for ( unsigned int i = 0; i < HWWVal::vtxs_sumpt().size(); ++i ){
          if (!isGoodVertex(i)) continue;
          if ( HWWVal::vtxs_sumpt().at(i) > sumPtMax ){
              iMax = i;
              sumPtMax = HWWVal::vtxs_sumpt().at(i);
          }
      }
      if (iMax<0) return 9999.;

      const LorentzVector& vtx = HWWVal::els_vertex_p4().at(index);
      const LorentzVector& p4 = HWWVal::els_trk_p4().at(index);
      const LorentzVector& pv = HWWVal::vtxs_position().at(iMax); 
      return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
  }

  double electron_d0PV_wwV1(unsigned int index){ 
      if ( HWWVal::vtxs_sumpt().empty() ) return 9999.;
      double sumPtMax = -1;
      int iMax = -1;
      for ( unsigned int i = 0; i < HWWVal::vtxs_sumpt().size(); ++i ){
          if (!isGoodVertex(i)) continue;
          if ( HWWVal::vtxs_sumpt().at(i) > sumPtMax ){
              iMax = i;
              sumPtMax = HWWVal::vtxs_sumpt().at(i);
          }
      }
      if (iMax<0) return 9999.;
      double dxyPV = HWWVal::els_d0().at(index)-
          HWWVal::vtxs_position().at(iMax).x()*sin(HWWVal::els_trk_p4().at(index).phi())+
          HWWVal::vtxs_position().at(iMax).y()*cos(HWWVal::els_trk_p4().at(index).phi());
      return dxyPV;
  }


  // this is now redundant
  float electronIsoValuePF2012_FastJetEffArea_HWW(int index){

      const float etaAbs = fabs(HWWVal::els_etaSC().at(index));
      const float pt     = HWWVal::els_p4().at(index).pt();

      // get effective area
      const float AEff = fastJetEffArea04_v1(etaAbs);

      // pf iso
      // calculate from the ntuple for now...
      const float pfiso_ch = HWWVal::els_iso04_pf2012_ch().at(index);
      const float pfiso_em = HWWVal::els_iso04_pf2012_em().at(index);
      const float pfiso_nh = HWWVal::els_iso04_pf2012_nh().at(index);

      // rho
      const float rhoPrime = std::max(HWWVal::evt_ww_rho(), 0.0f);
      const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);  
      const float pfiso = (pfiso_ch + pfiso_n) / pt;   

    // debug
    if(0) {
       LogDebug("electronSelections") << "AEff : " << AEff << " "
                                      << "rho : " << rhoPrime << " "
                                      << "pfiso_ch : " << pfiso_ch << " "
                                      << "pfiso_em : " << pfiso_em << " "
                                      << "pfiso_nh : " << pfiso_nh << " "
                                      << "pfiso_n : " << pfiso_n << " "
                                      << "pfiso : " << pfiso << " "
                                      << "pt : " << pt << " "
                                      << "etaAbs : " << etaAbs;
    }

      return pfiso;
  }

}
