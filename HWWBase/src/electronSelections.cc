/////////////////////////
// Electron Selections //
/////////////////////////

#include "Math/VectorUtil.h"

#include "HWWValidation/HWWBase/interface/electronSelections.h"
#include "HWWValidation/HWWBase/interface/eventSelections.h"
#include "HWWValidation/HWWBase/interface/MITConversionUtilities.h"
#include "HWWValidation/HWWBase/interface/muonSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"

using namespace std;
using namespace wp2012;

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
    if (fabs(hww.els_p4().at(index).eta()) < 1.479){
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
    if (electronId_smurf_v1ss(index)) cuts_passed |= (1ll<<ELEID_SMURFV1SS_EXTRA);
    if (electronId_smurf_v2ss(index)) cuts_passed |= (1ll<<ELEID_SMURFV2SS_EXTRA);

    // 2012 ID
    // electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, LOOSE);
    // if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) cuts_passed |= (1ll<<ELEID_WP2012_LOOSE_NOISO);
    // if ((answer_loose_2012 & PassWP2012CutsNoIsoNoIP) == PassWP2012CutsNoIsoNoIP) cuts_passed |= (1ll<<ELEID_WP2012_LOOSE_NOISO_NOIP);

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
    // VBTF70 (optimised in 35X)
    // answer_vbtf = electronId_VBTF(index, VBTF_35Xr2_70, applyAlignmentCorrection, removedEtaCutInEndcap);
    // if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_70);
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

    // CIC ID  
    // MEDIUM (V03 optimisation)
    electronIdComponent_t answer_cic = electronId_CIC(index, 3, CIC_MEDIUM, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_CIC_V03_MEDIUM);

    //////////////////////////
    // Conversion Rejection //
    //////////////////////////
    if (!isFromConversionPartnerTrack(index)) cuts_passed |= (1ll<<ELENOTCONV_DISTDCOT002);
    if (!isFromConversionPartnerTrack_v2(index)) cuts_passed |= (1ll<<ELENOTCONV_DISTDCOT002_OLD);
    if (!isFromConversionHitPattern(index)) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN);
    if (hww.els_exp_innerlayers().at(index) == 0) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN_0MHITS);
    if(!isFromConversionMIT(index)) cuts_passed |= (1ll<<ELENOTCONV_MIT);

    /////////////////
    // Fiduciality //
    /////////////////
    if ((hww.els_type().at(index) & (1ll<<ISECALDRIVEN))) cuts_passed |= (1ll<<ELESEED_ECAL);
    if (fabs(hww.els_p4().at(index).eta()) < 2.5) cuts_passed |= (1ll<<ELEETA_250);
    if (fabs(hww.els_p4().at(index).eta()) < 2.4) cuts_passed |= (1ll<<ELEETA_240);
    //if (electronId_noMuon(index)) cuts_passed |= (1ll<<ELENOMUON_010);
    //if (electronId_noMuon_SS(index)) cuts_passed |= (1ll<<ELENOMUON_010_SS);

    ////////
    // Pt //
    ////////
    if( hww.els_p4().at(index).pt() > 10.0 ) cuts_passed |= (1ll<<ELEPT_010);

    // Veto electron in transition region
    if( fabs(hww.els_etaSC().at(index)) < 1.4442 || fabs(hww.els_etaSC().at(index)) > 1.566 )  cuts_passed |= (1ll<<ELE_NOT_TRANSITION);

    /////////////////
    // Charge Flip //
    /////////////////
    if (!isChargeFlip3agree(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP3AGREE);

    // return which selections passed
    return cuts_passed;

}

///////////////
// Muon Veto //
//////////////

bool electronId_noMuon( const unsigned int index ) {
    if ( hww.els_closestMuon().at(index) != -1) return false;
    return true;
}

bool electronId_noMuon_SS( const unsigned int index ) {
    int idx = hww.els_closestMuon().at(index);
    if (idx < 0) return true;
    if (muonId(idx, NominalSSv5)) return false;
    return true;
}


////////////////////
// Identification //
////////////////////

// SMURF

bool electronId_smurf_v1(const unsigned int index)
{

  if (hww.els_p4().at(index).pt() > 20.0) return true;

  if (hww.els_fbrem().at(index) > 0.15) return true;
  
  if (fabs(hww.els_etaSC().at(index)) < 1.) {
    if (hww.els_eOverPIn().at(index) > 0.95 && fabs(hww.els_dPhiIn().at(index)*hww.els_charge().at(index)) < 0.006) return true; 
  }

  return false;
}

bool electronId_smurf_v2(const unsigned int index)	 
 {	 
 	 
   if (hww.els_p4().at(index).pt() > 20.0) return true;	 
 	 
   if (hww.els_fbrem().at(index) > 0.15) return true;	 
 	 
   if (fabs(hww.els_etaSC().at(index)) < 1.) {	 
     if (hww.els_eOverPIn().at(index) > 0.95) return true;	 
   }	 
 	 
   return false;	 
 }	 
 	 
 bool electronId_smurf_v3(const unsigned int index)	 
 {	 
   if (fabs(hww.els_etaSC().at(index)) > 1.479 && hww.els_hOverE().at(index) > 0.1) return false;
 	 
   if (hww.els_p4().at(index).pt() > 20.0) return true;	 
 	 
   electronIdComponent_t answer_vbtf = 0;	 
   answer_vbtf = electronId_VBTF(index, VBTF_70_NOHOEEND, false, false);	 
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {	 
 	 
     if (hww.els_fbrem().at(index) > 0.15) return true;	 
 	 
     if (fabs(hww.els_etaSC().at(index)) < 1.) {	 
       if (hww.els_eOverPIn().at(index) > 0.95) return true;	 
     }	 
 	 
   }	 
 	 
   return false;	 
 }

 bool electronId_smurf_v1ss(const unsigned int index)
 {

   if (hww.els_p4().at(index).pt() > 20.0) return true;

   electronIdComponent_t answer_vbtf = 0;
   answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, false, false);
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {

     if (hww.els_fbrem().at(index) > 0.15) return true;

     if (fabs(hww.els_etaSC().at(index)) < 1.) {
       if (hww.els_eOverPIn().at(index) > 0.95) return true;
     }

   }

   return false;
 }

 bool electronId_smurf_v2ss(const unsigned int index)
 {

   if (hww.els_p4().at(index).pt() > 20.0) return true;

   electronIdComponent_t answer_vbtf = 0;
   answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, false, false);
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {

       if (fabs(hww.els_etaSC().at(index)) > 1.479)
           if (hww.els_hOverE().at(index) < 0.15) return true;

     if (hww.els_fbrem().at(index) > 0.15) return true;

     if (fabs(hww.els_etaSC().at(index)) < 1.) {
       if (hww.els_eOverPIn().at(index) > 0.95) return true;
     }

   }

   return false;
 }


// CIC
electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness, bool applyAlignementCorrection, bool removedEtaCutInEndcap)
{

    // check that a valid version number was supplied
    if (version != 2 && version != 3 && version != 4 && version != 6) {
        std::cout << "[electronId_CIC] Error! Version must be 2, 3, 4 or 6 - fail" << std::endl;
        return 0;
    }

    //
    // get the variables that are going to be cut on
    //

    double scTheta = (2*atan(exp(-1*hww.els_etaSC().at(index))));
    double scEt = hww.els_eSC().at(index)*sin(scTheta); 
    double scEta = hww.els_etaSC().at(index);
    double eSeedOverPin = hww.els_eSeedOverPIn().at(index);
    double fBrem = hww.els_fbrem().at(index);
    double hOverE = hww.els_hOverE().at(index);
    double sigmaee = hww.els_sigmaIEtaIEta().at(index);
    int mishits = hww.els_exp_innerlayers().at(index);
    double dist = (hww.els_conv_dist().at(index) == -9999.? 9999:hww.els_conv_dist().at(index));
    double dcot = (hww.els_conv_dcot().at(index) == -9999.? 9999:hww.els_conv_dcot().at(index));
    float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);
    double tkIso = hww.els_tkIso().at(index);
    double ecalIso = hww.els_ecalIso04().at(index);
    double hcalIso = hww.els_hcalIso04().at(index);

    double ip = fabs(hww.els_d0().at(index));
    if (hww.vtxs_sumpt().size() > 0) {
        const float vx = hww.vtxs_position().at(0).x();
        const float vy = hww.vtxs_position().at(0).y();
        const float px = hww.els_trk_p4().at(index).px();
        const float py = hww.els_trk_p4().at(index).py();
        ip = fabs(-1*hww.els_d0().at(index) + (vx*py - vy*px) / hww.els_trk_p4().at(index).pt());
    }

    //
    // get corrected dEtaIn and dPhiIn
    //

    float deltaEtaIn = hww.els_dEtaIn().at(index);
    float deltaPhiIn = hww.els_dPhiIn().at(index);
    if (applyAlignementCorrection) electronCorrection_pos(index, deltaEtaIn, deltaPhiIn);

    // find the catagory for this electron
    unsigned int cat = eidClassify(version, index);

    // determine if in EB or EE
    int eb;
    if (hww.els_fiduciality().at(index) & (1<<ISEB)) eb = 0;
    else eb = 1;

    // Version V02
    if (version == 2) {

        // set the parameters for the chosen tightness
        std::vector<double> cutdeta;
        std::vector<double> cutdphi;
        std::vector<double> cuteopin;
        std::vector<double> cutet;
        std::vector<double> cuthoe;
        std::vector<double> cutip;
        std::vector<double> cutisoecal;
        std::vector<double> cutisohcal;
        std::vector<double> cutisotk;
        std::vector<double> cutmishits;
        std::vector<double> cutsee;

        eidGetCIC_V02(tightness, cutdeta, cutdphi, cuteopin, cutet, cuthoe, 
                cutip, cutisoecal, cutisohcal, cutisotk, cutmishits, cutsee);


        unsigned int result = 0;
        int bin = 0;
        if (scEt < 20.)
            bin = 2;
        else if (scEt > 30.)
            bin = 0;
        else
            bin = 1;

        if (fBrem > 0)
            eSeedOverPin = eSeedOverPin + fBrem;

        if (bin != 2) {     
            tkIso = tkIso*pow(40./scEt, 2); 
            ecalIso = ecalIso*pow(40./scEt, 2); 
            hcalIso = hcalIso*pow(40./scEt, 2); 
        }

        if ((tkIso < cutisotk[cat+3*eb+bin*6]) &&
                (ecalIso < cutisoecal[cat+3*eb+bin*6]) &&
                (hcalIso < cutisohcal[cat+3*eb+bin*6]))
            result |= (1<<ELEID_ISO);

        if (fBrem < -2)
            return result;

        bool passdEtaCut = fabs(deltaEtaIn) < cutdeta[cat+3*eb+bin*6];


        if(removedEtaCutInEndcap && eb == 1) passdEtaCut = true;
        if (hOverE < cuthoe[cat+3*eb+bin*6] &&
                sigmaee < cutsee[cat+3*eb+bin*6] &&
                fabs(deltaPhiIn) < cutdphi[cat+3*eb+bin*6] &&
                passdEtaCut &&
                eSeedOverPin > cuteopin[cat+3*eb+bin*6])
            result |= (1<<ELEID_ID);

        if (ip < cutip[cat+3*eb+bin*6])
            result |= (1<<ELEID_IP);

        if (mishits < cutmishits[cat+3*eb+bin*6])
            result |= (1<<ELEID_CONV);

        return result;
    }

    // version V03, V04 or V05
    if (version == 3 || version == 4 || version == 5) {

        //
        // set the parameters for the chosen tightness
        //
        std::vector<double> cutdcotdist;
        std::vector<double> cutdetain;
        std::vector<double> cutdphiin;
        std::vector<double> cuteseedopcor;
        std::vector<double> cutet;
        std::vector<double> cutfmishits;
        std::vector<double> cuthoe;
        std::vector<double> cutip_gsf;
        std::vector<double> cutiso_sum;
        std::vector<double> cutiso_sumoet;
        std::vector<double> cutsee;

        // V03 uses Et binning
        bool wantBinning = true;

        if (version == 3) {
            eidGetCIC_V03(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
        }

        if (version == 4) {
            eidGetCIC_V04(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
            // V04 does not use Et binning
            wantBinning = false;
        }

        // this is certainly true for V03
        // but not sure the meaning of V04,05
        unsigned int result = 0;
        int bin = 0;
        if (wantBinning) {
            if (scEt < 20.)
                bin = 2;
            else if (scEt > 30.)
                bin = 0;
            else
                bin = 1;
        }

        if (fBrem > 0)
            eSeedOverPin = eSeedOverPin + fBrem;

        float iso_sum = tkIso + ecalIso + hcalIso;
        float iso_sum_corrected = iso_sum*pow(40./scEt, 2);

        if ((iso_sum < cutiso_sum[cat+bin*9]) &&
                (iso_sum_corrected < cutiso_sumoet[cat+bin*9]))
            result |= (1<<ELEID_ISO);

        if (fBrem > -2) {
            if ((hOverE < cuthoe[cat+bin*9]) and
                    (sigmaee < cutsee[cat+bin*9]) and
                    (fabs(deltaPhiIn) < cutdphiin[cat+bin*9]) and
                    (fabs(deltaEtaIn) < cutdetain[cat+bin*9]) and
                    (eSeedOverPin > cuteseedopcor[cat+bin*9]) and
                    (scEt > cutet[cat+bin*9]))
                result |= (1<<ELEID_ID);
        }

        if (ip < cutip_gsf[cat+bin*9])
            result |= (1<<ELEID_IP);

        if ((mishits < cutfmishits[cat+bin*9]) and 
                (dcotdistcomb < cutdcotdist[cat+bin*9]))
            result |= (1<<ELEID_CONV);

        return result;
    }

    if (version == 6) {

        std::vector<double> cutIsoSum;
        std::vector<double> cutIsoSumCorr;
        std::vector<double> cuthoe;
        std::vector<double> cutsee;
        std::vector<double> cutdphi;
        std::vector<double> cutdeta;
        std::vector<double> cuteopin;
        std::vector<double> cutmishits;
        std::vector<double> cutdcotdist;
        std::vector<double> cutip;
        std::vector<double> cutIsoSumCorrl;
        std::vector<double> cuthoel;
        std::vector<double> cutseel;
        std::vector<double> cutdphil;
        std::vector<double> cutdetal;
        std::vector<double> cutipl;

        eidGetCIC_V06(tightness, cutIsoSum, cutIsoSumCorr, cuthoe, cutsee, cutdphi, 
                cutdeta, cuteopin, cutmishits, cutdcotdist, cutip, cutIsoSumCorrl, 
                cuthoel, cutseel, cutdphil, cutdetal, cutipl);

        int result = 0;
        const int ncuts = 10;
        std::vector<bool> cut_results(ncuts, false);

        float iso_sum = tkIso + ecalIso + hcalIso;
        if(fabs(scEta)>1.5) 
            iso_sum += (fabs(scEta)-1.5)*1.09;

        float iso_sumoet = iso_sum*(40./scEt);
        float eseedopincor = eSeedOverPin + fBrem;
        if(fBrem < 0) eseedopincor = eSeedOverPin;

        for (int cut=0; cut<ncuts; cut++) {

            switch (cut) {
                case 0:
                    cut_results[cut] = eidComputeCut(fabs(deltaEtaIn), scEt, cutdetal[cat], cutdeta[cat]);
                    break;
                case 1:
                    cut_results[cut] = eidComputeCut(fabs(deltaPhiIn), scEt, cutdphil[cat], cutdphi[cat]);
                    break;
                case 2:
                    cut_results[cut] = (eseedopincor > cuteopin[cat]);
                    break;
                case 3:
                    cut_results[cut] = eidComputeCut(hOverE, scEt, cuthoel[cat], cuthoe[cat]);
                    break;
                case 4:
                    cut_results[cut] = eidComputeCut(sigmaee, scEt, cutseel[cat], cutsee[cat]);
                    break;
                case 5:
                    cut_results[cut] = eidComputeCut(iso_sumoet, scEt, cutIsoSumCorrl[cat], cutIsoSumCorr[cat]);
                    break;
                case 6:
                    cut_results[cut] = (iso_sum < cutIsoSum[cat]);
                    break;
                case 7:
                    cut_results[cut] = eidComputeCut(fabs(ip), scEt, cutipl[cat], cutip[cat]);
                    break;
                case 8:
                    cut_results[cut] = (mishits < cutmishits[cat]);
                    break;
                case 9:
                    cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
                    break;
            }

        }

        // ID part
        if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4]) result |= (1<<ELEID_ID);   

        // ISO part
        if (cut_results[5] & cut_results[6]) result |= (1<<ELEID_ISO);

        // IP part
        if (cut_results[7]) result |= (1<<ELEID_IP);

        // Conversion part
        if (cut_results[8] and cut_results[9]) result |= (1<<ELEID_CONV);

        return result;

    }

    std::cout << "[electronId_CIC] Error! got to the end and didn't apply any cuts - fail" << std::endl;
    return 0;

}

unsigned int eidClassify(const unsigned int version, const unsigned int index) {

    // variables used for classification
    double eta = fabs(hww.els_etaSC().at(index));
    double eOverP = hww.els_eOverPIn().at(index);
    double fBrem = hww.els_fbrem().at(index);
    bool isEB = (hww.els_fiduciality().at(index) & (1<<ISEB));
    bool isEE = (hww.els_fiduciality().at(index) & (1<<ISEE));
    bool ecalDriven = (hww.els_type().at(index) & (1<<ISECALDRIVEN));
    bool trackerDriven = (hww.els_type().at(index) & (1<<ISTRACKERDRIVEN));

    int cat = -1;

    //
    // version V00 or V01
    //
    if (version == 0 || version == 1) {
        if((isEB && fBrem<0.06) || (isEE && fBrem<0.1)) 
            cat=1;
        else if (eOverP < 1.2 && eOverP > 0.8) 
            cat=0;
        else 
            cat=2;
        return cat;
    }

    // version V02
    if (version == 2) {
        if (isEB) {       // BARREL
            if(fBrem < 0.12)
                cat=1;
            else if (eOverP < 1.2 && eOverP > 0.9) 
                cat=0;
            else 
                cat=2;
        } else {                     // ENDCAP
            if(fBrem < 0.2)
                cat=1;
            else if (eOverP < 1.22 && eOverP > 0.82) 
                cat=0;
            else 
                cat=2;
        }
        return cat;
    }

    //
    // version V03, V04 or V05
    // took this from V00-03-07-02 of 
    // ElectronIdentification/src/CutBasedElectronID.cc
    // NOTE - DLE - now adding version 6
    // this means newCategories is switched on
    //


    if (version == 3 || version == 4 || version == 5 || version == 6) {

        // this is certainly true for V03 and V04
        // not sure about V05
        bool newCategories = false;
        if (version == 6) newCategories = true;

        if (isEB) {
            if ((fBrem >= 0.12) && (eOverP > 0.9) && (eOverP < 1.2))
                cat = 0;
            else if (((eta >  .445   && eta <  .45  ) ||
                        (eta >  .79    && eta <  .81  ) ||
                        (eta > 1.137   && eta < 1.157 ) ||
                        (eta > 1.47285 && eta < 1.4744)) && newCategories)
                cat = 6;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.12)
                cat = 1;
            else
                cat = 2;
        } else {
            if ((fBrem >= 0.2) && (eOverP > 0.82) && (eOverP < 1.22))
                cat = 3;
            else if (eta > 1.5 && eta <  1.58 && newCategories)
                cat = 7;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.2)
                cat = 4;
            else
                cat = 5;
        }

        return cat;
    }


    return cat;

}

bool eidComputeCut(double x, double et, double cut_min, double cut_max, bool gtn)
{

  float et_min = 10;
  float et_max = 40;
  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  //std::cout << x << " " << cut_min << " " << cut << " " << cut_max << " " << et << " " << accept << std::endl;
  return accept;
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
    unsigned int det = ((hww.els_fiduciality().at(index) & (1<<ISEB)) == (1<<ISEB)) ? 0 : 1;
    float etaAbs = fabs(hww.els_etaSC().at(index));
    float pt     = hww.els_p4().at(index).pt();

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
    float pfiso_ch = hww.els_iso03_pf2012_ch().at(index);
    float pfiso_em = hww.els_iso03_pf2012_em().at(index);
    float pfiso_nh = hww.els_iso03_pf2012_nh().at(index);

    // rho
    float rhoPrime = std::max(hww.evt_kt6pf_foregiso_rho(), float(0.0));
    float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
    float pfiso = (pfiso_ch + pfiso_n) / pt;

    // |1/E - 1/p|
    float ooemoop = fabs( (1.0/hww.els_ecalEnergy().at(index)) - (hww.els_eOverPIn().at(index)/hww.els_ecalEnergy().at(index)) );

    // MIT conversion vtx fit
    bool vtxFitConversion = isMITConversion(index, 0,   1e-6,   2.0,   true,  false);

    // d0
    float d0vtx = electron_d0PV_smurfV3(index);
    float dzvtx = electron_dzPV_smurfV3(index);
 
    // test cuts
    if (fabs(hww.els_dEtaIn().at(index)) < dEtaInThresholds[det])             mask |= wp2012::DETAIN;
    if (fabs(hww.els_dPhiIn().at(index)) < dPhiInThresholds[det])             mask |= wp2012::DPHIIN;
    if (hww.els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[det])     mask |= wp2012::SIGMAIETAIETA;
    if (hww.els_hOverE().at(index) < hoeThresholds[det])                      mask |= wp2012::HOE;
    if (ooemoop < ooemoopThresholds[det])                                   mask |= wp2012::OOEMOOP;
    if (fabs(d0vtx) < d0VtxThresholds[det])                                 mask |= wp2012::D0VTX;
    if (fabs(dzvtx) < dzVtxThresholds[det])                                 mask |= wp2012::DZVTX;
    if (!vtxFitThresholds[det] || !vtxFitConversion)                        mask |= wp2012::VTXFIT;
    if (hww.els_exp_innerlayers().at(index) <= mHitsThresholds[det])          mask |= wp2012::MHITS;
    if (pt >= 20.0 && pfiso < isoHiThresholds[det])                         mask |= wp2012::ISO;
    if (pt < 20.0 && pfiso < isoLoThresholds[det])                          mask |= wp2012::ISO;

    // return the mask
    return mask;

}


// do be used in the functions below for consize parameter
bool is_equal(float lhs, float rhs)
{
	static const float epsilon = 0.0001;
	return (fabs(lhs-rhs) < epsilon);
}

// calculate Effective area (updated to value from Egamma)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
// Topic revision: r12 - 28-Nov-2012
float fastJetEffArea03_v2(const float eta)
{
	// use absolute eta
    const float etaAbs = fabs(eta);

    // get effective area
    if      (etaAbs <= 1.0                  ) {return 0.13;}
    else if (etaAbs > 1.0 && etaAbs <= 1.479) {return 0.14;}
    else if (etaAbs > 1.479 && etaAbs <= 2.0) {return 0.07;}
    else if (etaAbs > 2.0 && etaAbs <= 2.2  ) {return 0.09;}
    else if (etaAbs > 2.2 && etaAbs <= 2.3  ) {return 0.11;}
    else if (etaAbs > 2.3 && etaAbs <= 2.4  ) {return 0.11;}
    else if (etaAbs > 2.4                   ) {return 0.14;}
    return -9999.0f;
}

float fastJetEffArea04_v2(const float eta)
{
	// use absolute eta
    const float etaAbs = fabs(eta);

    // get effective area
    if      (etaAbs <= 1.0                  ) {return 0.21;}
    else if (etaAbs > 1.0 && etaAbs <= 1.479) {return 0.21;}
    else if (etaAbs > 1.479 && etaAbs <= 2.0) {return 0.11;}
    else if (etaAbs > 2.0 && etaAbs <= 2.2  ) {return 0.14;}
    else if (etaAbs > 2.2 && etaAbs <= 2.3  ) {return 0.18;}
    else if (etaAbs > 2.3 && etaAbs <= 2.4  ) {return 0.19;}
    else if (etaAbs > 2.4                   ) {return 0.26;}
    return -9999.0f;
}

// calculate Effective area (updated to value from Egamma)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
// Topic revision: r10 - 11-Apr-2012
float fastJetEffArea03_v1(const float eta)
{
	// use absolute eta
    const float etaAbs = fabs(eta);

    // get effective area
    if      (etaAbs <= 1.0                  ) {return 0.10; }
    else if (etaAbs > 1.0 && etaAbs <= 1.479) {return 0.12; }
    else if (etaAbs > 1.479 && etaAbs <= 2.0) {return 0.085;}
    else if (etaAbs > 2.0 && etaAbs <= 2.2  ) {return 0.11; }
    else if (etaAbs > 2.2 && etaAbs <= 2.3  ) {return 0.12; }
    else if (etaAbs > 2.3 && etaAbs <= 2.4  ) {return 0.12; }
    else if (etaAbs > 2.4                   ) {return 0.13; }
    return -9999.0f;
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

    float dEtaIn = hww.els_dEtaIn().at(index);
    float dPhiIn = hww.els_dPhiIn().at(index);
    if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

    // barrel
    if (fabs(hww.els_etaSC().at(index)) < 1.479) {

        if (electronIsolation_rel(index, true) < relisoThresholds[0])
            answer |= (1<<ELEID_ISO);

        if (fabs(dEtaIn) < dEtaInThresholds[0] &&
                fabs(dPhiIn) < dPhiInThresholds[0] &&
                hww.els_hOverE().at(index) < hoeThresholds[0] &&
                hww.els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[0])
            answer |= (1<<ELEID_ID);
    }

    // endcap
    if (fabs(hww.els_etaSC().at(index)) > 1.479) {
        if (electronIsolation_rel(index, true) < relisoThresholds[1])
            answer |= (1<<ELEID_ISO);
        bool passdEtaCut = fabs(dEtaIn) < dEtaInThresholds[1];
        if(removedEtaCutInEndcap) passdEtaCut = true;
        if ( passdEtaCut &&
                fabs(dPhiIn) < dPhiInThresholds[1] &&
                hww.els_hOverE().at(index) < hoeThresholds[1] &&
                hww.els_sigmaIEtaIEta().at(index) < sigmaIEtaIEtaThresholds[1])
            answer |= (1<<ELEID_ID);
    }

    return answer;

}

electronIdComponent_t passLikelihoodId(unsigned int index, float lhValue, int workingPoint) {
  unsigned int answer = 0;
  float etaSC = hww.els_etaSC().at(index);
  unsigned int nbrem = hww.els_nSeed().at(index);
  if ( workingPoint == 95 ) {
    if (
      ( fabs(etaSC) < 1.479 && nbrem ==0 && lhValue > -4.274 ) ||
      ( fabs(etaSC) < 1.479 && nbrem >=1 && lhValue >- 3.773 ) ||
      ( fabs(etaSC) > 1.479 && nbrem ==0 && lhValue > -5.092 ) ||
      ( fabs(etaSC) > 1.479 && nbrem >=1 && lhValue > -2.796 )
    ) {
      answer |= (1<<ELEID_ID);
    }
  }
  else if ( workingPoint == 90 ) {
    if (
      ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > -1.497 ) ||
      ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > -1.521 ) ||
      ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > -2.571 ) ||
      ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > -0.657 )
    ) {
      answer |= (1<<ELEID_ID);
    }
  }
  else if ( workingPoint == 85 ) {
    if (
      ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > +0.163 ) ||
      ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > +0.065 ) ||
      ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > -0.683 ) ||
      ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > +1.564 )
    ) {
      answer |= (1<<ELEID_ID);
    }
  }
  else if ( workingPoint == 80 ) {
    if (
      ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > +1.193 ) ||
      ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > +1.345 ) ||
      ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > +0.810 ) ||
      ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > +3.021 )
    ) {
      answer |= (1<<ELEID_ID);
    }
  }
  else if ( workingPoint == 70 ) {
    if (
      ( fabs(etaSC) < 1.479 && nbrem == 0 && lhValue > +1.781 ) ||
      ( fabs(etaSC) < 1.479 && nbrem >= 1 && lhValue > +2.397 ) ||
      ( fabs(etaSC) > 1.479 && nbrem == 0 && lhValue > +2.361 ) ||
      ( fabs(etaSC) > 1.479 && nbrem >= 1 && lhValue > +4.052 )
    ) {
      answer |= (1<<ELEID_ID);
    }
  }
  else {
    cout << "Error! Likelihood WP not supported: " << workingPoint << ". Please choose 70, 80, 85, 90, 95" << endl;
  }
  return answer;
}

bool passLikelihoodId_v2(unsigned int index, float lhValue, int workingPoint)
{

    float etaSC = hww.els_etaSC().at(index);
    float pt = hww.els_p4().at(index).Pt();
    unsigned int nbrem = hww.els_nSeed().at(index);

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
        cout << "Error! Likelihood WP not supported: " 
             << workingPoint << ". Please choose 0 for Emanuele 8th September" << endl;
    }

    return false;
}




 

///////////////
// Isolation //
///////////////

// relative truncated
float electronIsolation_rel( const unsigned int index, bool use_calo_iso ) {
    float sum = hww.els_tkIso().at(index);
    if (use_calo_iso) {
        if (fabs(hww.els_etaSC().at(index)) > 1.479) sum += hww.els_ecalIso().at(index);
        if (fabs(hww.els_etaSC().at(index)) <= 1.479) sum += max(0., (hww.els_ecalIso().at(index) -1.));
        sum += hww.els_hcalIso().at(index);
    }
    double pt = hww.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

// relative truncated, fast-jet corrected
float electronIsolation_rel_FastJet( const unsigned int index, bool use_calo_iso ) {
    float sum = hww.els_tkIso().at(index);
    float offset = (hww.evt_rho()) * TMath::Pi() * pow( 0.3 , 2 );
    if (use_calo_iso) {
        float caloiso = 0.;

        if (fabs(hww.els_etaSC().at(index)) >  1.479) caloiso += hww.els_ecalIso().at(index);
        if (fabs(hww.els_etaSC().at(index)) <= 1.479) caloiso += max(0., (hww.els_ecalIso().at(index) -1.));
        caloiso += hww.els_hcalIso().at(index);
        
        caloiso -= offset;
        if( caloiso > 0 ) sum += caloiso;
    }
    double pt = hww.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

// Relative Isolation, Non-Truncated
float electronIsolation_rel_v1( const unsigned int index, bool use_calo_iso ) {
    float pt               = hww.els_p4().at(index).pt();          // Electron Pt
    float TRCK_sum_over_pt = hww.els_tkIso().at(index) / pt;       // Tracker Relative Isolation, Non-Truncated
    float ECAL_sum_over_pt = electronIsolation_ECAL_rel_v1(index);  // ECAL    Relative Isolation, Non-Truncated
    float HCAL_sum_over_pt = electronIsolation_HCAL_rel_v1(index);  // HCAL    Relative Isolation, Non-Truncated

    float sum_over_pt      = TRCK_sum_over_pt;                      // Combined Subdetector Relative Isolation, Non-Truncated
    if(use_calo_iso){
      sum_over_pt += ECAL_sum_over_pt;
      sum_over_pt += HCAL_sum_over_pt;
    }
    return sum_over_pt;
}

// corrected, relative isolation, non-truncated
float electronIsolation_cor_rel_v1(const unsigned int index, bool use_calo_iso) {
  float ntiso = electronIsolation_rel_v1(index, use_calo_iso);
  float pt = hww.els_p4().at(index).pt();
  int nvtxs = numberOfGoodVertices();
  float coriso = ntiso - ((TMath::Log(pt)*nvtxs)/(30*pt));
  return coriso;
}

// Relative Isolation, Non-Truncated, FastJet-corrected
float electronIsolation_rel_v1_FastJet( const unsigned int index, bool use_calo_iso ){
    float pt               = hww.els_p4().at(index).pt();          // Electron Pt
    float TRCK_sum_over_pt = hww.els_tkIso().at(index) / pt;       // Tracker Relative Isolation, Non-Truncated
    float ECAL_sum_over_pt = electronIsolation_ECAL_rel_v1(index);  // ECAL    Relative Isolation, Non-Truncated
    float HCAL_sum_over_pt = electronIsolation_HCAL_rel_v1(index);  // HCAL    Relative Isolation, Non-Truncated
    float offset           = el_fastjet_rel_offset(index);          // fastjet offset = pi X dR^2 X rho
    float sum_over_pt      = TRCK_sum_over_pt;                      // Combined Subdetector Relative Isolation, Non-Truncated
    if(use_calo_iso){
      float calo_iso = ECAL_sum_over_pt + HCAL_sum_over_pt - offset;
      if( calo_iso > 0 ) sum_over_pt += calo_iso;
    }
    return sum_over_pt;
}

// ECAL Relative Isolation, Non-Truncated
float electronIsolation_ECAL_rel_v1( const unsigned int index, bool useEBps ) {
  float pt               = hww.els_p4().at(index).pt();                                                                  // Electron Pt
  float ecal_sum_over_pt = 0.0;                                                                                           // ECAL Relative Isolation, NT
  if( fabs(hww.els_etaSC().at(index)) > 1.479  ) ecal_sum_over_pt += hww.els_ecalIso().at(index);                       // EE: Ecal Endcap  
  if( fabs(hww.els_etaSC().at(index)) <= 1.479 ) {
      if (useEBps)
          ecal_sum_over_pt += max( 0.0, ( hww.els_ecalIso().at(index) - 1.0 ) ); // EB: Ecal Barrel
      else
          ecal_sum_over_pt += hww.els_ecalIso().at(index); // EB: Ecal Barrel
  }
  ecal_sum_over_pt /= pt;
  return ecal_sum_over_pt;
}

// HCAL Relative Isolation, Non-Truncated
float electronIsolation_HCAL_rel_v1( const unsigned int index ){
  float pt               = hww.els_p4().at(index).pt();      // Electron Pt
  float hcal_sum_over_pt = hww.els_hcalIso().at(index) / pt; // HCAL Relative Isolation, NT
  return hcal_sum_over_pt;
}

float el_fastjet_rel_offset( const unsigned int index ){
  double pt     = hww.els_p4().at(index).pt();
  double offset = TMath::Pi() * pow( 0.3 , 2 ) * (hww.evt_rho());
  return offset / pt;
}

// ECAL Relative Isolation, Truncated
float electronIsolation_ECAL_rel( const unsigned int index ){
  float pt               = hww.els_p4().at(index).pt();                                                                  // Electron Pt
  float ecal_sum_over_pt = 0.0;                                                                                           // ECAL Relative Isolation
  if( fabs(hww.els_etaSC().at(index)) > 1.479  ) ecal_sum_over_pt += hww.els_ecalIso().at(index);                       // EE: Ecal Endcap
  if( fabs(hww.els_etaSC().at(index)) <= 1.479 ) ecal_sum_over_pt += max( 0.0, ( hww.els_ecalIso().at(index) - 1.0 ) ); // EB: Ecal Barrel
  ecal_sum_over_pt /= max(pt,(float)20.0);
  return ecal_sum_over_pt;
}

// HCAL Relative Isolation, Truncated
float electronIsolation_HCAL_rel( const unsigned int index ){
  float pt               = max( hww.els_p4().at(index).pt() , (float) 20.0 ) ;      // Electron Pt
  float hcal_sum_over_pt = hww.els_hcalIso().at(index) / pt;                        // HCAL Relative Isolation, NT
  return hcal_sum_over_pt;
}

// electron isolation definitions for WW analysis
float electronIsolation_rel_ww( const unsigned int index, bool use_calo_iso ) {
    float sum = hww.els_tkIso().at(index);
    if(use_calo_iso)
        sum += max(0., (hww.els_ecalIso().at(index) -1.));
    sum += hww.els_hcalIso().at(index);
    double pt = hww.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

float electronIsoValuePF( const unsigned int iel, unsigned int ivtx, float coner, float minptn, float dzcut, float footprintdr, float gammastripveto, float elestripveto, int filterId ) {

  int elgsftkid = hww.els_gsftrkidx().at(iel);
  int eltkid = hww.els_trkidx().at(iel);
  //take dz from gsf, and if it does not exist (should always exist) take it from ctf track
  float eldz = elgsftkid>=0 ? gsftrks_dz_pv( elgsftkid,ivtx ).first : trks_dz_pv(eltkid,ivtx).first;
  float eleta = hww.els_p4().at(iel).eta();

  float pfciso = 0.;
  float pfniso = 0.;
  float pffootprint = 0.;
  float pfjurveto = 0.;
  float pfjurvetoq = 0.;
  for (unsigned int ipf=0; ipf<hww.pfcands_p4().size(); ++ipf){

    float dR = ROOT::Math::VectorUtil::DeltaR( hww.pfcands_p4().at(ipf), hww.els_p4().at(iel) );

    if (dR>coner) continue;

    float pfpt = hww.pfcands_p4().at(ipf).pt();    
    float pfeta = hww.pfcands_p4().at(ipf).eta();    
    float deta = fabs(pfeta - eleta);
    int pfid = abs(hww.pfcands_particleId().at(ipf));

    if (filterId!=0 && filterId!=pfid) continue;

    if (hww.pfcands_charge().at(ipf)==0) {
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
      int pftkid = hww.pfcands_trkidx().at(ipf);
      if (eltkid>=0 && pftkid>=0 && eltkid==pftkid) continue;
      if (pfid==11 && hww.pfcands_pfelsidx().at(ipf)>=0 && hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))>=0) {
	int pfgsfid = hww.els_gsftrkidx().at(hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))); 
	if (elgsftkid>=0 && pfgsfid>=0 && elgsftkid==pfgsfid) continue;
      }
      //check electrons with gsf track
      if (pfid==11 && hww.pfcands_pfelsidx().at(ipf)>=0 && hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))>=0) {
	int gsfid = hww.els_gsftrkidx().at(hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))); 
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
	if(fabs( trks_dz_pv(hww.pfcands_trkidx().at(ipf),ivtx).first - eldz )<dzcut) {//dz cut
	  pfciso+=pfpt;
	  if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
	}
      }
    } 
  }

  return (pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq)/hww.els_p4().at(iel).pt();

}

//////////////////////////
// Conversion Rejection //
//////////////////////////

bool isFromConversionHitPattern( const unsigned int index ) {
    if(hww.els_exp_innerlayers().at(index) > 1) return true;
    return false;
}

bool isFromConversionPartnerTrack(const unsigned int index) {
    if( fabs(hww.els_conv_dist().at(index)) < 0.02 && fabs(hww.els_conv_dcot().at(index)) < 0.02 ) return true;
    return false;
}

bool isFromConversionPartnerTrack_v2(const unsigned int index) {
    if (fabs(hww.els_conv_old_dist().at(index)) < 0.02 && fabs(hww.els_conv_old_dcot().at(index)) < 0.02 ) return true;
    return false;
}

bool isFromConversionMIT(const unsigned int index){
  return isMITConversion(index, 0,   1e-6,   2.0,   true,  false);
}


/////////////////
// Charge Flip //
/////////////////

int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits) {
    if(hww.els_sccharge().at(elIdx)*hww.els_trk_charge().at(elIdx) > 0 || (hww.els_trkidx().at(elIdx) < 0 || hww.els_trkshFrac().at(elIdx) < minFracSharedHits))
        return hww.els_sccharge().at(elIdx);
    else 
        return  hww.trks_charge().at(hww.els_trkidx().at(elIdx));

}

bool isChargeFlip(int elIndex){
    if ((hww.els_trkidx().at(elIndex) >= 0) && (hww.els_trk_charge().at(elIndex) != hww.trks_charge().at(hww.els_trkidx().at(elIndex))) ) return true;
    if ((hww.els_trkidx().at(elIndex) < 0)  && (hww.els_trk_charge().at(elIndex) != hww.els_sccharge().at(elIndex))) return true;
    return false;
}

bool isChargeFlip3agree(int elIndex){
  if (hww.els_trkidx().at(elIndex) >= 0) {
  // false if 3 charge measurements agree
    if(
        (hww.els_trk_charge().at(elIndex)                        // gsf
        == hww.trks_charge().at(hww.els_trkidx().at(elIndex)))  // ctf 
        &&
        (hww.trks_charge().at(hww.els_trkidx().at(elIndex))     // ctf 
        == hww.els_sccharge().at(elIndex)) )                     // sc
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

    dEtaIn = hww.els_dEtaIn().at(index);
    dPhiIn = hww.els_dPhiIn().at(index);

    //
    // if configered not to apply correction
    // or in barrel or no valid super cluster
    // return uncorrected values
    //

    if (!(hww.els_fiduciality().at(index) & 1<<ISEE)) return;
    if (hww.els_scindex().at(index) == -1) return;

    //
    // set up correction parameters for EE+ and EE-
    // RecoEgamma/EgammaTools/python/correctedElectronsProducer_cfi.py?revision=1.2
    //

    //                                      X',     Y',     Z'
    float scPositionCorrectionEEP[3] = {   0.52,   -0.81,  0.81};
    float scPositionCorrectionEEM[3] = {    -0.02,  -0.81,  -0.94};

    LorentzVector initial_pos = hww.scs_pos_p4().at(hww.els_scindex().at(index));
    LorentzVector corrected_pos;

    //
    // work out corrected position
    //

    if (hww.els_etaSC().at(index) < 0) {
        corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEM[0],
                initial_pos.y() + scPositionCorrectionEEM[1],
                initial_pos.z() + scPositionCorrectionEEM[2], 0.0);
    }
    if (hww.els_etaSC().at(index) > 0) {
        corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEP[0],
                initial_pos.y() + scPositionCorrectionEEP[1],
                initial_pos.z() + scPositionCorrectionEEP[2], 0.0);
    }

    //
    // work out correction to dEtaIn and dPhiIn
    //

    float deta_sc = corrected_pos.Eta() - initial_pos.Eta();
    float dphi_sc = acos(cos(corrected_pos.Phi() - initial_pos.Phi()));
    dEtaIn = deta_sc + hww.els_dEtaIn().at(index);
    dPhiIn = acos(cos(dphi_sc + hww.els_dPhiIn().at(index)));

}

////////
// d0 //
////////

double electron_d0PV(unsigned int index){
    if ( hww.vtxs_sumpt().empty() ) return false;
    unsigned int iMax = 0;
    double sumPtMax = hww.vtxs_sumpt().at(0);
    for ( unsigned int i = iMax+1; i < hww.vtxs_sumpt().size(); ++i )
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    double dxyPV = hww.els_d0().at(index)-
        hww.vtxs_position().at(iMax).x()*sin(hww.els_trk_p4().at(index).phi())+
        hww.vtxs_position().at(iMax).y()*cos(hww.els_trk_p4().at(index).phi());
    return dxyPV;
}


double electron_d0PV_smurfV3(unsigned int index){
  int vtxIndex = 0;
  double dxyPV = hww.els_d0().at(index)-
    hww.vtxs_position().at(vtxIndex).x()*sin(hww.els_trk_p4().at(index).phi())+
    hww.vtxs_position().at(vtxIndex).y()*cos(hww.els_trk_p4().at(index).phi());
  return dxyPV;
}

double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}
double electron_dzPV_smurfV3(unsigned int index){
  int vtxIndex = 0;
  double dzpv = dzPV(hww.els_vertex_p4().at(index), hww.els_trk_p4().at(index), hww.vtxs_position().at(vtxIndex));
  return dzpv;
}

double electron_dzPV_wwV1(unsigned int index){ 
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;

    const LorentzVector& vtx = hww.els_vertex_p4().at(index);
    const LorentzVector& p4 = hww.els_trk_p4().at(index);
    const LorentzVector& pv = hww.vtxs_position().at(iMax); 
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
}

double electron_d0PV_wwV1(unsigned int index){ 
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    double dxyPV = hww.els_d0().at(index)-
        hww.vtxs_position().at(iMax).x()*sin(hww.els_trk_p4().at(index).phi())+
        hww.vtxs_position().at(iMax).y()*cos(hww.els_trk_p4().at(index).phi());
    return dxyPV;
}

double electron_d0PV_mindz(unsigned int index){ 
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    double minDz = 999.;
    int iMin = -1;
    for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
      
      const LorentzVector& vtx = hww.gsftrks_vertex_p4().at(hww.els_gsftrkidx().at(index));
      const LorentzVector& p4 = hww.gsftrks_p4().at(hww.els_gsftrkidx().at(index));
      const LorentzVector& pv = hww.vtxs_position().at(i); 
      float dz = fabs((vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt()); 
      
      if ( dz < minDz ){
	iMin = i;
	minDz = dz;
      }
    }
    if (iMin<0) return 9999.;
    double dxyPV = hww.els_d0().at(index)-
        hww.vtxs_position().at(iMin).x()*sin(hww.els_trk_p4().at(index).phi())+
        hww.vtxs_position().at(iMin).y()*cos(hww.els_trk_p4().at(index).phi());
    return dxyPV;
}


// this is now redundant
float electronIsoValuePF2012_FastJetEffArea_HWW(int index){

    const float etaAbs = fabs(hww.els_etaSC().at(index));
    const float pt     = hww.els_p4().at(index).pt();

    // get effective area
    const float AEff = fastJetEffArea04_v1(etaAbs);

    // pf iso
    // calculate from the ntuple for now...
    const float pfiso_ch = hww.els_iso04_pf2012_ch().at(index);
    const float pfiso_em = hww.els_iso04_pf2012_em().at(index);
    const float pfiso_nh = hww.els_iso04_pf2012_nh().at(index);

    // rho
    const float rhoPrime = std::max(hww.evt_ww_rho(), 0.0f);
    const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);  
    const float pfiso = (pfiso_ch + pfiso_n) / pt;   

	// debug
	if(0) {
	cout << "AEff : " << AEff << " "
		 << "rho : " << rhoPrime << " "
		 << "pfiso_ch : " << pfiso_ch << " "
		 << "pfiso_em : " << pfiso_em << " "
		 << "pfiso_nh : " << pfiso_nh << " "
		 << "pfiso_n : " << pfiso_n << " "
		 << "pfiso : " << pfiso << " "
		 << "pt : " << pt << " "
		 << "etaAbs : " << etaAbs 
		 << endl; 
	}

    return pfiso;
}
