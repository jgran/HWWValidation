#include <iostream>

#include "Math/VectorUtil.h"

#include "HWWValidation/HWWBase/interface/muonSelections.h"
#include "HWWValidation/HWWBase/interface/eventSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

using namespace std;

////////////////////
// Identification //
////////////////////

bool muonId(unsigned int index, SelectionType type){

    float isovalue;
    bool  truncated = true;

    switch(type) {

        ///////////////
        // Higgs, WW //
        ///////////////

        // WW
    case NominalWWV0:
    case NominalWWV1:
        isovalue = 0.15;
        break;
    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_ww:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_smurf_04:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0,0.3) < 0.40;
        break;
    case muonSelectionFO_mu_wwV1_iso10_d0:
    case muonSelectionFO_mu_wwV1_iso10:
    case muonSelectionFO_mu_ww_iso10:
        isovalue = 1.0;
        break;

        // SMURF
    case muonSelectionFO_mu_smurf_10:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0,0.3) < 1.0;
        break;
    case NominalSmurfV3:
        if (!muonIdNotIsolated( index, type )) return false;
        if (hww.mus_p4().at(index).pt()<20) 
            return muonIsoValue(index,false) < 0.1;
        else
            return muonIsoValue(index,false) < 0.15;
        break;
    case NominalSmurfV4:
        if (!muonIdNotIsolated( index, type )) return false;
        if (hww.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(hww.mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0) < 0.22;
            else return muonIsoValuePF(index,0) < 0.20;
        } else {
            return muonIsoValuePF(index,0) < 0.11;
        }
        break;
    case NominalSmurfV5:
    case NominalSmurfV6:
        if (!muonIdNotIsolated( index, type )) return false;
        if (hww.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(hww.mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.13;
            else return muonIsoValuePF(index,0,0.3) < 0.09;
        } else {
            if (TMath::Abs(hww.mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.06;
            else return muonIsoValuePF(index,0,0.3) < 0.05;
        }
        break;



        /////////////
        // Default //
        /////////////
    default:
        std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
        exit(1);
        return false;
    } 
    return 
        muonIdNotIsolated( index, type ) &&   // Id
        muonIsoValue(index,truncated) < isovalue;           // Isolation cut
}


bool isGoodStandardMuon( unsigned int index ){
    if ( TMath::Abs( hww.mus_p4().at(index).eta() ) > 2.4 )              return false;
    if ( hww.mus_gfit_chi2().at(index) / hww.mus_gfit_ndof().at(index) >= 50 )  return false;
    if ( ( ( hww.mus_type().at(index) ) & (1<<1) ) == 0 )                return false;
    if ( ( ( hww.mus_type().at(index) ) & (1<<2) ) == 0 )                return false;
    if ( hww.mus_validHits().at(index) < 11 )                            return false;
    if ( hww.mus_gfit_validSTAHits().at(index) == 0)                     return false;
    return true;
}

////////////////////
// Identification //
////////////////////

bool muonIdNotIsolated(unsigned int index, SelectionType type) {

    if ( hww.mus_p4().at(index).pt() < 5.0) {
        // std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
        return false;
    }

    int vtxidx = firstGoodVertex();
    int trkidx = hww.mus_trkidx().at(index);
    
    // Muon Selections that are standard for Analysis & Fake Selections
    bool standardMuon = true;
    if ( TMath::Abs( hww.mus_p4().at(index).eta() ) > 2.4 )              standardMuon = false;
    if ( hww.mus_gfit_chi2().at(index) / hww.mus_gfit_ndof().at(index) >= 50 )  standardMuon = false;
    if ( ( ( hww.mus_type().at(index) ) & (1<<1) ) == 0 )                standardMuon = false;
    if ( ( ( hww.mus_type().at(index) ) & (1<<2) ) == 0 )                standardMuon = false;
    if ( hww.mus_validHits().at(index) < 11 )                            standardMuon = false;
    if ( hww.mus_gfit_validSTAHits().at(index) == 0)                     standardMuon = false;
    if ( hww.mus_ptErr().at(index) / hww.mus_p4().at(index).pt() > 0.1 )        standardMuon = false;


    //
    switch (type) {

    case NominalWWV0:
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;
        break;

    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_wwV1_iso10:
    case NominalWWV1:
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.02)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1) return false;
        if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(index))==0) return false;
        if (hww.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_wwV1_iso10_d0: // same as muonSelectionFO_mu_wwV1_iso10 but with looser d0 cut
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.2)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1) return false;
        if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(index))==0) return false;
        if (hww.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_ww:
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;

    case muonSelectionFO_mu_ww_iso10:
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;

    case OSGeneric_v4:
        //baseline selector for 2011 OS analysis
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSGeneric_v3:
        //baseline selector for 2011 OS analysis
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSGeneric_v3_FO:
	    // Fakes for 2011: reliso < 0.4, d0 < 0.2, chisq/ndof < 50
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.2)                             return false; // d0(PV) < 0.2 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSZ_v4:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)                       return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case OSZ_v3:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)                       return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(hww.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case OSZ_v2:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)                       return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(hww.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case NominalSmurfV3:
    case NominalSmurfV4:
    case NominalSmurfV5:
        if (type == NominalSmurfV3 || type == NominalSmurfV4 || type == NominalSmurfV5){
            if (hww.mus_p4().at(index).pt()<20){
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
            } else {
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
            }
        } else {
            if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
        }
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (hww.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1) return false;
        if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(index))==0) return false;
        if (hww.mus_nmatches().at(index)<2) return false;
        return true;
        break;

        //baseline selector for 2011 SS analysis
    case NominalSSv4:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || hww.mus_trkidx().at(index) < 0) {
            if (fabs(hww.mus_d0corr().at(index)) > 0.02)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(hww.mus_trkidx().at(index), vtxidx).first) > 0.02)
                return false;            
        }
        return true;
        break;
        //baseline selector for 2011 SS analysis
    case muonSelectionFO_ssV4:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((hww.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (hww.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || hww.mus_trkidx().at(index) < 0) {
            if (fabs(hww.mus_d0corr().at(index)) > 0.2)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(hww.mus_trkidx().at(index), vtxidx).first) > 0.2)
                return false;            
        }
        return true;
        break;

        // muon POG tight muon requirements, with d0 cut tightened to 0.02 cm
	// see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId

    case ZMet2012_v1:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (hww.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
	if (hww.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (trkidx < 0)                                                          return false; // require a matching track
        if (vtxidx < 0 || trkidx < 0) {
	        std::cout << __FILE__ << " " << __LINE__ << std::endl;
	        std::cout << "WARNING: didn't find any good vertices, should never get here" << std::endl;

            if (fabs(hww.mus_d0corr().at(index)) > 0.02) return false;
            if (fabs(hww.mus_z0corr().at(index)) > 0.5)  return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.5)  return false;
        }
        else return false;     

        if (hww.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit
        if (hww.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits

        return true;
        break;

        //baseline selector for 2012 SS analysis
    case NominalSSv5:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (hww.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (hww.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        if (trkidx < 0)                                                          return false; // require a matching track
        if (hww.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits
        if (hww.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit

        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (hww.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || trkidx < 0) {
            if (fabs(hww.mus_d0corr().at(index)) > 0.02)
                return false;
            if (fabs(hww.mus_z0corr().at(index)) > 0.1)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02)
                return false;
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.1)
                return false;
        }
        else return false;     
        return true;
        break;
        //baseline FO selector for 2012 SS analysis
    case muonSelectionFO_ssV5:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((hww.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (hww.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (hww.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        if (trkidx < 0)                                                          return false; // require a matching track
        if (hww.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits
        if (hww.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit

        if (hww.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || trkidx < 0) {
            if (fabs(hww.mus_d0corr().at(index)) > 0.2)
                return false;
            if (fabs(hww.mus_z0corr().at(index)) > 0.1)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(hww.mus_trkidx().at(index), vtxidx).first) > 0.2)
                return false;            
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.1)
                return false;
        }
        else return false;
        return true;
        break;

    case muonSelectionFO_mu_smurf_04:
    case muonSelectionFO_mu_smurf_10:
    case NominalSmurfV6:
    {
        if (type == NominalSmurfV6){
            if (hww.mus_p4().at(index).pt()<20){
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
            } else {
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
            }
        } else {
            if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
        }
        if ( TMath::Abs(hww.mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
        if (hww.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt()>0.1) return false;
        if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(index))==0) return false;
        bool goodMuonGlobalMuon = false;
        if (((hww.mus_type().at(index)) & (1<<1)) != 0) { // global muon
            goodMuonGlobalMuon = true;
            if (hww.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
            if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
            if (hww.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false; // Glb fit must have hits in mu chambers
        } 
        bool goodMuonTrackerMuon = false;
        if (((hww.mus_type().at(index)) & (1<<2)) != 0) { // tracker muon
            goodMuonTrackerMuon = true;
            if (hww.mus_pid_TMLastStationTight().at(index) == 0 ) goodMuonTrackerMuon = false; // last station tight
        }
        return goodMuonGlobalMuon || goodMuonTrackerMuon;
        break;
    }
//jgran
/*
    case NominalTTZ_loose_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
		if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;
    case NominalTTZ_tight_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
        if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        if (hww.mus_iso_ecalvetoDep().at(index) > 4) return false; // ECalE < 4 
        if (hww.mus_iso_hcalvetoDep().at(index) > 6) return false; // HCalE < 6
        return true;
        break;
    case NominalTTZ_looseFO_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 50) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.20) return false;
		if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;
    case NominalTTZ_tightFO_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 50) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.20) return false;
        if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;        
        //baseline selector for 2012 OS RPV STOP analysis
    case NominalOSv1:
        if (fabs(hww.mus_p4().at(index).eta()) > 2.4)                           return false;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) return false;
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT))                 return false;
        if (hww.mus_ptErr().at(index)/hww.mus_p4().at(index).pt() > 0.1)       return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02)                       return false;
        if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.1)                        return false;
        return true;
        break;
*/
    default:
        std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
        return false;
    }
}




////////////////////////////
// Isolation Calculations //
////////////////////////////

double muonIsoValue(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) );
}
double muonIsoValue_FastJet(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + TMath::Max( muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) - mu_fastjet_rel_offset(index,truncated) , 0.0 ) );
}
double mu_fastjet_rel_offset(unsigned int index, bool truncated ){
    double pt        = hww.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    double offset = TMath::Pi() * pow( 0.3 , 2 ) * (hww.evt_rho());
    return offset / pt;
}
double muonIsoValue_TRK(unsigned int index, bool truncated ){
    double pt        = hww.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return hww.mus_iso03_sumPt().at(index) / pt;
}
double muonIsoValue_ECAL(unsigned int index, bool truncated ){
    double pt  = hww.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return hww.mus_iso03_emEt().at(index) / pt;
}
double muonIsoValue_HCAL(unsigned int index, bool truncated){
    double pt  = hww.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return hww.mus_iso03_hadEt().at(index) / pt;
}
double muonCorIsoValue (unsigned int index, bool truncated) {
    double ntiso  = muonIsoValue(index, truncated);
    double pt     = hww.mus_p4().at(index).pt();
    int nvtxs     = numberOfGoodVertices();
    double coriso = ntiso - ((TMath::Log(pt)*nvtxs)/(30*pt));
    return coriso;
}

#ifdef PFISOFROMNTUPLE
double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
    if (fabs(coner-0.3)<0.0001) {
        if (hww.mus_iso03_pf().at(imu)<-99.) return 9999.;
        return hww.mus_iso03_pf().at(imu)/hww.mus_p4().at(imu).pt();
    } else if (fabs(coner-0.4)<0.0001) {
        if (hww.mus_iso04_pf().at(imu)<-99.) return 9999.;
        return hww.mus_iso04_pf().at(imu)/hww.mus_p4().at(imu).pt();
    } else {
        std::cout << "muonIsoValuePF: CONE SIZE NOT SUPPORTED" << std::endl;
        return 9999.;
    }
}
#else
double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
    float pfciso = 0;
    float pfniso = 0;
    int mutkid = hww.mus_trkidx().at(imu);
    float mudz = mutkid>=0 ? trks_dz_pv(mutkid,ivtx).first : hww.mus_sta_z0corr().at(imu);
    for (unsigned int ipf=0; ipf<hww.pfcands_p4().size(); ++ipf){
        float dR = ROOT::Math::VectorUtil::DeltaR( hww.pfcands_p4().at(ipf), hww.mus_p4().at(imu) );
        if (dR>coner) continue;
        float pfpt = hww.pfcands_p4().at(ipf).pt();
        int pfid = abs(hww.pfcands_particleId().at(ipf));
        if (filterId!=0 && filterId!=pfid) continue;
        if (hww.pfcands_charge().at(ipf)==0) {
            //neutrals
            if (pfpt>minptn) pfniso+=pfpt;
        } else {
            //charged
            //avoid double counting of muon itself
            int pftkid = hww.pfcands_trkidx().at(ipf);
            if (mutkid>=0 && pftkid>=0 && mutkid==pftkid) continue;
            //first check electrons with gsf track
            if (abs(hww.pfcands_particleId().at(ipf))==11 && hww.pfcands_pfelsidx().at(ipf)>=0 && hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))>=0) {
                int gsfid = hww.els_gsftrkidx().at(hww.pfels_elsidx().at(hww.pfcands_pfelsidx().at(ipf))); 
                if (gsfid>=0) { 
                    if(fabs(gsftrks_dz_pv( gsfid,ivtx ).first - mudz )<dzcut) {//dz cut
                        pfciso+=pfpt;
                    }   
                    continue;//and avoid double counting
                }
            }
            //then check anything that has a ctf track
            if (hww.pfcands_trkidx().at(ipf)>=0) {//charged (with a ctf track)
                if(fabs( trks_dz_pv(hww.pfcands_trkidx().at(ipf),ivtx).first - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 
    return (pfciso+pfniso)/hww.mus_p4().at(imu).pt();
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove cosmics by looking for back-to-back muon-track pairs ( http://indico.cern.ch/contributionDisplay.py?contribId=2&confId=86834 ) //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isCosmics(unsigned int index){
    for (int itrk=0; itrk < int(hww.trks_trk_p4().size()); ++itrk) {
        const LorentzVector& mu_p4  = hww.mus_trk_p4().at(index);
        const LorentzVector& trk_p4 = hww.trks_trk_p4().at(itrk);
        double sprod = mu_p4.px()*trk_p4.px()+mu_p4.py()*trk_p4.py()+mu_p4.pz()*trk_p4.pz();
        if ( acos( -(sprod/trk_p4.P()/mu_p4.P()) ) < 0.01 &&
             fabs(trk_p4.pt()-mu_p4.pt())/mu_p4.pt() < 0.05 )
            return true;
    }
    return false;
}

/////////////////////////////
// Muon d0 corrected by PV //
////////////////////////////

double mud0PV(unsigned int index){
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    unsigned int iMax = 0;
    double sumPtMax = hww.vtxs_sumpt().at(0);
    for ( unsigned int i = iMax+1; i < hww.vtxs_sumpt().size(); ++i )
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    double dxyPV = hww.mus_d0().at(index)-
        hww.vtxs_position().at(iMax).x()*sin(hww.mus_trk_p4().at(index).phi())+
        hww.vtxs_position().at(iMax).y()*cos(hww.mus_trk_p4().at(index).phi());
    return dxyPV;
}

double mud0PV_wwV1(unsigned int index){
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (hww.vtxs_isFake().at(i)) continue;
        if (hww.vtxs_ndof().at(i) < 4.) continue;
        if (hww.vtxs_position().at(i).Rho() > 2.0) continue;
        if (fabs(hww.vtxs_position().at(i).Z()) > 24.0) continue;
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    double dxyPV = hww.mus_d0().at(index)-
        hww.vtxs_position().at(iMax).x()*sin(hww.mus_trk_p4().at(index).phi())+
        hww.vtxs_position().at(iMax).y()*cos(hww.mus_trk_p4().at(index).phi());
    return dxyPV;
}

double mud0PV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dxyPV = hww.mus_d0().at(index)-
        hww.vtxs_position().at(vtxIndex).x()*sin(hww.mus_trk_p4().at(index).phi())+
        hww.vtxs_position().at(vtxIndex).y()*cos(hww.mus_trk_p4().at(index).phi());
    return dxyPV;
}

double dzPV_mu(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

double mudzPV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dzpv = dzPV_mu(hww.mus_vertex_p4().at(index), hww.mus_trk_p4().at(index), hww.vtxs_position().at(vtxIndex));
    return dzpv;
}

double mudzPV_wwV1(unsigned int index){
    if ( hww.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (hww.vtxs_isFake().at(i)) continue;
        if (hww.vtxs_ndof().at(i) < 4.) continue;
        if (hww.vtxs_position().at(i).Rho() > 2.0) continue;
        if (fabs(hww.vtxs_position().at(i).Z()) > 24.0) continue;
        if ( hww.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = hww.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    // double dzpv = hww.mus_z0corr().at(index)-hww.vtxs_position()[iMax].z();
    const LorentzVector& vtx = hww.mus_vertex_p4().at(index);
    const LorentzVector& p4 = hww.mus_trk_p4().at(index);
    const LorentzVector& pv = hww.vtxs_position().at(iMax);
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
}

bool isPFMuon( int index , bool requireSamePt , float dpt_max ){

    int ipf = hww.mus_pfmusidx().at( index );

    //--------------------------
    // require matched pfmuon
    //--------------------------

    if( ipf >= int(hww.pfmus_p4().size()) || ipf < 0 ) return false;

    //----------------------------------------------------
    // require PFMuon pt = reco muon pt (within dpt_max)
    //----------------------------------------------------

    if( requireSamePt ){

        float pt_pf = hww.pfmus_p4().at(ipf).pt();
        float pt    = hww.mus_p4().at(index).pt();

        if( fabs( pt_pf - pt ) > dpt_max ) return false;

    }

    return true;

}

void muonIsoValuePF2012 (float &pfiso_ch, float &pfiso_em, float &pfiso_nh, const float R, const unsigned int imu, const int ivtx, float neutral_et_threshold)
{

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0; 
    pfiso_nh = 0.0;
       
    // loop on pfcandidates
    for (unsigned int ipf = 0; ipf < hww.pfcands_p4().size(); ++ipf) {
            
        // skip electrons and muons
        const int particleId = abs(hww.pfcands_particleId().at(ipf));
        if (particleId == 11)    continue;
        if (particleId == 13)    continue;
    
        // deltaR between electron and cadidate
        const float dR = ROOT::Math::VectorUtil::DeltaR(hww.pfcands_p4().at(ipf), hww.mus_p4().at(imu));
        if (dR > R)              continue;

        // charged hadrons closest vertex
        // should be the primary vertex
        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211) {
            if (hww.pfcands_vtxidx().at(ipf) != ivtx) continue;
            if (dR < 0.0001)
                continue;
        }
        if (particleId == 22 || particleId == 130 || particleId == 111 || particleId == 310 || particleId == 2112) {
            if (hww.pfcands_p4().at(ipf).pt() < neutral_et_threshold)
                continue;
            if (dR < 0.01)
                continue;
        }

        // add to isolation sum
        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211)      pfiso_ch += hww.pfcands_p4().at(ipf).pt();
        if (particleId == 22)                                                                          pfiso_em += hww.pfcands_p4().at(ipf).pt();
        if (particleId == 130 || particleId == 111 || particleId == 310 || particleId == 2112)         pfiso_nh += hww.pfcands_p4().at(ipf).pt();
    }
}

float muonIsoValuePF2012_FastJetEffArea(int index, float conesize, float effective_area, int ivtx)
{
    float pt     = hww.mus_p4().at(index).pt();

    // pf iso
    // calculate from the ntuple for now...
    float pfiso_ch = 0.0;
    float pfiso_em = 0.0;
    float pfiso_nh = 0.0;
    muonIsoValuePF2012(pfiso_ch, pfiso_em, pfiso_nh, conesize, index, ivtx);

    // rho
    float rhoPrime = std::max(hww.evt_rho(), float(0.0));
    float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * effective_area, float(0.0));
    float pfiso = (pfiso_ch + pfiso_n) / pt;   

    return pfiso;    
}

float muonRadialIsolation (unsigned int imu, float &chiso, float &nhiso, float &emiso, float neutral_et_threshold, float cone_size, bool verbose)
{
    float radial_iso = 0.;
    chiso = 0.;
    nhiso = 0.;
    emiso = 0.;

    int ivtx = firstGoodVertex();

    LorentzVector p4 = (hww.mus_trk_p4().at(imu).pt() > 0.01) ? hww.mus_trk_p4().at(imu) : hww.mus_sta_p4().at(imu);
    if (p4.pt() < 0.01)
        return -9999.;

    for (unsigned int ipf = 0; ipf < hww.pfcands_p4().size(); ipf++) {

        // skip electrons and muons
        const int particleId = abs(hww.pfcands_particleId().at(ipf));
        if (particleId == 11) {
            if (verbose)
                std::cout << "Skipping electron with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }
        if (particleId == 13) {
            if (verbose)
                std::cout << "Skipping muon with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }

        // in the event that the muon is not a PF muon, need to remove any other PF cand reconstructed using the same track as the muon
        if (!hww.mus_pid_PFMuon().at(imu) && hww.mus_trkidx().at(imu) >= 0 && hww.mus_trkidx().at(imu) == hww.pfcands_trkidx().at(ipf)) {
            if (verbose)
                std::cout << "Skipping PF cand with same track as muon with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }

        const float dr = ROOT::Math::VectorUtil::DeltaR(hww.pfcands_p4().at(ipf), hww.mus_p4().at(imu));
        if (dr > cone_size) {
            if (verbose)
                std::cout << "Skipping PF candidate outside of cone with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;            
            continue;
        }
        if (dr < 0.01) {
            if (verbose)
                std::cout << "Skipping PF candidate in veto cone with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;            
            continue;
        }

        // deal with charged
        if (hww.pfcands_charge().at(ipf) != 0) {
            if (hww.pfcands_vtxidx().at(ipf) != ivtx) {
                if (verbose)
                    std::cout << "Skipping PF candidate from other vertex  with id, pt, eta, ivtx = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " 
                              << hww.pfcands_p4().at(ipf).eta() << ", " << hww.pfcands_vtxidx().at(ipf) << std::endl;
                continue;
            }
            radial_iso += hww.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / hww.mus_p4().at(imu).pt();
            chiso += hww.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / hww.mus_p4().at(imu).pt();
            if (verbose)
                std::cout << "Summing CH with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;            
        }
        else if (hww.pfcands_p4().at(ipf).pt() > neutral_et_threshold) {
            radial_iso += hww.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / hww.mus_p4().at(imu).pt();
            if (particleId == 22) {
                emiso += hww.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / hww.mus_p4().at(imu).pt();
                if (verbose)
                    std::cout << "Summing EM with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;            
            }
            else {
                nhiso += hww.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / hww.mus_p4().at(imu).pt();
                if (verbose)
                    std::cout << "Summing NH with id, pt, eta = " << hww.pfcands_particleId().at(ipf) << ", " << hww.pfcands_p4().at(ipf).pt() << ", " << hww.pfcands_p4().at(ipf).eta() << std::endl;            
            }
        }
    } // loop over pfcands

    return radial_iso;
}

float muonIsoValuePF2012_deltaBeta(unsigned int imu)
{
    const float chiso = hww.mus_isoR03_pf_ChargedHadronPt().at(imu);
    const float nhiso = hww.mus_isoR03_pf_NeutralHadronEt().at(imu);
    const float emiso = hww.mus_isoR03_pf_PhotonEt().at(imu);
    const float deltaBeta = hww.mus_isoR03_pf_PUPt().at(imu);
    const float pt = hww.mus_p4().at(imu).pt();

    const float absiso = chiso + max(0.0, nhiso + emiso - 0.5 * deltaBeta);
    return (absiso / pt);
}

// muon POG selections:
bool passes_muid_wp2012(const unsigned int index, const mu2012_tightness::value_type tightness)
{

    const bool is_global  = ((hww.mus_type().at(index) & (1<<1)) != 0);
    const bool is_tracker = ((hww.mus_type().at(index) & (1<<2)) != 0);
    const bool is_pfmu    = ((hww.mus_type().at(index) & (1<<5)) != 0);

    switch (tightness) 
    {
        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
        case mu2012_tightness::LOOSE: 
            {
                if (not is_pfmu)              {return false;}
                if (is_global && !is_tracker) {return false;}     

                // if we're here, then it passes
                return true;
            }
            break;

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon
        case mu2012_tightness::TIGHT:
            {
                // event has a good vertex and muon has associated ctf track
                const int ctfidx = hww.mus_trkidx().at(index);
                const int vtxidx = firstGoodVertex();
                if (ctfidx < 0 || vtxidx < 0)  {return false;}

                const float d0       = trks_d0_pv(ctfidx, vtxidx).first;
                const float dz       = trks_dz_pv(ctfidx, vtxidx).first;
                const float chi2ndof = hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index);

                if (not is_global)                                {return false;} // The candidate is reconstructed as a Global Muon
                if (not is_pfmu)                                  {return false;} // Particle-Flow muon id 
                if (chi2ndof >= 10)                               {return false;} // χ2/ndof of the global-muon track fit < 10
                if (hww.mus_numberOfMatchedStations().at(index) <= 1) {return false;} // Muon segments in at least two muon stations
                if (hww.mus_gfit_validSTAHits().at(index) <= 0)       {return false;} // At least one muon chamber hit included in the global-muon track fit 
                if (fabs(d0) > 0.2)                               {return false;} // Its tracker track has transverse impact parameter dxy < 2 mm w.r.t. the primary vertex
                if (fabs(dz) > 0.5)                               {return false;} // The longitudinal distance of the tracker track wrt. the primary vertex is dz < 5 mm
                if (hww.trks_valid_pixelhits().at(ctfidx) <= 0)       {return false;} // Number of pixel hits > 0
                if (hww.trks_nlayers().at(ctfidx) <= 5)               {return false;} // Cut on number of tracker layers with hits > 5

                // if we're here, then it passes
                return true;
            }
            break;

        default: {/*do nothing*/}
    } // end switch block

    // if we're here, then return false
    return false;
}

