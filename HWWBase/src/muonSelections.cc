#include <iostream>

#include "Math/VectorUtil.h"

#include "HWWValidation/HWWBase/interface/muonSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

using namespace std;

namespace HWWFunctions {

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
          if (HWWVal::mus_p4().at(index).pt()<20) 
              return muonIsoValue(index,false) < 0.1;
          else
              return muonIsoValue(index,false) < 0.15;
          break;
      case NominalSmurfV4:
          if (!muonIdNotIsolated( index, type )) return false;
          if (HWWVal::mus_p4().at(index).pt()>20) {
              if (TMath::Abs(HWWVal::mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0) < 0.22;
              else return muonIsoValuePF(index,0) < 0.20;
          } else {
              return muonIsoValuePF(index,0) < 0.11;
          }
          break;
      case NominalSmurfV5:
      case NominalSmurfV6:
          if (!muonIdNotIsolated( index, type )) return false;
          if (HWWVal::mus_p4().at(index).pt()>20) {
              if (TMath::Abs(HWWVal::mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.13;
              else return muonIsoValuePF(index,0,0.3) < 0.09;
          } else {
              if (TMath::Abs(HWWVal::mus_p4().at(index).eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.06;
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


  ////////////////////
  // Identification //
  ////////////////////

  bool muonIdNotIsolated(unsigned int index, SelectionType type) {

      if ( HWWVal::mus_p4().at(index).pt() < 5.0) {
          return false;
      }

      switch (type) {

      case NominalWWV0:
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
          return true;
          break;

      case muonSelectionFO_mu_wwV1:
      case muonSelectionFO_mu_wwV1_iso10:
      case NominalWWV1:
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mud0PV_wwV1(index)) >= 0.02)         return false; // d0 from pvtx
          if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
          if (HWWVal::mus_ptErr().at(index)/HWWVal::mus_p4().at(index).pt()>0.1) return false;
          if (HWWVal::trks_valid_pixelhits().at(HWWVal::mus_trkidx().at(index))==0) return false;
          if (HWWVal::mus_nmatches().at(index)<2) return false;
          return true;
          break;

      case muonSelectionFO_mu_wwV1_iso10_d0: // same as muonSelectionFO_mu_wwV1_iso10 but with looser d0 cut
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mud0PV_wwV1(index)) >= 0.2)         return false; // d0 from pvtx
          if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
          if (HWWVal::mus_ptErr().at(index)/HWWVal::mus_p4().at(index).pt()>0.1) return false;
          if (HWWVal::trks_valid_pixelhits().at(HWWVal::mus_trkidx().at(index))==0) return false;
          if (HWWVal::mus_nmatches().at(index)<2) return false;
          return true;
          break;

      case muonSelectionFO_mu_ww:
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
          return true;

      case muonSelectionFO_mu_ww_iso10:
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
          return true;

      case NominalSmurfV3:
      case NominalSmurfV4:
      case NominalSmurfV5:
          if (type == NominalSmurfV3 || type == NominalSmurfV4 || type == NominalSmurfV5){
              if (HWWVal::mus_p4().at(index).pt()<20){
                  if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
              } else {
                  if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
              }
          } else {
              if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
          }
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
          if (((HWWVal::mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
          if (((HWWVal::mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (HWWVal::mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
          if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
          if (HWWVal::mus_ptErr().at(index)/HWWVal::mus_p4().at(index).pt()>0.1) return false;
          if (HWWVal::trks_valid_pixelhits().at(HWWVal::mus_trkidx().at(index))==0) return false;
          if (HWWVal::mus_nmatches().at(index)<2) return false;
          return true;
          break;


      case muonSelectionFO_mu_smurf_04:
      case muonSelectionFO_mu_smurf_10:
      case NominalSmurfV6:
      {
          if (type == NominalSmurfV6){
              if (HWWVal::mus_p4().at(index).pt()<20){
                  if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
              } else {
                  if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
              }
          } else {
              if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
          }
          if ( TMath::Abs(HWWVal::mus_p4().at(index).eta()) > 2.4)  return false; // eta cut
          if (HWWVal::mus_validHits().at(index) < 11)            return false; // # of tracker hits  
          if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
          if (HWWVal::mus_ptErr().at(index)/HWWVal::mus_p4().at(index).pt()>0.1) return false;
          if (HWWVal::trks_valid_pixelhits().at(HWWVal::mus_trkidx().at(index))==0) return false;
          bool goodMuonGlobalMuon = false;
          if (((HWWVal::mus_type().at(index)) & (1<<1)) != 0) { // global muon
              goodMuonGlobalMuon = true;
              if (HWWVal::mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
              if (HWWVal::mus_gfit_chi2().at(index)/HWWVal::mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
              if (HWWVal::mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false; // Glb fit must have hits in mu chambers
          } 
          bool goodMuonTrackerMuon = false;
          if (((HWWVal::mus_type().at(index)) & (1<<2)) != 0) { // tracker muon
              goodMuonTrackerMuon = true;
              if (HWWVal::mus_pid_TMLastStationTight().at(index) == 0 ) goodMuonTrackerMuon = false; // last station tight
          }
          return goodMuonGlobalMuon || goodMuonTrackerMuon;
          break;
      }
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
  double muonIsoValue_TRK(unsigned int index, bool truncated ){
      double pt        = HWWVal::mus_p4().at(index).pt();
      if(truncated) pt = max( pt, 20.0 );
      return HWWVal::mus_iso03_sumPt().at(index) / pt;
  }
  double muonIsoValue_ECAL(unsigned int index, bool truncated ){
      double pt  = HWWVal::mus_p4().at(index).pt();
      if(truncated) pt = max( pt, 20.0 );
      return HWWVal::mus_iso03_emEt().at(index) / pt;
  }
  double muonIsoValue_HCAL(unsigned int index, bool truncated){
      double pt  = HWWVal::mus_p4().at(index).pt();
      if(truncated) pt = max( pt, 20.0 );
      return HWWVal::mus_iso03_hadEt().at(index) / pt;
  }

  #ifdef PFISOFROMNTUPLE
  double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
      if (fabs(coner-0.3)<0.0001) {
          if (HWWVal::mus_iso03_pf().at(imu)<-99.) return 9999.;
          return HWWVal::mus_iso03_pf().at(imu)/HWWVal::mus_p4().at(imu).pt();
      } else if (fabs(coner-0.4)<0.0001) {
          if (HWWVal::mus_iso04_pf().at(imu)<-99.) return 9999.;
          return HWWVal::mus_iso04_pf().at(imu)/HWWVal::mus_p4().at(imu).pt();
      } else {
          std::cout << "muonIsoValuePF: CONE SIZE NOT SUPPORTED" << std::endl;
          return 9999.;
      }
  }
  #else
  double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
      float pfciso = 0;
      float pfniso = 0;
      int mutkid = HWWVal::mus_trkidx().at(imu);
      float mudz = mutkid>=0 ? trks_dz_pv(mutkid,ivtx).first : HWWVal::mus_sta_z0corr().at(imu);
      for (unsigned int ipf=0; ipf<HWWVal::pfcands_p4().size(); ++ipf){
          float dR = ROOT::Math::VectorUtil::DeltaR( HWWVal::pfcands_p4().at(ipf), HWWVal::mus_p4().at(imu) );
          if (dR>coner) continue;
          float pfpt = HWWVal::pfcands_p4().at(ipf).pt();
          int pfid = abs(HWWVal::pfcands_particleId().at(ipf));
          if (filterId!=0 && filterId!=pfid) continue;
          if (HWWVal::pfcands_charge().at(ipf)==0) {
              //neutrals
              if (pfpt>minptn) pfniso+=pfpt;
          } else {
              //charged
              //avoid double counting of muon itself
              int pftkid = HWWVal::pfcands_trkidx().at(ipf);
              if (mutkid>=0 && pftkid>=0 && mutkid==pftkid) continue;
              //first check electrons with gsf track
              if (abs(HWWVal::pfcands_particleId().at(ipf))==11 && HWWVal::pfcands_pfelsidx().at(ipf)>=0 && HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))>=0) {
                  int gsfid = HWWVal::els_gsftrkidx().at(HWWVal::pfels_elsidx().at(HWWVal::pfcands_pfelsidx().at(ipf))); 
                  if (gsfid>=0) { 
                      if(fabs(gsftrks_dz_pv( gsfid,ivtx ).first - mudz )<dzcut) {//dz cut
                          pfciso+=pfpt;
                      }   
                      continue;//and avoid double counting
                  }
              }
              //then check anything that has a ctf track
              if (HWWVal::pfcands_trkidx().at(ipf)>=0) {//charged (with a ctf track)
                  if(fabs( trks_dz_pv(HWWVal::pfcands_trkidx().at(ipf),ivtx).first - mudz )<dzcut) {//dz cut
                      pfciso+=pfpt;
                  }
              } 
          }
      } 
      return (pfciso+pfniso)/HWWVal::mus_p4().at(imu).pt();
  }
  #endif


  /////////////////////////////
  // Muon d0 corrected by PV //
  ////////////////////////////

  double mud0PV(unsigned int index){
      if ( HWWVal::vtxs_sumpt().empty() ) return 9999.;
      unsigned int iMax = 0;
      double sumPtMax = HWWVal::vtxs_sumpt().at(0);
      for ( unsigned int i = iMax+1; i < HWWVal::vtxs_sumpt().size(); ++i )
          if ( HWWVal::vtxs_sumpt().at(i) > sumPtMax ){
              iMax = i;
              sumPtMax = HWWVal::vtxs_sumpt().at(i);
          }
      double dxyPV = HWWVal::mus_d0().at(index)-
          HWWVal::vtxs_position().at(iMax).x()*sin(HWWVal::mus_trk_p4().at(index).phi())+
          HWWVal::vtxs_position().at(iMax).y()*cos(HWWVal::mus_trk_p4().at(index).phi());
      return dxyPV;
  }

  double mud0PV_wwV1(unsigned int index){
      if ( HWWVal::vtxs_sumpt().empty() ) return 9999.;
      double sumPtMax = -1;
      int iMax = -1;
      for ( unsigned int i = 0; i < HWWVal::vtxs_sumpt().size(); ++i ){
          // if (!isGoodVertex(i)) continue;
          // Copied from eventSelections.cc 
          if (HWWVal::vtxs_isFake().at(i)) continue;
          if (HWWVal::vtxs_ndof().at(i) < 4.) continue;
          if (HWWVal::vtxs_position().at(i).Rho() > 2.0) continue;
          if (fabs(HWWVal::vtxs_position().at(i).Z()) > 24.0) continue;
          if ( HWWVal::vtxs_sumpt().at(i) > sumPtMax ){
              iMax = i;
              sumPtMax = HWWVal::vtxs_sumpt().at(i);
          }
      }
      if (iMax<0) return 9999.;
      double dxyPV = HWWVal::mus_d0().at(index)-
          HWWVal::vtxs_position().at(iMax).x()*sin(HWWVal::mus_trk_p4().at(index).phi())+
          HWWVal::vtxs_position().at(iMax).y()*cos(HWWVal::mus_trk_p4().at(index).phi());
      return dxyPV;
  }

  double mud0PV_smurfV3(unsigned int index){
      int vtxIndex = 0;
      double dxyPV = HWWVal::mus_d0().at(index)-
          HWWVal::vtxs_position().at(vtxIndex).x()*sin(HWWVal::mus_trk_p4().at(index).phi())+
          HWWVal::vtxs_position().at(vtxIndex).y()*cos(HWWVal::mus_trk_p4().at(index).phi());
      return dxyPV;
  }

  double dzPV_mu(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
      return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
  }

  double mudzPV_smurfV3(unsigned int index){
      int vtxIndex = 0;
      double dzpv = dzPV_mu(HWWVal::mus_vertex_p4().at(index), HWWVal::mus_trk_p4().at(index), HWWVal::vtxs_position().at(vtxIndex));
      return dzpv;
  }

  double mudzPV_wwV1(unsigned int index){
      if ( HWWVal::vtxs_sumpt().empty() ) return 9999.;
      double sumPtMax = -1;
      int iMax = -1;
      for ( unsigned int i = 0; i < HWWVal::vtxs_sumpt().size(); ++i ){
          if (HWWVal::vtxs_isFake().at(i)) continue;
          if (HWWVal::vtxs_ndof().at(i) < 4.) continue;
          if (HWWVal::vtxs_position().at(i).Rho() > 2.0) continue;
          if (fabs(HWWVal::vtxs_position().at(i).Z()) > 24.0) continue;
          if ( HWWVal::vtxs_sumpt().at(i) > sumPtMax ){
              iMax = i;
              sumPtMax = HWWVal::vtxs_sumpt().at(i);
          }
      }
      if (iMax<0) return 9999.;
      const LorentzVector& vtx = HWWVal::mus_vertex_p4().at(index);
      const LorentzVector& p4 = HWWVal::mus_trk_p4().at(index);
      const LorentzVector& pv = HWWVal::vtxs_position().at(iMax);
      return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
  }

  bool isPFMuon( int index , bool requireSamePt , float dpt_max ){

      int ipf = HWWVal::mus_pfmusidx().at( index );

      //--------------------------
      // require matched pfmuon
      //--------------------------

      if( ipf >= int(HWWVal::pfmus_p4().size()) || ipf < 0 ) return false;

      //----------------------------------------------------
      // require PFMuon pt = reco muon pt (within dpt_max)
      //----------------------------------------------------

      if( requireSamePt ){

          float pt_pf = HWWVal::pfmus_p4().at(ipf).pt();
          float pt    = HWWVal::mus_p4().at(index).pt();

          if( fabs( pt_pf - pt ) > dpt_max ) return false;

      }

      return true;

  }

}
