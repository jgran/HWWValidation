#include "DataFormats/PatCandidates/interface/Electron.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "HWWValidation/HWWBase/interface/EgammaFiduciality.h"
#include "HWWValidation/HWWBase/interface/ElectronMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"


using namespace reco;
using namespace edm;
using namespace std;

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

ElectronMaker::ElectronMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector) {

  TrackCollection_          = iCollector.consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackInputTag"));
  GSFTrackCollection_       = iCollector.consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("gsftrksInputTag"));
  GSFElectron_              = iCollector.consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronsInputTag"));
  GSFElectronCollection_    = iCollector.consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsInputTag"));
  PFCandidateCollection_    = iCollector.consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandsInputTag"));
  thePVCollection_          = iCollector.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));
  BeamSpot_                 = iCollector.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  ConversionCollection_     = iCollector.consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("recoConversionInputTag"));

}


void ElectronMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_els_fiduciality();
  HWWVal::Load_els_type();
  HWWVal::Load_els_ecalEnergy();
  HWWVal::Load_els_trk_p4();
  HWWVal::Load_els_p4();
  HWWVal::Load_els_vertex_p4();
  HWWVal::Load_els_ecalIso();
  HWWVal::Load_els_hcalIso();
  HWWVal::Load_els_tkIso();
  HWWVal::Load_els_ecalIso04();
  HWWVal::Load_els_hcalIso04();
  HWWVal::Load_els_iso03_pf_ch();
  HWWVal::Load_els_iso03_pf_gamma05();
  HWWVal::Load_els_iso03_pf_nhad05();
  HWWVal::Load_els_iso04_pf_ch();
  HWWVal::Load_els_iso04_pf_gamma05();
  HWWVal::Load_els_iso04_pf_nhad05();
  HWWVal::Load_els_iso03_pf2012_ch();
  HWWVal::Load_els_iso03_pf2012_em();
  HWWVal::Load_els_iso03_pf2012_nh();
  HWWVal::Load_els_iso04_pf2012_ch();
  HWWVal::Load_els_iso04_pf2012_em();
  HWWVal::Load_els_iso04_pf2012_nh();
  HWWVal::Load_els_etaSC();
  HWWVal::Load_els_eSC();
  HWWVal::Load_els_eSCRaw();
  HWWVal::Load_els_eSCPresh();
  HWWVal::Load_els_nSeed();
  HWWVal::Load_els_e1x5();
  HWWVal::Load_els_e5x5();
  HWWVal::Load_els_sigmaIEtaIEta();
  HWWVal::Load_els_etaSCwidth();
  HWWVal::Load_els_phiSCwidth();
  HWWVal::Load_els_sigmaIPhiIPhi();
  HWWVal::Load_els_e3x3();
  HWWVal::Load_els_hOverE();
  HWWVal::Load_els_eOverPIn();
  HWWVal::Load_els_eSeedOverPOut();
  HWWVal::Load_els_eSeedOverPIn();
  HWWVal::Load_els_eOverPOut();
  HWWVal::Load_els_fbrem();
  HWWVal::Load_els_dEtaIn();
  HWWVal::Load_els_dEtaOut();
  HWWVal::Load_els_dPhiIn();
  HWWVal::Load_els_dPhiOut();
  HWWVal::Load_els_chi2();
  HWWVal::Load_els_ndof();
  HWWVal::Load_els_gsftrkidx();
  HWWVal::Load_els_charge();
  HWWVal::Load_els_trk_charge();
  HWWVal::Load_els_sccharge();
  HWWVal::Load_els_d0();
  HWWVal::Load_els_d0corr();
  HWWVal::Load_els_z0corr();
  HWWVal::Load_els_trkidx();
  HWWVal::Load_els_trkshFrac();
  HWWVal::Load_els_ip3d();
  HWWVal::Load_els_ip3derr();
  HWWVal::Load_els_exp_innerlayers();
  HWWVal::Load_els_conv_dist();
  HWWVal::Load_els_conv_dcot();
  HWWVal::Load_els_conv_old_dist();
  HWWVal::Load_els_conv_old_dcot();


  // access the tracker
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

  ////////////////
  // Get Tracks //
  ////////////////
 
  Handle<TrackCollection> tracks_h;
  iEvent.getByToken(TrackCollection_, tracks_h);


  ////////////////
  // GSF Tracks //
  ////////////////

  Handle<GsfTrackCollection> gsftracks_h;
  iEvent.getByToken(GSFTrackCollection_, gsftracks_h);


  /////////////
  // B Field //
  /////////////

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  float evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
 

  ///////////////
  // Electrons //
  ///////////////

  Handle<View<GsfElectron> > els_h;
  iEvent.getByToken(GSFElectron_, els_h);
  View<GsfElectron> gsfElColl = *(els_h.product());

  Handle<GsfElectronCollection> els_coll_h;
  iEvent.getByToken(GSFElectronCollection_, els_coll_h);    


  //////////////
  // PF Cands //
  //////////////

  iEvent.getByToken(PFCandidateCollection_, pfCand_h);


  ////////////
  // Vertex //
  ////////////

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(thePVCollection_, vertexHandle);


  ///////////////////////////
  // TransientTrackBuilder //
  ///////////////////////////

  ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);


  ////////////////////////////////////////////////
  // Get tools to get cluster shape information //
  ////////////////////////////////////////////////

  EcalClusterLazyTools* clusterTools_;
  clusterTools_ = new EcalClusterLazyTools( iEvent, iSetup, InputTag("reducedEcalRecHitsEB"), InputTag("reducedEcalRecHitsEE") );


  //////////////
  // Beamspot //
  //////////////

  Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByToken(BeamSpot_, beamspot_h);
  const reco::BeamSpot &beamSpotreco = *(beamspot_h.product());




  /////////////////////////
  // Loop Over Electrons //
  /////////////////////////

  double mass     = 0.000510998918;
  size_t elsIndex = 0;
  for( View<GsfElectron>::const_iterator el = els_h->begin(); el != els_h->end(); el++, elsIndex++ ) {

      ////////////////
      // References //
      ////////////////

      const Track*                 el_track         = (const Track*)(el->gsfTrack().get());
      const RefToBase<GsfElectron> gsfElRef         = els_h->refAt(elsIndex);    
      const TrackRef               ctfTkRef         = el->closestTrack();
      const GsfTrackRef            gsfTkRef         = el->gsfTrack();
      const VertexCollection*      vertexCollection = vertexHandle.product();

      ////////////
      // Vertex //
      ////////////
      VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
      int firstGoodVertexIdx = 0;
      for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx) {
          if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
              firstGoodVertex = vtx;
              break;
          }
      }

      //////////////////////
      // Fiduciality Mask //
      //////////////////////

      int fiducialityMask = 0;  // the enum is in interface/EgammaFiduciality.h
      if ( el->isEB()        ) fiducialityMask |= 1 << ISEB;
      if ( el->isEBEEGap()   ) fiducialityMask |= 1 << ISEBEEGAP;
      if ( el->isEE()        ) fiducialityMask |= 1 << ISEE;
      if ( el->isEEGap()     ) fiducialityMask |= 1 << ISEEGAP;
      if ( el->isEBEtaGap()  ) fiducialityMask |= 1 << ISEBETAGAP;
      if ( el->isEBPhiGap()  ) fiducialityMask |= 1 << ISEBPHIGAP;
      if ( el->isEEDeeGap()  ) fiducialityMask |= 1 << ISEEDEEGAP;
      if ( el->isEERingGap() ) fiducialityMask |= 1 << ISEERINGGAP;
      if ( el->isGap()       ) fiducialityMask |= 1 << ISGAP;


      ///////////////////////////
      // Corrections & Seeding //
      ///////////////////////////

      int electronTypeMask = 0;
      if ( el->isEcalEnergyCorrected()        ) electronTypeMask |= 1 << ISECALENERGYCORRECTED;
      if ( el->trackerDrivenSeed()            ) electronTypeMask |= 1 << ISTRACKERDRIVEN;
      if ( el->ecalDrivenSeed()               ) electronTypeMask |= 1 << ISECALDRIVEN;
      if ( el->passingCutBasedPreselection()  ) electronTypeMask |= 1 << ISCUTPRESELECTED;
      if ( el->passingMvaPreselection()       ) electronTypeMask |= 1 << ISMVAPRESELECTED;


      /////////////////////
      // Lorentz Vectors //
      /////////////////////

      LorentzVector    p4In;
      LorentzVector    p4Out;
      LorentzVector    trk_p4( el_track->px(), el_track->py(), el_track->pz(), el_track->p() );
      math::XYZVectorF p3In  = el->trackMomentumAtVtx();
      math::XYZVectorF p3Out = el->trackMomentumOut();
      p4In.SetXYZT (   p3In.x() , p3In.y() , p3In.z() , sqrt( mass*mass + p3In.R() *p3In.R()  ) );
      p4Out.SetXYZT(   p3Out.x(), p3Out.y(), p3Out.z(), sqrt( mass*mass + p3Out.R()*p3Out.R() ) );


      //////////////
      // Electron //
      //////////////

      HWWVal::els_fiduciality()        .push_back( fiducialityMask                                 );
      HWWVal::els_type()               .push_back( electronTypeMask                                );
      HWWVal::els_ecalEnergy()         .push_back( el->correctedEcalEnergy()                       ); 
      HWWVal::els_p4()                 .push_back( LorentzVector( el->p4() )                       );
      HWWVal::els_trk_p4()             .push_back( trk_p4                                          );
      HWWVal::els_vertex_p4()          .push_back( LorentzVector(el->vx(), el->vy(), el->vz(), 0.) );


      ///////////////
      // Isolation //
      ///////////////

      HWWVal::els_ecalIso()               .push_back( el->dr03EcalRecHitSumEt()                  );
      HWWVal::els_hcalIso()               .push_back( el->dr03HcalTowerSumEt()                   );
      HWWVal::els_tkIso()                 .push_back( el->dr03TkSumPt()                          );
      HWWVal::els_ecalIso04()             .push_back( el->dr04EcalRecHitSumEt()                  );
      HWWVal::els_hcalIso04()             .push_back( el->dr04HcalTowerSumEt()                   );


      //////////////////
      // PF Isolation //
      //////////////////

      if ( firstGoodVertex!=vertexCollection->end() ) {

          HWWVal::els_iso03_pf_ch()      .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 99999., 0.1, 0.07, 0.025, 0.025, 0  ) );
          HWWVal::els_iso03_pf_gamma05() .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 0.5   , 0.1, 0.07, 0.025, 0.025, 22 ) );
          HWWVal::els_iso03_pf_nhad05()  .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 0.5   , 0.1, 0.07, 0.025, 0.025, 130) );

          HWWVal::els_iso04_pf_ch()      .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4, 99999., 0.1, 0.07, 0.025, 0.025, 0  ) );
          HWWVal::els_iso04_pf_gamma05() .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4, 0.5   , 0.1, 0.07, 0.025, 0.025, 22 ) );
          HWWVal::els_iso04_pf_nhad05()  .push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4,  0.5  , 0.1, 0.07, 0.025, 0.025, 130) );


          // pf iso 2012
          float pfiso_ch = 0.0;
          float pfiso_em = 0.0;
          float pfiso_nh = 0.0;

          PFIsolation2012(*el, vertexCollection, firstGoodVertexIdx, 0.3, pfiso_ch, pfiso_em, pfiso_nh);
          HWWVal::els_iso03_pf2012_ch() .push_back( pfiso_ch );
          HWWVal::els_iso03_pf2012_em() .push_back( pfiso_em );
          HWWVal::els_iso03_pf2012_nh() .push_back( pfiso_nh );

          PFIsolation2012(*el, vertexCollection, firstGoodVertexIdx, 0.4, pfiso_ch, pfiso_em, pfiso_nh);
          HWWVal::els_iso04_pf2012_ch() .push_back( pfiso_ch );
          HWWVal::els_iso04_pf2012_em() .push_back( pfiso_em );
          HWWVal::els_iso04_pf2012_nh() .push_back( pfiso_nh );	    

      } else {

          HWWVal::els_iso03_pf_ch()      .push_back( -9999. );
          HWWVal::els_iso03_pf_gamma05() .push_back( -9999. );
          HWWVal::els_iso03_pf_nhad05()  .push_back( -9999. );

          HWWVal::els_iso04_pf_ch()      .push_back( -9999. );
          HWWVal::els_iso04_pf_gamma05() .push_back( -9999. );
          HWWVal::els_iso04_pf_nhad05()  .push_back( -9999. );

          HWWVal::els_iso03_pf2012_ch() .push_back( -9999. );
          HWWVal::els_iso03_pf2012_em() .push_back( -9999. );
          HWWVal::els_iso03_pf2012_nh() .push_back( -9999. );
          HWWVal::els_iso04_pf2012_ch() .push_back( -9999. );
          HWWVal::els_iso04_pf2012_em() .push_back( -9999. );
          HWWVal::els_iso04_pf2012_nh() .push_back( -9999. );
      }

      //////////////////
      // Supercluster //
      //////////////////

      HWWVal::els_etaSC()         .push_back( el->superCluster()->eta()             );
      HWWVal::els_eSC()           .push_back( el->superCluster()->energy()          );
      HWWVal::els_eSCRaw()        .push_back( el->superCluster()->rawEnergy()       );
      HWWVal::els_eSCPresh()      .push_back( el->superCluster()->preshowerEnergy() );
      HWWVal::els_nSeed()         .push_back( el->basicClustersSize() - 1           );
      HWWVal::els_e1x5()          .push_back( el->e1x5()                            );
      HWWVal::els_e5x5()          .push_back( el->e5x5()                            );
      HWWVal::els_sigmaIEtaIEta() .push_back( el->sigmaIetaIeta()                   );
      HWWVal::els_etaSCwidth()    .push_back( el->superCluster()->etaWidth()        );
      HWWVal::els_phiSCwidth()    .push_back( el->superCluster()->phiWidth()        );


      ///////////////////////////////////////////////////////
      // Get cluster info that is not stored in the object //
      ///////////////////////////////////////////////////////

      if( el->superCluster()->seed().isAvailable() ) { 

          const BasicCluster&  clRef              = *(el->superCluster()->seed());
          const vector<float>& lcovs              = clusterTools_->localCovariances(clRef);                    // get the local covariances computed in a 5x5 around the seed
          const vector<float>  localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));  // get the local covariances computed using all crystals in the SC

          HWWVal::els_sigmaIPhiIPhi()   .push_back( isfinite(lcovs[2])  ? lcovs[2] > 0  ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2])  : -9999. );
          HWWVal::els_e3x3()            .push_back( clusterTools_->e3x3(clRef) );
      } 
      else {

          HWWVal::els_sigmaIPhiIPhi()   .push_back(-9999.);
          HWWVal::els_e3x3()            .push_back(-9999.);

      } 

 
      ////////
      // ID //
      ////////

      HWWVal::els_hOverE()                        .push_back( el->hcalOverEcal()                   );
      HWWVal::els_eOverPIn()                      .push_back( el->eSuperClusterOverP()             );
      HWWVal::els_eSeedOverPOut()                 .push_back( el->eSeedClusterOverPout()           );
      HWWVal::els_eSeedOverPIn()                  .push_back( el->eSeedClusterOverP()              );
      HWWVal::els_eOverPOut()                     .push_back( el->eEleClusterOverPout()            );
      HWWVal::els_fbrem()                         .push_back( el->fbrem()                          );

      HWWVal::els_dEtaIn()                        .push_back( el->deltaEtaSuperClusterTrackAtVtx() );
      HWWVal::els_dEtaOut()                       .push_back( el->deltaEtaSeedClusterTrackAtCalo() );
      HWWVal::els_dPhiIn()                        .push_back( el->deltaPhiSuperClusterTrackAtVtx() );
      HWWVal::els_dPhiOut()                       .push_back( el->deltaPhiSeedClusterTrackAtCalo() );

      ////////////
      // Tracks //
      ////////////

      HWWVal::els_chi2()                  .push_back( el_track->chi2()                          );
      HWWVal::els_ndof()                  .push_back( el_track->ndof()                          );
      HWWVal::els_gsftrkidx()             .push_back( static_cast<int>((el->gsfTrack()).key())  );
      HWWVal::els_charge()                .push_back( el->charge()                              );
      HWWVal::els_trk_charge()            .push_back( el_track->charge()                        );
      HWWVal::els_sccharge()              .push_back( el->scPixCharge()                         );
      HWWVal::els_d0()                    .push_back( el_track->d0()                            );
      HWWVal::els_d0corr()                .push_back( -1*(el_track->dxy(beamSpotreco))              );
      HWWVal::els_z0corr()                .push_back( el_track->dz(beamSpotreco.position(el_track->vz()))                    );
 

      /////////
      // CTF //
      /////////

      if( ctfTkRef.isNonnull() ) {
          HWWVal::els_trkidx()    . push_back( static_cast<int>  ( ctfTkRef.key()        )                                  );
          HWWVal::els_trkshFrac() . push_back( static_cast<float>( el->ctfGsfOverlap() )                                    );
      } 
      else {
          HWWVal::els_trkidx()    . push_back(-9999.);
          HWWVal::els_trkshFrac() . push_back(-9999.);
      }

      
      ////////////////////
      // Regular Vertex //
      ////////////////////        
      TransientTrack tt = theTTBuilder->build(el->gsfTrack());
  
      if ( firstGoodVertex!=vertexCollection->end() ) {
          Measurement1D ip3D_regular = IPTools::absoluteImpactParameter3D(tt, *firstGoodVertex).second;

          HWWVal::els_ip3d()      . push_back( ip3D_regular.value() );
          HWWVal::els_ip3derr()   . push_back( ip3D_regular.error() );
      } else {

          HWWVal::els_ip3d()      . push_back( -999. );
          HWWVal::els_ip3derr()   . push_back( -999. );
      }


      /////////////////
      // Hit Pattern //
      /////////////////

      const HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 

      HWWVal::els_exp_innerlayers().push_back(p_inner.numberOfHits());


      /////////////////
      // Conversions //
      /////////////////

      ConversionFinder convFinder; //vector of conversion infos - all the candidate conversion tracks
      vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(el->core()), tracks_h, gsftracks_h, evt_bField);
  
      vector<int>           v_tkidx;
      vector<int>           v_gsftkidx;
      vector<int>           v_delmisshits;
      vector<int>           v_flag;
      vector<float>         v_dist;
      vector<float>         v_dcot;
      vector<float>         v_rad;
      vector<LorentzVector> v_pos_p4;

      for(unsigned int i_conv = 0; i_conv < v_convInfos.size(); i_conv++) {
    
          math::XYZPoint convPoint  = v_convInfos.at(i_conv).pointOfConversion();
          float          convPointX = isfinite(convPoint.x()) ? convPoint.x() : -9999.;
          float          convPointY = isfinite(convPoint.y()) ? convPoint.y() : -9999.;
          float          convPointZ = isfinite(convPoint.z()) ? convPoint.z() : -9999.;

          v_dist        .push_back( isfinite(v_convInfos.at(i_conv).dist()) ? v_convInfos.at(i_conv).dist() : -9999.  );
          v_dcot        .push_back( v_convInfos.at(i_conv).dcot()                                                     );
          v_rad         .push_back( v_convInfos.at(i_conv).radiusOfConversion()                                       );
          v_delmisshits .push_back( v_convInfos.at(i_conv).deltaMissingHits()                                         );
          v_flag        .push_back( v_convInfos.at(i_conv).flag()                                                     );
          v_pos_p4      .push_back( LorentzVector(convPointX, convPointY, convPointZ, 0)                              );

          if( v_convInfos.at(i_conv).conversionPartnerCtfTk().isNonnull() ) {
              v_tkidx.push_back(v_convInfos.at(i_conv).conversionPartnerCtfTk().key());
          }
          else {
              v_tkidx.push_back(-9999);
          }

          //
          if( v_convInfos.at(i_conv).conversionPartnerGsfTk().isNonnull() ) {
              v_gsftkidx.push_back(v_convInfos.at(i_conv).conversionPartnerGsfTk().key());
          }
          else { 
              v_gsftkidx.push_back(-9999);
          }

      } // end for loop


      ConversionInfo convInfo   = convFinder.getConversionInfo( *el, tracks_h, gsftracks_h, evt_bField );

      HWWVal::els_conv_dist().push_back( isfinite(convInfo.dist()) ? convInfo.dist() : -9999. );
      HWWVal::els_conv_dcot().push_back( convInfo.dcot()                                      );


      //////////////////////////////////////////////
      // Flag For Vertex Fit Conversion Rejection //
      //////////////////////////////////////////////

      Handle<ConversionCollection> convs_h;
      iEvent.getByToken(ConversionCollection_, convs_h);


      //////////////////////////////
      // Old Conversion Rejection //
      //////////////////////////////

      HWWVal::els_conv_old_dist()        . push_back( isfinite(el->convDist())   ? el->convDist()   : -9999. );
      HWWVal::els_conv_old_dcot()        . push_back( isfinite(el->convDcot())   ? el->convDcot()   : -9999. );


      //////////////////////
      // 2012 Electron ID //
      //////////////////////

      GsfElectronRef ele(els_coll_h, elsIndex);

  } // end Loop on Electrons

}


double ElectronMaker::electronIsoValuePF(const GsfElectron& el, const Vertex& vtx, float coner, float minptn, float dzcut,
                                         float footprintdr, float gammastripveto, float elestripveto, int filterId){

    float pfciso = 0.;
    float pfniso = 0.;
    float pffootprint = 0.;
    float pfjurveto = 0.;
    float pfjurvetoq = 0.;

    TrackRef siTrack     = el.closestTrack();
    GsfTrackRef gsfTrack = el.gsfTrack();

    if (gsfTrack.isNull() && siTrack.isNull()) return -9999.;

    float eldz = gsfTrack.isNonnull() ? gsfTrack->dz(vtx.position()) : siTrack->dz(vtx.position());
    float eleta = el.eta();

    for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

        float pfeta = pf->eta();    
        float dR = deltaR(pfeta, pf->phi(), eleta, el.phi());
        if (dR>coner) continue;

        float deta = fabs(pfeta - eleta);
        int pfid = abs(pf->pdgId());
        float pfpt = pf->pt();

        if (filterId!=0 && filterId!=pfid) continue;

        if (pf->charge()==0) {
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
            const TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //below pfid==1 is commented out: in some cases the pfCand has a gsf even if it is not an electron... this is to improve the sync with MIT
            if (/*pfid==11 &&*/ pf->gsfTrackRef().isNonnull()) {
                if (gsfTrack.isNonnull() && gsfTrack.key()==pf->gsfTrackRef().key()) continue;
            } 
            //check electrons with gsf track
            if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if(fabs( pfTrack->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
            }
        } 
    }
    return pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq;
}


void ElectronMaker::PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection,
        const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh)
{

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0;
    pfiso_nh = 0.0;

    // loop on pfcandidates
    reco::PFCandidateCollection::const_iterator pf = pfCand_h->begin();
    for (pf = pfCand_h->begin(); pf != pfCand_h->end(); ++pf) {

        // skip electrons and muons
        if (pf->particleId() == reco::PFCandidate::e)     continue;
        if (pf->particleId() == reco::PFCandidate::mu)    continue;

        // deltaR between electron and cadidate
        const float dR = deltaR(pf->eta(), pf->phi(), el.eta(), el.phi());
        if (dR > R)                             continue;

        PFPileUpAlgo *pfPileUpAlgo_ = new PFPileUpAlgo();

        if (pf->particleId() == reco::PFCandidate::h) {
            int pfVertexIndex = pfPileUpAlgo_->chargedHadronVertex(*vertexCollection, *pf);
            if (pfVertexIndex != vertexIndex) continue;
        }

        // endcap region
        if (!el.isEB()) {
            if (pf->particleId() == reco::PFCandidate::h      && dR <= 0.015)   continue;
            if (pf->particleId() == reco::PFCandidate::gamma  && dR <= 0.08)    continue;
        }

        // add to isolation sum
        if (pf->particleId() == reco::PFCandidate::h)       pfiso_ch += pf->pt();
        if (pf->particleId() == reco::PFCandidate::gamma)   pfiso_em += pf->pt();
        if (pf->particleId() == reco::PFCandidate::h0)      pfiso_nh += pf->pt();

    }

}

