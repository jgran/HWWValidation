#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "HWWValidation/HWWBase/interface/TrackMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;
using reco::Track;
using reco::TrackBase;

TrackMaker::TrackMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  TrackCollection_ = iCollector.consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackInputTag"));

}

void TrackMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;

  HWWVal::Load_trks_trk_p4();
  HWWVal::Load_trks_d0();
  HWWVal::Load_trks_chi2();
  HWWVal::Load_trks_ndof();
  HWWVal::Load_trks_z0();
  HWWVal::Load_trks_d0Err();
  HWWVal::Load_trks_z0Err();
  HWWVal::Load_trks_etaErr();
  HWWVal::Load_trks_phiErr();
  HWWVal::Load_trks_d0phiCov();
  HWWVal::Load_trks_charge();
  HWWVal::Load_trks_qualityMask();
  HWWVal::Load_trks_valid_pixelhits();
  HWWVal::Load_trks_nlayers();

  ////////////////
  // Get Tracks //
  ////////////////

  Handle<edm::View<reco::Track> > track_h;
  iEvent.getByToken(TrackCollection_, track_h); 
  if( !track_h.isValid() ) {
    LogInfo("OutputInfo") << " failed to retrieve track collection";
    LogInfo("OutputInfo") << " TrackMaker cannot continue...!";
    return;
  }

  /////////////////
  // Get B Field //
  /////////////////

  ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  if( !theMagField.isValid() ) {
    LogInfo("OutputInfo") << " failed to retrieve the magnetic field";
    LogInfo("OutputInfo") << " TrackMaker cannot continue...!";
    return;
  }
  const MagneticField* bf = theMagField.product();

  //////////////////////////
  // Get Tracker Geometry //
  //////////////////////////  

  ESHandle<TrackerGeometry> theG;
  iSetup.get<TrackerDigiGeometryRecord>().get(theG);
  
  //////////////////////
  // Loop over Tracks //
  //////////////////////

  View<reco::Track>::const_iterator tracks_end = track_h->end();
  unsigned int iTIndex=-1;
  for (View<reco::Track>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {

    iTIndex++;

    HWWVal::trks_trk_p4()       .push_back( LorentzVector( i->px(), i->py(), i->pz(), i->p() )       );
    HWWVal::trks_d0()           .push_back( i->d0()                                                  );
    HWWVal::trks_chi2()         .push_back( i->chi2()                                                );
    HWWVal::trks_ndof()         .push_back( i->ndof()                                                );
    HWWVal::trks_z0()           .push_back( i->dz()                                                  );
    HWWVal::trks_d0Err()        .push_back( i->d0Error()                                             );
    HWWVal::trks_z0Err()        .push_back( i->dzError()                                             );
    HWWVal::trks_etaErr()       .push_back( i->etaError()                                            );
    HWWVal::trks_phiErr()       .push_back( i->phiError()                                            );
    HWWVal::trks_d0phiCov()     .push_back( -i->covariance(TrackBase::i_phi, TrackBase::i_dxy)	      ); 
    HWWVal::trks_charge()       .push_back( i->charge()                                              );
    HWWVal::trks_qualityMask()  .push_back( i->qualityMask()                                         );

    GlobalPoint  tpVertex   ( i->vx(), i->vy(), i->vz() );
    GlobalVector tpMomentum ( i->px(), i->py(), i->pz() );
    int tpCharge ( i->charge() );
    
    FreeTrajectoryState fts ( tpVertex, tpMomentum, tpCharge, bf);

    const float zdist  = 314.;
    const float radius = 130.;
    const float corner = 1.479;

    Plane::PlanePointer lendcap      = Plane::build( Plane::PositionType (0, 0, -zdist), Plane::RotationType () );    
    Plane::PlanePointer rendcap      = Plane::build( Plane::PositionType (0, 0,  zdist), Plane::RotationType () );
    Cylinder::CylinderPointer barrel = Cylinder::build( Cylinder::PositionType (0, 0, 0), Cylinder::RotationType (), radius);
    AnalyticalPropagator myAP (bf, alongMomentum, 2*M_PI);
    TrajectoryStateOnSurface tsos;    
        
    if( i->eta() < -corner ) {
      tsos = myAP.propagate( fts, *lendcap);
    }
    else if( fabs(i->eta()) < corner ) {
      tsos = myAP.propagate( fts, *barrel);
    }
    else if( i->eta() > corner ) {
      tsos = myAP.propagate( fts, *rendcap);
    }
    
    const reco::HitPattern& pattern = i->hitPattern();    
    HWWVal::trks_valid_pixelhits() .push_back(pattern.numberOfValidPixelHits());
      
    if(i->extra().isAvailable()) {
      bool valid_hit      = false;
      uint32_t hit_pattern; 
      int i_layer       = 1;
      bool pixel_hit   = false;
      bool strip_hit   = false;


      typedef Ref<edmNew::DetSetVector<SiStripCluster>,SiStripCluster > ClusterRef;
      typedef Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;


      for(trackingRecHit_iterator ihit = i->recHitsBegin(); ihit != i->recHitsEnd(); ++ihit){
        if(i_layer > 1) break;
        int k = ihit-i->recHitsBegin();
        hit_pattern = pattern.getHitPattern(k);
        valid_hit = pattern.validHitFilter(hit_pattern);
        pixel_hit = pattern.pixelHitFilter(hit_pattern);
        strip_hit = pattern.stripHitFilter(hit_pattern);
        if(!valid_hit) continue;
        if(pixel_hit){
          const SiPixelRecHit *pixel_hit_cast = dynamic_cast<const SiPixelRecHit*>(&(**ihit));
          assert(pixel_hit_cast != 0);
          if(i_layer == 1){
            i_layer++;

          }
        }
        else if (strip_hit){
          const SiStripRecHit1D *strip_hit_cast = dynamic_cast<const SiStripRecHit1D*>(&(**ihit));
          const SiStripRecHit2D *strip2d_hit_cast = dynamic_cast<const SiStripRecHit2D*>(&(**ihit));
          ClusterRef cluster;
          if(strip_hit_cast == NULL)
            cluster = strip2d_hit_cast->cluster();
          else 
            cluster = strip_hit_cast->cluster();
          int cluster_size   = (int)cluster->amplitudes().size();
          int cluster_charge = 0;
          double   cluster_weight_size = 0.0;
          int max_strip_i = std::max_element(cluster->amplitudes().begin(),cluster->amplitudes().end())-cluster->amplitudes().begin();
          for(int istrip = 0; istrip < cluster_size; istrip++){
            cluster_charge += (int)cluster->amplitudes().at(istrip);
            cluster_weight_size += (istrip-max_strip_i)*(istrip-max_strip_i)*(cluster->amplitudes().at(istrip));
          }
          cluster_weight_size = sqrt(cluster_weight_size/cluster_charge);
          if(i_layer == 1){
            i_layer++;
          }
        }
      }
    }
    
    HWWVal::trks_nlayers()    .push_back( i->hitPattern().trackerLayersWithMeasurement() );

  } // End loop on tracks

}
