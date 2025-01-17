#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "HWWValidation/HWWBase/interface/GSFTrackMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

GSFTrackMaker::GSFTrackMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector) {

  GSFTrack_ = iCollector.consumes<edm::View<reco::GsfTrack> > (iConfig.getParameter<edm::InputTag>("gsftracksInputTag"));

}

void GSFTrackMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using std::vector;
  using reco::Track;
  using reco::TrackBase;

  HWWVal::Load_gsftrks_p4();
  HWWVal::Load_gsftrks_vertex_p4();
  HWWVal::Load_gsftrks_d0();
  HWWVal::Load_gsftrks_z0();
  HWWVal::Load_gsftrks_d0Err();
  HWWVal::Load_gsftrks_z0Err();
  HWWVal::Load_gsftrks_etaErr();
  HWWVal::Load_gsftrks_phiErr();
  HWWVal::Load_gsftrks_d0phiCov();

  Handle<edm::View<reco::GsfTrack> > track_h;
  iEvent.getByToken(GSFTrack_, track_h);

  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  
  edm::ESHandle<TrackerGeometry> theG;
  iSetup.get<TrackerDigiGeometryRecord>().get(theG);
  
  edm::View<reco::GsfTrack>::const_iterator tracks_end = track_h->end();
  
  for (edm::View<reco::GsfTrack>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {

    HWWVal::gsftrks_p4()           .push_back( LorentzVector( i->px(), i->py(), i->pz(), i->p() )       );
    HWWVal::gsftrks_vertex_p4()    .push_back( LorentzVector(i->vx(),i->vy(), i->vz(), 0.)              );
    HWWVal::gsftrks_d0()           .push_back( i->d0()                                                  );
    HWWVal::gsftrks_z0()           .push_back( i->dz()                                                  );
    HWWVal::gsftrks_d0Err()        .push_back( i->d0Error()                                             );
    HWWVal::gsftrks_z0Err()        .push_back( i->dzError()                                             );
    HWWVal::gsftrks_etaErr()       .push_back( i->etaError()                                            );
    HWWVal::gsftrks_phiErr()       .push_back( i->phiError()                                            );
    HWWVal::gsftrks_d0phiCov()     .push_back( -i->covariance(TrackBase::i_phi, TrackBase::i_dxy)       );

  }
}
