// user include files
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "HWWValidation/HWWBase/interface/GSFTrackMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using std::vector;
using reco::Track;
using reco::TrackBase;

//
// class decleration
//

//
// constructors and destructor
//
GSFTrackMaker::GSFTrackMaker(const edm::ParameterSet& iConfig) {
       
  produces<vector<LorentzVector> >	("gsftrksp4"		).setBranchAlias("gsftrks_p4"     );	// track p4						
  produces<vector<LorentzVector> >	("gsftrksvertexp4"	).setBranchAlias("gsftrks_vertex_p4"  );	// track p4
  produces<vector<float> >		("gsftrksd0"		).setBranchAlias("gsftrks_d0"         );	// impact parameter at the point of closest approach	
  produces<vector<float> >		("gsftrksz0"		).setBranchAlias("gsftrks_z0"         );	// z position of the point of closest approach		
  produces<vector<float> >		("gsftrksd0Err"		).setBranchAlias("gsftrks_d0Err"      );	// error on the impact parameter			
  produces<vector<float> >		("gsftrksz0Err"		).setBranchAlias("gsftrks_z0Err"      );	// error on z position of the point of closest approach	
  produces<vector<float> >		("gsftrksetaErr"	).setBranchAlias("gsftrks_etaErr"     );	// track eta error					
  produces<vector<float> >		("gsftrksphiErr"	).setBranchAlias("gsftrks_phiErr"     );	// track phi error					
  produces<vector<float> >		("gsftrksd0phiCov"	).setBranchAlias("gsftrks_d0phiCov"   ); // track cov(d0, phi) 

  gsftracksInputTag_ = iConfig.getParameter<edm::InputTag>("gsftracksInputTag");
  beamSpot_tag_     = iConfig.getParameter<edm::InputTag> ("beamSpotTag"     );


}

void GSFTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::auto_ptr<vector<LorentzVector> >	gsftrks_p4		(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> >	gsftrks_vertex_p4	(new vector<LorentzVector>	);
  std::auto_ptr<vector<float> >		gsftrks_d0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_z0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_d0Err		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_z0Err		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_etaErr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_phiErr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_d0phiCov		(new vector<float>		);      

  // get tracks
  Handle<edm::View<reco::GsfTrack> > track_h;
  iEvent.getByLabel(gsftracksInputTag_, track_h);

  if( !track_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve track collection";
    edm::LogInfo("OutputInfo") << " GSFTrackMaker cannot continue...!";
    return;
  }

  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);

  if( !theMagField.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve the magnetic field";
    edm::LogInfo("OutputInfo") << " GSFTrackMaker cannot continue...!";
    return;
  }
  
  //get tracker geometry   
  edm::ESHandle<TrackerGeometry> theG;
  iSetup.get<TrackerDigiGeometryRecord>().get(theG);
  
  edm::View<reco::GsfTrack>::const_iterator tracks_end = track_h->end();
  
  for (edm::View<reco::GsfTrack>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {

    gsftrks_p4           ->push_back( LorentzVector( i->px(), i->py(), i->pz(), i->p() )       );
    gsftrks_vertex_p4    ->push_back( LorentzVector(i->vx(),i->vy(), i->vz(), 0.)              );
    gsftrks_d0           ->push_back( i->d0()                                                  );
    gsftrks_z0           ->push_back( i->dz()                                                  );
    gsftrks_d0Err        ->push_back( i->d0Error()                                             );
    gsftrks_z0Err        ->push_back( i->dzError()                                             );
    gsftrks_etaErr       ->push_back( i->etaError()                                            );
    gsftrks_phiErr       ->push_back( i->phiError()                                            );
    gsftrks_d0phiCov     ->push_back( -i->covariance(TrackBase::i_phi, TrackBase::i_dxy)       );  // minus sign because we want cov(d0, phi) = cov(-dxy, phi)

  }

  // store vectors
  iEvent.put(gsftrks_p4			, "gsftrksp4"			);
  iEvent.put(gsftrks_vertex_p4		, "gsftrksvertexp4"		);
  iEvent.put(gsftrks_d0			, "gsftrksd0"			);
  iEvent.put(gsftrks_z0			, "gsftrksz0"			);
  iEvent.put(gsftrks_d0Err		, "gsftrksd0Err"		);
  iEvent.put(gsftrks_z0Err		, "gsftrksz0Err"		);
  iEvent.put(gsftrks_etaErr		, "gsftrksetaErr"		);
  iEvent.put(gsftrks_phiErr		, "gsftrksphiErr"		);
  iEvent.put(gsftrks_d0phiCov		, "gsftrksd0phiCov"		);

}

// ------------ method called once each job just before starting event loop  ------------
void GSFTrackMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GSFTrackMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GSFTrackMaker);
