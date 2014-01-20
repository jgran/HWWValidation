#include "FWCore/Framework/interface/MakerMacros.h"
#include "HWWValidation/HWWBase/interface/VertexMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;

VertexMaker::VertexMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  thePVCollection_       = iCollector.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));

}

void VertexMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  HWWVal::Load_vtxs_position();
  HWWVal::Load_vtxs_xError();
  HWWVal::Load_vtxs_yError();
  HWWVal::Load_vtxs_zError();
  HWWVal::Load_vtxs_ndof();
  HWWVal::Load_vtxs_isFake();
  HWWVal::Load_vtxs_sumpt();
  HWWVal::Load_vtxs_covMatrix();

  // get the primary vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(thePVCollection_, vertexHandle); 
  const reco::VertexCollection *vertexCollection = vertexHandle.product();

  unsigned int index = 0;
  const unsigned int covMatrix_dim = 3;

  for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++index) {
    HWWVal::vtxs_position()         .push_back( LorentzVector( vtx->position().x(), vtx->position().y(), vtx->position().z(), 0 ) );
    HWWVal::vtxs_xError()           .push_back( vtx->xError()            );
    HWWVal::vtxs_yError()           .push_back( vtx->yError()            );
    HWWVal::vtxs_zError()           .push_back( vtx->zError()            );
    HWWVal::vtxs_ndof()             .push_back( vtx->ndof()              );
    HWWVal::vtxs_isFake()           .push_back( vtx->isFake()            );
    double sumpt = 0;
    for (reco::Vertex::trackRef_iterator i = vtx->tracks_begin(); i != vtx->tracks_end(); ++i) sumpt += (*i)->pt();

    HWWVal::vtxs_sumpt()		 .push_back(sumpt);

    std::vector<float> temp_vec;
    temp_vec.clear();

    for( unsigned int i = 0; i < covMatrix_dim; i++ ) {
      for( unsigned int j = 0; j < covMatrix_dim; j++ ) {
        temp_vec.push_back( vtx->covariance(i, j) );
      }
    }

    HWWVal::vtxs_covMatrix().push_back( temp_vec );
    
  }

}
