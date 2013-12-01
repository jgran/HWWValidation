#ifndef TRACKMAKER_H
#define TRACKMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/Point3D.h"

typedef math::XYZPoint Point;

class TrackMaker : public edm::EDProducer {
public:
  explicit TrackMaker (const edm::ParameterSet&);
  double calculateTrkIsolation(const edm::View<reco::Track>*, const reco::Track&, const Point&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag trackInputTag;
  edm::InputTag primaryVertexInputTag_;
  
  float dRConeMin_;
  float dRConeMax_;
  float vtxDiffZMax_;
  float tkVtxDMax_;
  float ptMin_;
  int   nHits_;
  std::string aliasprefix_;
};

#endif
