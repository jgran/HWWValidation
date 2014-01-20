#ifndef TRACKMAKER_H
#define TRACKMAKER_H

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class TrackMaker {

  public:

    TrackMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<edm::View<reco::Track> > TrackCollection_;

};

#endif
