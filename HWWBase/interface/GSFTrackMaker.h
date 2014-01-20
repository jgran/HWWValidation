#ifndef GSFTRACKMAKER_H
#define GSFTRACKMAKER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

class GSFTrackMaker {

  public:

    GSFTrackMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<edm::View<reco::GsfTrack> >              GSFTrack_;

};
#endif
