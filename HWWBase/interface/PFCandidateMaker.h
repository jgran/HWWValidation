#ifndef PFCANDIDATEMAKER_H
#define PFCANDIDATEMAKER_H

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

class PFCandidateMaker {

  public:

    PFCandidateMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<reco::PFCandidateCollection>             PFCandidateCollection_;
    edm::EDGetTokenT<edm::ValueMap<reco::PFCandidatePtr> >    PFElectrons_;
    edm::EDGetTokenT<reco::TrackCollection>                   TrackCollection_;
    edm::EDGetTokenT<reco::VertexCollection>                  thePVCollection_;

};

#endif

