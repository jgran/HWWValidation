#ifndef BTAGMAKER_H
#define BTAGMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"

class BTagMaker {

  public:

    BTagMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<edm::View<reco::Jet> >                PFJets_;
    edm::EDGetTokenT<reco::JetFloatAssociation::Container> BJetTags_;

};

#endif

