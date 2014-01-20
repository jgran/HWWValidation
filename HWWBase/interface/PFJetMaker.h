#ifndef PFJETMAKER_H
#define PFJETMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/JetReco/interface/PFJet.h"

class PFJetMaker {

  public:

    PFJetMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<reco::PFJetCollection> PFJetCollection_;

};

#endif
