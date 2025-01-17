#ifndef RECOCONVERSIONMAKER_H
#define RECOCONVERSIONMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

class RecoConversionMaker {

  public:

    RecoConversionMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<edm::View<reco::Conversion> >     Conversion_;
    edm::EDGetTokenT<reco::BeamSpot>                   BeamSpot_;

};

#endif
