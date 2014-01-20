#ifndef MVAJETIDMAKER_H
#define MVAJETIDMAKER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "HWWValidation/HWWBase/interface/PileupJetIdAlgo.h"

class MVAJetIdMaker {

  public:

    MVAJetIdMaker(const edm::ParameterSet&, edm::ConsumesCollector);
    void SetVars(const edm::Event&, const edm::EventSetup&);

  private:

    edm::EDGetTokenT<reco::PFJetCollection>       PFJetCollection_;
    edm::EDGetTokenT<reco::PFJetCollection>       CorrPFJetCollection_;
    edm::EDGetTokenT<reco::VertexCollection>      thePVCollection_;
    PileupJetIdAlgo  *fPUJetIdAlgo;

};

#endif
