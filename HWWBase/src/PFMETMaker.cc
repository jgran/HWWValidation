#include "HWWValidation/HWWBase/interface/PFMETMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  PFMET_ = iCollector.consumes<edm::View<reco::PFMET> >(iConfig.getParameter<edm::InputTag>("pfmetInputTag"));

}

void PFMETMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    HWWVal::Load_evt_pfmet();
    HWWVal::Load_evt_pfmetPhi();
  
    edm::Handle<edm::View<reco::PFMET> > met_h;
    iEvent.getByToken(PFMET_, met_h);

    HWWVal::evt_pfmet()    = ( met_h->front() ).et();
    HWWVal::evt_pfmetPhi() = ( met_h->front() ).phi();
}
