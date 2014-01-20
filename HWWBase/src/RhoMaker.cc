#include "HWWValidation/HWWBase/interface/RhoMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

RhoMaker::RhoMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  Rho_         = iCollector.consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInputTag"));
  wwRho_       = iCollector.consumes<double>(iConfig.getParameter<edm::InputTag>("wwrhoInputTag"));
  wwRhoVor_    = iCollector.consumes<double>(iConfig.getParameter<edm::InputTag>("wwrhovorInputTag"));
  RhoForEGIso_ = iCollector.consumes<double>(iConfig.getParameter<edm::InputTag>("forEGIsoInputTag"));

}

void RhoMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_evt_rho();
  HWWVal::Load_evt_ww_rho();
  HWWVal::Load_evt_ww_rho_vor();
  HWWVal::Load_evt_kt6pf_foregiso_rho();

  edm::Handle<double> rhoH;
  edm::Handle<double> ww_rhoH;
  edm::Handle<double> ww_rho_vorH;
  edm::Handle<double> kt6pf_foregiso_rhoH;

  iEvent.getByToken( Rho_, rhoH);
  if(rhoH.isValid()){
    HWWVal::evt_rho() = *rhoH;
  }
  iEvent.getByToken( wwRho_, ww_rhoH);
  if(ww_rhoH.isValid()){
    HWWVal::evt_ww_rho() = *ww_rhoH;
  }
  iEvent.getByToken( wwRhoVor_ , ww_rho_vorH);
  if(ww_rho_vorH.isValid()){
    HWWVal::evt_ww_rho_vor() = *ww_rho_vorH;
  }
  iEvent.getByToken( RhoForEGIso_ , kt6pf_foregiso_rhoH);
  if(kt6pf_foregiso_rhoH.isValid()){
    HWWVal::evt_kt6pf_foregiso_rho() = *kt6pf_foregiso_rhoH;
  }
}
