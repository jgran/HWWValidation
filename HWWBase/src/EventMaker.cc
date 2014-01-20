#include "HWWValidation/HWWBase/interface/EventMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

EventMaker::EventMaker() {}

void EventMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     HWWVal::Load_evt_run();
     HWWVal::Load_evt_event();
     HWWVal::Load_evt_lumiBlock();
     HWWVal::Load_evt_isRealData();

     HWWVal::evt_run()                 = iEvent.id().run()       ;
     HWWVal::evt_event()               = iEvent.id().event()     ;
     HWWVal::evt_lumiBlock()           = iEvent.luminosityBlock();
     HWWVal::evt_isRealData()          = iEvent.isRealData()     ;

}
