#ifndef EVENTMAKER_H
#define EVENTMAKER_H

#include "FWCore/Framework/interface/Event.h"

class EventMaker {

  public:

    EventMaker();
    void SetVars(const edm::Event&, const edm::EventSetup&);

};

#endif
