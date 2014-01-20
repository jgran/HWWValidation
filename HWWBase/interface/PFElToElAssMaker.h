#ifndef PFELTOELASSMAKER_H
#define PFELTOELASSMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"

class PFElToElAssMaker {

  public:

    PFElToElAssMaker() {};
    void SetVars(const edm::Event&, const edm::EventSetup&);

};

#endif
