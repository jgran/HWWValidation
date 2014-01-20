#ifndef HYPDILEPMAKER_H
#define HYPDILEPMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"

class HypDilepMaker {

  public:

    HypDilepMaker(){};
    void SetVars(const edm::Event&, const edm::EventSetup&);

};

#endif
