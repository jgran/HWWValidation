#ifndef TRKMETMAKER_H
#define TRKMETMAKER_H

#include "FWCore/Framework/interface/MakerMacros.h"

class TrkMETMaker {

  public:

    TrkMETMaker() {};
    void SetVars(const edm::Event&, const edm::EventSetup&);

};

#endif

