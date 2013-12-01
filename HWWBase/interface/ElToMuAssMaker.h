#ifndef ELTOMUASSMAKER_H
#define ELTOMUASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class ElToMuAssMaker : public edm::EDProducer {
public:
     explicit ElToMuAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
  // ----------member data ---------------------------
  double m_minDR_;
  std::string aliasprefix_;
  edm::InputTag elsInputTag_;
  edm::InputTag musInputTag_;
};


#endif
