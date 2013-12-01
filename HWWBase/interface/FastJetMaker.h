#ifndef FASTJETMAKER_H
#define FASTJETMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class FastJetMaker : public edm::EDProducer {
public:
  explicit FastJetMaker (const edm::ParameterSet&);
  ~FastJetMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag rhoJEC_tag;
  edm::InputTag rhoIso_tag;
  std::string aliasprefix_;
};


#endif
