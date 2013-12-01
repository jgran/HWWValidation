#ifndef JETCOLLECTIONPRUNER_H
#define JETCOLLECTIONPRUNER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

//
// class declaration
//

class JetCollectionPruner : public edm::EDProducer {
public:
  explicit JetCollectionPruner (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag inputUncorrectedJPTJetCollection_;
  edm::InputTag inputUncorrectedPFJetCollection_;
  edm::InputTag inputUncorrectedTrkJetCollection_;

  edm::InputTag inputCaloJetIdCollection_;
  double uncorrectedJPTJetPtCut_;
  double uncorrectedPFJetPtCut_;
  double uncorrectedTrkJetPtCut_;
};

#endif
