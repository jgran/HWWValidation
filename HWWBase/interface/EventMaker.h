#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TString.h"
//
// class decleration
//

class EventMaker : public edm::EDProducer {
public:
     explicit EventMaker (const edm::ParameterSet&);
     ~EventMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     edm::InputTag dcsTag_;
	   std::string aliasprefix_;
     bool isData_;
     //std::string datasetName_;
     //std::string CMS2tag_;
};


#endif
