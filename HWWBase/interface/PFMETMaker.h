#ifndef PFMETMAKER_H
#define PFMETMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class PFMETMaker : public edm::EDProducer {
public:
    explicit PFMETMaker (const edm::ParameterSet&);
    ~PFMETMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::InputTag pfMetInputTag;
    edm::InputTag pfMetCorInputTag;
	  std::string aliasprefix_;
    
};


#endif
