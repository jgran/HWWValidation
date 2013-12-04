#ifndef PFELECTRONMAKER_H
#define PFELECTRONMAKER_H

//
// class decleration
//

class PFElectronMaker : public edm::EDProducer {
public:
     explicit PFElectronMaker (const edm::ParameterSet&);
     ~PFElectronMaker();

private:
//  virtual void beginJob() ;
     virtual void beginJob() ;
     virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
    
     // ----------member data ---------------------------
     edm::InputTag pfCandidatesTag_;
/*
     edm::InputTag isoc_vm_tag_;
     edm::InputTag ison_vm_tag_;
     edm::InputTag isop_vm_tag_;
     edm::InputTag isoc04_vm_tag_;
     edm::InputTag ison04_vm_tag_;
     edm::InputTag isop04_vm_tag_;
     edm::InputTag pfAllElectrons_tag_;  
*/
};

#endif

