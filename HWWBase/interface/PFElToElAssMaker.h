#ifndef PFELTOELASSMAKER_H
#define PFELTOELASSMAKER_H

//
// class declaration
//

class PFElToElAssMaker : public edm::EDProducer {
public:
     explicit PFElToElAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
     double m_minDR_;
     std::string aliasprefix_;
     edm::InputTag elsInputTag_;
     edm::InputTag pfelsInputTag_;
};


#endif
