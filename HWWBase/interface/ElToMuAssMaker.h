#ifndef ELTOMUASSMAKER_H
#define ELTOMUASSMAKER_H

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
