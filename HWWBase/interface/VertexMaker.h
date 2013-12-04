#ifndef CMS2_SCMAKER_H
#define CMS2_SCMAKER_H

//
// class declaration
//

class VertexMaker : public edm::EDProducer {
public:
  explicit VertexMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag primaryVertexInputTag_;

	std::string aliasprefix_;
};

#endif
