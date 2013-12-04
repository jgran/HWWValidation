#ifndef FASTJETMAKER_H
#define FASTJETMAKER_H

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
