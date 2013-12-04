#ifndef MVAJETIDMAKER_H
#define MVAJETIDMAKER_H

//
// class decleration
//

class MVAJetIdMaker : public edm::EDProducer {
public:
  explicit MVAJetIdMaker(const edm::ParameterSet&);
  ~MVAJetIdMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool passPFLooseId(const reco::PFJet *iJet);     
 
 // ----------member data ---------------------------
  edm::InputTag pfJetsInputTag_;
  edm::InputTag fVertexNameTag_;
  edm::InputTag fCorrJetNameData;
  edm::InputTag fCorrJetNameMC;
  edm::InputTag fUnCorrJetName;
   
  double            fJetPtMin; 
  PileupJetIdAlgo  *fPUJetIdAlgo;
  
  std::string aliasprefix_;
  std::string PFJetCorrectorL2L3_;
  std::string PFJetCorrectorL1L2L3_;
  std::string PFJetCorrectorL1FastL2L3_;
};
#endif
