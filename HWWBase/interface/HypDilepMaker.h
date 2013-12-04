#ifndef HYPDILEPMAKER_H
#define HYPDILEPMAKER_H

//
// class decleration
//

class HypDilepMaker : public edm::EDProducer {
public:
  
    

  explicit HypDilepMaker (const edm::ParameterSet&);
  ~HypDilepMaker();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool testJetForLeptons(const math::XYZTLorentzVectorF& jetP4, const math::XYZTLorentzVectorF& lepP4);
   
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;
  edm::InputTag jetsInputTag;
  double        hypJetMaxEtaCut;
  double        hypJetMinPtCut;
  double        tightptcut;
  double        looseptcut;
    
	std::string aliasprefix_;
};


#endif
