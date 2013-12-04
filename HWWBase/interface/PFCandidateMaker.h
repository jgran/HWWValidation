#ifndef PFCANDIDATEMAKER_H
#define PFCANDIDATEMAKER_H

//
// class decleration
//

class PFCandidateMaker : public edm::EDProducer {
public:
     explicit PFCandidateMaker (const edm::ParameterSet&);
     ~PFCandidateMaker();

private:
  //  virtual void beginJob() ;
  virtual void beginJob() ;
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  float getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins);
  
  // ----------member data ---------------------------
  double minDR_electron_;
  edm::InputTag pfElectronsTag_;
  edm::InputTag pfCandidatesTag_;
  edm::InputTag tracksInputTag_;
  edm::InputTag vertexInputTag_;

  const reco::PFCandidateCollection *pfCandidates;

    PFPileUpAlgo *pfPileUpAlgo_;

};

#endif

