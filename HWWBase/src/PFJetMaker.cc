#include "HWWValidation/HWWBase/interface/PFJetMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;

PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  PFJetCollection_       = iCollector.consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));

}


bool sortPFJetsByPt(reco::PFJet jet1, reco::PFJet jet2) {
  return jet1.pt() > jet2.pt();
}


void PFJetMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

  HWWVal::Load_pfjets_p4();
  HWWVal::Load_pfjets_area();

  edm::Handle<reco::PFJetCollection>    pfjetsHandle;
  iEvent.getByToken(PFJetCollection_,  pfjetsHandle	);

  std::auto_ptr<reco::PFJetCollection> outPFJetCollection (new reco::PFJetCollection        );
    
  for(reco::PFJetCollection::const_iterator jet_it = pfjetsHandle->begin();
      jet_it != pfjetsHandle->end(); jet_it++) {
       
    if(jet_it->pt() > 0.0 )
      outPFJetCollection->push_back(*jet_it);
  }
  
  std::sort(outPFJetCollection->begin(),  outPFJetCollection->end(), sortPFJetsByPt);

  for(reco::PFJetCollection::const_iterator pfjet_it = outPFJetCollection->begin(); pfjet_it != outPFJetCollection->end(); pfjet_it++) {

    HWWVal::pfjets_p4()     .push_back( LorentzVector( pfjet_it->p4() )      );
    HWWVal::pfjets_area()   .push_back(pfjet_it->jetArea()                   );

  }

}
