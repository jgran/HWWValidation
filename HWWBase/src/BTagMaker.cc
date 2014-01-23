#include "HWWValidation/HWWBase/interface/CommonUtils.h"
#include "HWWValidation/HWWBase/interface/BTagMaker.h" 
#include "HWWValidation/HWWBase/interface/HWW.h" 

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

BTagMaker::BTagMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector) {
  
  PFJets_   = iCollector.consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  BJetTags_ = iCollector.consumes<reco::JetFloatAssociation::Container>(iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"));

}

edm::RefToBase<reco::Jet> getReferenceJetRef(const edm::View<reco::Jet>* refJets, const reco::Jet* jet) {

  double mindR = 0.01;
  edm::RefToBase<reco::Jet> retRef = edm::RefToBase<reco::Jet>();

  for(edm::View<reco::Jet>::const_iterator it = refJets->begin(); it!= refJets->end(); it++) {

    double dR = ROOT::Math::VectorUtil::DeltaR(it->p4(), jet->p4());

    if(dR < mindR) {
      mindR = dR;
      unsigned int idx = it - refJets->begin();
      retRef = refJets->refAt(idx);
    }

  }

  if (mindR == 0.01) std::cout << "\n didn't find a match!\n";

  if(!retRef.isNonnull())
    throw cms::Exception("Reference jet not found in BTagMaker");
  return retRef;
					  
}

void BTagMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_pfjets_trackCountingHighEffBJetTag();

  Handle<View<reco::Jet> > CaloJetsHandle;
  iEvent.getByToken(PFJets_, CaloJetsHandle);
  const View<reco::Jet> *CaloJets = CaloJetsHandle.product();
  
  Handle<View<reco::Jet> > referenceCaloJetsHandle;
  iEvent.getByToken(PFJets_, referenceCaloJetsHandle);
  const View<reco::Jet> *referenceCaloJets = referenceCaloJetsHandle.product();

  edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighEffBJetTags;
  iEvent.getByToken(BJetTags_, trackCountingHighEffBJetTags);
  
  for( edm::View<reco::Jet>::const_iterator it =  CaloJets->begin(); it != CaloJets->end(); it++ ) {

    edm::RefToBase<reco::Jet> jetRef   = getReferenceJetRef(referenceCaloJets, &(*it));

    HWWVal::pfjets_trackCountingHighEffBJetTag()  .push_back( CommonUtils::isinf((*trackCountingHighEffBJetTags)[jetRef]) 
							  ? -9999.  : (*trackCountingHighEffBJetTags)[jetRef]		); 
  }
}
