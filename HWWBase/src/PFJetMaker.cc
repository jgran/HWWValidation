#include "HWWValidation/HWWBase/interface/PFJetMaker.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;

PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  PFJetCollection_  = iCollector.consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  BJetTags_         = iCollector.consumes<reco::JetFloatAssociation::Container>(iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"));

}

void PFJetMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  HWWVal::Load_pfjets_p4();
  HWWVal::Load_pfjets_area();
  HWWVal::Load_pfjets_trackCountingHighEffBJetTag();

  edm::Handle<edm::View<reco::Jet> >    pfjetsHandle;
  iEvent.getByToken(PFJetCollection_,  pfjetsHandle	);

  edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighEffBJetTags;
  iEvent.getByToken(BJetTags_, trackCountingHighEffBJetTags);
    
  for(edm::View<reco::Jet>::const_iterator jet_it = pfjetsHandle->begin(); jet_it != pfjetsHandle->end(); jet_it++) {
       
    if(jet_it->pt() <= 0.0 ) continue;

    HWWVal::pfjets_p4()     .push_back( LorentzVector( jet_it->p4() )      );
    HWWVal::pfjets_area()   .push_back(jet_it->jetArea()                   );

    unsigned int idx = jet_it-pfjetsHandle->begin();
    edm::RefToBase<reco::Jet> jetRef = pfjetsHandle->refAt(idx);

    HWWVal::pfjets_trackCountingHighEffBJetTag().push_back( (*trackCountingHighEffBJetTags)[jetRef]	); 

  }

}
