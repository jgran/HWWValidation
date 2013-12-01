// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"

#include "HWWValidation/HWWBase/interface/CommonUtils.h"
#include "HWWValidation/HWWBase/interface/MCUtilities.h"
#include "HWWValidation/HWWBase/interface/BTagMaker.h" 

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

BTagMaker::BTagMaker(const edm::ParameterSet& iConfig) {
  
  // parameters from configuration

  CaloJetsTag_                       = iConfig.getParameter<edm::InputTag>("CaloJetsTag"                      );
  referenceCaloJetsTag_                  = iConfig.getParameter<edm::InputTag>("referenceCaloJetsTag"                 );
  aliasprefix_                           = iConfig.getParameter<std::string>  ("AliasPrefix"                          );
  trackCountingHighEffBJetTags_          = iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"         );

  //btagging info
  produces<vector<float> >   (aliasprefix_+"trackCountingHighEffBJetTag"         ).setBranchAlias(aliasprefix_+"_trackCountingHighEffBJetTag"	      );
  
}


BTagMaker::~BTagMaker() {}

void  BTagMaker::beginJob() {
}

void BTagMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BTagMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  Handle<View<reco::Jet> > CaloJetsHandle;
  iEvent.getByLabel(CaloJetsTag_, CaloJetsHandle);
  const View<reco::Jet> *CaloJets = CaloJetsHandle.product();
  
  Handle<View<reco::Jet> > referenceCaloJetsHandle;
  iEvent.getByLabel(referenceCaloJetsTag_, referenceCaloJetsHandle);
  const View<reco::Jet> *referenceCaloJets = referenceCaloJetsHandle.product();

  auto_ptr<vector<float> >     jets_trackCountingHighEffBJetTag          (new vector<float>  );

    edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighEffBJetTags;
  iEvent.getByLabel(trackCountingHighEffBJetTags_, trackCountingHighEffBJetTags);
  
  
  for( edm::View<reco::Jet>::const_iterator it =  CaloJets->begin();
	   it != CaloJets->end(); it++ ) {
    edm::RefToBase<reco::Jet> jetRef   = getReferenceJetRef(referenceCaloJets, &(*it));

    jets_trackCountingHighEffBJetTag         ->push_back( CommonUtils::isinf((*trackCountingHighEffBJetTags)[jetRef]) 
							  ? -9999.  : (*trackCountingHighEffBJetTags)[jetRef]		); 
}

  iEvent.put(jets_trackCountingHighEffBJetTag          ,aliasprefix_+"trackCountingHighEffBJetTag"         );	  

}

//---------------------------------------------------------------------------------------
edm::RefToBase<reco::Jet> BTagMaker::getReferenceJetRef(const edm::View<reco::Jet>* refJets, const reco::Jet* jet) {

  double mindR = 0.01;
  edm::RefToBase<reco::Jet> retRef = edm::RefToBase<reco::Jet>();
  for(edm::View<reco::Jet>::const_iterator it = refJets->begin();  
      it!= refJets->end(); it++) {

    double dR = ROOT::Math::VectorUtil::DeltaR(it->p4(), jet->p4());
    if(dR < mindR) {
      mindR = dR;
      unsigned int idx = it - refJets->begin();
      retRef = refJets->refAt(idx);
    }
  }

  if (mindR == 0.01)
       std::cout << "\n didn't find a match!\n";

  if(!retRef.isNonnull())
    throw cms::Exception("Reference jet not found in BTagMaker");
  return retRef;
					  
}
//define this as a plug-in
DEFINE_FWK_MODULE(BTagMaker);

