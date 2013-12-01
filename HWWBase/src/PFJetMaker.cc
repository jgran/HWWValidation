#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HWWValidation/HWWBase/interface/PFJetMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetsp4"                               ).setBranchAlias( "pfjets_p4"                               );
  produces<vector<float> >         ( "pfjetschargedHadronE"                   ).setBranchAlias( "pfjets_chargedHadronE"                   );
  produces<vector<float> >         ( "pfjetsneutralHadronE"                   ).setBranchAlias( "pfjets_neutralHadronE"                   );
  produces<vector<float> >         ( "pfjetschargedEmE"                       ).setBranchAlias( "pfjets_chargedEmE"                       );
  produces<vector<float> >         ( "pfjetsneutralEmE"                       ).setBranchAlias( "pfjets_neutralEmE"                       );
  produces<vector<float> >         ( "pfjetsphotonE"                          ).setBranchAlias( "pfjets_photonE"                          );
  produces<vector<float> >         ( "pfjetselectronE"                        ).setBranchAlias( "pfjets_electronE"                        );
  produces<vector<float> >         ( "pfjetsmuonE"                            ).setBranchAlias( "pfjets_muonE"                            );
  produces<vector<float> >         ( "pfjetshfHadronE"                        ).setBranchAlias( "pfjets_hfHadronE"                        );
  produces<vector<float> >         ( "pfjetshfEmE"                            ).setBranchAlias( "pfjets_hfEmE"                            );
  produces<vector<int> >           ( "pfjetschargedHadronMultiplicity"        ).setBranchAlias( "pfjets_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsneutralHadronMultiplicity"        ).setBranchAlias( "pfjets_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsphotonMultiplicity"               ).setBranchAlias( "pfjets_photonMultiplicity"               );
  produces<vector<int> >           ( "pfjetselectronMultiplicity"             ).setBranchAlias( "pfjets_electronMultiplicity"             );
  produces<vector<int> >           ( "pfjetsmuonMultiplicity"                 ).setBranchAlias( "pfjets_muonMultiplicity"                 );
  produces<vector<int> >           ( "pfjetshfHadronMultiplicity"             ).setBranchAlias( "pfjets_hfHadronMultiplicity"             );
  produces<vector<int> >           ( "pfjetshfEmMultiplicity"                 ).setBranchAlias( "pfjets_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( "pfjetschargedMultiplicity"              ).setBranchAlias( "pfjets_chargedMultiplicity"              );
  produces<vector<int>   >         ( "pfjetsneutralMultiplicity"              ).setBranchAlias( "pfjets_neutralMultiplicity"              );
  produces<vector<float> >         ( "pfjetsarea"                             ).setBranchAlias( "pfjets_area"                             );

  //
  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfCandidatesTag_                  = iConfig.getParameter<InputTag>   ( "pfCandidatesTag"                  );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
}

// Destructor
PFJetMaker::~PFJetMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  auto_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  

  Handle<View<PFJet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  //get pfcandidates
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);

  for(View<PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )      );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );

    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < PFJet > > ( pfJetsHandle , idx ) );

  }

  iEvent.put(pfjets_p4                        , "pfjetsp4"                        );
  iEvent.put(pfjets_area                      , "pfjetsarea"                      );
}


//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
