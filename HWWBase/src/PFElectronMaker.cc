// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "HWWValidation/HWWBase/interface/PFElectronMaker.h"

//
typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
typedef edm::ValueMap<reco::PFCandidatePtr> PFCandMap;

//
using namespace reco;
using namespace edm;
using namespace std;

//
PFElectronMaker::PFElectronMaker(const ParameterSet& iConfig) {
      
  //
  pfCandidatesTag_    = iConfig.getParameter<InputTag> ("pfCandidatesTag");

  //
  produces<vector<LorentzVector> >  ("pfelsp4"                  ).setBranchAlias("pfels_p4"                 );

} //

PFElectronMaker::~PFElectronMaker() {}
void  PFElectronMaker::beginRun( Run&, const EventSetup& es ) {}
void PFElectronMaker::beginJob() {}
void PFElectronMaker::endJob() {}

//
void PFElectronMaker::produce( Event& iEvent, const EventSetup& iSetup ) {

  auto_ptr<vector<LorentzVector> > pfels_p4                  (new vector<LorentzVector> );

  //
  Handle<PFCandMap > pfCandidatesHandle;
  iEvent.getByLabel( pfCandidatesTag_, pfCandidatesHandle );
  const ValueMap<reco::PFCandidatePtr> *pfCandidates  = pfCandidatesHandle.product();
  
  //
  PFCandMap::const_iterator pf_pit = pfCandidates->begin();
  unsigned int nC = pf_pit.size();
  for( unsigned int iC = 0; iC < nC; ++iC ) {

   const PFCandidatePtr& pf_it = pf_pit[iC];
   if ( pf_it.isNull() ) continue;
   int pfflags = 0;

   for( unsigned int i = 0; i < 17; i++ ) {
     if(pf_it->flag((PFCandidate::Flags)i)) pfflags |= (1<<i);
   }
  
   pfels_p4               ->push_back(LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p()) );
    
  } //
  

  //
  iEvent.put( pfels_p4                  , "pfelsp4"                 );
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFElectronMaker);

