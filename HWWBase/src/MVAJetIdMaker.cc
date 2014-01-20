#include "HWWValidation/HWWBase/interface/MVAJetIdMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;

MVAJetIdMaker::MVAJetIdMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  PFJetCollection_     = iCollector.consumes<reco::PFJetCollection> (iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  CorrPFJetCollection_ = iCollector.consumes<reco::PFJetCollection> (iConfig.getParameter<edm::InputTag>("corrPFJetsInputTag"));
  thePVCollection_     = iCollector.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));

  fPUJetIdAlgo = new PileupJetIdAlgo(iConfig);

}

bool passPFLooseId(const reco::PFJet *iJet) {
	if(iJet->energy()== 0)                                  return false;
	if(iJet->neutralHadronEnergy()/iJet->energy() > 0.99)   return false;
	if(iJet->neutralEmEnergy()/iJet->energy()     > 0.99)   return false;
	if(iJet->nConstituents() <  2)                          return false;
	if(iJet->chargedHadronEnergy()/iJet->energy() <= 0 && fabs(iJet->eta()) < 2.4 ) return false;
	if(iJet->chargedEmEnergy()/iJet->energy() >  0.99  && fabs(iJet->eta()) < 2.4 ) return false;
	if(iJet->chargedMultiplicity()            < 1      && fabs(iJet->eta()) < 2.4 ) return false;
	return true;
}
          
void MVAJetIdMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

  HWWVal::Load_pfjets_corr_p4();
  HWWVal::Load_pfjets_mvavalue();
  HWWVal::Load_pfjets_JEC();

  //Uncorrected Jets
  Handle<PFJetCollection>       lHUCJets;
  iEvent.getByToken(PFJetCollection_, lHUCJets);
  PFJetCollection               lUCJets = *lHUCJets;

  //Corrected Jets
  Handle<PFJetCollection>       lHCJets;
  iEvent.getByToken(CorrPFJetCollection_  , lHCJets);
  PFJetCollection               lCJets = *lHCJets;

  // vertices    
  Handle<reco::VertexCollection> lHVertices;
  iEvent.getByToken(thePVCollection_, lHVertices); 
  VertexCollection lVertices = *lHVertices;

  //store corrected pfjets
  for(unsigned int ijet=0; ijet<lCJets.size(); ijet++){

		const PFJet     *pCJet  = &(lCJets.at(ijet));
    HWWVal::pfjets_corr_p4() .push_back( LorentzVector( pCJet->p4() ) );

  }


  // select good vertices 
  // make new collection to put into computeIdVariables(...)
  VertexCollection lGoodVertices;
  for(int ivtx    = 0; ivtx < (int)lVertices.size(); ivtx++)
  {
	  const Vertex       *vtx = &(lVertices.at(ivtx));
	  if( vtx->isFake()               		)  continue;
	  if( vtx->ndof()<=4              		)  continue;
	  if( vtx->position().Rho()>2.0   		)  continue;
	  if( fabs(vtx->position().Z())>24.0    )  continue;
	  lGoodVertices.push_back(*vtx);
  }

  // loop over jets 
  for(int i0   = 0; i0 < (int) lUCJets.size(); i0++) {   // uncorrected jets collection                                           
	  const PFJet       *pUCJet = &(lUCJets.at(i0));
	  for(int i1 = 0; i1 < (int) lCJets.size(); i1++) {   // corrected jets collection                                         
		  const PFJet     *pCJet  = &(lCJets.at(i1));
		  if( pUCJet->jetArea() != pCJet->jetArea()                  	) continue;
		  if( fabs(pUCJet->eta() - pCJet->eta())         > 0.001         ) continue;
      	  if( pUCJet->pt()                               < 0.0    ) continue;
		  double lJec = pCJet ->pt()/pUCJet->pt();

		  // calculate mva value only when there are good vertices 
		  // otherwise store -999
		  if( lGoodVertices.size()>0 ) 
		  {
		  	PileupJetIdentifier lPUJetId =  fPUJetIdAlgo->computeIdVariables(pCJet,lJec,&lGoodVertices[0],lGoodVertices,true);
		   	HWWVal::pfjets_mvavalue() .push_back( lPUJetId.mva()              );
        HWWVal::pfjets_JEC() .push_back( lJec ); 
		  
			// print out MVA inputs 
			if(0)
			{
				std::cout << setprecision(5)
					<< "Debug Jet MVA : "
				 	<< iEvent.id() 			<< " : "
					<< lPUJetId.nvtx()      << " "
					<< pCJet->pt()         	<< " "
					<< lPUJetId.jetEta()    << " "
					<< lPUJetId.jetPhi()    << " "
					<< lPUJetId.d0()        << " "
					<< lPUJetId.dZ()        << " "
					<< lPUJetId.beta()      << " "
					<< lPUJetId.betaStar()  << " "
					<< lPUJetId.nCharged()  << " "
					<< lPUJetId.nNeutrals() << " "
					<< lPUJetId.dRMean()    << " "
					<< lPUJetId.frac01()    << " "
					<< lPUJetId.frac02()    << " "
					<< lPUJetId.frac03()    << " "
					<< lPUJetId.frac04()    << " "
					<< lPUJetId.frac05()
					<< " === : === "; 
				cout << lPUJetId.mva() << endl;
			}
		  }
		  else             
		  	HWWVal::pfjets_mvavalue() .push_back( -999. );

		  break;

	  } // lCJets
  } // lUCJets
}
