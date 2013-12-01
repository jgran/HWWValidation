#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "HWWValidation/HWWBase/interface/MVAJetIdMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;



// Constructor
MVAJetIdMaker::MVAJetIdMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetscorrp4"   ).setBranchAlias( "pfjets_corr_p4"  );
  produces<vector<float> >         ( "pfjetsJEC"      ).setBranchAlias( "pfjets_JEC"	    );
  produces<vector<float> >         ( "pfjetsmvavalue" ).setBranchAlias( "pfjets_mvavalue"	);
  
  //
  fVertexNameTag_   = iConfig.getParameter<InputTag>	( "VertexName" 		);
  fCorrJetNameData  = iConfig.getParameter<InputTag>	( "CorrJetNameData"		);
  fCorrJetNameMC    = iConfig.getParameter<InputTag>	( "CorrJetNameMC"		);
  fUnCorrJetName  	= iConfig.getParameter<InputTag>	( "JetName"			);
  fJetPtMin       	= iConfig.getParameter<double>      ("JetPtMin");
  // 
  fPUJetIdAlgo    	= new PileupJetIdAlgo(iConfig); 

}

// Destructor
MVAJetIdMaker::~MVAJetIdMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void MVAJetIdMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void MVAJetIdMaker::endJob() {}

// ------------ method called to produce the data  ------------
void MVAJetIdMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){



  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  auto_ptr<vector<LorentzVector> > pfjets_corr_p4                 (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_JEC                     (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_mvavalue                (new vector<float>          );  


  //Uncorrected Jets
  Handle<PFJetCollection>       lHUCJets;
  iEvent.getByLabel(fUnCorrJetName, lHUCJets);
  PFJetCollection               lUCJets = *lHUCJets;

  //Corrected Jets
  Handle<PFJetCollection>       lHCJets;
  if(iEvent.isRealData()) iEvent.getByLabel(fCorrJetNameData  , lHCJets);
  else iEvent.getByLabel(fCorrJetNameMC  , lHCJets);
  PFJetCollection               lCJets = *lHCJets;

  // vertices    
  Handle<reco::VertexCollection> lHVertices;
  iEvent.getByLabel(fVertexNameTag_      , lHVertices); 
  VertexCollection lVertices = *lHVertices;


  //store corrected pfjets
  for(unsigned int ijet=0; ijet<lCJets.size(); ijet++){

		const PFJet     *pCJet  = &(lCJets.at(ijet));
    pfjets_corr_p4 ->push_back( LorentzVector( pCJet->p4() ) );

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
      	  if( pUCJet->pt()                               < fJetPtMin    ) continue;
		  double lJec = pCJet ->pt()/pUCJet->pt();

		  // calculate mva value only when there are good vertices 
		  // otherwise store -999
		  if( lGoodVertices.size()>0 ) 
		  {
		  	PileupJetIdentifier lPUJetId =  fPUJetIdAlgo->computeIdVariables(pCJet,lJec,&lGoodVertices[0],lGoodVertices,true);
		   	pfjets_mvavalue              	->push_back( lPUJetId.mva()              );
        pfjets_JEC                    ->push_back( lJec ); 
		  
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
		  	pfjets_mvavalue                 ->push_back( -999.              );

		  break;

	  } // lCJets
  } // lUCJets

  // 
  iEvent.put(pfjets_corr_p4,                 "pfjetscorrp4"                   );
  iEvent.put(pfjets_JEC,                  	 "pfjetsJEC"                      );
  iEvent.put(pfjets_mvavalue,             	 "pfjetsmvavalue"                 );

}

bool MVAJetIdMaker::passPFLooseId(const reco::PFJet *iJet) {
	if(iJet->energy()== 0)                                  return false;
	if(iJet->neutralHadronEnergy()/iJet->energy() > 0.99)   return false;
	if(iJet->neutralEmEnergy()/iJet->energy()     > 0.99)   return false;
	if(iJet->nConstituents() <  2)                          return false;
	if(iJet->chargedHadronEnergy()/iJet->energy() <= 0 && fabs(iJet->eta()) < 2.4 ) return false;
	if(iJet->chargedEmEnergy()/iJet->energy() >  0.99  && fabs(iJet->eta()) < 2.4 ) return false;
	if(iJet->chargedMultiplicity()            < 1      && fabs(iJet->eta()) < 2.4 ) return false;
	return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MVAJetIdMaker);
