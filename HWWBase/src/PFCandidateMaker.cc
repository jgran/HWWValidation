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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "HWWValidation/HWWBase/interface/PFCandidateMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//
PFCandidateMaker::PFCandidateMaker(const edm::ParameterSet& iConfig) {

     pfCandidatesTag_		= iConfig.getParameter<InputTag>	("pfCandidatesTag");
     pfElectronsTag_		= iConfig.getParameter<InputTag>	("pfElectronsTag");
     tracksInputTag_            = iConfig.getParameter<InputTag>        ("tracksInputTag");
     vertexInputTag_            = iConfig.getParameter<InputTag>        ("vertexInputTag");
     minDR_electron_            = iConfig.getParameter<double>          ("minDRelectron");

     produces<vector<LorentzVector>	> ("pfcandsp4"              ).setBranchAlias("pfcands_p4"			          );
     produces<vector<int>	>           ("pfcandscharge"		      ).setBranchAlias("pfcands_charge"		        );
     produces<vector<int> >           ("pfcandsparticleId"		  ).setBranchAlias("pfcands_particleId"		    );
     produces<vector<int>	>           ("pfcandspfelsidx"		    ).setBranchAlias("pfcands_pfelsidx"		      );
     produces<vector<int>	>           ("pfcandstrkidx"		      ).setBranchAlias("pfcands_trkidx"		        );
     produces<vector<int> >          ("pfcandsvtxidx"             ).setBranchAlias("pfcands_vtxidx"       );

    // for matching to vertices using the "PFNoPileup" method
    // hint: it is just track vertex association 
    pfPileUpAlgo_ = new PFPileUpAlgo();

}

PFCandidateMaker::~PFCandidateMaker() 
{
    if (pfPileUpAlgo_ != 0) delete pfPileUpAlgo_;
}

void  PFCandidateMaker::beginRun(edm::Run&, const edm::EventSetup& es) {}
void PFCandidateMaker::beginJob() {}
void PFCandidateMaker::endJob()   {}

// ------------ method called to produce the data  ------------
void PFCandidateMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<vector<LorentzVector> >	pfcands_p4		            (new vector<LorentzVector>  );
     auto_ptr<vector<int> >		        pfcands_charge		        (new vector<int>		        );
     auto_ptr<vector<int> >		        pfcands_particleId	      (new vector<int>		        );
     auto_ptr<vector<int> >	          pfcands_trkidx		        (new vector<int>        	  );
     auto_ptr<vector<int> >	          pfcands_pfelsidx	        (new vector<int>        	  );
     auto_ptr<vector<int> >           pfcands_vtxidx            (new vector<int>              );
    
     //get pfcandidates
     Handle<PFCandidateCollection> pfCandidatesHandle;
     iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
     pfCandidates  = pfCandidatesHandle.product();

     //get pfelectrons
     typedef edm::ValueMap<reco::PFCandidatePtr> PFCandMap;
     Handle<PFCandMap> pfElectronsHandle;
     iEvent.getByLabel(pfElectronsTag_, pfElectronsHandle);
     const PFCandMap *pfElectrons  = pfElectronsHandle.product();

     // get tracks
     Handle<reco::TrackCollection>  track_h;
     iEvent.getByLabel(tracksInputTag_, track_h);

     // get vertices
     Handle<reco::VertexCollection> vertex_h;
     iEvent.getByLabel(vertexInputTag_, vertex_h);
     const reco::VertexCollection *vertices = vertex_h.product();

      int iCand = 0;
      for( PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
  
        //
	      int pfflags = 0;
	      for( unsigned int i = 0; i < 17; i++ ) {
          if(pf_it->flag((PFCandidate::Flags)i)) pfflags |= (1<<i);
	      }

        //
        iCand++;

  	    pfcands_p4			          ->push_back( LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p())	  );
  	    pfcands_charge		        ->push_back( pf_it->charge()						                                        );
  	    pfcands_particleId		    ->push_back( pf_it->translateTypeToPdgId(pf_it->particleId())	                  );

          //for charged pfcandidates, find corresponding track index
          //here we take the track directly from PFCandidate::trackRef()
          if( pf_it->charge() != 0 ){

            // match pf cand to a vertex
            // using the no pileup algo
            pfcands_vtxidx->push_back(pfPileUpAlgo_->chargedHadronVertex(*vertices, *pf_it));

            reco::TrackRef pftrack = pf_it->trackRef();

            int trkidx = 0;
            bool foundTrack = false;

            reco::TrackCollection::const_iterator tracks_end = track_h->end();

            for (reco::TrackCollection::const_iterator itrk = track_h->begin(); itrk != tracks_end; ++itrk) {
              
              reco::TrackRef trkref( track_h , itrk - track_h->begin() );

              if( pftrack.key() == trkref.key() ){
              
                //sanity check
                float dpt = pftrack->pt() - trkref->pt();
                if( fabs( dpt ) > 0.1 ){
                  cout << "Warning: pfcandidate track pt - matched track pt = " << dpt << ", possible mismatch" << endl;
                }

                //found corresponding track
                pfcands_trkidx->push_back( trkidx );
                foundTrack = true;
                break;
              }
              
              ++trkidx;
            }

            if( !foundTrack ){
              //no matched track found, set trkidx to -1
              pfcands_trkidx->push_back(-1);
            }
            
          }else{
            //neutral particle, set trkidx to -2
            pfcands_trkidx->push_back(-2);
            // no vertex match
            pfcands_vtxidx->push_back(-2);
          }

          //find corresponding PFMuon index

          //find corresponding PFElectron index
          if (pf_it->particleId() == PFCandidate::e){
	    int index = -1; 
	    if (pf_it->gsfTrackRef().isNonnull()) {
	      int pfGsfTkId = pf_it->gsfTrackRef().key();
	      unsigned int elsIndex = 0;
	      PFCandMap::const_iterator el_pit = pfElectrons->begin();
	      unsigned int nE = el_pit.size();
	      for(unsigned int iE = 0; iE < nE; ++iE){
		const PFCandidatePtr& el_it = el_pit[iE];
		if (el_it.isNull()) continue;
		int elGsfTkId = -1;
		if (el_it->gsfTrackRef().isNonnull()) elGsfTkId = el_it->gsfTrackRef().key();
		if (elGsfTkId==pfGsfTkId) {
		  index = elsIndex;
		  break;
		}
		elsIndex++;
	      }            
	    }
	    pfcands_pfelsidx->push_back(index);
	  } else {
            pfcands_pfelsidx->push_back(-2);
          }
 
          
     }//loop over candidate collection

     iEvent.put(pfcands_p4,			            "pfcandsp4"		          );
     iEvent.put(pfcands_charge,			        "pfcandscharge"		      );
     iEvent.put(pfcands_particleId,		      "pfcandsparticleId"	    );
     iEvent.put(pfcands_trkidx,			        "pfcandstrkidx"		      );
     iEvent.put(pfcands_vtxidx,               "pfcandsvtxidx"         );
     iEvent.put(pfcands_pfelsidx,		        "pfcandspfelsidx"	      );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateMaker);

