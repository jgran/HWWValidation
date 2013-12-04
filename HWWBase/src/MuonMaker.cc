#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "HWWValidation/HWWBase/interface/MuonMaker.h"


//////////////
// typedefs //
//////////////

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;


////////////////
// namespaces //
////////////////

using namespace std;
using namespace reco;
using namespace edm;


/////////////////
// Constructor //
/////////////////

MuonMaker::MuonMaker( const ParameterSet& iConfig ) {

  /////////////////////////////
  // Branch & Alias prefixes //
  /////////////////////////////

  aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
  branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


  //////////////////////
  // Input Parameters //
  //////////////////////

  muonsInputTag    = iConfig.getParameter<InputTag> ("muonsInputTag"   );
  beamSpot_tag_    = iConfig.getParameter<InputTag> ("beamSpotTag"     );
  pfCandsInputTag  = iConfig.getParameter<InputTag> ("pfCandsInputTag" );
  vtxInputTag      = iConfig.getParameter<InputTag> ("vtxInputTag"     );
  tevMuonsName     = iConfig.getParameter<string>   ("tevMuonsName"    );
  src_             = iConfig.getParameter<InputTag> ("cosmicCompat"    ); 
  showerTag_       = iConfig.getParameter<InputTag> ("muonShower"      ); 

  produces<vector<float> >          ( branchprefix_ + "stad0"                     ).setBranchAlias( aliasprefix_ + "_sta_d0"            ); // d0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0"                     ).setBranchAlias( aliasprefix_ + "_sta_z0"            ); // z0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "stad0Err"                  ).setBranchAlias( aliasprefix_ + "_sta_d0Err"         ); // d0Err from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0Err"                  ).setBranchAlias( aliasprefix_ + "_sta_z0Err"         ); // z0Err from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "stad0corr"                 ).setBranchAlias( aliasprefix_ + "_sta_d0corr"        ); // Beamspot corrected d0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0corr"                 ).setBranchAlias( aliasprefix_ + "_sta_z0corr"        ); // Beamspot corrected z0 from STA fit, if it exists
  produces<vector<LorentzVector> >  ( branchprefix_ + "stap4"                     ).setBranchAlias( aliasprefix_ + "_sta_p4"            ); // 




  produces<vector<float> >          ( branchprefix_ + "gfitchi2"                  ).setBranchAlias( aliasprefix_ + "_gfit_chi2"          ); // chi2 of the global muon fit 
  produces<vector<float> >          ( branchprefix_ + "gfitndof"                  ).setBranchAlias( aliasprefix_ + "_gfit_ndof"          ); // number of degree of freedom of the global muon fit 
  produces<vector<int> >            ( branchprefix_ + "gfitvalidSTAHits"          ).setBranchAlias( aliasprefix_ + "_gfit_validSTAHits"  ); // number of hits in the stand alone fit that made it into the gfit

  /////////////
  // Quality //
  /////////////

  produces<vector<float> >          ( branchprefix_ + "trkKink"                   ).setBranchAlias( aliasprefix_ + "_trkKink"             );  // Muon Quality - trkKink

  //////////
  // Muon //
  //////////

  produces<vector<int> >            ( branchprefix_ + "type"                      ).setBranchAlias( aliasprefix_ + "_type"                   ); // type
  produces<vector<int> >            ( branchprefix_ + "goodmask"                  ).setBranchAlias( aliasprefix_ + "_goodmask"               ); // good mask
  produces<vector<int> >            ( branchprefix_ + "charge"                    ).setBranchAlias( aliasprefix_ + "_charge"                 ); // charge from muon object             
  produces<vector<int> >            ( branchprefix_ + "nmatches"                  ).setBranchAlias( aliasprefix_ + "_nmatches"               ); // number of stations with matched segments 
  produces<vector<float> >          ( branchprefix_ + "caloCompatibility"         ).setBranchAlias( aliasprefix_ + "_caloCompatibility"      ); // calo compatibility variable
  produces<vector<float> >          ( branchprefix_ + "segmCompatibility"         ).setBranchAlias( aliasprefix_ + "_segmCompatibility"      );
  produces<vector<LorentzVector> >  ( branchprefix_ + "p4"                        ).setBranchAlias( aliasprefix_ + "_p4"                     ); // candidate p4->this can either be gfit p4, tracker p4 or STA p4 (only for STA muoons)     
  produces<vector<int> >            ( branchprefix_ + "numberOfMatchedStations"   ).setBranchAlias( aliasprefix_ + "_numberOfMatchedStations"); // number of muon stations with muon segements used in the fit

  ////////
  // ID //
  ////////

  produces<vector<int> >            ( branchprefix_ + "pidTMLastStationTight"     ).setBranchAlias( aliasprefix_ + "_pid_TMLastStationTight"     ); // tight tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >            ( branchprefix_ + "pidPFMuon" ).setBranchAlias( aliasprefix_ + "_pid_PFMuon" ); // is particle flow muon

  ////////////
  // Energy //
  ////////////

  produces<vector<float> >          ( branchprefix_ + "eem"                       ).setBranchAlias( aliasprefix_ + "_e_em"                ); // energy in crossed ECAL crystalls 
  produces<vector<float> >          ( branchprefix_ + "ehad"                      ).setBranchAlias( aliasprefix_ + "_e_had"               ); // energy in crossed HCAL towers 
  produces<vector<float> >          ( branchprefix_ + "eho"                       ).setBranchAlias( aliasprefix_ + "_e_ho"                ); // energy in crossed HO towers 
  produces<vector<float> >          ( branchprefix_ + "eemS9"                     ).setBranchAlias( aliasprefix_ + "_e_emS9"              ); // energy in 3x3 ECAL crystall shape 
  produces<vector<float> >          ( branchprefix_ + "ehadS9"                    ).setBranchAlias( aliasprefix_ + "_e_hadS9"             ); // energy in 3x3 HCAL towers 
  produces<vector<float> >          ( branchprefix_ + "ehoS9"                     ).setBranchAlias( aliasprefix_ + "_e_hoS9"              ); // energy in 3x3 HO towers 


  ///////////////
  // Isolation //
  ///////////////
                                    
  produces<vector<float> >          ( branchprefix_ + "isoecalvetoDep"            ).setBranchAlias( aliasprefix_ + "_iso_ecalvetoDep"     ); // sumEt in the veto cone, ecal
  produces<vector<float> >          ( branchprefix_ + "isohcalvetoDep"            ).setBranchAlias( aliasprefix_ + "_iso_hcalvetoDep"     ); // sumPt in the veto cone, hcal
  produces<vector<float> >          ( branchprefix_ + "iso03sumPt"                ).setBranchAlias( aliasprefix_ + "_iso03_sumPt"         ); // sum of track Pt for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso03emEt"                 ).setBranchAlias( aliasprefix_ + "_iso03_emEt"          ); // sum of ecal Et for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso03hadEt"                ).setBranchAlias( aliasprefix_ + "_iso03_hadEt"         ); // sum of hcal Et for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso05sumPt"                ).setBranchAlias( aliasprefix_ + "_iso05_sumPt"         ); // sum of track Pt for cone of 0.5 
  produces<vector<float> >          ( branchprefix_ + "iso05emEt"                 ).setBranchAlias( aliasprefix_ + "_iso05_emEt"          ); // sum of ecal Et for cone of 0.5 
  produces<vector<float> >          ( branchprefix_ + "iso05hadEt"                ).setBranchAlias( aliasprefix_ + "_iso05_hadEt"         ); // sum of hcal Et for cone of 0.5 

  ////////////
  // Tracks //
  ////////////

  produces<vector<LorentzVector> >  ( branchprefix_ + "trkp4"                     ).setBranchAlias( aliasprefix_ + "_trk_p4"              ); // track p4            
  produces<vector<LorentzVector> >  ( branchprefix_ + "vertexp4"                  ).setBranchAlias( aliasprefix_ + "_vertex_p4"           ); // from the silicon fit
  produces<vector<int>   >          ( branchprefix_ + "trkidx"                    ).setBranchAlias( aliasprefix_ + "_trkidx"              ); // track index matched to muon
  produces<vector<float> >          ( branchprefix_ + "d0"                        ).setBranchAlias( aliasprefix_ + "_d0"                  ); // impact parameter at the point of closest approach  using the tracker fit
  produces<vector<float> >          ( branchprefix_ + "z0"                        ).setBranchAlias( aliasprefix_ + "_z0"                  ); // z position of the point of closest approach. From the si track    
  produces<vector<float> >          ( branchprefix_ + "d0corr"                    ).setBranchAlias( aliasprefix_ + "_d0corr"              ); // corrected impact parameter at the point of closest approach. From si track  
  produces<vector<float> >          ( branchprefix_ + "z0corr"                    ).setBranchAlias( aliasprefix_ + "_z0corr"              ); // corrected z position of the point of closest approach. From si track    
  produces<vector<float> >          ( branchprefix_ + "chi2"                      ).setBranchAlias( aliasprefix_ + "_chi2"                ); // chi2 of the silicon tracker fit      
  produces<vector<float> >          ( branchprefix_ + "ndof"                      ).setBranchAlias( aliasprefix_ + "_ndof"                ); // number of degrees of freedom of the si tracker fit    
  produces<vector<int> >            ( branchprefix_ + "validHits"                 ).setBranchAlias( aliasprefix_ + "_validHits"           ); // number of used hits in the sitracker fit      
  produces<vector<float> >          ( branchprefix_ + "d0Err"                     ).setBranchAlias( aliasprefix_ + "_d0Err"               ); // error on the impact parameter, si track fit      
  produces<vector<float> >          ( branchprefix_ + "z0Err"                     ).setBranchAlias( aliasprefix_ + "_z0Err"               ); // error on z position of the point of closest approach, si track fit  
  produces<vector<float> >          ( branchprefix_ + "ptErr"                     ).setBranchAlias( aliasprefix_ + "_ptErr"               ); // si track Pt error          
  produces<vector<float> >          ( branchprefix_ + "etaErr"                    ).setBranchAlias( aliasprefix_ + "_etaErr"              ); // si track eta error          
  produces<vector<float> >          ( branchprefix_ + "phiErr"                    ).setBranchAlias( aliasprefix_ + "_phiErr"              ); // si track phi error          


  ////////
  // PF //
  ////////

  produces<vector<float> >          ( branchprefix_ + "isoR03pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoR03_pf_ChargedHadronPt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoR03_pf_ChargedParticlePt" );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoR03_pf_NeutralHadronEt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PhotonEt"          );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoR03_pf_NeutralHadronEtHighThreshold");
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PhotonEtHighThreshold"       );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PUPt"              );

  ///////////
  // IP 3D //
  ///////////
  
  produces<vector<float> >          ( branchprefix_ + "ip3d"                      ).setBranchAlias( aliasprefix_ + "_ip3d"                ); // Ip3d from standard vertex
  produces<vector<float> >          ( branchprefix_ + "ip3derr"                   ).setBranchAlias( aliasprefix_ + "_ip3derr"             ); // Ip3d error from standard vertex

} // end Constructor

void MuonMaker::beginJob () {}  // method called once each job just before starting event loop
void MuonMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void MuonMaker::produce(Event& iEvent, const EventSetup& iSetup) {

  /////////
  // STA //
  /////////

  auto_ptr<vector<float> >         vector_mus_sta_d0                      ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0                      ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_d0Err                   ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0Err                   ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_d0corr                  ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0corr                  ( new vector<float>   );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_p4                      ( new vector<LorentzVector>  );

  ////////////  
  // Global //
  ////////////

  auto_ptr<vector<float> >         vector_mus_gfit_chi2                   ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_ndof                   ( new vector<float>         );
  auto_ptr<vector<int> >           vector_mus_gfit_validSTAHits           ( new vector<int>           );

  /////////////
  // Quality //
  /////////////

  auto_ptr<vector<float> >  vector_mus_trkKink             ( new vector<float>  );

  ///////////
  // Muons //
  ///////////

  auto_ptr<vector<int> >           vector_mus_type                    ( new vector<int>           );        
  auto_ptr<vector<int> >           vector_mus_goodmask                ( new vector<int>           );        
  auto_ptr<vector<int> >           vector_mus_nmatches                ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_charge                  ( new vector<int>           );        
  auto_ptr<vector<float> >         vector_mus_caloCompatibility       ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_segmCompatibility       ( new vector<float>         );
  auto_ptr<vector<LorentzVector> > vector_mus_p4                      ( new vector<LorentzVector> );
  auto_ptr<vector<int> >           vector_mus_numberOfMatchedStations ( new vector<int>           );

  ////////
  // ID //
  ////////

  auto_ptr<vector<int> >           vector_mus_pid_TMLastStationTight      ( new vector<int>     );
  auto_ptr<vector<int> >           vector_mus_pid_PFMuon                  ( new vector<int>     );

  ////////////
  // Energy //
  ////////////

  auto_ptr<vector<float> >         vector_mus_e_em                ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_had               ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_ho                ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_emS9              ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_hadS9             ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_hoS9              ( new vector<float>          );

  ///////////////
  // Isolation //
  ///////////////

  auto_ptr<vector<float> >         vector_mus_iso_ecalvetoDep     ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso_hcalvetoDep     ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_sumPt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_emEt          ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_hadEt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_sumPt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_emEt          ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_hadEt         ( new vector<float>          );

  ////////////
  // Tracks //
  ////////////

  auto_ptr<vector<LorentzVector> > vector_mus_trk_p4              ( new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> > vector_mus_vertex_p4           ( new vector<LorentzVector> );
  auto_ptr<vector<int>   >         vector_mus_trkidx              ( new vector<int>            );
  auto_ptr<vector<float> >         vector_mus_d0                  ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0                  ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_d0corr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0corr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_chi2                ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_ndof                ( new vector<float>          );      
  auto_ptr<vector<int> >           vector_mus_validHits           ( new vector<int>            );        
  auto_ptr<vector<float> >         vector_mus_d0Err               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0Err               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_ptErr               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_etaErr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_phiErr              ( new vector<float>          );      

  ////////
  // PF //
  ////////

  auto_ptr< vector<float> >         vector_mus_isoR03_pf_ChargedHadronPt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_ChargedParticlePt               ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_NeutralHadronEt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_PhotonEt                        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_PUPt                            ( new vector<float>   );

  ///////////
  // IP 3D //
  ///////////

  auto_ptr<vector<float> >         vector_mus_ip3d                        ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_ip3derr                     ( new vector<float>   );



////////////////////////////
// --- Fill Muon Data --- //
////////////////////////////

  ///////////////
  // Get Muons //
  ///////////////

  Handle<View<Muon> > muon_h;
  iEvent.getByLabel( muonsInputTag , muon_h );

  
  /////////////////////////////////
  // Get Muon Shower Information //
  /////////////////////////////////

  Handle<ValueMap<MuonShower> > showerMap;
  iEvent.getByLabel( showerTag_ , showerMap );


  ////////////////////////////////
  // Get pf-Muon from reco-muon //
  ////////////////////////////////

  Handle<ValueMap<PFCandidatePtr> > pfMap;
  iEvent.getByLabel( "particleFlow"  , muonsInputTag.label(), pfMap );


  //////////////////
  // Get Vertices //
  //////////////////

  iEvent.getByLabel( vtxInputTag , vertexHandle );  


  ///////////////////////
  // Get PF Candidates //
  ///////////////////////

  iEvent.getByLabel( pfCandsInputTag , pfCand_h );


  //////////////
  // Beamspot //
  //////////////

  Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByLabel(beamSpot_tag_, beamspot_h);
  const reco::BeamSpot &beamSpotreco = *(beamspot_h.product());
  //const Point beamSpot = beamSpotreco.isValid() ? Point(beamSpotreco.x0(), beamSpotreco.y0(), beamSpotreco.z0()) : Point(0,0,0);
  const Point beamSpot = Point(beamSpotreco.x0(), beamSpotreco.y0(), beamSpotreco.z0());

  
  ////////////////////////// 
  // Cosmic Compatibility //
  //////////////////////////

  Handle<MuonCollection> muons;
  iEvent.getByLabel( "muons", muons );
  Handle<ValueMap<MuonCosmicCompatibility> > CosmicMap;
  iEvent.getByLabel( src_, CosmicMap );


  ///////////
  // Muons // 
  ///////////
  
  unsigned int muonIndex = 0;
  View<Muon>::const_iterator muons_end = muon_h->end();  // Iterator
  for ( View<Muon>::const_iterator muon = muon_h->begin(); muon != muons_end; ++muon ) {

    // References
    const RefToBase<Muon>         muonRef                 = muon_h->refAt(muonIndex); 
    const TrackRef                globalTrack             = muon->globalTrack();
    const TrackRef                siTrack                 = muon->innerTrack();
    const TrackRef                staTrack                = muon->outerTrack();
    const MuonQuality             quality                 = muon->combinedQuality();
    const VertexCollection*       vertexCollection        = vertexHandle.product();

    // Iterators
    VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();


    /////////
    // STA //
    /////////

    vector_mus_sta_d0            -> push_back( staTrack.isNonnull()  ? staTrack->d0()                   :  -9999.        );
    vector_mus_sta_z0            -> push_back( staTrack.isNonnull()  ? staTrack->dz()                   :  -9999.        );
    vector_mus_sta_d0Err         -> push_back( staTrack.isNonnull()  ? staTrack->d0Error()              :  -9999.        );
    vector_mus_sta_z0Err         -> push_back( staTrack.isNonnull()  ? staTrack->dzError()              :  -9999.        );
    vector_mus_sta_d0corr        -> push_back( staTrack.isNonnull()  ? -1*(staTrack->dxy(beamSpot))     :  -9999.        );
    vector_mus_sta_z0corr        -> push_back( staTrack.isNonnull()  ? staTrack->dz(beamSpot)           :  -9999.        );
    vector_mus_sta_p4            -> push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->px() , staTrack->py() , staTrack->pz() , staTrack->p() ) : LorentzVector(0, 0, 0, 0) );


    ////////////
    // Global //
    ////////////

    vector_mus_gfit_chi2         -> push_back( globalTrack.isNonnull() ? globalTrack->chi2()               :  -9999.        );
    vector_mus_gfit_ndof         -> push_back( globalTrack.isNonnull() ? globalTrack->ndof()               :  -9999         );
    vector_mus_gfit_validSTAHits -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidMuonHits()    : -9999         );

    //////////////////
    // Muon Quality //
    //////////////////

    vector_mus_trkKink             -> push_back( quality.trkKink             );


    //////////////////////////
    // Cosmic Compatibility //
    //////////////////////////


    //////////
    // Muon //
    //////////

    // Calculate Overlaps
    int mus_overlap0 = -1, muInd = -1, mus_nOverlaps = 0;
    for ( View<Muon>::const_iterator muonJ = muon_h->begin(); muonJ != muons_end; ++muonJ ) {
      muInd++;
      if ( muonJ != muon ){
        if ( muon::overlap( *muon, *muonJ ) ) {
          if ( mus_overlap0 == -1) mus_overlap0 = muInd;
          mus_nOverlaps++;
        }
      }
    }

    // Calculate Muon position at ECAL
    math::XYZPoint ecal_p( -9999.0, -9999.0, -9999.0 );
    if( muon->isEnergyValid() ) ecal_p = muon->calEnergy().ecal_position;

    // Calculate Mask
    int goodMask = 0;
    for ( int iG = 0; iG < 24; ++iG ) { //overkill here
      if( isGoodMuon( *muon, (muon::SelectionType)iG ) ) goodMask |= (1<<iG);
    }

    
    /////////////////////
    // Muon Quantities //
    /////////////////////

    vector_mus_type                    -> push_back( muon->type()                                              );
    vector_mus_goodmask                -> push_back( goodMask                                                  );
    vector_mus_charge                  -> push_back( muon->charge()                                            );
    vector_mus_nmatches                -> push_back( muon->isMatchesValid() ? muon->numberOfMatches() :  -9999 );
    vector_mus_caloCompatibility       -> push_back( muon->caloCompatibility()                                 );
    vector_mus_segmCompatibility       -> push_back( muon::segmentCompatibility(*muon)                         );
    vector_mus_p4                      -> push_back( LorentzVector( muon->p4()                              )  );
    vector_mus_numberOfMatchedStations ->push_back( muon->numberOfMatchedStations()                            );


    /////////////////////////////
    // Muon Shower Information //
    /////////////////////////////

    const MuonShower muShower = showerMap.isValid() ? (*showerMap)[muonRef] : MuonShower();


    ////////
    // ID //
    ////////

    bool matchIsValid = muon->isMatchesValid();

    vector_mus_pid_TMLastStationTight     -> push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationTight     ) : -9999  );
    vector_mus_pid_PFMuon                 -> push_back( muon->isPFMuon() );

    ////////////
    // Energy //
    ////////////

    bool energyIsValid = muon->isEnergyValid();

    vector_mus_e_em               -> push_back( energyIsValid ? muon->calEnergy().em                                 :  -9999.       );
    vector_mus_e_had              -> push_back( energyIsValid ? muon->calEnergy().had                                :  -9999.       );
    vector_mus_e_ho               -> push_back( energyIsValid ? muon->calEnergy().ho                                 :  -9999.       );
    vector_mus_e_emS9             -> push_back( energyIsValid ? muon->calEnergy().emS9                               :  -9999.       );
    vector_mus_e_hadS9            -> push_back( energyIsValid ? muon->calEnergy().hadS9                              :  -9999.       );
    vector_mus_e_hoS9             -> push_back( energyIsValid ? muon->calEnergy().hoS9                               :  -9999.       );


    ///////////////
    // Isolation //
    ///////////////

    vector_mus_iso_ecalvetoDep    -> push_back( muon->isEnergyValid()    ? muon->isolationR03().emVetoEt       : -9999.        );
    vector_mus_iso_hcalvetoDep    -> push_back( muon->isEnergyValid()    ? muon->isolationR03().hadVetoEt      : -9999.        );
    vector_mus_iso03_sumPt        -> push_back( muon->isIsolationValid() ? muon->isolationR03().sumPt          : -9999.        );
    vector_mus_iso03_emEt         -> push_back( muon->isIsolationValid() ? muon->isolationR03().emEt           : -9999.        );
    vector_mus_iso03_hadEt        -> push_back( muon->isIsolationValid() ? muon->isolationR03().hadEt          : -9999.        );
    vector_mus_iso05_sumPt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().sumPt          : -9999.        );
    vector_mus_iso05_emEt         -> push_back( muon->isIsolationValid() ? muon->isolationR05().emEt           : -9999.        );
    vector_mus_iso05_hadEt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().hadEt          : -9999.        );

    ////////////
    // Tracks //
    ////////////

    vector_mus_trk_p4             -> push_back( siTrack.isNonnull()     ? LorentzVector( siTrack.get()->px() , siTrack.get()->py() , siTrack.get()->pz() , siTrack.get()->p() ) : LorentzVector(     0.0,     0.0,     0.0,     0.0) );
    vector_mus_vertex_p4          -> push_back( siTrack.isNonnull()     ? LorentzVector( siTrack->vx()       , siTrack->vy()       , siTrack->vz()       , 0.0                ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );
    vector_mus_trkidx             -> push_back( siTrack.isNonnull()     ? static_cast<int>(siTrack.key())                      : -9999         );
    vector_mus_d0                 -> push_back( siTrack.isNonnull()     ? siTrack->d0()                                        : -9999.        );
    vector_mus_z0                 -> push_back( siTrack.isNonnull()     ? siTrack->dz()                                        : -9999.        );
    //vector_mus_d0corr             -> push_back( siTrack.isNonnull()     ? -1*(siTrack->dxy(beamSpotreco))                          : -9999.        );
    //vector_mus_z0corr             -> push_back( siTrack.isNonnull()     ? siTrack->dz(beamSpotreco)                                : -9999.        );
    vector_mus_chi2               -> push_back( siTrack.isNonnull()     ? siTrack->chi2()                                      : -9999.        );
    vector_mus_ndof               -> push_back( siTrack.isNonnull()     ? siTrack->ndof()                                      : -9999.        );
    vector_mus_validHits          -> push_back( siTrack.isNonnull()     ? siTrack->numberOfValidHits()                         : -9999         );
    vector_mus_d0Err              -> push_back( siTrack.isNonnull()     ? siTrack->d0Error()                                   :  -9999.       );
    vector_mus_z0Err              -> push_back( siTrack.isNonnull()     ? siTrack->dzError()                                   :  -9999.       );
    vector_mus_ptErr              -> push_back( siTrack.isNonnull()     ? siTrack->ptError()                                   :  -9999.       );
    vector_mus_etaErr             -> push_back( siTrack.isNonnull()     ? siTrack->etaError()                                  :  -9999.       );
    vector_mus_phiErr             -> push_back( siTrack.isNonnull()     ? siTrack->phiError()                                  :  -9999.       );


    ////////
    // PF //
    ////////

    MuonPFIsolation pfStructR03 = muon->pfIsolationR03();

    vector_mus_isoR03_pf_ChargedHadronPt                 -> push_back( pfStructR03.sumChargedHadronPt              );
    vector_mus_isoR03_pf_ChargedParticlePt               -> push_back( pfStructR03.sumChargedParticlePt            );
    vector_mus_isoR03_pf_NeutralHadronEt                 -> push_back( pfStructR03.sumNeutralHadronEt              );
    vector_mus_isoR03_pf_PhotonEt                        -> push_back( pfStructR03.sumPhotonEt                     );
    vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructR03.sumNeutralHadronEtHighThreshold );
    vector_mus_isoR03_pf_sumPhotonEtHighThreshold        -> push_back( pfStructR03.sumPhotonEtHighThreshold        );
    vector_mus_isoR03_pf_PUPt                            -> push_back( pfStructR03.sumPUPt                         );


    ///////////
    // IP 3D //
    ///////////

    if ( siTrack.isNonnull() && firstGoodVertex != vertexCollection->end() ) {
      
      //TransientTrack tt       = theTTBuilder->build( siTrack );
      //Measurement1D ip3D      = IPTools::absoluteImpactParameter3D( tt, *firstGoodVertex ).second;
      //vector_mus_ip3d         -> push_back( ip3D.value()                    );
      //vector_mus_ip3derr      -> push_back( ip3D.error()                    );

    } else {

      vector_mus_ip3d         -> push_back( -9999. );
      vector_mus_ip3derr      -> push_back( -9999. );

    } 

    muonIndex++;

  } // end loop on muons


  /////////
  // STA //
  /////////

  iEvent.put( vector_mus_sta_d0                       , branchprefix_ + "stad0"              );
  iEvent.put( vector_mus_sta_z0                       , branchprefix_ + "staz0"              );
  iEvent.put( vector_mus_sta_d0Err                    , branchprefix_ + "stad0Err"           );
  iEvent.put( vector_mus_sta_z0Err                    , branchprefix_ + "staz0Err"           );
  iEvent.put( vector_mus_sta_d0corr                   , branchprefix_ + "stad0corr"          );
  iEvent.put( vector_mus_sta_z0corr                   , branchprefix_ + "staz0corr"          );
  iEvent.put( vector_mus_sta_p4                       , branchprefix_ + "stap4"              );
  

  ////////////                                                                       
  // Global //
  /////////////

  iEvent.put( vector_mus_gfit_chi2                    , branchprefix_ + "gfitchi2"           );
  iEvent.put( vector_mus_gfit_ndof                    , branchprefix_ + "gfitndof"           );
  iEvent.put( vector_mus_gfit_validSTAHits            , branchprefix_ + "gfitvalidSTAHits"   );


  //////////////////
  // Muon Quality //
  //////////////////

  iEvent.put( vector_mus_trkKink            , branchprefix_ + "trkKink"            );


  ///////////
  // Muons //
  ///////////

  iEvent.put( vector_mus_type                    , branchprefix_ + "type"                    );
  iEvent.put( vector_mus_goodmask                , branchprefix_ + "goodmask"                );
  iEvent.put( vector_mus_charge                  , branchprefix_ + "charge"                  );
  iEvent.put( vector_mus_nmatches                , branchprefix_ + "nmatches"                );
  iEvent.put( vector_mus_caloCompatibility       , branchprefix_ + "caloCompatibility"       );
  iEvent.put( vector_mus_segmCompatibility       , branchprefix_ + "segmCompatibility"       );
  iEvent.put( vector_mus_p4                      , branchprefix_ + "p4"                      );
  iEvent.put( vector_mus_numberOfMatchedStations , branchprefix_ + "numberOfMatchedStations" );


  ////////
  // ID //
  ////////

  iEvent.put( vector_mus_pid_TMLastStationTight       , branchprefix_ + "pidTMLastStationTight"    );
  iEvent.put( vector_mus_pid_PFMuon                   , branchprefix_ + "pidPFMuon");


  ////////////
  // Energy //
  ////////////

  iEvent.put( vector_mus_e_em               , branchprefix_ + "eem"                );
  iEvent.put( vector_mus_e_had              , branchprefix_ + "ehad"               );
  iEvent.put( vector_mus_e_ho               , branchprefix_ + "eho"                );
  iEvent.put( vector_mus_e_emS9             , branchprefix_ + "eemS9"              );
  iEvent.put( vector_mus_e_hadS9            , branchprefix_ + "ehadS9"             );
  iEvent.put( vector_mus_e_hoS9             , branchprefix_ + "ehoS9"              );

  
  ///////////////
  // Isolation //
  ///////////////

  iEvent.put( vector_mus_iso_ecalvetoDep    , branchprefix_ + "isoecalvetoDep"     );
  iEvent.put( vector_mus_iso_hcalvetoDep    , branchprefix_ + "isohcalvetoDep"     );
  iEvent.put( vector_mus_iso03_sumPt        , branchprefix_ + "iso03sumPt"         );
  iEvent.put( vector_mus_iso03_emEt         , branchprefix_ + "iso03emEt"          );
  iEvent.put( vector_mus_iso03_hadEt        , branchprefix_ + "iso03hadEt"         );
  iEvent.put( vector_mus_iso05_sumPt        , branchprefix_ + "iso05sumPt"         );
  iEvent.put( vector_mus_iso05_emEt         , branchprefix_ + "iso05emEt"          );
  iEvent.put( vector_mus_iso05_hadEt        , branchprefix_ + "iso05hadEt"         );


  ////////////
  // Tracks //
  ////////////

  iEvent.put( vector_mus_trk_p4             , branchprefix_ + "trkp4"              );
  iEvent.put( vector_mus_vertex_p4          , branchprefix_ + "vertexp4"           );
  iEvent.put( vector_mus_trkidx             , branchprefix_ + "trkidx"             );
  iEvent.put( vector_mus_d0                 , branchprefix_ + "d0"                 );
  iEvent.put( vector_mus_z0                 , branchprefix_ + "z0"                 );
  iEvent.put( vector_mus_d0corr             , branchprefix_ + "d0corr"             );
  iEvent.put( vector_mus_z0corr             , branchprefix_ + "z0corr"             );
  iEvent.put( vector_mus_chi2               , branchprefix_ + "chi2"               );
  iEvent.put( vector_mus_ndof               , branchprefix_ + "ndof"               );
  iEvent.put( vector_mus_validHits          , branchprefix_ + "validHits"          );
  iEvent.put( vector_mus_d0Err              , branchprefix_ + "d0Err"              );
  iEvent.put( vector_mus_z0Err              , branchprefix_ + "z0Err"              );
  iEvent.put( vector_mus_ptErr              , branchprefix_ + "ptErr"              );
  iEvent.put( vector_mus_etaErr             , branchprefix_ + "etaErr"             );
  iEvent.put( vector_mus_phiErr             , branchprefix_ + "phiErr"             );

  ////////                                  
  // PF //
  ////////

  iEvent.put( vector_mus_isoR03_pf_ChargedHadronPt                , branchprefix_ + "isoR03pfChargedHadronPt"             );
  iEvent.put( vector_mus_isoR03_pf_ChargedParticlePt              , branchprefix_ + "isoR03pfChargedParticlePt"           );
  iEvent.put( vector_mus_isoR03_pf_NeutralHadronEt                , branchprefix_ + "isoR03pfNeutralHadronEt"             );
  iEvent.put( vector_mus_isoR03_pf_PhotonEt                       , branchprefix_ + "isoR03pfPhotonEt"                    );
  iEvent.put( vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold, branchprefix_ + "isoR03pfNeutralHadronEtHighThreshold");
  iEvent.put( vector_mus_isoR03_pf_sumPhotonEtHighThreshold       , branchprefix_ + "isoR03pfPhotonEtHighThreshold"       );
  iEvent.put( vector_mus_isoR03_pf_PUPt                           , branchprefix_ + "isoR03pfPUPt"                        );


  ///////////
  // IP 3D //
  ///////////

  iEvent.put( vector_mus_ip3d                         , branchprefix_ + "ip3d"               );
  iEvent.put( vector_mus_ip3derr                      , branchprefix_ + "ip3derr"            );



} 

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
