#include "HWWValidation/HWWBase/interface/MuonMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

using namespace std;
using namespace reco;
using namespace edm;


MuonMaker::MuonMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector) {

  Muon_                     = iCollector.consumes<edm::View<reco::Muon> > (iConfig.getParameter<edm::InputTag>("muonsInputTag"));
  MuonShower_               = iCollector.consumes<edm::ValueMap<reco::MuonShower> > (iConfig.getParameter<edm::InputTag>("muonShower"));
  thePVCollection_          = iCollector.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));
  PFCandidateCollection_    = iCollector.consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandsInputTag"));
  BeamSpot_                 = iCollector.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  MuonCollection_           = iCollector.consumes<reco::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonsInputTag"));

}


void MuonMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_mus_sta_d0();
  HWWVal::Load_mus_sta_z0corr();
  HWWVal::Load_mus_sta_p4();
  HWWVal::Load_mus_gfit_chi2();
  HWWVal::Load_mus_gfit_ndof();
  HWWVal::Load_mus_gfit_validSTAHits();
  HWWVal::Load_mus_trkKink();
  HWWVal::Load_mus_type();
  HWWVal::Load_mus_goodmask();
  HWWVal::Load_mus_charge();
  HWWVal::Load_mus_nmatches();
  HWWVal::Load_mus_caloCompatibility();
  HWWVal::Load_mus_segmCompatibility();
  HWWVal::Load_mus_p4();
  HWWVal::Load_mus_numberOfMatchedStations();
  HWWVal::Load_mus_pid_TMLastStationTight();
  HWWVal::Load_mus_pid_PFMuon();
  HWWVal::Load_mus_e_em();
  HWWVal::Load_mus_e_had();
  HWWVal::Load_mus_e_ho();
  HWWVal::Load_mus_e_emS9();
  HWWVal::Load_mus_e_hadS9();
  HWWVal::Load_mus_e_hoS9();
  HWWVal::Load_mus_iso_ecalvetoDep();
  HWWVal::Load_mus_iso_hcalvetoDep();
  HWWVal::Load_mus_iso03_sumPt();
  HWWVal::Load_mus_iso03_emEt();
  HWWVal::Load_mus_iso03_hadEt();
  HWWVal::Load_mus_iso05_sumPt();
  HWWVal::Load_mus_iso05_emEt();
  HWWVal::Load_mus_iso05_hadEt();
  HWWVal::Load_mus_trk_p4();
  HWWVal::Load_mus_vertex_p4();
  HWWVal::Load_mus_trkidx();
  HWWVal::Load_mus_d0();
  HWWVal::Load_mus_chi2();
  HWWVal::Load_mus_ndof();
  HWWVal::Load_mus_validHits();
  HWWVal::Load_mus_ptErr();
  HWWVal::Load_mus_isoR03_pf_ChargedHadronPt();
  HWWVal::Load_mus_isoR03_pf_NeutralHadronEt();
  HWWVal::Load_mus_isoR03_pf_PhotonEt();
  HWWVal::Load_mus_isoR03_pf_PUPt();
  HWWVal::Load_mus_ip3d();
  HWWVal::Load_mus_ip3derr();



  ///////////////
  // Get Muons //
  ///////////////

  Handle<View<Muon> > muon_h;
  iEvent.getByToken( Muon_ , muon_h );

  
  /////////////////////////////////
  // Get Muon Shower Information //
  /////////////////////////////////

  Handle<ValueMap<MuonShower> > showerMap;
  iEvent.getByToken( MuonShower_ , showerMap );


  //////////////////
  // Get Vertices //
  //////////////////

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken( thePVCollection_ , vertexHandle );  


  ///////////////////////
  // Get PF Candidates //
  ///////////////////////

  edm::Handle<reco::PFCandidateCollection> pfCand_h;
  iEvent.getByToken( PFCandidateCollection_ , pfCand_h );


  //////////////
  // Beamspot //
  //////////////

  Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByToken(BeamSpot_, beamspot_h);
  const reco::BeamSpot &beamSpotreco = *(beamspot_h.product());
  const Point beamSpot = Point(beamSpotreco.x0(), beamSpotreco.y0(), beamSpotreco.z0());

  
  ////////////////////////// 
  // Cosmic Compatibility //
  //////////////////////////

  Handle<MuonCollection> muons;
  iEvent.getByToken( MuonCollection_, muons );

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

    HWWVal::mus_sta_d0()            .push_back( staTrack.isNonnull()  ? staTrack->d0()                   :  -9999.        );
    HWWVal::mus_sta_z0corr()        .push_back( staTrack.isNonnull()  ? staTrack->dz(beamSpot)           :  -9999.        );
    HWWVal::mus_sta_p4()            .push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->px() , staTrack->py() , staTrack->pz() , staTrack->p() ) : LorentzVector(0, 0, 0, 0) );


    ////////////
    // Global //
    ////////////

    HWWVal::mus_gfit_chi2()         .push_back( globalTrack.isNonnull() ? globalTrack->chi2()               :  -9999.        );
    HWWVal::mus_gfit_ndof()         .push_back( globalTrack.isNonnull() ? globalTrack->ndof()               :  -9999         );
    HWWVal::mus_gfit_validSTAHits() .push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidMuonHits()    : -9999         );

    //////////////////
    // Muon Quality //
    //////////////////

    HWWVal::mus_trkKink()             .push_back( quality.trkKink             );



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

    HWWVal::mus_type()                    .push_back( muon->type()                                              );
    HWWVal::mus_goodmask()                .push_back( goodMask                                                  );
    HWWVal::mus_charge()                  .push_back( muon->charge()                                            );
    HWWVal::mus_nmatches()                .push_back( muon->isMatchesValid() ? muon->numberOfMatches() :  -9999 );
    HWWVal::mus_caloCompatibility()       .push_back( muon->caloCompatibility()                                 );
    HWWVal::mus_segmCompatibility()       .push_back( muon::segmentCompatibility(*muon)                         );
    HWWVal::mus_p4()                      .push_back( LorentzVector( muon->p4()                              )  );
    HWWVal::mus_numberOfMatchedStations() .push_back( muon->numberOfMatchedStations()                            );


    /////////////////////////////
    // Muon Shower Information //
    /////////////////////////////

    const MuonShower muShower = showerMap.isValid() ? (*showerMap)[muonRef] : MuonShower();


    ////////
    // ID //
    ////////

    bool matchIsValid = muon->isMatchesValid();

    HWWVal::mus_pid_TMLastStationTight()     .push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationTight     ) : -9999  );
    HWWVal::mus_pid_PFMuon()                 .push_back( muon->isPFMuon() );

    ////////////
    // Energy //
    ////////////

    bool energyIsValid = muon->isEnergyValid();

    HWWVal::mus_e_em()               .push_back( energyIsValid ? muon->calEnergy().em                                 :  -9999.       );
    HWWVal::mus_e_had()              .push_back( energyIsValid ? muon->calEnergy().had                                :  -9999.       );
    HWWVal::mus_e_ho()               .push_back( energyIsValid ? muon->calEnergy().ho                                 :  -9999.       );
    HWWVal::mus_e_emS9()             .push_back( energyIsValid ? muon->calEnergy().emS9                               :  -9999.       );
    HWWVal::mus_e_hadS9()            .push_back( energyIsValid ? muon->calEnergy().hadS9                              :  -9999.       );
    HWWVal::mus_e_hoS9()             .push_back( energyIsValid ? muon->calEnergy().hoS9                               :  -9999.       );


    ///////////////
    // Isolation //
    ///////////////

    HWWVal::mus_iso_ecalvetoDep()    .push_back( muon->isEnergyValid()    ? muon->isolationR03().emVetoEt       : -9999.        );
    HWWVal::mus_iso_hcalvetoDep()    .push_back( muon->isEnergyValid()    ? muon->isolationR03().hadVetoEt      : -9999.        );
    HWWVal::mus_iso03_sumPt()        .push_back( muon->isIsolationValid() ? muon->isolationR03().sumPt          : -9999.        );
    HWWVal::mus_iso03_emEt()         .push_back( muon->isIsolationValid() ? muon->isolationR03().emEt           : -9999.        );
    HWWVal::mus_iso03_hadEt()        .push_back( muon->isIsolationValid() ? muon->isolationR03().hadEt          : -9999.        );
    HWWVal::mus_iso05_sumPt()        .push_back( muon->isIsolationValid() ? muon->isolationR05().sumPt          : -9999.        );
    HWWVal::mus_iso05_emEt()         .push_back( muon->isIsolationValid() ? muon->isolationR05().emEt           : -9999.        );
    HWWVal::mus_iso05_hadEt()        .push_back( muon->isIsolationValid() ? muon->isolationR05().hadEt          : -9999.        );

    ////////////
    // Tracks //
    ////////////

    HWWVal::mus_trk_p4()             .push_back( siTrack.isNonnull()     ? LorentzVector( siTrack.get()->px() , siTrack.get()->py() , siTrack.get()->pz() , siTrack.get()->p() ) : LorentzVector(     0.0,     0.0,     0.0,     0.0) );
    HWWVal::mus_vertex_p4()          .push_back( siTrack.isNonnull()     ? LorentzVector( siTrack->vx()       , siTrack->vy()       , siTrack->vz()       , 0.0                ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );
    HWWVal::mus_trkidx()             .push_back( siTrack.isNonnull()     ? static_cast<int>(siTrack.key())                      : -9999         );
    HWWVal::mus_d0()                 .push_back( siTrack.isNonnull()     ? siTrack->d0()                                        : -9999.        );
    HWWVal::mus_chi2()               .push_back( siTrack.isNonnull()     ? siTrack->chi2()                                      : -9999.        );
    HWWVal::mus_ndof()               .push_back( siTrack.isNonnull()     ? siTrack->ndof()                                      : -9999.        );
    HWWVal::mus_validHits()          .push_back( siTrack.isNonnull()     ? siTrack->numberOfValidHits()                         : -9999         );
    HWWVal::mus_ptErr()              .push_back( siTrack.isNonnull()     ? siTrack->ptError()                                   :  -9999.       );


    ////////
    // PF //
    ////////

    MuonPFIsolation pfStructR03 = muon->pfIsolationR03();

    HWWVal::mus_isoR03_pf_ChargedHadronPt()                 .push_back( pfStructR03.sumChargedHadronPt              );
    HWWVal::mus_isoR03_pf_NeutralHadronEt()                 .push_back( pfStructR03.sumNeutralHadronEt              );
    HWWVal::mus_isoR03_pf_PhotonEt()                        .push_back( pfStructR03.sumPhotonEt                     );
    HWWVal::mus_isoR03_pf_PUPt()                            .push_back( pfStructR03.sumPUPt                         );


    ///////////
    // IP 3D //
    ///////////

    if ( siTrack.isNonnull() && firstGoodVertex != vertexCollection->end() ) {
      
      //TransientTrack tt       = theTTBuilder->build( siTrack );
      //Measurement1D ip3D      = IPTools::absoluteImpactParameter3D( tt, *firstGoodVertex ).second;
      //HWWVal::mus_ip3d         -> push_back( ip3D.value()                    );
      //HWWVal::mus_ip3derr      -> push_back( ip3D.error()                    );

    } else {

      HWWVal::mus_ip3d()         .push_back( -9999. );
      HWWVal::mus_ip3derr()      .push_back( -9999. );

    } 

    muonIndex++;

  } // end loop on muons

} 
