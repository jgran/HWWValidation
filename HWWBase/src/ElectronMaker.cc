
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"


#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "HWWValidation/HWWBase/interface/EgammaFiduciality.h"
#include "HWWValidation/HWWBase/interface/ElectronMaker.h"


using namespace reco;
using namespace edm;
using namespace std;

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

// constructors and destructor
ElectronMaker::ElectronMaker(const ParameterSet& iConfig) {

    //get setup parameters
    electronsInputTag_        = iConfig.getParameter<edm::InputTag>   ("electronsInputTag"        );
    trksInputTag_             = iConfig.getParameter<edm::InputTag>   ("trksInputTag"             );
    gsftracksInputTag_        = iConfig.getParameter<edm::InputTag>   ("gsftracksInputTag"        );
    pfCandsInputTag           = iConfig.getParameter<edm::InputTag>   ("pfCandsInputTag"          );
    vtxInputTag               = iConfig.getParameter<edm::InputTag>   ("vtxInputTag"              );
    recoConversionInputTag_   = iConfig.getParameter<edm::InputTag>   ("recoConversionInputTag"   );
    beamSpot_tag_             = iConfig.getParameter<edm::InputTag>   ("beamSpotTag"              );
    minAbsDist_               = iConfig.getParameter<double>          ("minAbsDist"              );
    minAbsDcot_               = iConfig.getParameter<double>          ("minAbsDcot"              );
    minSharedFractionOfHits_  = iConfig.getParameter<double>          ("minSharedFractionOfHits" );
    aliasprefix_              = iConfig.getUntrackedParameter<string> ("aliasPrefix"             );

    mtsTransform_ = 0;
    clusterTools_ = 0;

    // ECAL related (superCluster) variables
    produces<vector<int> >       ("elsnSeed"                   ).setBranchAlias("els_nSeed"                  );
    produces<vector<float> >     ("elseSC"                     ).setBranchAlias("els_eSC"                    );
    produces<vector<float> >     ("elsetaSC"                   ).setBranchAlias("els_etaSC"                  );
    produces<vector<float> >     ("elsphiSC"                   ).setBranchAlias("els_phiSC"                  );
    produces<vector<float> >     ("elseSCRaw"                  ).setBranchAlias("els_eSCRaw"                 );
    produces<vector<float> >     ("elseSCPresh"                ).setBranchAlias("els_eSCPresh"               );
    produces<vector<int> >       ("elsscindex"                 ).setBranchAlias("els_scindex"                );
    produces<vector<float> >     ("elsetaSCwidth"              ).setBranchAlias("els_etaSCwidth"             );
    produces<vector<float> >     ("elsphiSCwidth"              ).setBranchAlias("els_phiSCwidth"             );
    produces<vector<int> >       ("elsfiduciality"             ).setBranchAlias("els_fiduciality"            );
    produces<vector<int> >       ("elstype"                    ).setBranchAlias("els_type"                   );

    // Corrections and uncertainties
    produces<vector<float> >     ("elsecalEnergy"              ).setBranchAlias("els_ecalEnergy"             );
    produces<vector<float> >     ("elsecalEnergyError"         ).setBranchAlias("els_ecalEnergyError"        );

    // ID variables
    produces<vector<float> >     ("elslh"                      ).setBranchAlias("els_lh"                     );
    produces<vector<float> >     ("elsdEtaIn"                  ).setBranchAlias("els_dEtaIn"                 );
    produces<vector<float> >     ("elsdEtaOut"                 ).setBranchAlias("els_dEtaOut"                );
    produces<vector<float> >     ("elsdPhiIn"                  ).setBranchAlias("els_dPhiIn"                 );
    produces<vector<float> >     ("elsdPhiOut"                 ).setBranchAlias("els_dPhiOut"                );
    produces<vector<float> >     ("elsdPhiInPhiOut"            ).setBranchAlias("els_dPhiInPhiOut"           );
    produces<vector<float> >     ("elsfbrem"                   ).setBranchAlias("els_fbrem"                  );
    produces<vector<float> >     ("elseOverPIn"                ).setBranchAlias("els_eOverPIn"               );
    produces<vector<float> >     ("elseSeedOverPOut"           ).setBranchAlias("els_eSeedOverPOut"          );
    produces<vector<float> >     ("elseSeedOverPIn"            ).setBranchAlias("els_eSeedOverPIn"           );
    produces<vector<float> >     ("elseOverPOut"               ).setBranchAlias("els_eOverPOut"              );
    produces<vector<float> >     ("elshOverE"                  ).setBranchAlias("els_hOverE"                 );
    produces<vector<float> >     ("elssigmaIPhiIPhi"           ).setBranchAlias("els_sigmaIPhiIPhi"          );
    produces<vector<float> >     ("elssigmaIEtaIEta"           ).setBranchAlias("els_sigmaIEtaIEta"          );
    produces<vector<float> >     ("else1x5"                    ).setBranchAlias("els_e1x5"                   );
    produces<vector<float> >     ("else5x5"                    ).setBranchAlias("els_e5x5"                   );
    produces<vector<float> >     ("else3x3"                    ).setBranchAlias("els_e3x3"                   );

    // isolation variables
    produces<vector<float> >     ("elstkIso"                  ).setBranchAlias("els_tkIso"                  );
    produces<vector<float> >     ("elsecalIso"                ).setBranchAlias("els_ecalIso"                );
    produces<vector<float> >     ("elshcalIso"                ).setBranchAlias("els_hcalIso"                );
    produces<vector<float> >     ("elsecalIso04"              ).setBranchAlias("els_ecalIso04"              );
    produces<vector<float> >     ("elshcalIso04"              ).setBranchAlias("els_hcalIso04"              );
    produces<vector<float> >     ("elsiso03pfch"              ).setBranchAlias("els_iso03_pf_ch"            ); // pf isolation in cone of 0.3, charged only
    produces<vector<float> >     ("elsiso03pfgamma05"         ).setBranchAlias("els_iso03_pf_gamma05"       ); // pf isolation in cone of 0.3, photons only with threshold 0.5 GeV
    produces<vector<float> >     ("elsiso03pfnhad05"          ).setBranchAlias("els_iso03_pf_nhad05"        ); // pf isolation in cone of 0.3, neutral hadrons only with threshold 0.5 GeV
    produces<vector<float> >     ("elsiso04pfch"              ).setBranchAlias("els_iso04_pf_ch"            ); // pf isolation in cone of 0.3, charged only
    produces<vector<float> >     ("elsiso04pfgamma05"         ).setBranchAlias("els_iso04_pf_gamma05"       ); // pf isolation in cone of 0.3, photons only with threshold 0.5 GeV
    produces<vector<float> >     ("elsiso04pfnhad05"          ).setBranchAlias("els_iso04_pf_nhad05"        ); // pf isolation in cone of 0.3, neutral hadrons only with threshold 0.5 GeV

    // 2012 Electron Particle Flow Isolation
    produces<vector<float> >     ("elsiso03pf2012ch"             ).setBranchAlias("els_iso03_pf2012_ch"    );
    produces<vector<float> >     ("elsiso03pf2012em"             ).setBranchAlias("els_iso03_pf2012_em"    );
    produces<vector<float> >     ("elsiso03pf2012nh"             ).setBranchAlias("els_iso03_pf2012_nh"    );
    produces<vector<float> >     ("elsiso04pf2012ch"             ).setBranchAlias("els_iso04_pf2012_ch"    );
    produces<vector<float> >     ("elsiso04pf2012em"             ).setBranchAlias("els_iso04_pf2012_em"    );
    produces<vector<float> >     ("elsiso04pf2012nh"             ).setBranchAlias("els_iso04_pf2012_nh"    );

    // track variables
    produces<vector<int> >       ("elscharge"        ).setBranchAlias("els_charge"         ); //candidate charge
    produces<vector<int> >       ("elssccharge"      ).setBranchAlias("els_sccharge"       );
    produces<vector<int> >       ("elstrkcharge"     ).setBranchAlias("els_trk_charge"     );
    produces<vector<float> >     ("elschi2"          ).setBranchAlias("els_chi2"           );
    produces<vector<float> >     ("elsndof"          ).setBranchAlias("els_ndof"           );
    produces<vector<float> >     ("elsd0"            ).setBranchAlias("els_d0"             );
    produces<vector<float> >     ("elsz0"            ).setBranchAlias("els_z0"             );
    produces<vector<float> >     ("elsd0Err"         ).setBranchAlias("els_d0Err"          );
    produces<vector<float> >     ("elsz0Err"         ).setBranchAlias("els_z0Err"          );
    produces<vector<float> >     ("elsd0corr"        ).setBranchAlias("els_d0corr"         );
    produces<vector<float> >     ("elsz0corr"        ).setBranchAlias("els_z0corr"         );
    produces<vector<float> >     ("elsptErr"         ).setBranchAlias("els_ptErr"          );
    produces<vector<float> >     ("elsetaErr"        ).setBranchAlias("els_etaErr"         );
    produces<vector<float> >     ("elsphiErr"        ).setBranchAlias("els_phiErr"         );
    produces<vector<int> >       ("elsgsftrkidx"     ).setBranchAlias("els_gsftrkidx"      );
    produces<vector<float> >     ("elsip3d"          ).setBranchAlias("els_ip3d"           ); // Ip3d from normal vertex
    produces<vector<float> >     ("elsip3derr"       ).setBranchAlias("els_ip3derr"        ); // Ip3d error from normal vertex

    // LorentzVectors
    produces<vector<LorentzVector> >  ("elsp4"    ).setBranchAlias("els_p4"     );
    produces<vector<LorentzVector> >  ("elstrkp4" ).setBranchAlias("els_trk_p4" );

    // Vertex
    produces<vector<LorentzVector> >  ("elsvertexp4").setBranchAlias("els_vertex_p4");

    //Hit Pattern information
    produces<vector<int> >            ("elsexpinnerlayers" ).setBranchAlias("els_exp_innerlayers" );
    produces<vector<int> >            ("elsexpouterlayers" ).setBranchAlias("els_exp_outerlayers" );   

    //CTF track matching stuff
    produces<vector<int>    >    ("elstrkidx"    ).setBranchAlias("els_trkidx"    );// track index matched to electron
    produces<vector<float>  >    ("elstrkshFrac" ).setBranchAlias("els_trkshFrac" );

    //conversion stuff
    produces<vector<float>    >       ("elsconvdist"        ).setBranchAlias("els_conv_dist"        );
    produces<vector<float>    >       ("elsconvdcot"        ).setBranchAlias("els_conv_dcot"        );
    produces<vector<float>    >       ("elsconvolddist"        ).setBranchAlias("els_conv_old_dist"        );
    produces<vector<float>    >       ("elsconvolddcot"        ).setBranchAlias("els_conv_old_dcot"        );

    produces<vector<float> >     ("elsiso03pf2012extch"             ).setBranchAlias("els_iso03_pf2012ext_ch"    );
    produces<vector<float> >     ("elsiso03pf2012extem"             ).setBranchAlias("els_iso03_pf2012ext_em"    );
    produces<vector<float> >     ("elsiso03pf2012extnh"             ).setBranchAlias("els_iso03_pf2012ext_nh"    );
    produces<vector<float> >     ("elsiso04pf2012extch"             ).setBranchAlias("els_iso04_pf2012ext_ch"    );
    produces<vector<float> >     ("elsiso04pf2012extem"             ).setBranchAlias("els_iso04_pf2012ext_em"    );
    produces<vector<float> >     ("elsiso04pf2012extnh"             ).setBranchAlias("els_iso04_pf2012ext_nh"    );

    pfPileUpAlgo_ = new PFPileUpAlgo();
}

ElectronMaker::~ElectronMaker()
{
  if (pfPileUpAlgo_) delete pfPileUpAlgo_;
  if (clusterTools_) delete clusterTools_;
  if (mtsTransform_) delete mtsTransform_;
}

void ElectronMaker::beginJob() {
}

void ElectronMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void ElectronMaker::produce(Event& iEvent, const EventSetup& iSetup) {

    // access the tracker
    edm::ESHandle<TrackerGeometry> theTrackerGeometry;
    iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

    // Define vectors to be filled
    auto_ptr<unsigned int>   evt_nels(new unsigned int) ;

    // ECAL related (superCluster) variables
    auto_ptr<vector<int> >   els_nSeed       (new vector<int>   );
    auto_ptr<vector<float> > els_etaSC       (new vector<float> );
    auto_ptr<vector<float> > els_phiSC       (new vector<float> );
    auto_ptr<vector<float> > els_eSC         (new vector<float> );
    auto_ptr<vector<float> > els_eSCRaw      (new vector<float> );
    auto_ptr<vector<float> > els_eSCPresh    (new vector<float> );
    auto_ptr<vector<int> >   els_fiduciality (new vector<int>   );
    auto_ptr<vector<int> >   els_type        (new vector<int>   );
    auto_ptr<vector<int> >   els_scindex     (new vector<int>   ); 
    auto_ptr<vector<float> > els_etaSCwidth  (new vector<float> );
    auto_ptr<vector<float> > els_phiSCwidth  (new vector<float> );
    auto_ptr<vector<float> > els_ecalEnergy            (new vector<float>);
    auto_ptr<vector<float> > els_ecalEnergyError       (new vector<float>);
  
    // ID variables
    auto_ptr<vector<float> > els_lh                            (new vector<float> );
    auto_ptr<vector<float> > els_mva                           (new vector<float> );
    auto_ptr<vector<float> > els_dEtaIn                        (new vector<float> );
    auto_ptr<vector<float> > els_dEtaOut                       (new vector<float> );
    auto_ptr<vector<float> > els_dPhiIn                        (new vector<float> );
    auto_ptr<vector<float> > els_dPhiOut                       (new vector<float> );
    auto_ptr<vector<float> > els_dPhiInPhiOut                  (new vector<float> );
    auto_ptr<vector<float> > els_fbrem                         (new vector<float> );
    auto_ptr<vector<float> > els_eOverPIn                      (new vector<float> );
    auto_ptr<vector<float> > els_eSeedOverPOut                 (new vector<float> );
    auto_ptr<vector<float> > els_eSeedOverPIn                  (new vector<float> );
    auto_ptr<vector<float> > els_eOverPOut                     (new vector<float> );
    auto_ptr<vector<float> > els_hOverE                        (new vector<float> );
    auto_ptr<vector<float> > els_sigmaIPhiIPhi                 (new vector<float> );
    auto_ptr<vector<float> > els_sigmaIEtaIEta                 (new vector<float> );
    auto_ptr<vector<float> > els_e1x5                          (new vector<float> );
    auto_ptr<vector<float> > els_e5x5                          (new vector<float> );
    auto_ptr<vector<float> > els_e3x3                          (new vector<float> );

    // isolation variables
    auto_ptr<vector<float> > els_tkIso                  (new vector<float> );
    auto_ptr<vector<float> > els_ecalIso                (new vector<float> );
    auto_ptr<vector<float> > els_hcalIso                (new vector<float> );
    auto_ptr<vector<float> > els_ecalIso04              (new vector<float> );
    auto_ptr<vector<float> > els_hcalIso04              (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf_ch            (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf_gamma05       (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf_nhad05        (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf_ch            (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf_gamma05       (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf_nhad05        (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf2012_ch        (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf2012_em        (new vector<float> );
    auto_ptr<vector<float> > els_iso03_pf2012_nh        (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf2012_ch        (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf2012_em        (new vector<float> );
    auto_ptr<vector<float> > els_iso04_pf2012_nh        (new vector<float> );


    // track variables
    auto_ptr<vector<int> >   els_charge     (new vector<int>   );
    auto_ptr<vector<int> >   els_trk_charge (new vector<int>   );
    auto_ptr<vector<int> >   els_sccharge   (new vector<int>   );
    auto_ptr<vector<float> > els_chi2       (new vector<float> );
    auto_ptr<vector<float> > els_ndof       (new vector<float> );
    auto_ptr<vector<float> > els_d0         (new vector<float> );
    auto_ptr<vector<float> > els_z0         (new vector<float> );
    auto_ptr<vector<float> > els_d0Err      (new vector<float> );
    auto_ptr<vector<float> > els_z0Err      (new vector<float> );
    auto_ptr<vector<float> > els_d0corr     (new vector<float> );
    auto_ptr<vector<float> > els_z0corr     (new vector<float> );
    auto_ptr<vector<float> > els_ptErr      (new vector<float> );
    auto_ptr<vector<float> > els_etaErr     (new vector<float> );
    auto_ptr<vector<float> > els_phiErr     (new vector<float> );
    auto_ptr<vector<int>   > els_gsftrkidx  (new vector<int>   );
    auto_ptr<vector<float> > els_ip3d       (new vector<float> );
    auto_ptr<vector<float> > els_ip3derr    (new vector<float> );
  
    // LorentzVectors
    auto_ptr<vector<LorentzVector> > els_p4     (new vector<LorentzVector>);
    auto_ptr<vector<LorentzVector> > els_trk_p4 (new vector<LorentzVector>);

    // Vertex
    auto_ptr<vector<LorentzVector> > els_vertex_p4 (new vector<LorentzVector>);

    //HitPattern information
    auto_ptr<vector<int> >                    els_exp_innerlayers      (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_exp_outerlayers      (new vector<int>           ); 

    auto_ptr<vector<float> >                  els_trkshFrac            (new vector<float>         );
    auto_ptr<vector<int> >                    els_trkidx               (new vector<int>           );

    //conversion stuff
    auto_ptr<vector<float> >                  els_conv_dist            (new vector<float>          );
    auto_ptr<vector<float> >                  els_conv_dcot            (new vector<float>          );
    auto_ptr<vector<float> >                  els_conv_old_dist        (new vector<float>          );
    auto_ptr<vector<float> >                  els_conv_old_dcot        (new vector<float>          );


    // --- Get Input Collections --- //

    ////////////////
    // Get Tracks //
    ////////////////
   
    Handle<TrackCollection> tracks_h;
    iEvent.getByLabel(trksInputTag_, tracks_h);

  
    ////////////////
    // GSF Tracks //
    ////////////////

    Handle<GsfTrackCollection> gsftracks_h;
    iEvent.getByLabel(gsftracksInputTag_, gsftracks_h);


    /////////////
    // B Field //
    /////////////

    ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    float evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();


    ///////////////
    // Electrons //
    ///////////////

    Handle<View<GsfElectron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    View<GsfElectron> gsfElColl = *(els_h.product());

    Handle<GsfElectronCollection> els_coll_h;
    iEvent.getByLabel(electronsInputTag_, els_coll_h);    


    //////////////
    // PF Cands //
    //////////////

    iEvent.getByLabel(pfCandsInputTag, pfCand_h);

  
    ////////////
    // Vertex //
    ////////////

    iEvent.getByLabel(vtxInputTag, vertexHandle);


    ///////////////////////////
    // TransientTrackBuilder //
    ///////////////////////////
    ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);


    ////////////////////////////////////////////////
    // Get tools to get cluster shape information //
    ////////////////////////////////////////////////

    if ( clusterTools_ ) delete clusterTools_;
    clusterTools_ = new EcalClusterLazyTools( iEvent, iSetup, InputTag("reducedEcalRecHitsEB"), InputTag("reducedEcalRecHitsEE") );


    //////////////
    // Beamspot //
    //////////////

    Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpot_tag_, beamspot_h);
    const reco::BeamSpot &beamSpotreco = *(beamspot_h.product());

    // --- Fill --- //

    /////////////////////////
    // Loop Over Electrons //
    /////////////////////////

    //remove *evt_nels       = els_h->size();
    double mass     = 0.000510998918;
    size_t elsIndex = 0;
    for( View<GsfElectron>::const_iterator el = els_h->begin(); el != els_h->end(); el++, elsIndex++ ) {

        ////////////////
        // References //
        ////////////////

        const Track*                 el_track         = (const Track*)(el->gsfTrack().get());
        const RefToBase<GsfElectron> gsfElRef         = els_h->refAt(elsIndex);    

        //const TrackRef               ctfTkRef         = el->closestCtfTrackRef();
        const TrackRef               ctfTkRef         = el->closestTrack();
        const GsfTrackRef            gsfTkRef         = el->gsfTrack();
        const VertexCollection*      vertexCollection = vertexHandle.product();

        ////////////
        // Vertex //
        ////////////
        VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
        int firstGoodVertexIdx = 0;
        for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx) {
            if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
                firstGoodVertex = vtx;
                break;
            }
        }

        //////////////////////
        // Fiduciality Mask //
        //////////////////////

        int fiducialityMask = 0;  // the enum is in interface/EgammaFiduciality.h
        if ( el->isEB()        ) fiducialityMask |= 1 << ISEB;
        if ( el->isEBEEGap()   ) fiducialityMask |= 1 << ISEBEEGAP;
        if ( el->isEE()        ) fiducialityMask |= 1 << ISEE;
        if ( el->isEEGap()     ) fiducialityMask |= 1 << ISEEGAP;
        if ( el->isEBEtaGap()  ) fiducialityMask |= 1 << ISEBETAGAP;
        if ( el->isEBPhiGap()  ) fiducialityMask |= 1 << ISEBPHIGAP;
        if ( el->isEEDeeGap()  ) fiducialityMask |= 1 << ISEEDEEGAP;
        if ( el->isEERingGap() ) fiducialityMask |= 1 << ISEERINGGAP;
        if ( el->isGap()       ) fiducialityMask |= 1 << ISGAP;

        ///////////////////////////
        // Corrections & Seeding //
        ///////////////////////////

        int electronTypeMask = 0;
        if ( el->isEcalEnergyCorrected()        ) electronTypeMask |= 1 << ISECALENERGYCORRECTED;
        if ( el->trackerDrivenSeed()            ) electronTypeMask |= 1 << ISTRACKERDRIVEN;
        if ( el->ecalDrivenSeed()               ) electronTypeMask |= 1 << ISECALDRIVEN;
        if ( el->passingCutBasedPreselection()  ) electronTypeMask |= 1 << ISCUTPRESELECTED;
        if ( el->passingMvaPreselection()       ) electronTypeMask |= 1 << ISMVAPRESELECTED;

        /////////////////////
        // Lorentz Vectors //
        /////////////////////

        LorentzVector    p4In;
        LorentzVector    p4Out;
        LorentzVector    trk_p4( el_track->px(), el_track->py(), el_track->pz(), el_track->p() );
        math::XYZVectorF p3In  = el->trackMomentumAtVtx();
        math::XYZVectorF p3Out = el->trackMomentumOut();
        p4In.SetXYZT (   p3In.x() , p3In.y() , p3In.z() , sqrt( mass*mass + p3In.R() *p3In.R()  ) );
        p4Out.SetXYZT(   p3Out.x(), p3Out.y(), p3Out.z(), sqrt( mass*mass + p3Out.R()*p3Out.R() ) );

        //////////////
        // Electron //
        //////////////

        els_fiduciality        ->push_back( fiducialityMask                                 );
        els_type               ->push_back( electronTypeMask                                );
        els_ecalEnergy         ->push_back( el->correctedEcalEnergy()                       );  // energy corrections and uncertainties
        els_ecalEnergyError    ->push_back( el->correctedEcalEnergyError()                  );
        els_p4                 ->push_back( LorentzVector( el->p4() )                       );
        els_trk_p4             ->push_back( trk_p4                                          );
        els_vertex_p4          ->push_back( LorentzVector(el->vx(), el->vy(), el->vz(), 0.) );

        ///////////////
        // Isolation //
        ///////////////

        els_ecalIso               ->push_back( el->dr03EcalRecHitSumEt()                  );
        els_hcalIso               ->push_back( el->dr03HcalTowerSumEt()                   );
        els_tkIso                 ->push_back( el->dr03TkSumPt()                          );
        els_ecalIso04             ->push_back( el->dr04EcalRecHitSumEt()                  );
        els_hcalIso04             ->push_back( el->dr04HcalTowerSumEt()                   );

        //////////////////
        // PF Isolation //
        //////////////////

        if ( firstGoodVertex!=vertexCollection->end() ) {

            els_iso03_pf_ch      -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 99999., 0.1, 0.07, 0.025, 0.025, 0  ) );
            els_iso03_pf_gamma05 -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 0.5   , 0.1, 0.07, 0.025, 0.025, 22 ) );
            els_iso03_pf_nhad05  -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.3, 0.5   , 0.1, 0.07, 0.025, 0.025, 130) );

            els_iso04_pf_ch      -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4, 99999., 0.1, 0.07, 0.025, 0.025, 0  ) );
            els_iso04_pf_gamma05 -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4, 0.5   , 0.1, 0.07, 0.025, 0.025, 22 ) );
            els_iso04_pf_nhad05  -> push_back( electronIsoValuePF( *el, *firstGoodVertex, 0.4,  0.5  , 0.1, 0.07, 0.025, 0.025, 130) );


            // pf iso 2012
            float pfiso_ch = 0.0;
            float pfiso_em = 0.0;
            float pfiso_nh = 0.0;
            PFIsolation2012(*el, vertexCollection, firstGoodVertexIdx, 0.3, pfiso_ch, pfiso_em, pfiso_nh);
            els_iso03_pf2012_ch ->push_back( pfiso_ch );
            els_iso03_pf2012_em ->push_back( pfiso_em );
            els_iso03_pf2012_nh ->push_back( pfiso_nh );

            PFIsolation2012(*el, vertexCollection, firstGoodVertexIdx, 0.4, pfiso_ch, pfiso_em, pfiso_nh);
            els_iso04_pf2012_ch ->push_back( pfiso_ch );
            els_iso04_pf2012_em ->push_back( pfiso_em );
            els_iso04_pf2012_nh ->push_back( pfiso_nh );	    

        } else {

            els_iso03_pf_ch      -> push_back( -9999. );
            els_iso03_pf_gamma05 -> push_back( -9999. );
            els_iso03_pf_nhad05  -> push_back( -9999. );

            els_iso04_pf_ch      -> push_back( -9999. );
            els_iso04_pf_gamma05 -> push_back( -9999. );
            els_iso04_pf_nhad05  -> push_back( -9999. );

            els_iso03_pf2012_ch ->push_back( -9999. );
            els_iso03_pf2012_em ->push_back( -9999. );
            els_iso03_pf2012_nh ->push_back( -9999. );
            els_iso04_pf2012_ch ->push_back( -9999. );
            els_iso04_pf2012_em ->push_back( -9999. );
            els_iso04_pf2012_nh ->push_back( -9999. );
        }

        //////////////////
        // Supercluster //
        //////////////////

        els_etaSC         ->push_back( el->superCluster()->eta()             );
        els_phiSC         ->push_back( el->superCluster()->phi()             );
        els_eSC           ->push_back( el->superCluster()->energy()          );
        els_eSCRaw        ->push_back( el->superCluster()->rawEnergy()       );
        els_eSCPresh      ->push_back( el->superCluster()->preshowerEnergy() );
        els_nSeed         ->push_back( el->basicClustersSize() - 1           );
        els_e1x5          ->push_back( el->e1x5()                            );
        els_e5x5          ->push_back( el->e5x5()                            );
        els_sigmaIEtaIEta ->push_back( el->sigmaIetaIeta()                   );
        els_etaSCwidth    ->push_back( el->superCluster()->etaWidth()        );
        els_phiSCwidth    ->push_back( el->superCluster()->phiWidth()        );


        ///////////////////////////////////////////////////////
        // Get cluster info that is not stored in the object //
        ///////////////////////////////////////////////////////

        if( el->superCluster()->seed().isAvailable() ) { 

            //
            const BasicCluster&  clRef              = *(el->superCluster()->seed());
            const vector<float>& lcovs              = clusterTools_->localCovariances(clRef);                    // get the local covariances computed in a 5x5 around the seed
            const vector<float>  localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));  // get the local covariances computed using all crystals in the SC

            //
            els_sigmaIPhiIPhi   ->push_back( isfinite(lcovs[2])              ? lcovs[2] > 0               ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2])                             : -9999. );

            //
            els_e3x3            ->push_back( clusterTools_->e3x3(clRef) );
        } 
        else {

            //
            els_sigmaIPhiIPhi   ->push_back(-9999.);

            //
            els_e3x3            ->push_back(-9999.);

        } //
 
   
        ////////
        // ID //
        ////////

        els_hOverE                        ->push_back( el->hcalOverEcal()                   );
        els_eOverPIn                      ->push_back( el->eSuperClusterOverP()             );
        els_eSeedOverPOut                 ->push_back( el->eSeedClusterOverPout()           );
        els_eSeedOverPIn                  ->push_back( el->eSeedClusterOverP()              );
        els_eOverPOut                     ->push_back( el->eEleClusterOverPout()            );
        els_fbrem                         ->push_back( el->fbrem()                          );

        els_dEtaIn                        ->push_back( el->deltaEtaSuperClusterTrackAtVtx() );
        els_dEtaOut                       ->push_back( el->deltaEtaSeedClusterTrackAtCalo() );
        els_dPhiIn                        ->push_back( el->deltaPhiSuperClusterTrackAtVtx() );
        els_dPhiOut                       ->push_back( el->deltaPhiSeedClusterTrackAtCalo() );

        ////////////
        // Tracks //
        ////////////

        float pt       = el_track->pt();
        float p        = el_track->p();
        float q        = el_track->charge();
        float pz       = el_track->pz();
        float trkpterr = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))+2*pt*p/q*pz*(el_track->covariance(0,1))+ pz*pz*(el_track->covariance(1,1) ) ) : -9999.;
            
        els_chi2                  ->push_back( el_track->chi2()                          );
        els_ndof                  ->push_back( el_track->ndof()                          );
        els_d0Err                 ->push_back( el_track->d0Error()                       );
        els_z0Err                 ->push_back( el_track->dzError()                       );
        els_ptErr                 ->push_back( trkpterr                                  );
        els_etaErr                ->push_back( el_track->etaError()                      );
        els_phiErr                ->push_back( el_track->phiError()                      );  
        els_gsftrkidx             ->push_back( static_cast<int>((el->gsfTrack()).key())  );
        els_charge                ->push_back( el->charge()                              );
        els_trk_charge            ->push_back( el_track->charge()                        );
        els_sccharge              ->push_back( el->scPixCharge()                         );
        els_d0                    ->push_back( el_track->d0()                            );
        els_z0                    ->push_back( el_track->dz()                            );
        els_d0corr                ->push_back( -1*(el_track->dxy(beamSpotreco))              );
        els_z0corr                ->push_back( el_track->dz(beamSpotreco.position(el_track->vz()))                    );
   

        /////////
        // CTF //
        /////////

        if( ctfTkRef.isNonnull() ) {
            els_trkidx    -> push_back( static_cast<int>  ( ctfTkRef.key()        )                                  );
            els_trkshFrac -> push_back( static_cast<float>( el->ctfGsfOverlap() )                                    );
        } 
        else {
            els_trkidx    -> push_back(-9999.);
            els_trkshFrac -> push_back(-9999.);
        }

        
        ////////////////////
        // Regular Vertex //
        ////////////////////        
        TransientTrack tt = theTTBuilder->build(el->gsfTrack());
    
        if ( firstGoodVertex!=vertexCollection->end() ) {
            Measurement1D ip3D_regular = IPTools::absoluteImpactParameter3D(tt, *firstGoodVertex).second;
            //
            els_ip3d      -> push_back( ip3D_regular.value() );
            els_ip3derr   -> push_back( ip3D_regular.error() );
        } else {
            //
            els_ip3d      -> push_back( -999. );
            els_ip3derr   -> push_back( -999. );
        }


        /////////////////
        // Hit Pattern //
        /////////////////

        if( el_track->extra().isAvailable() ) {
        } else {
        }
    
        const HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
        const HitPattern& p_outer = el_track->trackerExpectedHitsOuter();

        els_exp_innerlayers -> push_back(p_inner.numberOfHits());
        els_exp_outerlayers -> push_back(p_outer.numberOfHits());



        /////////////////
        // Conversions //
        /////////////////

        ConversionFinder convFinder; //vector of conversion infos - all the candidate conversion tracks
        vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(el->core()), tracks_h, gsftracks_h, evt_bField);
    
        vector<int>           v_tkidx;
        vector<int>           v_gsftkidx;
        vector<int>           v_delmisshits;
        vector<int>           v_flag;
        vector<float>         v_dist;
        vector<float>         v_dcot;
        vector<float>         v_rad;
        vector<LorentzVector> v_pos_p4;

        for(unsigned int i_conv = 0; i_conv < v_convInfos.size(); i_conv++) {
      
            //
            math::XYZPoint convPoint  = v_convInfos.at(i_conv).pointOfConversion();
            float          convPointX = isfinite(convPoint.x()) ? convPoint.x() : -9999.;
            float          convPointY = isfinite(convPoint.y()) ? convPoint.y() : -9999.;
            float          convPointZ = isfinite(convPoint.z()) ? convPoint.z() : -9999.;

            //
            v_dist        .push_back( isfinite(v_convInfos.at(i_conv).dist()) ? v_convInfos.at(i_conv).dist() : -9999.  );
            v_dcot        .push_back( v_convInfos.at(i_conv).dcot()                                                     );
            v_rad         .push_back( v_convInfos.at(i_conv).radiusOfConversion()                                       );
            v_delmisshits .push_back( v_convInfos.at(i_conv).deltaMissingHits()                                         );
            v_flag        .push_back( v_convInfos.at(i_conv).flag()                                                     );
            v_pos_p4      .push_back( LorentzVector(convPointX, convPointY, convPointZ, 0)                              );

            //
            if( v_convInfos.at(i_conv).conversionPartnerCtfTk().isNonnull() ) {
                v_tkidx.push_back(v_convInfos.at(i_conv).conversionPartnerCtfTk().key());
            }
            else {
                v_tkidx.push_back(-9999);
            }

            //
            if( v_convInfos.at(i_conv).conversionPartnerGsfTk().isNonnull() ) {
                v_gsftkidx.push_back(v_convInfos.at(i_conv).conversionPartnerGsfTk().key());
            }
            else { 
                v_gsftkidx.push_back(-9999);
            }

        } // end for loop


        //
        ConversionInfo convInfo   = convFinder.getConversionInfo( *el, tracks_h, gsftracks_h, evt_bField );

        //
        els_conv_dist        -> push_back( isfinite(convInfo.dist()) ? convInfo.dist() : -9999. );
        els_conv_dcot        -> push_back( convInfo.dcot()                                      );


        //////////////////////////////
        // Flag For Vertex Fit Conversion Rejection //
        //////////////////////////////

        Handle<ConversionCollection> convs_h;
        iEvent.getByLabel(recoConversionInputTag_, convs_h);

        //////////////////////////////
        // Old Conversion Rejection //
        //////////////////////////////


        els_conv_old_dist        -> push_back( isfinite(el->convDist())   ? el->convDist()   : -9999. );
        els_conv_old_dcot        -> push_back( isfinite(el->convDcot())   ? el->convDcot()   : -9999. );

        //////////////////////
        // 2012 Electron ID //
        //////////////////////
        GsfElectronRef ele(els_coll_h, elsIndex);

    } // end Loop on Electrons
  

    // Put the results into the event

    // Track parameters
    //
    iEvent.put(els_d0         , "elsd0"        );
    iEvent.put(els_z0         , "elsz0"        );
    iEvent.put(els_d0corr     , "elsd0corr"    );
    iEvent.put(els_z0corr     , "elsz0corr"    );
    iEvent.put(els_chi2       , "elschi2"      );
    iEvent.put(els_ndof       , "elsndof"      );
    iEvent.put(els_d0Err      , "elsd0Err"     );
    iEvent.put(els_z0Err      , "elsz0Err"     );
    iEvent.put(els_ptErr      , "elsptErr"     );
    iEvent.put(els_etaErr     , "elsetaErr"    );
    iEvent.put(els_phiErr     , "elsphiErr"    );
    iEvent.put(els_gsftrkidx  , "elsgsftrkidx" );
    iEvent.put(els_ip3d       , "elsip3d"      );
    iEvent.put(els_ip3derr    , "elsip3derr"   );
  
    iEvent.put(els_charge     , "elscharge"    );
    iEvent.put(els_sccharge   , "elssccharge"  );
    iEvent.put(els_trk_charge , "elstrkcharge" );

    // Supercluster parameters
    //
    iEvent.put(els_nSeed       , "elsnSeed"       );
    iEvent.put(els_etaSC       , "elsetaSC"       );
    iEvent.put(els_phiSC       , "elsphiSC"       );
    iEvent.put(els_eSC         , "elseSC"         );
    iEvent.put(els_eSCRaw      , "elseSCRaw"      );
    iEvent.put(els_eSCPresh    , "elseSCPresh"    );
    iEvent.put(els_e1x5        , "else1x5"        );
    iEvent.put(els_e3x3        , "else3x3"        );
    iEvent.put(els_e5x5        , "else5x5"        );
    iEvent.put(els_fiduciality , "elsfiduciality" );
    iEvent.put(els_type        , "elstype"        );
    iEvent.put(els_scindex     , "elsscindex"     );
    iEvent.put(els_etaSCwidth  , "elsetaSCwidth"  );
    iEvent.put(els_phiSCwidth  , "elsphiSCwidth"  );


    // Electron ID
    //
    iEvent.put(els_sigmaIPhiIPhi      , "elssigmaIPhiIPhi"      );
    iEvent.put(els_sigmaIEtaIEta      , "elssigmaIEtaIEta"      );
    iEvent.put(els_hOverE             , "elshOverE"             );

    iEvent.put(els_eOverPIn                      , "elseOverPIn"                      );
    iEvent.put(els_eSeedOverPOut                 , "elseSeedOverPOut"                 );
    iEvent.put(els_eSeedOverPIn                  , "elseSeedOverPIn"                  );
    iEvent.put(els_eOverPOut                     , "elseOverPOut"                     );
    iEvent.put(els_fbrem                         , "elsfbrem"                         );
    iEvent.put(els_lh                            , "elslh"                            );
    iEvent.put(els_dEtaIn                        , "elsdEtaIn"                        );
    iEvent.put(els_dEtaOut                       , "elsdEtaOut"                       );
    iEvent.put(els_dPhiIn                        , "elsdPhiIn"                        );
    iEvent.put(els_dPhiOut                       , "elsdPhiOut"                       );

    // Lorentz vectors
    //
    iEvent.put(els_p4     , "elsp4"    );
    iEvent.put(els_trk_p4 , "elstrkp4" );

    // Vertex
    //
    iEvent.put(els_vertex_p4, "elsvertexp4");
    iEvent.put(els_ecalEnergy         , "elsecalEnergy"         );

    // Isolation
    //
    iEvent.put(els_tkIso                , "elstkIso"                );
    iEvent.put(els_ecalIso              , "elsecalIso"              );
    iEvent.put(els_hcalIso              , "elshcalIso"              );
    iEvent.put(els_ecalIso04              , "elsecalIso04"              );
    iEvent.put(els_hcalIso04              , "elshcalIso04"              );

    iEvent.put(els_iso03_pf_ch      , "elsiso03pfch"      );
    iEvent.put(els_iso03_pf_gamma05 , "elsiso03pfgamma05" );
    iEvent.put(els_iso03_pf_nhad05  , "elsiso03pfnhad05"  );
    iEvent.put(els_iso04_pf_ch      , "elsiso04pfch"      );
    iEvent.put(els_iso04_pf_gamma05 , "elsiso04pfgamma05" );
    iEvent.put(els_iso04_pf_nhad05  , "elsiso04pfnhad05"  );

    iEvent.put(els_iso03_pf2012_ch , "elsiso03pf2012ch" );
    iEvent.put(els_iso03_pf2012_em , "elsiso03pf2012em" );
    iEvent.put(els_iso03_pf2012_nh , "elsiso03pf2012nh" );
    iEvent.put(els_iso04_pf2012_ch , "elsiso04pf2012ch" );
    iEvent.put(els_iso04_pf2012_em , "elsiso04pf2012em" );
    iEvent.put(els_iso04_pf2012_nh , "elsiso04pf2012nh" );

    //Hit Pattern Information
    iEvent.put(els_exp_innerlayers , "elsexpinnerlayers" );
    iEvent.put(els_exp_outerlayers , "elsexpouterlayers" );

    //CTF track info
    //
    iEvent.put(els_trkidx          , "elstrkidx"        );
    iEvent.put(els_trkshFrac       , "elstrkshFrac"     );

    //conversion
    iEvent.put(els_conv_dist        , "elsconvdist"        );
    iEvent.put(els_conv_dcot        , "elsconvdcot"        );
    iEvent.put(els_conv_old_dist        , "elsconvolddist"        );
    iEvent.put(els_conv_old_dcot        , "elsconvolddcot"        );
}

double ElectronMaker::electronIsoValuePF(const GsfElectron& el, const Vertex& vtx, float coner, float minptn, float dzcut,
                                         float footprintdr, float gammastripveto, float elestripveto, int filterId){

    float pfciso = 0.;
    float pfniso = 0.;
    float pffootprint = 0.;
    float pfjurveto = 0.;
    float pfjurvetoq = 0.;

    //TrackRef siTrack     = el.closestCtfTrackRef();
    TrackRef siTrack     = el.closestTrack();
    GsfTrackRef gsfTrack = el.gsfTrack();

    if (gsfTrack.isNull() && siTrack.isNull()) return -9999.;

    float eldz = gsfTrack.isNonnull() ? gsfTrack->dz(vtx.position()) : siTrack->dz(vtx.position());
    float eleta = el.eta();

    for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

        float pfeta = pf->eta();    
        float dR = deltaR(pfeta, pf->phi(), eleta, el.phi());
        if (dR>coner) continue;

        float deta = fabs(pfeta - eleta);
        int pfid = abs(pf->pdgId());
        float pfpt = pf->pt();

        if (filterId!=0 && filterId!=pfid) continue;

        if (pf->charge()==0) {
            //neutrals
            if (pfpt>minptn) {
                pfniso+=pfpt;
                if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
                if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
            }
        } else {
            //charged  
            //avoid double counting of electron itself
            //if either the gsf or the ctf track are shared with the candidate, skip it
            const TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //below pfid==1 is commented out: in some cases the pfCand has a gsf even if it is not an electron... this is to improve the sync with MIT
            if (/*pfid==11 &&*/ pf->gsfTrackRef().isNonnull()) {
                if (gsfTrack.isNonnull() && gsfTrack.key()==pf->gsfTrackRef().key()) continue;
            } 
            //check electrons with gsf track
            if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if(fabs( pfTrack->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
            }
        } 
    }
    return pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq;
}


void ElectronMaker::PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection,
        const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh)
{

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0;
    pfiso_nh = 0.0;

    // loop on pfcandidates
    reco::PFCandidateCollection::const_iterator pf = pfCand_h->begin();
    for (pf = pfCand_h->begin(); pf != pfCand_h->end(); ++pf) {

        // skip electrons and muons
        if (pf->particleId() == reco::PFCandidate::e)     continue;
        if (pf->particleId() == reco::PFCandidate::mu)    continue;

        // deltaR between electron and cadidate
        const float dR = deltaR(pf->eta(), pf->phi(), el.eta(), el.phi());
        if (dR > R)                             continue;

        // charged hadrons closest vertex
        // should be the primary vertex
        if (pf->particleId() == reco::PFCandidate::h) {
            int pfVertexIndex = pfPileUpAlgo_->chargedHadronVertex(*vertexCollection, *pf);
            if (pfVertexIndex != vertexIndex) continue;
        }

        // endcap region
        if (!el.isEB()) {
            if (pf->particleId() == reco::PFCandidate::h      && dR <= 0.015)   continue;
            if (pf->particleId() == reco::PFCandidate::gamma  && dR <= 0.08)    continue;
        }

        // add to isolation sum
        if (pf->particleId() == reco::PFCandidate::h)       pfiso_ch += pf->pt();
        if (pf->particleId() == reco::PFCandidate::gamma)   pfiso_em += pf->pt();
        if (pf->particleId() == reco::PFCandidate::h0)      pfiso_nh += pf->pt();

    }

}


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
