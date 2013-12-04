#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HWWValidation/HWWBase/interface/SCMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

//
// constructors and destructor
//
SCMaker::SCMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    produces<std::vector<LorentzVector> >(branchprefix+"posp4").setBranchAlias(aliasprefix_+"_pos_p4");
    produces<std::vector<float> >(branchprefix+"eMax").setBranchAlias(aliasprefix_+"_eMax");
    produces<std::vector<float> >(branchprefix+"e1x3").setBranchAlias(aliasprefix_+"_e1x3");
    produces<std::vector<float> >(branchprefix+"e3x1").setBranchAlias(aliasprefix_+"_e3x1"); 
    produces<std::vector<float> >(branchprefix+"sigmaIEtaIPhi").setBranchAlias(aliasprefix_+"_sigmaIEtaIPhi");

    // add superclusters to the ntuple if they have ET > scEtMin_
    scEtMin_ = iConfig.getParameter<double>("scEtMin");

    // input tags for superclusters
    scInputTag_EE_ = iConfig.getParameter<edm::InputTag>("scInputTag_EE");
    scInputTag_EB_ = iConfig.getParameter<edm::InputTag>("scInputTag_EB");
    scInputTags_.clear();
    scInputTags_.push_back(scInputTag_EE_);
    scInputTags_.push_back(scInputTag_EB_);

    hitInputTags_.clear();
    ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
    ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
    hitInputTags_.push_back(ecalRecHitsInputTag_EE_);
    hitInputTags_.push_back(ecalRecHitsInputTag_EB_);

    // other input tags
    hcalRecHitsInputTag_HBHE_ = iConfig.getParameter<edm::InputTag>("hcalRecHitsInputTag_HBHE");
    primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    electronsInputTag_ = iConfig.getParameter<edm::InputTag>("electronsInputTag");

    MCTruthCollection_ = iConfig.getParameter<edm::InputTag>("MCTruthCollection");

    // initialise this
    cachedCaloGeometryID_ = 0;
    clusterTools_ = 0;

}

void SCMaker::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{

    iSetup.get<EcalLaserDbRecord>().get(laser_);

}

void SCMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // get the calo geometry
    if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
        cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
        iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
    }

    // get hcal rechits
    edm::Handle<HBHERecHitCollection> hcalRecHitsHandle;
    try {
        iEvent.getByLabel(hcalRecHitsInputTag_HBHE_, hcalRecHitsHandle);
    }
    catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") << "Error! can't get the HCAL Hits";
    }

    // get hcal rechit metacollection 
    HBHERecHitMetaCollection *mhbhe = 0;
    if (!hcalRecHitsHandle.failedToGet()) {
        mhbhe =  new HBHERecHitMetaCollection(*hcalRecHitsHandle);
    }

    // get the primary vertices
    edm::Handle<reco::VertexCollection> vertexHandle;
    try {
        iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
    }
    catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") << "Error! can't get the primary vertex";
    }
    const reco::VertexCollection *vertexCollection = vertexHandle.product();
    Point pv(0.0, 0.0, 0.0);
    if (vertexCollection->size() > 0) {
        pv = vertexCollection->at(0).position();
    }

    // get the electrons (for matching)
    edm::Handle<reco::GsfElectronCollection> electronsHandle;
    try {
        iEvent.getByLabel(electronsInputTag_, electronsHandle);
    }
    catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") << "Error! can't get the electrons";
    }

    // get hoe variable
    HoECalculator hoeCalc(caloGeometry_);



    std::auto_ptr<std::vector<LorentzVector> > vector_scs_pos_p4 (new std::vector<LorentzVector>);
    std::auto_ptr<std::vector<float> > vector_scs_eMax (new std::vector<float>);
    std::auto_ptr<std::vector<float> > vector_scs_e1x3 (new std::vector<float>);
    std::auto_ptr<std::vector<float> > vector_scs_e3x1 (new std::vector<float>);
    std::auto_ptr<std::vector<float> > vector_scs_sigmaIEtaIPhi(new std::vector<float>);

    edm::Handle<EcalRecHitCollection> tempH;
    bool haveHits = true;
    for (unsigned int i = 0; i < scInputTags_.size(); ++i) {
        iEvent.getByLabel(hitInputTags_[i], tempH);
        if(tempH.failedToGet()) {
            haveHits = false;
            break;
        }
    }

    if(clusterTools_) delete clusterTools_; 
    if(haveHits) 
        clusterTools_ = new EcalClusterLazyTools(iEvent, iSetup,
                ecalRecHitsInputTag_EB_, ecalRecHitsInputTag_EE_);

    // there are multiple supercluster collections. In the ntuple
    // these will become concatonated
    for (unsigned int i = 0; i < scInputTags_.size(); ++i)
    {

        // get superclusters
        edm::Handle<reco::SuperClusterCollection> scHandle;
        try {
            iEvent.getByLabel(scInputTags_[i], scHandle);
        }
        catch ( cms::Exception& ex ) {
            edm::LogError("SCMakerError") << "Error! can't get the SuperClusters";
        }
        const reco::SuperClusterCollection *scCollection = scHandle.product();

        // get hits
        edm::Handle<EcalRecHitCollection> rhcHandle;
        iEvent.getByLabel(hitInputTags_[i], rhcHandle);
        const EcalRecHitCollection *recHits;
        if(haveHits) //has been determined beforehand
            recHits = rhcHandle.product();
        else 
            recHits = NULL;


        size_t scIndex = 0;
        for (reco::SuperClusterCollection::const_iterator sc = scCollection->begin();
                sc != scCollection->end(); ++sc, ++scIndex) {

            // do ET cut
            if ( (sc->energy()/cosh(sc->eta())) < scEtMin_) continue;

            vector_scs_pos_p4->push_back( LorentzVector(sc->position().x(), sc->position().y(), sc->position().z(), 0.) );

            const std::vector<std::pair<DetId, float > > detIds = sc->hitsAndFractions() ;
            if(haveHits) {

                EcalRecHitCollection::const_iterator maxHit;
                float eMax = 0.0;
                float eSumNoLaser = 0.0;
                for (unsigned int i = 0; i < detIds.size(); ++i) {
                    EcalRecHitCollection::const_iterator hit = recHits->find(detIds[i].first);
                    if (hit != recHits->end()) {
                        eSumNoLaser += hit->energy() / laser_->getLaserCorrection(hit->id(), iEvent.time());
                        if (hit->energy() > eMax) {
                            eMax = hit->energy();
                            maxHit = hit;
                        }
                    }
                }


            }

            if(haveHits) {

                vector_scs_eMax->push_back( clusterTools_->eMax(*(sc->seed())) );
                vector_scs_e1x3->push_back( clusterTools_->e1x3(*(sc->seed())) );
                vector_scs_e3x1->push_back( clusterTools_->e3x1(*(sc->seed())) );

                // get the covariances computed in 5x5 around the seed
                std::vector<float> covariances = clusterTools_->covariances(*(sc->seed()));
                // get the local covariances computed in a 5x5 around the seed
                const std::vector<float> localCovariances = clusterTools_->localCovariances(*(sc->seed()));
                // get the local covariances computed using all crystals in the SC
                const std::vector<float> localCovariancesSC = clusterTools_->scLocalCovariances(*sc);

                // if seed basic cluster is in the endcap then correct sigma eta eta
                // according to the super cluster eta
                if(fabs(sc->seed()->eta()) > 1.479) 
                    covariances[0] -= 0.02*(fabs(sc->eta()) - 2.3);

                vector_scs_sigmaIEtaIPhi->push_back	 ( std::isfinite(localCovariances[1]  ) ? (localCovariances[1]   > 0 ? sqrt(localCovariances[1]  ) : -1 * sqrt(-1 * localCovariances[1]  ) ) : -9999.);
            } else {

                //replace everything in this section with -9999

                vector_scs_eMax->push_back( -9999 );
                vector_scs_e1x3->push_back( -9999 );
                vector_scs_e3x1->push_back( -9999 );
                vector_scs_sigmaIEtaIPhi->push_back	 ( -9999. );
            }

        } // end loop on scs

    } // end loop on sc input tags

    //*evt_nscs = vector_scs_p4->size();

    // put results into the event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(vector_scs_pos_p4		, branchprefix+"posp4"			);
    iEvent.put(vector_scs_eMax		, branchprefix+"eMax"			);
    iEvent.put(vector_scs_e1x3		, branchprefix+"e1x3"			);
    iEvent.put(vector_scs_e3x1		, branchprefix+"e3x1"			);
    iEvent.put(vector_scs_sigmaIEtaIPhi	, branchprefix+"sigmaIEtaIPhi"		);

    delete mhbhe;

}

math::XYZTLorentzVectorF SCMaker::initP4(const math::XYZPoint &pvPos, 
        const reco::SuperCluster &sc) {

    math::XYZVector scPos(sc.x(), sc.y(), sc.z());
    math::XYZVector pvPosVec(pvPos.x(), pvPos.y(), pvPos.z());
    math::XYZVector objPosition = scPos - pvPosVec;
    double scale = sc.energy() / objPosition.R();
    return math::XYZTLorentzVectorF(objPosition.x() * scale, 
            objPosition.y() * scale, 
            objPosition.z() * scale, 
            sc.energy());
}



// ------------ method called once each job just before starting event loop  ------------
    void 
SCMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
    void 
SCMaker::endJob() 
{
}

//
// Closest MC Particle
//
void SCMaker::closestMCParticle(const HepMC::GenEvent *genEvent, const reco::SuperCluster &sc,
        double &dRClosest, double &energyClosest)
{

    // SuperCluster eta, phi
    double scEta = sc.eta();
    double scPhi = sc.phi();

    // initialize dRClosest to a large number
    dRClosest = 999.9;

    // loop over the MC truth particles to find the
    // closest to the superCluster in dR space
    for(HepMC::GenEvent::particle_const_iterator currentParticle = genEvent->particles_begin();
            currentParticle != genEvent->particles_end(); currentParticle++ )
    {
        if((*currentParticle)->status() == 3 && abs((*currentParticle)->pdg_id()) == 11)
        {
            // need GenParticle in ECAL co-ordinates
            HepMC::FourVector vtx = (*currentParticle)->production_vertex()->position();
            double phiTrue = (*currentParticle)->momentum().phi();
            double etaTrue = ecalEta((*currentParticle)->momentum().eta(), vtx.z()/10., vtx.perp()/10.);

            double dPhi = reco::deltaPhi(phiTrue, scPhi);
            double dEta = scEta - etaTrue;
            double deltaR = std::sqrt(dPhi*dPhi + dEta*dEta);

            if(deltaR < dRClosest)
            {
                dRClosest = deltaR;
                energyClosest = (*currentParticle)->momentum().e();
            }

        } // end if stable particle     

    } // end loop on get particles

}

//
// Compute Eta in the ECAL co-ordinate system
//
float SCMaker::ecalEta(float EtaParticle , float Zvertex, float plane_Radius)
{
    const float R_ECAL           = 136.5;
    const float Z_Endcap         = 328.0;
    const float etaBarrelEndcap  = 1.479;

    if(EtaParticle != 0.)
    {
        float Theta = 0.0  ;
        float ZEcal = (R_ECAL-plane_Radius)*sinh(EtaParticle)+Zvertex;

        if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
        if(Theta<0.0) Theta = Theta+Geom::pi() ;

        float ETA = - log(tan(0.5*Theta));

        if( fabs(ETA) > etaBarrelEndcap )
        {
            float Zend = Z_Endcap ;
            if(EtaParticle<0.0 )  Zend = -Zend ;
            float Zlen = Zend - Zvertex ;
            float RR = Zlen/sinh(EtaParticle);
            Theta = atan((RR+plane_Radius)/Zend);
            if(Theta<0.0) Theta = Theta+Geom::pi() ;
            ETA = - log(tan(0.5*Theta));
        }

        return ETA;
    }
    else
    {
        edm::LogWarning("")  << "[EgammaSuperClusters::ecalEta] Warning: Eta equals to zero, not correcting" ;
        return EtaParticle;
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCMaker);
