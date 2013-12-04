#ifndef NTUPLEMAKER_ELECTRONMAKER_H
#define NTUPLEMAKER_ELECTRONMAKER_H


//
// class decleration
//

class ElectronMaker : public edm::EDProducer {
public:
    explicit ElectronMaker (const edm::ParameterSet&);
    ~ElectronMaker();

private:
//  virtual void beginJob() ;
    virtual void beginJob() ;
    //virtual void beginJob(const edm::EventSetup&) ;
    //virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    double electronIsoValuePF(const reco::GsfElectron& el, const reco::Vertex& vtx, float coner, float minptn, float dzcut,
                              float footprintdr, float gammastripveto, float elestripveto, int filterId);
  
    int classify(const edm::RefToBase<reco::GsfElectron> &);
    template<typename T> const edm::ValueMap<T>& getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag);
 
    // for 2012 pf isolation
    void PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection, 
                         const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh);
 
    // ----------member data ---------------------------
    edm::InputTag electronsInputTag_;
    //edm::InputTag beamSpotInputTag_;
    edm::InputTag trksInputTag_;
    edm::InputTag gsftracksInputTag_;
    edm::InputTag eidLHTag_;
    edm::InputTag cms2scsseeddetidInputTag_;
    edm::InputTag pfCandsInputTag;
    edm::InputTag vtxInputTag;

  edm::InputTag pfIsoCharged03InputTag;
  edm::InputTag pfIsoGamma03InputTag;
  edm::InputTag pfIsoNeutral03InputTag;
  edm::InputTag pfIsoCharged04InputTag;
  edm::InputTag pfIsoGamma04InputTag;
  edm::InputTag pfIsoNeutral04InputTag;

    edm::InputTag recoConversionInputTag_;

    EcalClusterLazyTools* clusterTools_;
    MultiTrajectoryStateTransform *mtsTransform_;

    double minAbsDist_;
    double minAbsDcot_;
    double minSharedFractionOfHits_;
    std::string aliasprefix_;

    edm::Handle<reco::PFCandidateCollection> pfCand_h;
    edm::Handle<reco::VertexCollection> vertexHandle;

    edm::InputTag rhoInputTag_;
    edm::InputTag beamSpot_tag_;

    PFPileUpAlgo *pfPileUpAlgo_;
};

#endif

