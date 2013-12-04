#ifndef BTAGMAKER_H
#define BTAGMAKER_H

//
// class decleration
//

class BTagMaker : public edm::EDProducer {
public:
     explicit BTagMaker (const edm::ParameterSet&);
      ~BTagMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  edm::RefToBase<reco::Jet> getReferenceJetRef(const edm::View<reco::Jet>*, const reco::Jet*);
  // ----------member data ---------------------------
  edm::InputTag CaloJetsTag_                          ;
  edm::InputTag referenceCaloJetsTag_                 ;
  std::string   aliasprefix_                          ;
  edm::InputTag combinedSecondaryVertexBJetTags_      ;
  edm::InputTag combinedSecondaryVertexMVABJetTags_   ;
  edm::InputTag ghostTrackBJetTags_                   ;
  edm::InputTag jetBProbabilityBJetTags_              ;
  edm::InputTag jetProbabilityBJetTags_               ;
  edm::InputTag simpleSecondaryVertexBJetTags_        ;
  edm::InputTag simpleSecondaryVertexHighEffBJetTags_ ;
  edm::InputTag simpleSecondaryVertexHighPurBJetTags_ ;
  edm::InputTag softElectronTags_                     ;
  edm::InputTag softElectronByIP3dBJetTags_           ;
  edm::InputTag softElectronByPtBJetTags_             ;
  edm::InputTag softMuonBJetTags_                     ;
  edm::InputTag softMuonByIP3dBJetTags_               ;
  edm::InputTag softMuonByPtBJetTags_                 ;
  edm::InputTag trackCountingHighEffBJetTags_         ;
  edm::InputTag trackCountingHighPurBJetTags_         ; 
  
  
};

#endif

