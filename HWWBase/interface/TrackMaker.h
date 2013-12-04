#ifndef TRACKMAKER_H
#define TRACKMAKER_H

typedef math::XYZPoint Point;

class TrackMaker : public edm::EDProducer {
public:
  explicit TrackMaker (const edm::ParameterSet&);
  double calculateTrkIsolation(const edm::View<reco::Track>*, const reco::Track&, const Point&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag trackInputTag;
  edm::InputTag primaryVertexInputTag_;
  
  float dRConeMin_;
  float dRConeMax_;
  float vtxDiffZMax_;
  float tkVtxDMax_;
  float ptMin_;
  int   nHits_;
  std::string aliasprefix_;
};

#endif
