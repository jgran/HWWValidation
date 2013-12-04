#ifndef GSFTRACKMAKER_H
#define GSFTRACKMAKER_H

typedef math::XYZPoint Point;

//
// class declaration
//

class GSFTrackMaker : public edm::EDProducer {
public:
  explicit GSFTrackMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag gsftracksInputTag_;
  edm::InputTag beamSpot_tag_;

  
};

#endif
