#ifndef RECOCONVERSIONMAKER_H
#define RECOCONVERSIONMAKER_H

//
// class declaration
//

class RecoConversionMaker : public edm::EDProducer {
public:
  explicit RecoConversionMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  double  lxy(const math::XYZPoint& myBeamSpot, const reco::Conversion&) const;
      
  // ----------member data ---------------------------
  std::string   aliasprefix_;
  edm::InputTag recoConversionInputTag_;
  edm::InputTag beamSpotInputTag_;
  


};

#endif
