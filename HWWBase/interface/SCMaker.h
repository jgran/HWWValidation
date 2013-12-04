#ifndef SCMAKER_H
#define SCMAKER_H

//
// class declaration
//

class SCMaker : public edm::EDProducer {
public:
  explicit SCMaker (const edm::ParameterSet&);

private:
  virtual void beginRun(edm::Run&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // used by mc debugging
  void closestMCParticle(const HepMC::GenEvent *genEvent, const reco::SuperCluster &sc,
			 double &dRClosest, double &energyClosest);
  float ecalEta(float EtaParticle , float Zvertex, float plane_Radius);
  //

  math::XYZTLorentzVectorF initP4(const math::XYZPoint &pvPos,                 
				  const reco::SuperCluster &sc);

  // ----------member data ---------------------------

  // mc debugging
  bool debug_;
  edm::InputTag MCTruthCollection_;

  // preselection cuts
  double scEtMin_;

  // supercluster input collections
  edm::InputTag scInputTag_EE_;
  edm::InputTag scInputTag_EB_;
  std::vector<edm::InputTag> scInputTags_;
  std::vector<edm::InputTag> hitInputTags_;

  // rechit input collections
  edm::InputTag hcalRecHitsInputTag_HBHE_;
  edm::InputTag ecalRecHitsInputTag_EE_;
  edm::InputTag ecalRecHitsInputTag_EB_;

  // primary vertex collection
  edm::InputTag primaryVertexInputTag_;

  // electrons
  edm::InputTag electronsInputTag_;

  // access to geometry
  unsigned long long cachedCaloGeometryID_;
  edm::ESHandle<CaloGeometry> caloGeometry_;
  const   EcalChannelStatus *channelStatus_;

  // access the laser
  edm::ESHandle<EcalLaserDbService>      laser_;

  EcalClusterLazyTools* clusterTools_;


  std::string aliasprefix_;
};

#endif

