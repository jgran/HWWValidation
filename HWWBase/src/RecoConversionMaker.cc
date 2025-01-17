#include "HWWValidation/HWWBase/interface/RecoConversionMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

using namespace std;
using namespace reco;
using namespace edm;
typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPointF Point;

RecoConversionMaker::RecoConversionMaker(const edm::ParameterSet& iConfig, edm::ConsumesCollector iCollector){

  Conversion_  = iCollector.consumes<edm::View<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("recoConversionInputTag"));
  BeamSpot_    = iCollector.consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));

}

double lxy(const math::XYZPoint& myBeamSpot, const Conversion& conv) {

  const reco::Vertex &vtx = conv.conversionVertex();
  if (!vtx.isValid()) return -9999.;

  math::XYZVectorF mom = conv.refittedPairMomentum();
  
  double dbsx = vtx.x() - myBeamSpot.x();
  double dbsy = vtx.y() - myBeamSpot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
  return lxy;  
  
}

void RecoConversionMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_convs_isConverted();
  HWWVal::Load_convs_quality();
  HWWVal::Load_convs_tkidx();
  HWWVal::Load_convs_tkalgo();
  HWWVal::Load_convs_nHitsBeforeVtx();
  HWWVal::Load_convs_ndof();
  HWWVal::Load_convs_chi2();
  HWWVal::Load_convs_dl();

  // get reco Conversions
  Handle<View<Conversion> > convs_h;
  iEvent.getByToken(Conversion_, convs_h);
  
  Handle<BeamSpot> beamSpotH;
  iEvent.getByToken(BeamSpot_, beamSpotH);

  for(View<Conversion>::const_iterator it = convs_h->begin();
      it != convs_h->end(); it++) { 

    HWWVal::convs_isConverted().push_back(it->isConverted());
    //quality 
    int qualityMask = 0;
    for(int iM = 0; iM < 32; ++iM) {
      if(it->quality((Conversion::ConversionQuality)iM)) qualityMask |= 1 << iM;
    }

    HWWVal::convs_quality().push_back(qualityMask);
    
    vector<edm::RefToBase<reco::Track> > v_temp_trks = it->tracks();
    vector<int> v_temp_out;
    vector<int> v_temp_outalgo;
    for(unsigned int i = 0; i < v_temp_trks.size(); i++) {
      v_temp_out.push_back(v_temp_trks.at(i).key());
      v_temp_outalgo.push_back(v_temp_trks.at(i)->algo());
    }
			       

    HWWVal::convs_tkidx().push_back(v_temp_out);
    HWWVal::convs_tkalgo().push_back(v_temp_outalgo);

    v_temp_out.clear();
    v_temp_outalgo.clear();
    vector<uint8_t> v_temp_nhits = it->nHitsBeforeVtx();

    for(unsigned int i = 0; i < v_temp_nhits.size(); i++) v_temp_out.push_back(v_temp_nhits.at(i));

    HWWVal::convs_nHitsBeforeVtx().push_back(v_temp_out);
    HWWVal::convs_ndof().push_back(it->conversionVertex().ndof() );
    HWWVal::convs_chi2().push_back(it->conversionVertex().chi2() );
    HWWVal::convs_dl().push_back(lxy(beamSpotH->position(), *it));

  }//reco conversion loop
}

