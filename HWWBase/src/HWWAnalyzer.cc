#include "DataFormats/Math/interface/LorentzVector.h"
#include "HWWValidation/HWWBase/interface/HWWAnalyzer.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace edm;
using namespace std;

HWWAnalyzer::HWWAnalyzer(const edm::ParameterSet& iConfig){

  eventInputTag         = iConfig.getParameter<InputTag>("eventInputTag");
  vertexInputTag        = iConfig.getParameter<InputTag>("vertexInputTag");
  trackInputTag         = iConfig.getParameter<InputTag>("trackInputTag");
  electronInputTag      = iConfig.getParameter<InputTag>("electronInputTag");
  muonInputTag          = iConfig.getParameter<InputTag>("muonInputTag");
  hypInputTag           = iConfig.getParameter<InputTag>("hypInputTag");
  conversionInputTag    = iConfig.getParameter<InputTag>("conversionInputTag");
  scInputTag            = iConfig.getParameter<InputTag>("scInputTag");
  pfCandidateInputTag   = iConfig.getParameter<InputTag>("pfCandidateInputTag");
  gsfTrackInputTag      = iConfig.getParameter<InputTag>("gsfTrackInputTag");
  kt6PFInputTag         = iConfig.getParameter<InputTag>("kt6PFInputTag");
  elToMuAssInputTag     = iConfig.getParameter<InputTag>("elToMuAssInputTag");
  pfElToElAssInputTag   = iConfig.getParameter<InputTag>("pfElToElAssInputTag");
  fastJetInputTag       = iConfig.getParameter<InputTag>("fastJetInputTag");
  wwRhoDefaultInputTag  = iConfig.getParameter<InputTag>("wwRhoDefaultInputTag");
  wwRhoVoronoiInputTag  = iConfig.getParameter<InputTag>("wwRhoVoronoiInputTag");
  pfmetInputTag         = iConfig.getParameter<InputTag>("pfmetInputTag");
  trkMetInputTag        = iConfig.getParameter<InputTag>("trkMetInputTag");
  pfJetInputTag         = iConfig.getParameter<InputTag>("pfJetInputTag");
  bTagPFJetInputTag     = iConfig.getParameter<InputTag>("bTagPFJetInputTag");
  mvaJetIdInputTag      = iConfig.getParameter<InputTag>("mvaJetIdInputTag");

  egammaMvaEleEstimator = 0;
  muonMVAEstimator = 0;

  // --------------- EGamma Id MVA  --------------------------
  vector<std::string> egammaweights;
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat1.weights.xml"); 
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat2.weights.xml"); 
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat3.weights.xml"); 
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat4.weights.xml"); 
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat5.weights.xml"); 
  egammaweights.push_back("../files/Electrons_BDTG_TrigV0_Cat6.weights.xml"); 
  egammaMvaEleEstimator = new EGammaMvaEleEstimator();
  egammaMvaEleEstimator->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, egammaweights );

  // --------------- Muon RingIso MVA  --------------------------
  vector<std::string> muonisoweights;
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muonisoweights.push_back("../files/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
  muonMVAEstimator = new MuonMVAEstimator();
  muonMVAEstimator->initialize( "MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muonisoweights );

/*
  edm::Service<TFileService> fs;
  cutflow_mm = fs->make<TH1F>("cutflow_mm" , "cutflow_mm" , 20 , 0 , 20);
  cutflow_ee = fs->make<TH1F>("cutflow_ee" , "cutflow_ee" , 20 , 0 , 20);
  cutflow_em = fs->make<TH1F>("cutflow_em" , "cutflow_em" , 20 , 0 , 20);
  cutflow_me = fs->make<TH1F>("cutflow_me" , "cutflow_me" , 20 , 0 , 20);
*/

}


HWWAnalyzer::~HWWAnalyzer(){ 
}

// ------------ method called for each event  ------------
void HWWAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //get variables needed for passFirstCuts
  HWWAnalyzer::GetFirstVariables(iEvent, iSetup);

  //check some basic event requirements
  std::vector<int> goodHyps;
  for(unsigned int i=0; i < hww.hyp_p4().size(); i++){
    if(!passFirstCuts(i)) continue;
    goodHyps.push_back(i);
  }
  
  //no need to continue if event didn't pass basic requirements
  if(goodHyps.size() > 0){

    //get variables needed for remaining selections
    HWWAnalyzer::GetVariables(iEvent, iSetup);

    //to hold indices of candidate dilepton pairs
    vector<int> candidates;

    //get lepton pairs that pass baseline selection
    for(unsigned int i=0; i < goodHyps.size(); i++){
      if(!passBaseline(goodHyps.at(i), egammaMvaEleEstimator, muonMVAEstimator)) continue;
      candidates.push_back(i);      
    }

    if(candidates.size()>0){

      //find best lepton pair
      int bestHyp = bestHypothesis(candidates);

      //create cutflow table
      doCutFlow(bestHyp, monitor, egammaMvaEleEstimator, muonMVAEstimator);

    }

  }  
}//end analyze


// ------------ method called once each job just before starting event loop  ------------
void 
HWWAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HWWAnalyzer::endJob() 
{

monitor.print();
monitor.makeHistograms();

}

// ------------ method called when starting to processes a run  ------------
void 
HWWAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HWWAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HWWAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HWWAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HWWAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HWWAnalyzer);
