#include "DataFormats/Math/interface/LorentzVector.h"
#include "HWWValidation/HWWBase/interface/HWWAnalyzer.h"
#include "HWWValidation/HWWBase/interface/HWW.h"


HWWAnalyzer::HWWAnalyzer(const edm::ParameterSet& iConfig)
            :vertexMaker        (iConfig, consumesCollector()),
             trackMaker         (iConfig, consumesCollector()),
             electronMaker      (iConfig, consumesCollector()),
             muonMaker          (iConfig, consumesCollector()),
             pfJetMaker         (iConfig, consumesCollector()),
             pfCandidateMaker   (iConfig, consumesCollector()),
             pfElectronMaker    (iConfig, consumesCollector()),
             gsfTrackMaker      (iConfig, consumesCollector()),
             recoConversionMaker(iConfig, consumesCollector()),
             rhoMaker           (iConfig, consumesCollector()),
             pfMETMaker         (iConfig, consumesCollector()),
             bTagMaker          (iConfig, consumesCollector()),
             mvaJetIdMaker      (iConfig, consumesCollector())
{

  egammaMvaEleEstimator = 0;
  muonMVAEstimator = 0;

  // --------------- EGamma Id MVA  --------------------------
  std::vector<std::string> egammaweights = iConfig.getParameter<std::vector<std::string> > ("InputEGammaWeights");
  egammaMvaEleEstimator = new EGammaMvaEleEstimator();
  egammaMvaEleEstimator->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, egammaweights );

  // --------------- Muon RingIso MVA  --------------------------
  std::vector<std::string> muonisoweights = iConfig.getParameter<std::vector<std::string> > ("InputMuonIsoWeights");
  muonMVAEstimator = new MuonMVAEstimator();
  muonMVAEstimator->initialize( "MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muonisoweights );

  //DQMStore
  //dbe_ = edm::Service<DQMStore>().operator->();

}


HWWAnalyzer::~HWWAnalyzer(){ 
}

// --- method called for each event  ------------
void HWWAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  using namespace HWWFunctions;

  HWWVal::Reset();

  //get variables
  eventMaker    .SetVars(iEvent, iSetup);
  vertexMaker   .SetVars(iEvent, iSetup);
  trackMaker    .SetVars(iEvent, iSetup);
  electronMaker .SetVars(iEvent, iSetup);
  muonMaker     .SetVars(iEvent, iSetup);
  pfJetMaker    .SetVars(iEvent, iSetup);
  hypDilepMaker .SetVars(iEvent, iSetup);

  //check some basic event requirements
  std::vector<int> goodHyps;
  for(unsigned int i=0; i < HWWVal::hyp_p4().size(); i++){
    if(!passFirstCuts(i)) continue;
    goodHyps.push_back(i);
  }
  
  //no need to continue if event didn't pass basic requirements
  if(goodHyps.size() > 0){

    //get variables
    pfCandidateMaker    .SetVars(iEvent, iSetup);
    pfElectronMaker     .SetVars(iEvent, iSetup);
    pfElToElAssMaker    .SetVars(iEvent, iSetup);
    gsfTrackMaker       .SetVars(iEvent, iSetup);
    recoConversionMaker .SetVars(iEvent, iSetup);
    rhoMaker            .SetVars(iEvent, iSetup);

    //to hold indices of candidate lepton pairs
    std::vector<int> candidates;

    //get lepton pairs that pass baseline selection
    for(unsigned int i=0; i < goodHyps.size(); i++){
      if(!passBaseline(goodHyps.at(i), egammaMvaEleEstimator, muonMVAEstimator)) continue;
      candidates.push_back(i);      
    }

    if(candidates.size()>0){

      //get variables
      pfMETMaker            .SetVars(iEvent, iSetup);
      trkMETMaker           .SetVars(iEvent, iSetup);
      bTagMaker             .SetVars(iEvent, iSetup);
      mvaJetIdMaker         .SetVars(iEvent, iSetup);

      //find best lepton pair
      int bestHyp = bestHypothesis(candidates);

      //perform remaining selections
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

//monitor.print();
monitor.makeHistograms();

}

DEFINE_FWK_MODULE(HWWAnalyzer);
