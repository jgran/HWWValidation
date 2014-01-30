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

}


HWWAnalyzer::~HWWAnalyzer(){ 
}


void HWWAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  using namespace HWWFunctions;

  //reset variables
  HWWVal::Reset();

  //reset counters to 0
  monitor.counters.clear();

  //count total events
  monitor.count(MM, "total events", 1.0);
  monitor.count(EE, "total events", 1.0);
  monitor.count(EM, "total events", 1.0);
  monitor.count(ME, "total events", 1.0);

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

  FillHistograms();

}//end analyze


void HWWAnalyzer::bookHistograms(DQMStore::IBooker & ibooker,edm::Run const &, edm::EventSetup const &){

  ibooker.setCurrentFolder("Physics_HWW");

  cutflow_mm = ibooker.book1D("cutflow_mm", "Cutflow mm", 20, 0, 20);	
  cutflow_ee = ibooker.book1D("cutflow_ee", "Cutflow ee", 20, 0, 20);	
  cutflow_em = ibooker.book1D("cutflow_em", "Cutflow em", 20, 0, 20);	
  cutflow_me = ibooker.book1D("cutflow_me", "Cutflow me", 20, 0, 20);	
  
}

void HWWAnalyzer::FillHistograms(){

  MonitorElement* hist[4];
  hist[0] = cutflow_mm;
  hist[1] = cutflow_ee;
  hist[2] = cutflow_em;
  hist[3] = cutflow_me;

  for (unsigned int i=0; i<4; i++){
    for (unsigned int j=0; j<monitor.counters.size(); ++j){
      int content = hist[i]->getBinContent(j+1) + monitor.counters[j].nevt[i];
      float error = sqrt(content);
      hist[i]->setBinContent(j+1, content);
      hist[i]->setBinError(j+1, error);
      hist[i]->setBinLabel(j+1, monitor.counters[j].name.c_str(), 1);
    }
  }

}

DEFINE_FWK_MODULE(HWWAnalyzer);
