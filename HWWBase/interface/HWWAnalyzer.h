#ifndef HWWANALYZER_H
#define HWWANALYZER_H  

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "HWWValidation/HWWBase/interface/analysisSelections.h"
#include "HWWValidation/HWWBase/interface/EGammaMvaEleEstimator.h"
#include "HWWValidation/HWWBase/interface/MuonMVAEstimator.h"
#include "HWWValidation/HWWBase/interface/monitor.h"


//typedef math::XYZTLorentzVectorF LorentzVector;

//
// class declaration
//

class HWWAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HWWAnalyzer(const edm::ParameterSet&);
      ~HWWAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void GetFirstVariables(const edm::Event&, const edm::EventSetup&);
      virtual void GetVariables(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


      // ----------member data ---------------------------

      edm::InputTag eventInputTag;
      edm::InputTag vertexInputTag;
      edm::InputTag trackInputTag;
      edm::InputTag electronInputTag;
      edm::InputTag muonInputTag;
      edm::InputTag hypInputTag;
      edm::InputTag conversionInputTag;
      edm::InputTag scInputTag;
      edm::InputTag pfCandidateInputTag;
      edm::InputTag gsfTrackInputTag;
      edm::InputTag kt6PFInputTag;
      edm::InputTag elToMuAssInputTag;
      edm::InputTag pfElToElAssInputTag;
      edm::InputTag fastJetInputTag;
      edm::InputTag wwRhoDefaultInputTag;
      edm::InputTag wwRhoVoronoiInputTag;
      edm::InputTag pfmetInputTag;
      edm::InputTag trkMetInputTag;
      edm::InputTag pfJetInputTag;
      edm::InputTag bTagPFJetInputTag;
      edm::InputTag mvaJetIdInputTag;

      int baselineCounter;
      int chargeCounter;
      int FullLepCounter;
      int ExtraLepCounter;
      int metCounter;
      int dilepMassCounter;
      int zvetoCounter;
      int minmetCounter;
      int minmet40Counter;
      int dphiCounter;
      int softMuonCounter;
      int topVetoCounter;
      int dilepPtCounter;
  
      EGammaMvaEleEstimator* egammaMvaEleEstimator;
      MuonMVAEstimator* muonMVAEstimator;
  
      hypo_monitor monitor;

/*
      TH1F *cutflow_mm;
      TH1F *cutflow_ee;
      TH1F *cutflow_em;
      TH1F *cutflow_me;
*/

};

#endif
