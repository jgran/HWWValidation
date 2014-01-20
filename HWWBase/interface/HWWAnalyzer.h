#ifndef HWWANALYZER_H
#define HWWANALYZER_H  

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "HWWValidation/HWWBase/interface/analysisSelections.h"
#include "HWWValidation/HWWBase/interface/EGammaMvaEleEstimator.h"
#include "HWWValidation/HWWBase/interface/MuonMVAEstimator.h"
#include "HWWValidation/HWWBase/interface/monitor.h"

#include "HWWValidation/HWWBase/interface/VertexMaker.h"

//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"


//
// class declaration
//

class HWWAnalyzer : public edm::EDAnalyzer {
//class HWWAnalyzer : public DQMEDAnalyzer {
   public:
      explicit HWWAnalyzer(const edm::ParameterSet&);
      ~HWWAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      EGammaMvaEleEstimator* egammaMvaEleEstimator;
      MuonMVAEstimator* muonMVAEstimator;
  
      EventMonitor::hypo_monitor monitor;

      VertexMaker         vertexMaker;
      EventMaker          eventMaker;
      TrackMaker          trackMaker;
      ElectronMaker       electronMaker;
      MuonMaker           muonMaker;
      PFJetMaker          pfJetMaker;
      HypDilepMaker       hypDilepMaker;
      PFCandidateMaker    pfCandidateMaker;
      PFElectronMaker     pfElectronMaker;
      PFElToElAssMaker    pfElToElAssMaker;
      GSFTrackMaker       gsfTrackMaker;
      RecoConversionMaker recoConversionMaker;
      RhoMaker            rhoMaker;
      PFMETMaker          pfMETMaker;
      TrkMETMaker         trkMETMaker;
      BTagMaker           bTagMaker;
      MVAJetIdMaker       mvaJetIdMaker;

};

#endif
