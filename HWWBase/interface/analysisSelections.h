#ifndef ANALYSIS_SELECTIONS_H
#define ANALYSIS_SELECTIONS_H

#include "HWWValidation/HWWBase/interface/EGammaMvaEleEstimator.h"
#include "HWWValidation/HWWBase/interface/monitor.h"
#include "HWWValidation/HWWBase/interface/ElectronIDMVA.h"
#include "HWWValidation/HWWBase/interface/MuonIDMVA.h"
#include "HWWValidation/HWWBase/interface/MuonMVAEstimator.h"
#include "HWWValidation/HWWBase/interface/wwtypes.h"
#include "HWWValidation/HWWBase/interface/analysisEnums.h"

namespace HWWFunctions {

  void doCutFlow(int, EventMonitor::hypo_monitor&, EGammaMvaEleEstimator*, MuonMVAEstimator*);
  bool passFirstCuts(int);
  bool passCharge(int);
  bool passBaseline(int, EGammaMvaEleEstimator*, MuonMVAEstimator*);
  bool passFullLep(int, EGammaMvaEleEstimator*, MuonMVAEstimator*);
  bool passExtraLeptonVeto(int, EGammaMvaEleEstimator*, MuonMVAEstimator*);
  bool passZVeto(int);
  bool passMinMet(int);
  bool passMinMet40(int);
  bool passDPhiDiLepJet(int);
  bool passSoftMuonVeto(int);
  bool passTopVeto(int);

  int  bestHypothesis(const std::vector<int>&);
  bool isGoodVertex(int);
  unsigned int nGoodVertex();


  std::vector<JetPair> getJets(int, double, double, bool, bool);
  std::vector<JetPair> getDefaultJets(unsigned int, bool);
  unsigned int numberOfJets(unsigned int);
  bool defaultBTag(unsigned int, float);
  Bool_t comparePt(JetPair lv1, JetPair lv2);
  bool passMVAJetId(double, double, double, unsigned int);

  //
  // Electron Id
  //

  bool   ww_elBase(unsigned int i);
  bool   ww_elId(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator);
  bool   ww_eld0(unsigned int i);
  bool   ww_eld0PV(unsigned int i);
  bool   ww_eldZPV(unsigned int i);
  bool   ww_elIso(unsigned int i);
  double ww_elIsoVal(unsigned int i);

  // combined analysis selectors
  bool goodElectronTMVA(EGammaMvaEleEstimator* egammaMvaEleEstimator, int useMVAeleId, unsigned int i); 
  bool goodElectronWithoutIsolation(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator);
  bool goodElectronIsolated(unsigned int i, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool lockToCoreSelectors);
  bool fakableElectron(unsigned int i,EleFOTypes);
  bool ElectronFOV4(unsigned int i);

  //
  // Muon Id
  //

  bool   ww_muBase(unsigned int i);
  bool   ww_muId(unsigned int i, bool useMVAmuId, MuonIDMVA *mva);
  bool   ww_muIso(unsigned int i);
  bool   ww_muIso(unsigned int i, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  double ww_muIsoVal(unsigned int i);
  bool   ww_mud0(unsigned int i);
  bool   ww_mud0PV(unsigned int i);
  bool   ww_mudZPV(unsigned int i, float cut=0.1);

  // combined analysis selectors
  bool goodMuonTMVA(MuonIDMVA *mva, unsigned int i);
  bool goodMuonWithoutIsolation(unsigned int i, bool useMVAmuId, MuonIDMVA *mva,
                  MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  bool goodMuonIsolated( unsigned int i, bool lockToCoreSelectors, bool useMVAmuId, MuonIDMVA *mva, 
               MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle );
  bool fakableMuon(unsigned int i, MuFOTypes, MuonMVAEstimator* muonMVAEstimator,
                  std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  bool passMuonRingsMVA(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  bool passMuonRingsMVAFO(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);


  //
  // leptons
  //
  std::vector<LeptonPair> getExtraLeptons(int i_hyp, double minPt,  bool useLHeleId, int useMVAeleId, 
                      EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
                      MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  unsigned int numberOfExtraLeptons(int i_hyp, double minPt, bool useLHeleId, int useMVAeleId, 
                    EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
                    MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle);
  unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated);

  double nearestDeltaPhi(double Phi, int i_hyp);
  double projectedMet(unsigned int i_hyp, double met, double phi);


  bool toptag(int i_hyp, double minPt, std::vector<JetPair> ignoreJets=std::vector<JetPair>());

}
#endif
