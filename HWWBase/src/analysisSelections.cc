#include "HWWValidation/HWWBase/interface/muonSelections.h"
#include "HWWValidation/HWWBase/interface/electronSelections.h"
#include "HWWValidation/HWWBase/interface/eventSelections.h"
#include "HWWValidation/HWWBase/interface/pfjetMVAtools.h"
#include "HWWValidation/HWWBase/interface/analysisSelections.h"

void doCutFlow(int i_hyp, hypo_monitor& monitor, EGammaMvaEleEstimator* egammaMvaEleEstimator, MuonMVAEstimator* muonMVAEstimator){

  HypothesisType type = getHypothesisTypeNew(i_hyp);
  double weight = 1.0;
  
  monitor.count(hww, type, "baseline", weight);
  if(!passCharge(i_hyp)) return;
  monitor.count(hww, type, "opposite sign", weight);
  if(!passFullLep(i_hyp, egammaMvaEleEstimator, muonMVAEstimator)) return;
  monitor.count(hww, type, "full lepton selection", weight);
  if(!passExtraLeptonVeto(i_hyp, egammaMvaEleEstimator, muonMVAEstimator)) return;
  monitor.count(hww, type, "extra lepton veto", weight);
  if(!(hww.evt_pfmet() > 20.0)) return;
  monitor.count(hww, type, "met > 20 GeV", weight);
  if(!(hww.hyp_p4().at(i_hyp).mass() > 12.0)) return;
  monitor.count(hww, type, "mll > 12 GeV", weight);
  if(!passZVeto(i_hyp)) return;
  monitor.count(hww, type, "|mll - mZ| > 15 GeV", weight);
  if(!passMinMet(i_hyp)) return;
  monitor.count(hww, type, "minMET > 20 GeV", weight);
  if(!passMinMet40(i_hyp)) return;
  monitor.count(hww, type, "minMET > 40 GeV for ee/mm", weight);
  if(!passDPhiDiLepJet(i_hyp)) return;
  monitor.count(hww, type, "dPhiDiLepJet < 165 dg for ee/mm", weight);
  if(!passSoftMuonVeto(i_hyp)) return;
  monitor.count(hww, type, "SoftMuons==0", weight);
  if(!passTopVeto(i_hyp)) return;
  monitor.count(hww, type, "top veto", weight);
  if(hww.hyp_p4().at(i_hyp).pt() <= 45.0) return;
  monitor.count(hww, type, "ptll > 45 GeV", weight);

  int njets = numberOfJets(i_hyp);
  std::vector<JetPair> sortedJets = getJets(jetType(), i_hyp, 30.0, 4.7, true, false);

  LorentzVector jet1;
  LorentzVector jet2;
  LorentzVector jet3;
  

  if(njets==0){
    monitor.count(hww, type, "njets == 0", weight);
    if(max(hww.hyp_ll_p4().at(i_hyp).pt(), hww.hyp_lt_p4().at(i_hyp).pt()) < 30) return; 
    monitor.count(hww, type, "max(lep1.pt(),lep2.pt())>30", weight);
    if(min(hww.hyp_ll_p4().at(i_hyp).pt(), hww.hyp_lt_p4().at(i_hyp).pt()) < 25) return; 
    monitor.count(hww, type, "min(lep1.pt(),lep2.pt())>25", weight);
  }
  if(njets==1){
    monitor.count(hww, type, "njets == 1", weight);
  }
  if ( njets==2 || njets==3 ) {
    monitor.count(hww, type, "njets == 2 or 3", weight);
    if (fabs(sortedJets[0].first.eta())>=4.7 || fabs(sortedJets[1].first.eta())>=4.7) return;
    monitor.count(hww,type,"abs(jet1.eta())<4.7 && abs(jet2.eta())<4.7",weight);
    if (njets==3 && sortedJets[2].first.pt()>30 && ((sortedJets[0].first.eta()-sortedJets[2].first.eta() > 0 && sortedJets[1].first.eta()-sortedJets[2].first.eta() < 0) ||
          (sortedJets[1].first.eta()-sortedJets[2].first.eta() > 0 && sortedJets[0].first.eta()-sortedJets[2].first.eta() < 0)) ) return;
    monitor.count(hww, type, "no central jets", weight);
  }

  
  return;

}

int bestHypothesis(const std::vector<int>& candidates){
  int best = -1;
  for( unsigned int i = 0; i < candidates.size(); ++i ) {
    unsigned int i_hyp = candidates.at(i);
    if (best<0){
      best = i_hyp;
      continue;
    }
    if ( std::max(hww.hyp_lt_p4().at(i_hyp).pt(), hww.hyp_ll_p4().at(i_hyp).pt()) >= //GC add = in case the lepton is the same 
	 std::max(hww.hyp_lt_p4().at(best).pt(),  hww.hyp_ll_p4().at(best).pt()) &&
	 std::min(hww.hyp_lt_p4().at(i_hyp).pt(), hww.hyp_ll_p4().at(i_hyp).pt()) >= 
	 std::min(hww.hyp_lt_p4().at(best).pt(),  hww.hyp_ll_p4().at(best).pt()) )
      best = i_hyp;
  }
  return best;
}

bool passFirstCuts(int i_hyp){

  if ( std::min(hww.hyp_lt_p4().at(i_hyp).pt(),hww.hyp_ll_p4().at(i_hyp).pt())<10 ) return false;
  if ( std::max(hww.hyp_lt_p4().at(i_hyp).pt(),hww.hyp_ll_p4().at(i_hyp).pt())<20 ) return false;
  if (hww.trks_d0().size()==0) return false;
  if (!isGoodVertex(0)) return false;
  if (nGoodVertex()<1) return false;
  if (hww.hyp_p4().at(i_hyp).mass2()<0) return false;

  return true;
}

bool passCharge(int i_hyp){

  if (hww.hyp_lt_id().at(i_hyp)*hww.hyp_ll_id().at(i_hyp)>0) return false;
  return true;
}

bool passBaseline(int i_hyp, EGammaMvaEleEstimator* egammaMvaEleEstimator, MuonMVAEstimator* muonMVAEstimator){

  if (abs(hww.hyp_lt_id().at(i_hyp)) == 13 && !ww_muBase(hww.hyp_lt_index().at(i_hyp)) ) return false;
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 13 && !ww_muBase(hww.hyp_ll_index().at(i_hyp)) ) return false;
  if (abs(hww.hyp_lt_id().at(i_hyp)) == 11 && !ww_elBase(hww.hyp_lt_index().at(i_hyp)) ) return false;
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 11 && !ww_elBase(hww.hyp_ll_index().at(i_hyp)) ) return false;

  bool lockToCoreSelectors = false;
  bool useLHeleId = false;
  int useMVAeleId = 1;//zero means off, otherwise it's the mva version
  bool useMVAmuId = false;

  MuonIDMVA* muonIdMVA = 0;
  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons  

  bool PASSED_LT_FINAL    = false;
  bool PASSED_LT_FO_MU2   = false;
  bool PASSED_LT_FO_ELEV4 = false; 
  bool PASSED_LL_FINAL    = false;
  bool PASSED_LL_FO_MU2   = false;
  bool PASSED_LL_FO_ELEV4 = false; 

  if (abs(hww.hyp_lt_id().at(i_hyp)) == 13){
    unsigned int index = hww.hyp_lt_index().at(i_hyp);
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) PASSED_LT_FINAL = true;
    if ( fakableMuon(index,MuFOV2, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) PASSED_LT_FO_MU2 = true;
  }
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 13){
    unsigned int index = hww.hyp_ll_index().at(i_hyp);
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) PASSED_LL_FINAL = true;
    if ( fakableMuon(index,MuFOV2, muonMVAEstimator,  nullMu, nullEle) && !goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle)) PASSED_LL_FO_MU2 = true;
  }
  if (abs(hww.hyp_lt_id().at(i_hyp)) == 11){
    unsigned int index = hww.hyp_lt_index().at(i_hyp);
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) PASSED_LT_FINAL = true;
    if ( fakableElectron(index,EleFOV4) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors)) PASSED_LT_FO_ELEV4 = true;
  }
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 11){
    unsigned int index = hww.hyp_ll_index().at(i_hyp);
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) PASSED_LL_FINAL = true;
    if ( fakableElectron(index,EleFOV4) && !goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors)) PASSED_LL_FO_ELEV4 = true;
  }
  
  if( PASSED_LT_FINAL && PASSED_LL_FINAL    ) return true;
  if( PASSED_LT_FINAL && PASSED_LL_FO_MU2   ) return true;
  if( PASSED_LT_FINAL && PASSED_LL_FO_ELEV4 ) return true;
  if( PASSED_LL_FINAL && PASSED_LT_FINAL    ) return true;
  if( PASSED_LL_FINAL && PASSED_LT_FO_MU2   ) return true;
  if( PASSED_LL_FINAL && PASSED_LT_FO_ELEV4 ) return true;

  return false;
}

bool passFullLep(int i_hyp, EGammaMvaEleEstimator* egammaMvaEleEstimator, MuonMVAEstimator* muonMVAEstimator){

  bool lockToCoreSelectors = false;
  bool useLHeleId = false;
  int useMVAeleId = 1;//zero means off, otherwise it's the mva version
  bool useMVAmuId = false;

  MuonIDMVA* muonIdMVA = 0;
  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons  

  bool PASSED_LT_FINAL    = false;
  bool PASSED_LL_FINAL    = false;

  if (abs(hww.hyp_lt_id().at(i_hyp)) == 13){
    unsigned int index = hww.hyp_lt_index().at(i_hyp);
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) PASSED_LT_FINAL = true;
  }
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 13){
    unsigned int index = hww.hyp_ll_index().at(i_hyp);
    if ( goodMuonIsolated(index, lockToCoreSelectors, useMVAmuId, muonIdMVA, muonMVAEstimator,  nullMu, nullEle) ) PASSED_LL_FINAL = true;
  }
  if (abs(hww.hyp_lt_id().at(i_hyp)) == 11){
    unsigned int index = hww.hyp_lt_index().at(i_hyp);
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) PASSED_LT_FINAL = true;
  }
  if (abs(hww.hyp_ll_id().at(i_hyp)) == 11){
    unsigned int index = hww.hyp_ll_index().at(i_hyp);
    if ( goodElectronIsolated(index, useLHeleId, useMVAeleId, egammaMvaEleEstimator, lockToCoreSelectors) ) PASSED_LL_FINAL = true;
  }
  
  if(PASSED_LT_FINAL && PASSED_LL_FINAL) return true;

  return false;
}


bool passExtraLeptonVeto(int i_hyp, EGammaMvaEleEstimator* egammaMvaEleEstimator, MuonMVAEstimator* muonMVAEstimator){

  bool useLHeleId = false;
  int useMVAeleId = 1;//zero means off, otherwise it's the mva version
  bool useMVAmuId = false;

  MuonIDMVA* muonIdMVA = 0;
  std::vector<Int_t> nullMu; // null identified muons 
  std::vector<Int_t> nullEle; // null identified electrons  

  return ( numberOfExtraLeptons(i_hyp,10, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, muonIdMVA,muonMVAEstimator,  nullMu, nullEle) == 0);

}

bool passZVeto(int i_hyp){

  if(hww.hyp_type().at(i_hyp) == 1 || hww.hyp_type().at(i_hyp) == 2) return true; //em or me dileptons have no cut on the z mass
  if(fabs(hww.hyp_p4().at(i_hyp).mass() - 91.1876) > 15.0) return true;
  return false; 
  
}


bool passMinMet(int i_hyp){

  double pmet      = projectedMet(i_hyp, hww.evt_pfmet(), hww.evt_pfmetPhi());
  double pTrackMet = projectedMet(i_hyp, hww.trk_met().at(i_hyp),   hww.trk_metPhi().at(i_hyp)  );
  return(min(pmet,pTrackMet)>20);

}

bool passMinMet40(int i_hyp){

  //require minMet > 40 GeV for ee and mm hypotheses

  if(hww.hyp_type().at(i_hyp) == 1) return true;
  if(hww.hyp_type().at(i_hyp) == 2) return true;
  double pmet      = projectedMet(i_hyp, hww.evt_pfmet(), hww.evt_pfmetPhi());
  double pTrackMet = projectedMet(i_hyp, hww.trk_met().at(i_hyp),   hww.trk_metPhi().at(i_hyp)  );
  return(min(pmet,pTrackMet)>40);

}


bool passDPhiDiLepJet(int i_hyp){

  //pass if em or me hypothesis
  if(hww.hyp_type().at(i_hyp)==1) return true;
  if(hww.hyp_type().at(i_hyp)==2) return true;

  int njets = numberOfJets(i_hyp);
  std::vector<JetPair> sortedJets = getJets(jetType(), i_hyp, 15.0, 4.7, true, false);

  if (sortedJets.size()>0) {
      if (njets<2 && fabs(ROOT::Math::VectorUtil::DeltaPhi(sortedJets[0].first,hww.hyp_p4().at(i_hyp))) >= (165.*TMath::Pi()/180.)) return false;
      else if ( (njets>= 2) && (fabs(ROOT::Math::VectorUtil::DeltaPhi( (sortedJets[0].first+sortedJets[1].first),hww.hyp_p4().at(i_hyp))) >= (165.*TMath::Pi()/180.)) ) return false;
      else return true;
  }
  return true;

}


bool passSoftMuonVeto(int i_hyp){

  return(numberOfSoftMuons(i_hyp,true)==0);

}


bool passTopVeto(int i_hyp){

  return(!toptag(jetType(), i_hyp, 10));

}



WWJetType jetType(){
    return pfJet;
}



std::vector<JetPair> getJets(WWJetType type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag){

    std::vector<JetPair> jets;
    const double vetoCone = 0.3;
    // bug fix for mva jet id
    vector <float> fixedpfjetmva_analobj; getGoodMVAs(fixedpfjetmva_analobj, "mvavalue"); 

    switch ( type ){

        case pfJet:
            for ( unsigned int i=0; i < hww.pfjets_p4().size(); ++i) {
                double jec = hww.pfjets_JEC().at(i);;

                if ( (hww.pfjets_p4().at(i).pt() * jec) < etThreshold ) continue;
                if ( btag && !defaultBTag(type,i, jec) ) continue;
                if ( TMath::Abs(hww.pfjets_p4().at(i).eta()) > etaMax ) continue;
                if ( (hww.hyp_lt_p4().at(i_hyp).pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_lt_p4().at(i_hyp), hww.pfjets_p4().at(i))) < vetoCone) ||
                     (hww.hyp_ll_p4().at(i_hyp).pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_ll_p4().at(i_hyp), hww.pfjets_p4().at(i))) < vetoCone) ) continue;
                if ( !passMVAJetId( hww.pfjets_p4().at(i).pt() * jec, hww.pfjets_p4().at(i).eta(), fixedpfjetmva_analobj[i], 2) ) continue;


                jets.push_back(JetPair(hww.pfjets_p4().at(i) * jec,i));
            }
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    }
    if ( sortJets ) std::sort(jets.begin(), jets.end(), comparePt);
    return jets;
}



std::vector<JetPair> getDefaultJets(unsigned int i_hyp, bool btagged){
    return getJets(jetType(), i_hyp, 30, 4.7, false, btagged); // V1 // new cut
}



unsigned int numberOfJets(unsigned int i_hyp){
    return getDefaultJets(i_hyp, false).size();
}



bool defaultBTag(WWJetType type, unsigned int iJet, float jec) {
    
	switch ( type ) {
        case pfJet:
            if ( hww.pfjets_trackCountingHighEffBJetTag().at(iJet) > 2.1) return true;
            break;
        default:
            std::cout << "ERROR: Please use PFJets : " << type << " FixIt!" << std::endl;
            assert(0);
    }
    return 0;
}

Bool_t comparePt(JetPair lv1, JetPair lv2) {
    return lv1.first.pt() > lv2.first.pt();
}



// tightness : 2=loose 1=medium 0=tight
bool passMVAJetId(double corjetpt, double jeteta, double mvavalue, unsigned int tightness)         
{
	if(tightness>2)
	{
		cout << "ERROR : tightness should be 0, 1, or 2. " << endl;
		return false;
	}

	double fMVACut[3][4][4];

	// These cuts are for 42X but used for 52X jet Id
	//Tight Id
	fMVACut[0][0][0] =  0.5; fMVACut[0][0][1] = 0.6; fMVACut[0][0][2] = 0.6; fMVACut[0][0][3] = 0.9;
	fMVACut[0][1][0] = -0.2; fMVACut[0][1][1] = 0.2; fMVACut[0][1][2] = 0.2; fMVACut[0][1][3] = 0.6;
	fMVACut[0][2][0] =  0.3; fMVACut[0][2][1] = 0.4; fMVACut[0][2][2] = 0.7; fMVACut[0][2][3] = 0.8;
	fMVACut[0][3][0] =  0.5; fMVACut[0][3][1] = 0.4; fMVACut[0][3][2] = 0.8; fMVACut[0][3][3] = 0.9;
	//Medium id
	fMVACut[1][0][0] =  0.2; fMVACut[1][0][1] = 0.4; fMVACut[1][0][2] = 0.2; fMVACut[1][0][3] = 0.6;
	fMVACut[1][1][0] = -0.3; fMVACut[1][1][1] = 0. ; fMVACut[1][1][2] = 0. ; fMVACut[1][1][3] = 0.5;
	fMVACut[1][2][0] =  0.2; fMVACut[1][2][1] = 0.2; fMVACut[1][2][2] = 0.5; fMVACut[1][2][3] = 0.7;
	fMVACut[1][3][0] =  0.3; fMVACut[1][3][1] = 0.2; fMVACut[1][3][2] = 0.7; fMVACut[1][3][3] = 0.8;
	//Loose Id 
	fMVACut[2][0][0] = -0.2; fMVACut[2][0][1] =  0. ; fMVACut[2][0][2] =  0.2; fMVACut[2][0][3] = 0.5;
	fMVACut[2][1][0] = -0.4; fMVACut[2][1][1] = -0.4; fMVACut[2][1][2] = -0.4; fMVACut[2][1][3] = 0.4;
	fMVACut[2][2][0] =  0. ; fMVACut[2][2][1] =  0. ; fMVACut[2][2][2] =  0.2; fMVACut[2][2][3] = 0.6;
	fMVACut[2][3][0] =  0. ; fMVACut[2][3][1] =  0. ; fMVACut[2][3][2] =  0.6; fMVACut[2][3][3] = 0.2;

	// pT categorization
	int ptId = 0;
	if( corjetpt > 10 && corjetpt < 20 ) ptId = 1;
	if( corjetpt > 20 && corjetpt < 30 ) ptId = 2;
	if( corjetpt > 30                  ) ptId = 3;

	// eta categorization
	int etaId = 0;
	if( fabs(jeteta) > 2.5  && fabs(jeteta) < 2.75 ) etaId = 1;
	if( fabs(jeteta) > 2.75 && fabs(jeteta) < 3.0  ) etaId = 2;
	if( fabs(jeteta) > 3.0  && fabs(jeteta) < 5.0  ) etaId = 3;

	// return  
	if( mvavalue > fMVACut[tightness][ptId][etaId] ) return true;
	return false;
}


 
unsigned int nGoodVertex() {
  unsigned int nVtx = 0;
  for ( unsigned int i = 0; i < hww.vtxs_sumpt().size(); ++i ){
    if (!isGoodVertex(i)) continue;
    nVtx++;
  }
  return nVtx;
}   

int primaryVertex(){
  return 0;
}

//
// Electron ID
//

bool goodElectronWithoutIsolation(unsigned int i,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator){
    return ww_elBase(i) && ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool goodElectronIsolated(unsigned int i,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool lockToCoreSelectors){
    bool ptcut = hww.els_p4().at(i).pt() >= 10.0;
    bool core = ptcut && pass_electronSelection( i, electronSelection_smurfV6);
    bool internal = ww_elBase(i) && ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_eldZPV(i) && ww_elIso(i);
    assert(!lockToCoreSelectors || core==internal); 
    return internal;
}

bool ElectronFOIdV4(unsigned int i) {

	float pt = hww.els_p4().at(i).pt();
	float etaSC = hww.els_etaSC().at(i);

	if (fabs(etaSC)<1.479) {
		if (hww.els_sigmaIEtaIEta().at(i)>0.01		||
			fabs(hww.els_dEtaIn().at(i))>0.007 	||
			fabs(hww.els_dPhiIn().at(i))>0.15 		||
			hww.els_hOverE().at(i)>0.12 			||
			hww.els_tkIso().at(i)/pt>0.2 			||
			(hww.els_ecalIso().at(i) - 1.0)/pt>0.2 ||
			hww.els_hcalIso().at(i)/pt>0.2 ) return false;
	} else {
		if (hww.els_sigmaIEtaIEta().at(i)>0.03		|| 
			fabs(hww.els_dEtaIn().at(i))>0.009 	||
			fabs(hww.els_dPhiIn().at(i))>0.10 		|| 
			hww.els_hOverE().at(i)>0.10 			||
			hww.els_tkIso().at(i)/pt>0.2 			||
			hww.els_ecalIso().at(i)/pt>0.2 		||
			hww.els_hcalIso().at(i)/pt>0.2 ) return false;
	}

    // MIT conversion
	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( hww.els_exp_innerlayers().at(i) > 0 ) return false;
	
	return true;
} 

bool ElectronFOV4(unsigned int i){
    return ww_elBase(i) && ElectronFOIdV4(i) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool fakableElectron(unsigned int i, EleFOTypes type){
    if ( hww.els_p4().at(i).pt() < 10.0 ) return false;
    switch (type){
        case EleFOV1: return pass_electronSelection( i, electronSelectionFO_el_smurf_v1);
        case EleFOV2: return pass_electronSelection( i, electronSelectionFO_el_smurf_v2);
        case EleFOV3: return pass_electronSelection( i, electronSelectionFO_el_smurf_v3);
        //case EleFOV4: return pass_electronSelection( i, electronSelectionFO_el_smurf_v4);
        case EleFOV4: return ElectronFOV4(i);
    }
    return false;
}

//
// Muon ID
//

bool goodMuonTMVA(MuonIDMVA* mva, unsigned int i) {
  //Find MVA Bin
  int subdet = 0;
  if (fabs(hww.mus_p4().at(i).eta()) < 1.479) subdet = 0;
  else subdet = 1;
  int ptBin = 0;
  if (hww.mus_p4().at(i).pt() > 14.5) ptBin = 1;
  if (hww.mus_p4().at(i).pt() > 20.0) ptBin = 2;

  int MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  double MVACut = -999.;
  //same signal eff as cut-based (using V10 - Detector Based Iso)
  if (MVABin == 0) MVACut = -0.5618;
  if (MVABin == 1) MVACut = -0.3002;
  if (MVABin == 2) MVACut = -0.4642;
  if (MVABin == 3) MVACut = -0.2478;
  if (MVABin == 4) MVACut = 0.1706;
  if (MVABin == 5) MVACut = 0.8146;

  double mvaValue=mva->MVAValue(i, 0);

  //Isolation
  double iso03 = 0;
  iso03 = muonIsoValuePF(i,0,0.3);

  //Explicitly Apply M2 Denominator Cuts
  bool pass = true;
  if (hww.mus_p4().at(i).pt() < 10) pass = false;
  if (fabs(hww.mus_p4().at(i).eta()) >= 2.4) pass = false;

  if (! ( (0==0)
         &&
          (
           (((hww.mus_type().at(i)) & (1<<1)) == (1<<1) 
            && hww.mus_gfit_chi2().at(i)/hww.mus_gfit_ndof().at(i) < 10.0
            && (hww.mus_gfit_validSTAHits().at(i) > 0)
            && (hww.mus_nmatches().at(i) > 1 )
            )
           || 
           ( ((hww.mus_type().at(i)) & (1<<2)) == (1<<2)   
             && hww.mus_pid_TMLastStationTight().at(i) == 1
             )
           )
          && ((hww.mus_type().at(i)) & (1<<2)) == (1<<2)
          && hww.mus_validHits().at(i) > 10
          && (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(i)) > 0)          
          //&& fabs(mu->d0) < 0.2
          //&& fabs(mu->dz) < 0.1
          && iso03 < 0.4
          && ( hww.mus_ptErr().at(i)/hww.mus_p4().at(i).pt() < 0.1)
          && hww.mus_trkKink().at(i) < 20.
          && mvaValue > MVACut
          )
      ) {
    pass = false;
  }
  return pass;
}

bool goodMuonWithoutIsolation(unsigned int i, bool useMVAmuId, MuonIDMVA *mva, 
								MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
  return ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i, useMVAmuId, mva) && passMuonRingsMVAFO(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
}

bool goodMuonIsolated(	unsigned int i, bool lockToCoreSelectors, bool useMVAmuId, MuonIDMVA *mva, 
						MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle ){
    bool ptcut = hww.mus_p4().at(i).pt() >= 10.0;
    bool core = ptcut && muonId(i, NominalSmurfV6);
    bool internal = ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i, useMVAmuId, mva) && ww_muIso(i, muonMVAEstimator, IdentifiedMu,  IdentifiedEle); 
    assert(!lockToCoreSelectors || core==internal);
    return internal;
}


///////////////////////////
/////// Electron ID ///////
///////////////////////////

bool ww_elBase(unsigned int index){
    if (hww.els_p4().at(index).pt() < 10.0) return false;
    if (fabs(hww.els_p4().at(index).eta()) > 2.5) return false;
    return true;
}


bool ww_elId(unsigned int index, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator) {

    if (useLHeleId) {
        if (hww.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,hww.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
        if (hww.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,hww.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
    }
    if (useMVAeleId>0){
      if (!goodElectronTMVA(egammaMvaEleEstimator, useMVAeleId, index)) return false;
    } else {
        if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
    }

    // MIT conversion
    if ( isFromConversionMIT(index) ) return false;
    // conversion rejection - hit based
    if ( hww.els_exp_innerlayers().at(index) > 0 ) return false;

    return true;
}

bool ww_eld0(unsigned int index){
    return fabs(hww.els_d0corr().at(index)) < 0.02;
}


bool ww_eld0PV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    double dxyPV = hww.els_d0().at(index)-
        hww.vtxs_position().at(vtxIndex).x()*sin(hww.els_trk_p4().at(index).phi())+
        hww.vtxs_position().at(vtxIndex).y()*cos(hww.els_trk_p4().at(index).phi());
    return fabs(dxyPV) < 0.02;
}

bool ww_eldZPV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    double dzpv = dzPV(hww.els_vertex_p4().at(index), hww.els_trk_p4().at(index), hww.vtxs_position().at(vtxIndex));
    return fabs(dzpv)<0.1;
}

double ww_elIsoVal(unsigned int index){
	return electronIsoValuePF2012_FastJetEffArea_HWW( index );
}

bool ww_elIso(unsigned int index){
	float pfiso = ww_elIsoVal( index ); 
	return pfiso<0.15;
}


///////////////////////////
///////// Muon ID /////////
///////////////////////////

bool ww_muBase(unsigned int index){
    if (hww.mus_p4().at(index).pt() < 10.0) return false;
    if (fabs(hww.mus_p4().at(index).eta()) > 2.4) return false;
    if (hww.mus_type().at(index) == 8) return false; // not STA
    return true;
}

bool ww_mud0(unsigned int index){
    return fabs(hww.mus_d0corr().at(index)) < 0.02;
}

double ww_mud0ValuePV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return 9999;
    double dxyPV = hww.mus_d0().at(index)-
        hww.vtxs_position().at(vtxIndex).x()*sin(hww.mus_trk_p4().at(index).phi())+
        hww.vtxs_position().at(vtxIndex).y()*cos(hww.mus_trk_p4().at(index).phi());
    return fabs(dxyPV);
}

bool ww_mud0PV(unsigned int index){
    if ( hww.mus_p4().at(index).pt() < 20. ) return ww_mud0ValuePV(index) < 0.01;
    return ww_mud0ValuePV(index) < 0.02;
}

bool ww_mudZPV(unsigned int index, float cut){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    double dzpv = dzPV(hww.mus_vertex_p4().at(index), hww.mus_trk_p4().at(index), hww.vtxs_position().at(vtxIndex));
    return fabs(dzpv)<cut;
}

bool ww_muId(unsigned int index, bool useMVAmuId, MuonIDMVA *mva){ 
    if (useMVAmuId){
      if (!goodMuonTMVA(mva,index)) return false;
      return true;
    }
    
	if (((hww.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (hww.trks_nlayers().at(hww.mus_trkidx().at(index)) < 6) return false; // # of tracker hits 
    if (hww.mus_ptErr().at(index)/hww.mus_trk_p4().at(index).pt()>0.1) return false; // Does pt come from track?
    if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(index))==0) return false;
    if (hww.mus_trkKink().at(index) > 20.) return false; //kink finder
    if (!hww.mus_pid_PFMuon().at(index)) return false; // should be a pfmuon
    // global muon
    bool goodMuonGlobalMuon = false;
    if (((hww.mus_type().at(index)) & (1<<1)) == (1<<1)){
        goodMuonGlobalMuon = true;
        if (hww.mus_gfit_chi2().at(index)/hww.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
        if (hww.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false;
        if (hww.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
    }
    return goodMuonGlobalMuon || 
        hww.mus_pid_TMLastStationTight().at(index) == 1; // TM id
}

double ww_muIsoVal(unsigned int index){
    double sum =  hww.mus_iso03_sumPt().at(index) +
        hww.mus_iso03_emEt().at(index)  +
        hww.mus_iso03_hadEt().at(index);
    double pt  = hww.mus_p4().at(index).pt();
    return sum/pt;
}

bool ww_muIso(unsigned int index){
    if (hww.mus_p4().at(index).pt()>20) {
        if (TMath::Abs(hww.mus_p4().at(index).eta())<1.479) 
            return muonIsoValuePF(index,0,0.3) < 0.13;
        else 
            return muonIsoValuePF(index,0,0.3) < 0.09;
    } else {
        if (TMath::Abs(hww.mus_p4().at(index).eta())<1.479) 
            return muonIsoValuePF(index,0,0.3) < 0.06;
        else 
            return muonIsoValuePF(index,0,0.3) < 0.05;
    }
}

bool ww_muIso(unsigned int index, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
	return passMuonRingsMVA(index, muonMVAEstimator, IdentifiedMu, IdentifiedEle); 
}



bool goodElectronTMVA(EGammaMvaEleEstimator* egammaMvaEleEstimator, int useMVAeleId, unsigned int i) 
{  

	float pt = hww.els_p4().at(i).pt();
  	float etaSC = hww.els_etaSC().at(i);

	//preselection
	if (fabs(etaSC)<1.479) {
		if (hww.els_sigmaIEtaIEta().at(i)>0.01 || 
				fabs(hww.els_dEtaIn().at(i))>0.007 ||
				fabs(hww.els_dPhiIn().at(i))>0.15 ||
				hww.els_hOverE().at(i)>0.12 ||
				hww.els_tkIso().at(i)/pt>0.2 ||
				TMath::Max(hww.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
				hww.els_hcalIso().at(i)/pt>0.20 ) return 0;
	} else {
		if (hww.els_sigmaIEtaIEta().at(i)>0.03 || 
				fabs(hww.els_dEtaIn().at(i))>0.009 ||
				fabs(hww.els_dPhiIn().at(i))>0.10 ||
				hww.els_hOverE().at(i)>0.10 ||
				hww.els_tkIso().at(i)/pt>0.2 ||
				hww.els_ecalIso().at(i)/pt>0.20 ||
				hww.els_hcalIso().at(i)/pt>0.20 ) return 0;
	}

	// MIT conversion
	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( hww.els_exp_innerlayers().at(i) > 0 ) return false;

	double mvavalue =  egammaMvaEleEstimator->mvaValue(i,false);

	if( pt > 20 ) {
		if( fabs(etaSC)>=1.479 && mvavalue>0.92)  return true;
		if( fabs(etaSC)>=0.8 && fabs(etaSC)<1.479 && mvavalue>0.85)  return true;
		if( fabs(etaSC)<0.8 && mvavalue>0.94)  return true;
		return false;
	}
	else {
		if( fabs(etaSC)>=1.479 && mvavalue>0.62)  return true;
		if( fabs(etaSC)>=0.8 && fabs(etaSC)<1.479 && mvavalue>0.1)  return true;
		if( fabs(etaSC)<0.8 && mvavalue>0.0)  return true;
		return false;
	}

	cout << "Something is wrong. You should not see this! " << endl; 
}

bool passMuonRingsMVA(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle)
{
	double mvavalue = muonMVAEstimator->mvaValueIso( mu, hww.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC,
	                                              	 IdentifiedEle, IdentifiedMu, false );

	double pt 	= hww.mus_trk_p4().at(mu).pt();
	double eta 	= hww.mus_trk_p4().at(mu).eta();

	if( pt>20 ) {
		if( fabs(eta)>=1.479 && fabs(eta)<2.4 && mvavalue>0.86 )  return true;
		if( fabs(eta)<1.479 && mvavalue>0.82 )  return true;
		return false;
	}
	else {
		if( fabs(eta)>=1.479 && fabs(eta)<2.4 && mvavalue>0.82 )  return true;
		if( fabs(eta)<1.479 && mvavalue>0.86 )  return true;
		return false;
	}	
	
	cout << "Something is wrong. You should not see this! " << endl; 
}

bool passMuonRingsMVAFO(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle)
{
	double mvavalue = muonMVAEstimator->mvaValueIso( mu, hww.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC,
	                                               	 IdentifiedEle, IdentifiedMu, false );

	if( mvavalue>-0.6 )  return true;
	return false;
}

bool MuonFOV2(	unsigned int i, MuonMVAEstimator* muonMVAEstimator, 
				std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){

	if (((hww.mus_type().at(i)) & (1<<2)) == 0)    return false; // tracker muon
    if (hww.trks_nlayers().at(hww.mus_trkidx().at(i)) < 6) return false; // # of tracker hits 
    if (hww.mus_ptErr().at(i)/hww.mus_trk_p4().at(i).pt()>0.1) return false; // Does pt come from track?
    if (hww.trks_valid_pixelhits().at(hww.mus_trkidx().at(i))==0) return false;
    if (hww.mus_trkKink().at(i) > 20.) return false; //kink finder
    if (!hww.mus_pid_PFMuon().at(i)) return false; // should be a pfmuon
    // global muon
    bool goodMuonGlobalMuon = false;
    if (((hww.mus_type().at(i)) & (1<<1)) == (1<<1)) {
        goodMuonGlobalMuon = true;
        if (hww.mus_gfit_chi2().at(i)/hww.mus_gfit_ndof().at(i) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
        if (hww.mus_gfit_validSTAHits().at(i)==0 ) goodMuonGlobalMuon = false;
        if (hww.mus_nmatches().at(i)<2) goodMuonGlobalMuon = false;
    }

    return 	(goodMuonGlobalMuon || hww.mus_pid_TMLastStationTight().at(i) == 1) 	&& // ---> Id
			ww_muBase(i) 															&& 
			ww_mud0ValuePV(i)<0.2 													&& 
			ww_mudZPV(i) 															&& 
  			passMuonRingsMVAFO(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
}


bool fakableMuon(unsigned int i, MuFOTypes type,MuonMVAEstimator* muonMVAEstimator,
                std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
    if ( hww.mus_p4().at(i).pt() < 10.0 ) return false;
    switch (type){
        case MuFOV1: return muonId(i, muonSelectionFO_mu_smurf_10);
        //case MuFOV2: return muonId(i, muonSelectionFO_mu_smurf_04);
        case MuFOV2: return MuonFOV2(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
    }
    return false;
}


std::vector<LeptonPair> getExtraLeptons(int i_hyp, double minPt,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
										MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle ){

	std::vector<LeptonPair> leptons;
    for (int i=0; i < int(hww.mus_charge().size()); ++i) {
        if ( hww.mus_p4().at(i).pt() < minPt ) continue;
        if ( TMath::Abs(hww.hyp_lt_id().at(i_hyp)) == 13 && hww.hyp_lt_index().at(i_hyp) == i ) continue;
        if ( TMath::Abs(hww.hyp_ll_id().at(i_hyp)) == 13 && hww.hyp_ll_index().at(i_hyp) == i ) continue;
        //if ( ! (ww_mud0PV(i) && ww_muId(i, useMVAmuId, mumva) && ww_muIso(i)&&
        if ( ! (ww_mud0PV(i) && ww_muId(i, useMVAmuId, mumva) && ww_muIso(i, muonMVAEstimator, IdentifiedMu,  IdentifiedEle) &&
                    fabs(hww.mus_p4().at(i).eta()) <2.4) ) continue;
        leptons.push_back(LeptonPair(true,i));
    }
    for (int i=0; i < int(hww.els_charge().size()); ++i) {
        if ( hww.els_p4().at(i).pt() < minPt ) continue;
        if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_lt_p4().at(i_hyp),hww.els_p4().at(i)) <0.1) ) continue;
        if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_ll_p4().at(i_hyp),hww.els_p4().at(i)) <0.1) ) continue;
        if ( !(ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_elIso(i) &&
                    fabs(hww.els_p4().at(i).eta()) < 2.5) ) continue;
        leptons.push_back(LeptonPair(false,i));
    }
    return leptons;
}

unsigned int numberOfExtraLeptons(	int i_hyp, double minPt, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
									MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
	return getExtraLeptons(i_hyp, minPt, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, mumva,muonMVAEstimator,IdentifiedMu,IdentifiedEle).size();
}


unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated)
{

const std::vector<JetPair> vetojets = std::vector<JetPair>();  //empty, so there is no jet veto 

    unsigned int nMuons = 0;
    for (int imu=0; imu < int(hww.mus_charge().size()); ++imu) {
        // quality cuts
        if (  ((hww.mus_goodmask().at(imu)) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
        if ( hww.mus_p4().at(imu).pt() < 3 ) continue;
        if ( ww_mud0ValuePV(imu) > 0.2) continue;
        if ( ! ww_mudZPV(imu,0.2) ) continue; //newcuts, was 0.1
        if (hww.trks_nlayers().at(hww.mus_trkidx().at(imu)) < 6) return false; // # of tracker hits 
        if ( TMath::Abs(hww.hyp_lt_id().at(i_hyp)) == 13 && hww.hyp_lt_index().at(i_hyp) == imu ) continue;
        if ( TMath::Abs(hww.hyp_ll_id().at(i_hyp)) == 13 && hww.hyp_ll_index().at(i_hyp) == imu ) continue;
        if ( nonisolated && ww_muIsoVal(imu)<0.1 && hww.mus_p4().at(imu).pt()>20 ) continue;
        bool skip = false;
        for ( std::vector<JetPair>::const_iterator ijet = vetojets.begin(); ijet != vetojets.end(); ++ijet ){
            if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,hww.mus_p4().at(imu))) < 0.3 ) skip=true;
        }
        if ( skip ) continue;
        ++nMuons;
    }   
    return nMuons;
}


double nearestDeltaPhi(double Phi, int i_hyp)
{
    double tightDPhi = fabs(hww.hyp_lt_p4().at(i_hyp).Phi() - Phi);
    tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
    double looseDPhi = fabs(hww.hyp_ll_p4().at(i_hyp).Phi() - Phi);
    looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
    return TMath::Min(tightDPhi, looseDPhi);
}


double projectedMet(unsigned int i_hyp, double met, double phi)
{
    double DeltaPhi = nearestDeltaPhi(phi,i_hyp);
    if (DeltaPhi < TMath::Pi()/2) return met*TMath::Sin(DeltaPhi);
    return met;
}



bool toptag(WWJetType type, int i_hyp, double minPt, std::vector<JetPair> ignoreJets)
{
    const double vetoCone    = 0.3;
    // bug fix for mva jet id
    vector <float> fixedpfjetmva_analsel; getGoodMVAs(fixedpfjetmva_analsel, "mvavalue"); 
    
    switch ( type ){
        case pfJet: 
            for ( unsigned int i=0; i < hww.pfjets_p4().size(); ++i) {
                
                if ( hww.pfjets_p4().at(i).pt() < minPt ) continue;
                bool ignoreJet = false;
                for ( std::vector<JetPair>::const_iterator ijet = ignoreJets.begin();
                        ijet != ignoreJets.end(); ++ijet )
                    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,hww.pfjets_p4().at(i))) < vetoCone ) ignoreJet=true;
                if ( ignoreJet ) continue;

                if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_lt_p4().at(i_hyp),hww.pfjets_p4().at(i))) < vetoCone ||
						TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hww.hyp_ll_p4().at(i_hyp),hww.pfjets_p4().at(i))) < vetoCone ) continue;

                double jec = hww.pfjets_JEC().at(i);
               
				if ( !passMVAJetId( hww.pfjets_p4().at(i).pt() * jec, hww.pfjets_p4().at(i).eta(), fixedpfjetmva_analsel[i], 2) ) continue;

				if ( !defaultBTag(type,i, jec) ) continue;
				
                return true;
            }
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    }
    return false;
}
