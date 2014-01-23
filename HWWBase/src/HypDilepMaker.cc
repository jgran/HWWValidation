#include "DataFormats/JetReco/interface/CaloJet.h"
#include "HWWValidation/HWWBase/interface/HypDilepMaker.h"
#include "HWWValidation/HWWBase/interface/HWW.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

bool testJetForLeptons(const LorentzVector& jetP4, const LorentzVector& lepp4) {
  
  
  bool matched = false;
  float lepphi  = lepp4.Phi();
  float jetphi = jetP4.Phi();
   
  float lepeta  = lepp4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = lepphi - jetphi;
  float deta = lepeta - jeteta;
  if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}

void HypDilepMaker::SetVars(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  HWWVal::Load_hyp_jets_p4();
  HWWVal::Load_hyp_type();
  HWWVal::Load_hyp_p4();
  HWWVal::Load_hyp_lt_charge();
  HWWVal::Load_hyp_lt_index();
  HWWVal::Load_hyp_lt_id();
  HWWVal::Load_hyp_lt_p4();
  HWWVal::Load_hyp_ll_charge();
  HWWVal::Load_hyp_ll_index();
  HWWVal::Load_hyp_ll_id();
  HWWVal::Load_hyp_ll_p4();

  double looseptcut = 10.0;
  double tightptcut = 20.0;
  double hypJetMaxEtaCut = 5.0;
  double hypJetMinPtCut = 30.0;

  // muon charge
  vector<int> *mus_charge = new vector<int>;
  *mus_charge = HWWVal::mus_charge(); 

  //muon p4
  vector<LorentzVector> *mus_p4 = new vector<LorentzVector>;
  *mus_p4 = HWWVal::mus_p4();

  //muon type
  vector<int> *mus_type = new vector<int>;
  *mus_type = HWWVal::mus_type();

  //-----------------------------------------------------------
  // electron variables
  //-----------------------------------------------------------
  vector<int> *els_charge = new vector<int>;
  *els_charge = HWWVal::els_charge();

  // electron p4
  vector<LorentzVector> *els_p4 = new vector<LorentzVector>;
  *els_p4 = HWWVal::els_p4();

  unsigned int nmus = mus_p4->size();
  unsigned int nels = els_p4->size();


  vector<LorentzVector> *jets_p4 = new vector<LorentzVector>;
  *jets_p4 = HWWVal::pfjets_p4();

  //------------------------------------------------------------
  // loop over the muons
  //------------------------------------------------------------
  //get the candidates and make hypotheses 
  for(unsigned int mus_index_1 = 0; mus_index_1 < nmus; mus_index_1++) {//first muon loop
    for(unsigned int mus_index_2 = 0; mus_index_2 < nmus; mus_index_2++) {//second muon loop

      if(mus_index_1 == mus_index_2) continue;
      if(mus_index_2 < mus_index_1)  continue;  //avoid double counting

      //don't look at standalone muons
      if(mus_type->at(mus_index_1) == 8) continue;
      if(mus_type->at(mus_index_2) == 8) continue;
      
      float mu_pt1 = mus_p4->at(mus_index_1).Pt();
      float mu_pt2 = mus_p4->at(mus_index_2).Pt();
      
      //if either fail the loose cut, go to the next muon
      if(mu_pt1 < looseptcut || mu_pt2 < looseptcut) continue;
      
      //if neither one passes the tight cut, go to the next muon
      if(mu_pt1 < tightptcut && mu_pt2 < tightptcut) continue;

      int tight_index = mus_index_1;
      int loose_index = mus_index_2;

      /*
	figure out which one should be tight and which should
	be loose in case one passes the tight cut and the other 
	does not
      */
      if(mu_pt1 < tightptcut && mu_pt2 > tightptcut) {
        tight_index = mus_index_2;
        loose_index = mus_index_1;
      }
      if(mu_pt2 < tightptcut && mu_pt1 > tightptcut) {
        tight_index = mus_index_1;
        loose_index = mus_index_2;
      }


      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<LorentzVector>  temp_jets_p4;      

	
      for(unsigned int i = 0; i<jets_p4->size(); i++) {
	
        // we don't want jets that overlap with electrons
        bool overlapsWithLepton = false;
        if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(loose_index))) 
          overlapsWithLepton = true;
        if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(tight_index))) 
          overlapsWithLepton = true;

        double jet_eta = jets_p4->at(i).eta();
        double jet_pt  = jets_p4->at(i).Pt();
        
        if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { //hyp jetas
          temp_jets_idx.push_back(i);
          temp_jets_p4                     .push_back(jets_p4              ->at(i));
        }
      }

      HWWVal::hyp_jets_p4()       .push_back(temp_jets_p4                          );
      HWWVal::hyp_type()          .push_back(0                                     );
      HWWVal::hyp_p4()            .push_back(mus_p4->at(tight_index)+mus_p4->at(loose_index));
      HWWVal::hyp_lt_charge()     .push_back(mus_charge       ->at(tight_index)  );
      HWWVal::hyp_lt_index()      .push_back(tight_index                         );
      HWWVal::hyp_lt_id()         .push_back(-13*(mus_charge   ->at(tight_index)));
      HWWVal::hyp_lt_p4()         .push_back(mus_p4           ->at(tight_index)  );
      HWWVal::hyp_ll_charge()     .push_back(mus_charge       ->at(loose_index)  );
      HWWVal::hyp_ll_index()      .push_back(loose_index                         );
      HWWVal::hyp_ll_id()         .push_back(-13*(mus_charge   ->at(loose_index)));
      HWWVal::hyp_ll_p4()         .push_back(mus_p4           ->at(loose_index)  );
    }
  }  

  //------------------------------------------------------------
  // loop over the elecrons
  //------------------------------------------------------------
  //get the candidates and make hypotheses 
  for(unsigned int els_index_1 = 0; els_index_1 < nels; els_index_1++) {
    for(unsigned int els_index_2 = 0; els_index_2 < nels; els_index_2++) {
      
      if(els_index_1 == els_index_2) continue;
      if(els_index_2 < els_index_1)  continue;  //avoid double counting
      
      float el_pt1 = els_p4->at(els_index_1).Pt();
      float el_pt2 = els_p4->at(els_index_2).Pt();
      
      //if either fail the loose cut, go to the next muon
      if(el_pt1 < looseptcut || el_pt2 < looseptcut) continue;
      
      //if neither one passes the tight cut, continue
      if(el_pt1 < tightptcut && el_pt2 < tightptcut) continue;
      
      int tight_index = els_index_1;
      int loose_index = els_index_2;
      
      /*
	figure out which one should be tight and which should
	be loose in case one passes the tight cut and the other 
	does not
      */
      if(el_pt1 < tightptcut && el_pt2 > tightptcut) {
        tight_index = els_index_2;
        loose_index = els_index_1;
      }
      if(el_pt2 < tightptcut && el_pt1 > tightptcut) {
        tight_index = els_index_1;
        loose_index = els_index_2;
      }


      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<LorentzVector>  temp_jets_p4;      

	
      for(unsigned int i = 0; i<jets_p4->size(); i++) {
	
        // we don't want jets that overlap with electrons
        bool overlapsWithLepton = false;
        if(!testJetForLeptons(jets_p4->at(i), els_p4->at(loose_index))) 
          overlapsWithLepton = true;
        if(!testJetForLeptons(jets_p4->at(i), els_p4->at(tight_index))) 
          overlapsWithLepton = true;

        double jet_eta = jets_p4->at(i).eta();
        double jet_pt  = jets_p4->at(i).Pt();
        
        if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { //hyp jetas
          temp_jets_idx.push_back(i);
          temp_jets_p4                     .push_back(jets_p4              ->at(i));
        }
      }

      HWWVal::hyp_jets_p4()       .push_back(temp_jets_p4                          );
      HWWVal::hyp_type()          .push_back(3);
      HWWVal::hyp_p4()            .push_back(els_p4->at(tight_index)+els_p4->at(loose_index));
      HWWVal::hyp_lt_charge()     .push_back(els_charge       ->at(tight_index)  );
      HWWVal::hyp_lt_index()      .push_back(tight_index                         );
      HWWVal::hyp_lt_id()         .push_back(-11*(els_charge   ->at(tight_index)));
      HWWVal::hyp_lt_p4()         .push_back(els_p4           ->at(tight_index)  );
      HWWVal::hyp_ll_charge()     .push_back(els_charge       ->at(loose_index)  );
      HWWVal::hyp_ll_index()      .push_back(loose_index                         );
      HWWVal::hyp_ll_id()         .push_back(-11*(els_charge   ->at(loose_index)));
      HWWVal::hyp_ll_p4()         .push_back(els_p4           ->at(loose_index)  );
    }
  }  
  
  /*------------------------------------------------------------
    The EMu, MuE cases
    To avoid double counting, only make MuE if Mu is tight and E is loose
  */

  for(unsigned int els_index = 0; els_index < nels; els_index++) {
    for(unsigned int mus_index = 0; mus_index < nmus; mus_index++) {

      if(mus_type->at(mus_index) == 8) continue;

      float el_pt = els_p4->at(els_index).Pt();
      float mu_pt = mus_p4->at(mus_index).Pt();

      //if either fail the loose cut, go to the next muon
      if(el_pt < looseptcut || mu_pt < looseptcut) continue;

      //if both fail the tight cut, continue
      if(el_pt < tightptcut && mu_pt < tightptcut) continue;
      
      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<LorentzVector>  temp_jets_p4;      

	
      for(unsigned int i = 0; i<jets_p4->size(); i++) {
	
        // we don't want jets that overlap with electrons
        bool overlapsWithLepton = false;
        if(!testJetForLeptons(jets_p4->at(i), els_p4->at(els_index))) 
          overlapsWithLepton = true;
        if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(mus_index))) 
          overlapsWithLepton = true;

        double jet_eta = jets_p4->at(i).eta();
        double jet_pt  = jets_p4->at(i).Pt();
        
        if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { //hyp jetas
          temp_jets_idx.push_back(i);
          temp_jets_p4                     .push_back(jets_p4              ->at(i));
        }
      }

      HWWVal::hyp_jets_p4()       .push_back(temp_jets_p4                          );
      HWWVal::hyp_p4()            .push_back(mus_p4->at(mus_index)+els_p4->at(els_index)                 );
	
      if(el_pt < tightptcut && mu_pt > tightptcut) {
        HWWVal::hyp_type()            .push_back(1);
        HWWVal::hyp_lt_charge()       .push_back(mus_charge       ->at(mus_index)  );
        HWWVal::hyp_lt_index()        .push_back(mus_index                         );
        HWWVal::hyp_lt_id()           .push_back(-13*(mus_charge   ->at(mus_index)));
        HWWVal::hyp_lt_p4()           .push_back(mus_p4           ->at(mus_index)  );
        HWWVal::hyp_ll_charge()       .push_back(els_charge       ->at(els_index)  );
        HWWVal::hyp_ll_index()        .push_back(els_index                         );
        HWWVal::hyp_ll_id()           .push_back(-11*(els_charge   ->at(els_index)));
        HWWVal::hyp_ll_p4()           .push_back(els_p4           ->at(els_index)  );
	
	  
      } else {
        HWWVal::hyp_type()            .push_back(2);
        HWWVal::hyp_lt_charge()       .push_back(els_charge       ->at(els_index)  );
        HWWVal::hyp_lt_index()        .push_back(els_index                         );
        HWWVal::hyp_lt_id()           .push_back(-11*(els_charge   ->at(els_index)));
        HWWVal::hyp_lt_p4()           .push_back(els_p4           ->at(els_index)  );
        HWWVal::hyp_ll_charge()       .push_back(mus_charge       ->at(mus_index)  );
        HWWVal::hyp_ll_index()        .push_back(mus_index                         );
        HWWVal::hyp_ll_id()           .push_back(-13*(mus_charge   ->at(mus_index)));
        HWWVal::hyp_ll_p4()           .push_back(mus_p4           ->at(mus_index)  );
	
      }
    }
  }
}

