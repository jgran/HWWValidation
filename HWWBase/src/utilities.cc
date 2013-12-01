// Header
#include "HWWValidation/HWWBase/interface/utilities.h"

// C++ includes
#include <vector>

// ROOT includes
#include "Math/VectorUtil.h"
#include "HWWValidation/HWWBase/interface/HWW.h"
using std::vector;

//return true if one of the leptons is the same in both hyps
bool hypsOverlap(int idxa, int idxb){
  int idlta   = hww.hyp_lt_id().at(idxa);
  int idlla   = hww.hyp_ll_id().at(idxa);
  int ilta    = hww.hyp_lt_index().at(idxa);
  int illa    = hww.hyp_ll_index().at(idxa);
  int idltb   = hww.hyp_lt_id().at(idxb);
  int idllb   = hww.hyp_ll_id().at(idxb);
  int iltb    = hww.hyp_lt_index().at(idxb);
  int illb    = hww.hyp_ll_index().at(idxb);
  int matches = (idlta == idltb && ilta == iltb) + (idlla == idllb && illa == illb) + (idlta == idllb && ilta == illb) + (idlla == idltb && illa == iltb);
  return matches>0;
}

int match4vector(const LorentzVector &lvec, const vector<LorentzVector> &vec, double cut=10.0 ){

  if( vec.size() == 0 ) return -1;
  //cout << "size of vec = " << vec.size() << endl;
  double dR = cut; 
  double x;
  int iret = -1;
  for ( unsigned int i=0; i < vec.size();++i) {
    x = ROOT::Math::VectorUtil::DeltaR( lvec, vec[i] );
    if (x < dR ) {dR = x; iret = i;}
  }
  return iret;
}

std::vector<LorentzVector> p4sInCone(const LorentzVector &refvec, const std::vector<LorentzVector> &invec, double coneSize=0.5 ) {
  vector<LorentzVector> result;
  if ( invec.size() == 0 ) return result;
  double dR = coneSize; 
  double x = 0.0;
  for ( unsigned int i=0; i < invec.size();++i) {
    x = ROOT::Math::VectorUtil::DeltaR( refvec, invec[i] );
    if (x < dR ) {result.push_back(invec[i]);}
  }
  return result;
}

std::vector<unsigned int> idxInCone(const LorentzVector &refvec, const std::vector<LorentzVector> &invec, double coneSize=0.5 ) {
  vector<unsigned int > result;
  if ( invec.size() == 0 ) return result;
  double dR = coneSize; 
  for ( unsigned int i=0; i < invec.size();++i) {
    if ( ROOT::Math::VectorUtil::DeltaR( refvec, invec[i] ) < dR ) {result.push_back(i);}
  }
  return result;
}
