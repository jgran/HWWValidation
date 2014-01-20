#include "HWWValidation/HWWBase/interface/MITConversionUtilities.h"
#include "HWWValidation/HWWBase/interface/HWW.h"
#include "TMath.h"


bool isMITConversion(unsigned int elidx, 
		     int nWrongHitsMax, 
		     float probMin,
		     float dlMin,
		     bool matchCTF,
		     bool requireArbitratedMerged) {

  unsigned int nconvs = HWWVal::convs_isConverted().size();
  if(nconvs == 0) 
    return false;
  bool isGoodConversion = false;

  for(unsigned int iconv = 0; iconv < nconvs; iconv++) {
    
    bool conversionMatchFound = false;
    for(unsigned int itk = 0; itk < HWWVal::convs_tkidx().at(iconv).size(); itk++) {

      if(HWWVal::convs_tkalgo().at(iconv)[itk] == 29 && HWWVal::convs_tkidx().at(iconv)[itk] == HWWVal::els_gsftrkidx().at(elidx))
	conversionMatchFound = true;
      if(matchCTF) {
	if(HWWVal::convs_tkalgo().at(iconv)[itk] > 3 && HWWVal::convs_tkalgo().at(iconv)[itk] < 14 && HWWVal::convs_tkalgo().at(iconv)[itk] != 12 && HWWVal::convs_tkidx().at(iconv)[itk] == HWWVal::els_trkidx().at(elidx))
	  conversionMatchFound = true;
      }
    
      if(conversionMatchFound)
	break;
    }
    
    
    if(conversionMatchFound==false)
      continue;
    
    if( TMath::Prob( HWWVal::convs_chi2().at(iconv), (Int_t)HWWVal::convs_ndof().at(iconv) )  > probMin && HWWVal::convs_dl().at(iconv) > dlMin ) isGoodConversion = true;
    if(requireArbitratedMerged) {
      if(HWWVal::convs_quality().at(iconv) & 4)
	isGoodConversion = true;
      else 
	isGoodConversion = false;
    }

    for(unsigned int j = 0; j < HWWVal::convs_nHitsBeforeVtx().at(iconv).size(); j++) {
      if(HWWVal::convs_nHitsBeforeVtx().at(iconv)[j] > nWrongHitsMax)
	isGoodConversion = false;
    }
      
    if(isGoodConversion)
      break;
      
      
  }//loop over convserions


  return isGoodConversion;
}
