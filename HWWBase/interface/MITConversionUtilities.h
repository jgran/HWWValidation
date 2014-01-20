#ifndef MITCONVERSIONUTILITIES_H
#define MITCONVERSIONUTILITIES_H

#include <stdint.h>
#include <vector>

namespace HWWFunctions {

bool isMITConversion(unsigned int elidx, 
		     int nWrongHitsMax, 
		     float probMin,
		     float dlMin,
		     bool matchCTF,
		     bool requireArbitratedMerged);

}

#endif
