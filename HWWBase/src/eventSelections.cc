#include "Math/LorentzVector.h"

#include "HWWValidation/HWWBase/interface/eventSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"
#include "HWWValidation/HWWBase/interface/HWW.h"

//----------------------------------------------------------------
// A ridicolusly simple function, but since the Z veto is used 
// in two places, might as well centralize it to keep consistency
//----------------------------------------------------------------

bool inZmassWindow (float mass) {
    if (mass > 76. && mass < 106.) return true;
    return false;
}

//
// function to select a good vertex
// 
bool isGoodVertex(size_t ivtx) {

  if (hww.vtxs_isFake().at(ivtx)) return false;
  if (hww.vtxs_ndof().at(ivtx) <= 4.) return false;
  if (hww.vtxs_position().at(ivtx).Rho() > 2.0) return false;
  if (fabs(hww.vtxs_position().at(ivtx).Z()) > 24.0) return false;
  return true;

}

//---------------------------------------------------------
//
// Find first good vertex
//
//---------------------------------------------------------
int firstGoodVertex () {
    for (unsigned int vidx = 0; vidx < hww.vtxs_position().size(); vidx++) {
        if (isGoodVertex(vidx))
            return vidx;
    }

    return -1;
}


//*****************************************************************************************/
// number of good vertices in the event
//*****************************************************************************************/
int numberOfGoodVertices(void) {
  int ngv = 0;
  for (unsigned int vidx = 0; vidx < hww.vtxs_position().size(); vidx++) {
    if (isGoodVertex(vidx)) ++ngv;
  }
  return ngv;
}
