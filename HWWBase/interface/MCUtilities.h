#ifndef MCUTILITIES_H
#define MCUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <vector>

typedef math::XYZTLorentzVectorF LorentzVector;

class MCUtilities {
public:
  MCUtilities();
  ~MCUtilities();

  static const reco::GenParticle* motherID(const reco::GenParticle& gp);
  static void writeDaughter(const reco::GenParticle& gp, int idx, std::vector<int>& genps_ld_id,
			    std::vector<int>& genps_ld_idx, std::vector<LorentzVector>& genps_ld_p4 );

private:

};

#endif
