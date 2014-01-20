#include <algorithm>
#include <utility>

#include "Math/VectorUtil.h"

#include "HWWValidation/HWWBase/interface/jetSelections.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"
#include "HWWValidation/jetcorr/JetCorrectorParameters.icc"
#include "HWWValidation/jetcorr/FactorizedJetCorrector.icc"
#include "HWWValidation/jetcorr/SimpleJetCorrector.icc"
#include "HWWValidation/jetcorr/JetCorrectionUncertainty.icc"
#include "HWWValidation/jetcorr/SimpleJetCorrectionUncertainty.icc"

using std::vector;
using std::pair;

// define this type for speed: allows us to get a vector of selected
// jets, potentially with a correction factor, without having to make
// copies
typedef vector<pair<const LorentzVector *, double> > jets_with_corr_t;


class FactorizedJetCorrector *makeJetCorrector (const char *l2corr, 
                                                const char *l3corr, 
                                                const char *l2l3_residual_corr)
{
    std::vector<std::string> corrs;
    corrs.reserve(3);
    corrs.push_back(l2corr);
    corrs.push_back(l3corr);
    corrs.push_back(l2l3_residual_corr);
    return makeJetCorrector(corrs);
}

class FactorizedJetCorrector *makeJetCorrector (const std::vector<std::string> &corrs)
{
    vector<JetCorrectorParameters> vParam;
    for (std::vector<std::string>::const_iterator i = corrs.begin(), i_end = corrs.end();
         i != i_end; ++i) {
        // do some rigmarole to evaluate env variables in the strings
        const std::string cmd = "echo ";
        FILE *f = popen((cmd + *i).c_str(), "r");
        if (!f) {
            perror((std::string("Opening pipe to execute ") + cmd + *i).c_str());
            return 0;
        }
        char corr_name[1024];
        int s = fscanf(f, " %1024s\n", corr_name);
        if (s != 1) {
            perror("reading file list");
        }
        assert(s == 1);
        JetCorrectorParameters JetCorPar(corr_name);
        // printf("%s\n", corr_name);
        vParam.push_back(JetCorrectorParameters(corr_name));
    }
    return new FactorizedJetCorrector(vParam);
}
