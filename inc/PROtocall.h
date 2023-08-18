#ifndef PRO_TO_CALL_H
#define PRO_TO_CALL_H

// C++ include 
#include <algorithm>
#include <unordered_map>
#include <string>

// PROfit include 
#include "PROlog.h"
#include "PROconfig.h"
#include "PROsc.h"
#include "PROpeller.h"
#include "PROspec.h"
#include "PROsyst.h"

namespace PROfit{


    /* Function: given a value for reconstructed variable, figure out which local bin in the histogram it belongs to
     * Note: bin index start from 0, not 1
     * Note: return value of -1 means the reco value is out of range
     */
    int FindLocalBin(const PROconfig &inconfig, double reco_value, int channel_index);
    /* Function: given a value for reconstructed variable, figure out which global bin it belongs to
     * Note: bin index start from 0, not 1
     * Note: if  the reco value is out of range, then return value of -1
     */
    int FindGlobalBin(const PROconfig &inconfig, double reco_value, int subchannel_index);
    int FindGlobalBin(const PROconfig &inconfig, double reco_value, const std::string& subchannel_fullname);

    /* Function: given a value for true variable, figure out which local bin in the histogram it belongs to
     * Note: bin index start from 0, not 1
     * Note: return value of -1 means the true value is out of range
     */
    int FindLocalTrueBin(const PROconfig &inconfig, double true_value, int channel_index);

    /* Function: given a value for truenstructed variable, figure out which global bin it belongs to
     * Note: bin index start from 0, not 1
     * Note: if  the true value is out of range, then return value of -1
     */
    int FindGlobalTrueBin(const PROconfig &inconfig, double true_value, int subchannel_index);
    int FindGlobalTrueBin(const PROconfig &inconfig, double true_value, const std::string& subchannel_fullname);

    /* Function: this is the master weighting function that combines all weights before filling into spectrum.
     * Given an event */
    float FillRecoSpectra(const PROconfig &infonfig, const PROpeller &inprop, const PROsyst &insyst, PROsc &inosc, std::vector<string> &insystname, std::vector<float> &inshifts, std::vector<float> physparams);
    float GetOscWeight(int ev_idx, const PROpeller &inprop, PROsc &inosc, std::vector<float> &inphysparams);
    float GetSystWeight(int ev_idx, const PROpeller &inprop, const PROsyst &insyst, float inshift, std::string insystname);

};

#endif
