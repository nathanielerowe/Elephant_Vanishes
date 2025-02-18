#ifndef PROCESS_H
#define PROCESS_H

#include <Eigen/Eigen>

// PROfit include 
#include "PROconfig.h"
#include "PROmodel.h"
#include "PROpeller.h"
#include "PROspec.h"
#include "PROsyst.h"

namespace PROfit{

    /* Function: 
     *  The master weighting function that combines all weights and fills into spectrum PROspec, event-by-event
     */

    PROspec FillCVSpectrum(const PROconfig &inconfig, const PROpeller &inprop, bool binned = false);
    PROspec FillRecoSpectra(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst, const PROmodel &inmodel, const Eigen::VectorXf &params, bool binned = true);
    PROspec FillSystRandomThrow(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst);
};

#endif
