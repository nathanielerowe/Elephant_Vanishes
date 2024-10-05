#ifndef PROSURF_H
#define PROSURF_H

#include "PROconfig.h"
#include "PROsc.h"
#include "PROspec.h"
#include "PROpeller.h"
#include "PROsyst.h"
#include "PROchi.h"

#include "LBFGSB.h"

#include <Eigen/Eigen>

#include <random>
#include <thread>
#include <future>

#include "TGraph.h"
#include "TMarker.h"
#include "TLine.h"

namespace PROfit {

struct surfOut{
    std::vector<int> grid_index;
    std::vector<float> grid_val;
    double chi;
};

int PROfile(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename);

class PROsurf {
public:
    size_t nbinsx, nbinsy;
    Eigen::VectorXd edges_x, edges_y;
    Eigen::MatrixXd surface;

    enum LogLin {
        LinAxis,
        LogAxis,
    };

    PROsurf(size_t nbinsx, const Eigen::VectorXd &edges_x, size_t nbinsy, const Eigen::VectorXd &edges_y) : nbinsx(nbinsx), nbinsy(nbinsy), edges_x(edges_x), edges_y(edges_y), surface(nbinsx, nbinsy) { }

    PROsurf(size_t nbinsx, LogLin llx, double x_lo, double x_hi, size_t nbinsy, LogLin lly, double y_lo, double y_hi);

    std::vector<surfOut> PointHelper(const PROconfig *config, const PROpeller *prop, const PROsyst *systs, const PROsc *osc, const PROspec *data, std::vector<surfOut> multi_physics_params, PROchi::EvalStrategy strat, bool binned_weighting, int start, int end);

    void FillSurfaceSimple(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, bool binned_weighting);
    void FillSurface(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, bool binned_weighting, int nthreads = 1);

};

}

#endif

