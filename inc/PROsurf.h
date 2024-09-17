#ifndef PROSURF_H
#define PROSURF_H

#include "PROconfig.h"
#include "PROsc.h"
#include "PROpeller.h"
#include "PROspec.h"
#include "PROsyst.h"

#include <Eigen/Eigen>

namespace PROfit {

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

    void FillSurfaceFast(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, int nthreads = 1);
    void FillSurface(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, int nthreads = 1);

};

}

#endif

