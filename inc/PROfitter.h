#ifndef PROFITTER_H
#define PROFITTER_H

#include "PROmetric.h"

#include <Eigen/Eigen>
#include "LBFGSB.h"

namespace PROfit {

class PROfitter {
public:
    Eigen::VectorXd ub, lb, best_fit;
    LBFGSpp::LBFGSBParam<double> param;
    LBFGSpp::LBFGSBSolver<double> solver;
    int n_multistart = 100, n_localfit = 5;

    PROfitter(const Eigen::VectorXd ub, const Eigen::VectorXd lb, const LBFGSpp::LBFGSBParam<double> &param = {})
        : ub(ub), lb(lb), param(param), solver(param) {}

    double Fit(PROmetric &metric);

    Eigen::VectorXd FinalGradient() const {return solver.final_grad();}
    double FinalGradientNorm() const {return solver.final_grad_norm();}
    Eigen::MatrixXd Hessian() const {return solver.final_approx_hessian();}
    Eigen::MatrixXd InverseHessian() const {return solver.final_approx_inverse_hessian();}
    Eigen::MatrixXd Covariance() const {return InverseHessian();}

    // If you don't belive the uncertainties on the parameters, you can use the final fit value to estimate the variance
    Eigen::MatrixXd ScaledCovariance(double chi2, int n_datapoint) const {return Covariance()*chi2/double(n_datapoint-best_fit.size());}

};

}

#endif

