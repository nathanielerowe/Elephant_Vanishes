#include "PROsurf.h"
#include "PROconfig.h"
#include "PROsc.h"
#include "PROspec.h"
#include "PROpeller.h"
#include "PROsyst.h"
#include "PROchi.h"

#include "LBFGSB.h"

#include <Eigen/Eigen>

#include <random>

using namespace PROfit;

PROsurf::PROsurf(size_t nbinsx, LogLin llx, double x_lo, double x_hi, size_t nbinsy, LogLin lly, double y_lo, double y_hi) : nbinsx(nbinsx), nbinsy(nbinsy), edges_x(Eigen::VectorXd::Constant(nbinsx + 1, 0)), edges_y(Eigen::VectorXd::Constant(nbinsy + 1, 0)), surface(nbinsx, nbinsy) {
    if(llx == LogAxis) {
        x_lo = std::log10(x_lo);
        x_hi = std::log10(x_hi);
    }
    if(lly == LogAxis) {
        y_lo = std::log10(y_lo);
        y_hi = std::log10(y_hi);
    }
    for(size_t i = 0; i < nbinsx + 1; i++)
        edges_x(i) = x_lo + i * (x_hi - x_lo) / nbinsx;
    for(size_t i = 0; i < nbinsy + 1; i++)
        edges_y(i) = y_lo + i * (y_hi - y_lo) / nbinsy;
}

void PROsurf::FillSurface(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, int nthreads) {
    std::random_device rd{};
    std::mt19937 rng{rd()};
    std::normal_distribution<float> d;
    std::ofstream chi_file;

    if(!filename.empty()){
        chi_file.open(filename);
    }

    for(size_t i = 0; i < nbinsx; i++) {
        for(size_t j = 0; j < nbinsy; j++) {
            std::cout << "Filling point " << i << " " << j << std::endl;
            LBFGSpp::LBFGSBParam<double> param;
            param.epsilon = 1e-6;
            param.max_iterations = 100;
            param.max_linesearch = 50;
            param.delta = 1e-6;
            LBFGSpp::LBFGSBSolver<double> solver(param);
            int nparams = systs.GetNSplines();
            std::vector<float> physics_params = {(float)edges_y(j), (float)edges_x(i)};//deltam^2, sin^22thetamumu
            PROchi chi("3plus1",&config,&prop,&systs,&osc, data, nparams, systs.GetNSplines(), physics_params);
            Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
            Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
            Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);

            double fx;
            int niter;
            std::vector<double> chi2s;
            int nfit = 0;
            do {
                nfit++;
                for(size_t i = 0; i < nparams; ++i)
                    x(i) = 0.3*d(rng);
                // x will be overwritten to be the best point found
                try {
                    niter = solver.minimize(chi, x, fx, lb, ub);
                } catch(std::runtime_error &except) {
                    log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
                    continue;
                }
                chi2s.push_back(fx);
            } while(chi2s.size() < 10 && nfit < 100);
            fx = *std::min_element(chi2s.begin(), chi2s.end());
            surface(i, j) = fx;
            if(!filename.empty()){
                chi_file<<" "<<edges_x(i)<<" "<<edges_y(j)<<" "<<fx;
            }
        }
    }
}

