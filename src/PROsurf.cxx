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

std::vector<std::vector<double>> latin_hypercube_sampling(size_t num_samples, size_t dimensions, std::uniform_real_distribution<float>&dis, std::mt19937 &gen) {
    std::vector<std::vector<double>> samples(num_samples, std::vector<double>(dimensions));

    for (size_t d = 0; d < dimensions; ++d) {

        std::vector<double> perm(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            perm[i] = (i + dis(gen)) / num_samples;  
        }
        std::shuffle(perm.begin(), perm.end(), gen);  
        for (size_t i = 0; i < num_samples; ++i) {
            samples[i][d] = perm[i]; 
        }
    }

    return samples;
}

std::vector<int> sorted_indices(const std::vector<double>& vec) {
    std::vector<int> indices(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&vec](int i1, int i2) { return vec[i1] < vec[i2]; });
    return indices;
}



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

void PROsurf::FillSurface(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, int nthreads) {
    std::random_device rd{};
    std::mt19937 rng{rd()};
    std::normal_distribution<float> d;
    std::uniform_real_distribution<float> d_uni(-2.0, 2.0);

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
            std::vector<float> physics_params = {edges_y(i),edges_x(i)};//deltam^2, sin^22thetamumu
            PROchi chi("3plus1",&config,&prop,&systs,&osc, data, nparams, systs.GetNSplines(), physics_params);
            Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
            Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
            Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);
            Eigen::VectorXd grad = Eigen::VectorXd::Constant(nparams, 0.0);


            //First do 100 simple function calls. 
            double fx;
            int niter;
            int N_multistart = 100;
            std::vector<double> chi2s_multistart;
            std::vector<std::vector<double>> latin_samples = latin_hypercube_sampling(N_multistart, nparams,d_uni,rng);
            for(int s=0; s<N_multistart; s++){
                x = Eigen::Map<Eigen::VectorXd>( latin_samples[s].data(), latin_samples[s].size());
                fx =  chi(x,grad);
                chi2s_multistart.push_back(fx);
            }
            std::vector<int> best_multistart = sorted_indices(chi2s_multistart);    
    

            int N_localfits = 5;
            std::vector<double> chi2s_localfits;
            int nfit = 0;
            for(int s=0; s<N_localfits; s++){
              x = Eigen::Map<Eigen::VectorXd>( latin_samples[best_multistart[s]].data(), latin_samples[best_multistart[s]].size());   
              try {
                    niter = solver.minimize(chi, x, fx, lb, ub);
                } catch(std::runtime_error &except) {
                    log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
                    continue;
                }
                chi2s_localfits.push_back(fx);
            }

            fx = *std::min_element(chi2s_localfits.begin(), chi2s_localfits.end());
            surface(i, j) = fx;
        }
    }
}

