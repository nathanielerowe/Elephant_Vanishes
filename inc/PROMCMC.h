#ifndef PROMCMC_H
#define PROMCMC_H

#include "PROmetric.h"
#include "PROsurf.h"
#include "TGraphAsymmErrors.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <random>

namespace PROfit {

template<class Proposal_FN, class Target_FN>
class Metropolis {
private:
    std::mt19937 rng;
    std::uniform_real_distribution<float> uniform;

public:
    Target_FN target;
    Proposal_FN proposal;
    Eigen::VectorXf current;
    PROmetric &metric;

    Metropolis(Target_FN target, Proposal_FN proposal, Eigen::VectorXf &initial, PROmetric &metric) 
        : target(target), proposal(proposal), current(initial), metric(metric) {
        std::random_device rd{};
        rng.seed(rd());
    }

    bool step() {
        Eigen::VectorXf p = proposal(current);
        float acceptance = std::min(1.0f, target(p)/target(current) * proposal.P(current, p)/proposal.P(p, current));
        float u = uniform(rng);
        if(u <= acceptance) {
            current = p;
            return true;
        }
        return false;
    }

    void run(size_t burnin, size_t steps, void (*action)(Eigen::VectorXf &)) {
        for(size_t i = 0; i < burnin; i++) {
            step();
        }
        for(size_t i = 0; i < steps; i++) {
            if(step() && action) action(current);
        }
    }

    struct simple_target {
        PROmetric &metric;

        float operator()(Eigen::VectorXf &value) {
            return std::exp(-0.5f*metric(value, value, false));
        }
    };

    struct simple_proposal {
        PROmetric &metric;
        float width;
        std::random_device rd;
        std::mt19937 rng;

        simple_proposal(PROmetric &metric, float width = 0.2) 
            : metric(metric), width(width) {
            rng.seed(rd());
        }

        Eigen::VectorXf operator()(Eigen::VectorXf &current) {
            Eigen::VectorXf ret = current;
            int nparams = metric.GetModel().nparams;
            for(int i = 0; i < ret.size(); ++i) {
                if(i < nparams) {
                    std::uniform_real_distribution<float> ud(metric.GetModel().lb(i), metric.GetModel().ub(i));
                    ret(i) = ud(rng);
                } else {
                    std::normal_distribution<float> nd(current(i), width);
                    float proposed_value = nd(rng);
                    ret(i) = std::clamp(proposed_value, metric.GetSysts().spline_lo[i-nparams], metric.GetSysts().spline_hi[i-nparams]);
                }
            }
            return ret;
        }

        float P(Eigen::VectorXf &value, Eigen::VectorXf &given) {
            float prob = 1.0;
            for(int i = 0; i < value.size(); ++i) {
                if(i < metric.GetModel().nparams) {
                    prob *= 1.0f / (metric.GetModel().ub(i) - metric.GetModel().lb(i));
                } else {
                    prob *= (1.0f / std::sqrt(2 * M_PI * width * width))
                            * std::exp(-(value(i) - given(i))*(value(i) - given(i))/(2 * width * width));
                }
            }
            return prob;
        }
    };
};

}

#endif

