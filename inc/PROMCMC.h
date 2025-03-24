#ifndef PROMCMC_H
#define PROMCMC_H

#include "PROmetric.h"
#include "PROsurf.h"
#include "TGraphAsymmErrors.h"
#include <Eigen/Eigen>
#include <random>

namespace PROfit {

template<class Proposal_FN, class Target_FN>
class Metropolis {
private:
    std::mt19937 rng;
    std::uniform_real_distribution<float> uniform;

public:
    //using target_fn = float (*)(Eigen::VectorXf &, PROmetric &);
    //using proposal_fn = Eigen::VectorXf (*)(Eigen::VectorXf &);

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
        float acceptance = std::min(1.0f, target(p, metric)/target(current, metric));
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

    static float simple_target(Eigen::VectorXf &value, PROmetric &metric) {
        return std::exp(-0.5f*metric(value, value, false));
    }

    static Eigen::VectorXf simple_proposal(Eigen::VectorXf &current) {
        static std::random_device rd{};
        static std::mt19937 rng(rd());
        std::normal_distribution<float> nd;
        Eigen::VectorXf ret = current;
        for(int i = 0; i < ret.size(); ++i) {
            // Gaussian proposal here
            // Uniform for physics parameters?
        }
        return ret;
    }

    static auto proposal_from_profile(const PROfile &prof) { 
        return [prof](const Eigen::VectorXf &current) {
            static random_device rd{};
            static mt19937 rng(rd());
            std::normal_distribution<float> normal;
            Eigen::VectorXf ret;
            const TGraphAsymmErrors &p = prof.onesig;
            for(int i = 0; i < p.GetN(); ++i) {
                float val = normal(rng);
                ret(i) = val > 0 ? p.GetPointY(i) + p.GetErrorYhigh(i) * val
                                 : p.GetPointY(i) + p.GetErrorYlow(i) * val;

            }
            return ret;
        };
    }

};

}

#endif

