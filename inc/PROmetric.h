#ifndef PROMETRIC_H
#define PROMETRIC_H

#include "PROsyst.h"

#include <Eigen/Eigen>

namespace PROfit {

class PROmetric {
public:
    enum EvalStrategy {
        EventByEvent,
        BinnedGrad,
        BinnedChi2
    };

    virtual int nParams() const = 0;
    virtual void set_physics_param_fixed(const std::vector<float> &physics_param) = 0;
    virtual void override_systs(const PROsyst &new_syst) = 0;
    virtual float operator()(const Eigen::VectorXf &param, Eigen::VectorXf &gradient) = 0;
    virtual float operator()(const Eigen::VectorXf &param, Eigen::VectorXf &gradient, bool nograd) = 0;
    virtual void reset() = 0;
    virtual PROmetric *Clone() const = 0;
    virtual ~PROmetric() {}
};

};

#endif

