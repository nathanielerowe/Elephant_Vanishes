#ifndef PROMODEL_H
#define PROMODEL_H

#include <Eigen/Eigen>

#include <functional>
#include <string>
#include <vector>

namespace PROfit {

class PROmodel {
public:
    size_t nparams;
    std::vector<std::string> param_names;
    Eigen::VectorXf lb, ub;
    std::vector<std::function<float(Eigen::VectorXf, float)>> model_functions;
    std::vector<Eigen::MatrixXf> hists;
};

}

#endif

