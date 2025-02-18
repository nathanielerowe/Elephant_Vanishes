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
    std::vector<float> lb, ub;
    std::vector<std::function<float(std::vector<float>, float)>> model_rules;
    std::vector<Eigen::MatrixXf> hists;
};

}

#endif

