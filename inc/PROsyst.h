#ifndef PROSYST_H_
#define PROSYST_H_

#include "PROcreate.h"

#include <vector>

namespace PROfit {

class PROsyst {
public:
    using Spline = std::vector<std::vector<std::array<float, 4>>>;

    PROsyst(const std::vector<SystStruct>& systs);

private:
    std::map<std::string, Spline> splines;
    std::map<std::string, Eigen::MatrixXd> covmats;
    //std::<std:string, MFA> mfa;
};

};

#endif
