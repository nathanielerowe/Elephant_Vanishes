#ifndef PROSYST_H_
#define PROSYST_H_

//C++ include 
#include <string>
#include <vector>
#include <unordered_map>

// Our include
#include "PROcreate.h"
#include "PROlog.h"

namespace PROfit {

class PROsyst {
public:
    using Spline = std::vector<std::vector<std::array<float, 4>>>;
   

    //constructor
    PROsyst(){}
    PROsyst(const std::vector<SystStruct>& systs);

    /*Function: add a covariance matrix to covmat_map */
    void AddMatrix(const SystStruct& syst);
    
private:
    std::unordered_map<std::string, Spline> splines;
    std::unordered_map<std::string, Eigen::MatrixXd> covmat_map;
    std::unordered_map<std::string, Eigen::MatrixXd> corrmat_map;
    //std::<std:string, MFA> mfa;
};

};

#endif
