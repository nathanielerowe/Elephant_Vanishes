#ifndef PROCESS_WEIGHTS_H_
#define PROCESS_WEIGHTS_H_

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>

// TINYXML2
#include "tinyxml2.h"

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"

namespace PROfit{

    struct SystStruct {

        SystStruct(const std::string& in_systname, int in_n_univ, const std::string& in_mode, const std::string& in_formula): systname(in_systname), n_univ(in_n_univ), mode(in_mode), formula(in_formula){}
        std::string systname;
        int n_univ;
        std::string mode;
        std::string formula;
        //map
        //hist

        std::vector<std::vector<float>> GetCovVec();
        std::vector<float> GetKnobs(int index, std::string variation);

    };


    PROcess(const PROconfig &inconfig);

    };
#endif
