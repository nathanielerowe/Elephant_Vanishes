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

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"
#include "PROconfig.h"

namespace PROfit{

    struct SystStruct {

        SystStruct(const std::string& in_systname, const int in_n_univ, const std::string& in_mode, const std::string& in_formula, const std::vector<double>& in_knobval, const &int in_index): systname(in_systname), n_univ(in_n_univ), mode(in_mode), formula(in_formula), knobval(in_knobval), index(in_index){}
        std::string systname;
        int n_univ;
        std::string mode;
        std::string formula;
        int index;
        std::vector<double> knobval;
        //map
        //hist

        std::vector<std::vector<float>> GetCovVec();
        std::vector<float> GetKnobs(int index, std::string variation);

    };


    int PROcess(const PROconfig &inconfig);
    int PROcess_CAFana(const PROconfig &inconfig);
    void ProcessEvent(const PROconfig &inconfig,
        const std::map<std::string, 
        std::vector<eweight_type> >& thisfWeight,
        size_t fileid,
        int entryid);
    };
#endif
