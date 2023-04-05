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

        SystStruct(const std::string& in_systname, const int in_n_univ): SystStruct(in_systname, in_n_univ, "multisim", "1"){}
        SystStruct(const std::string& in_systname, const int in_n_univ, const std::string& in_mode, const std::string& in_weight_formula): systname(in_systname), n_univ(in_n_univ), mode(in_mode), weight_formula(in_weight_formula){}
        std::string systname;
        int n_univ;
        std::string mode;
        std::string weight_formula;
        //map
        //hist
	std::vector<std::vector<eweight_type>> multi_vecspec;

	inline
	void SetMode(const std::string& in_mode){mode = in_mode; return;}

	inline
	void SetWeightFormula(const std::string& in_formula){weight_formula = in_formula; return;}

        std::vector<std::vector<eweight_type>> GetCovVec();
        std::vector<eweight_type> GetKnobs(int index, std::string variation);

	void SanityCheck() const;
	void CleanSpecs();
	void SetSpecDimension(int row, int col);
    };


    int PROcess(const PROconfig &inconfig);
    void ProcessEvent(const PROconfig &inconfig,
        const std::map<std::string, 
        std::vector<eweight_type> >& thisfWeight,
        size_t fileid,
        int entryid);
    };
#endif
