#ifndef PROCREATE_H_
#define PROCREATE_H_

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>
#include <ctime>
#include <cmath>

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"
#include "PROconfig.h"
#include "PROspec.h"
#include "PROtocall.h"

//CAFana
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/SRWeightPSet.h"

namespace PROfit{

    typedef Eigen::Matrix<eweight_type, Eigen::Dynamic, Eigen::Dynamic> multi_spec;

    struct SystStruct {

	//members
        std::string systname;
        int n_univ;
        std::string mode;
        std::string weight_formula;
        std::vector<float> knobval;
        int index;

        std::unique_ptr<multi_spec> p_multi_spec;

	// functions 
        SystStruct(const std::string& in_systname, const int in_n_univ): SystStruct(in_systname, in_n_univ, "multisim", "1",{},0){}
        SystStruct(const std::string& in_systname, const int in_n_univ, const std::string& in_mode, const std::string& in_formula, const std::vector<float>& in_knobval, const int in_index): systname(in_systname), n_univ(in_n_univ), mode(in_mode), weight_formula(in_formula), knobval(in_knobval), index(in_index){}


        inline
            void SetMode(const std::string& in_mode){mode = in_mode; return;}

        inline
            void SetWeightFormula(const std::string& in_formula){weight_formula = in_formula; return;}

        inline
            int GetNUniverse() const {return n_univ;}

        inline 
            const std::string& GetSysName() const {return systname;}

        inline 
            const std::string& GetWeightFormula() const {return weight_formula;}

        std::vector<std::vector<eweight_type>> GetCovVec();
        std::vector<eweight_type> GetKnobs(int index, std::string variation);

        //function might not needed
        void SanityCheck() const;
        void CleanSpecs();
        void CreateSpecs(int row, int col);
    };

    int PROcess_SBNfit(const PROconfig &inconfig);
    int PROcess_CAFana(const PROconfig &inconfig);
    void ProcessEvent(const PROconfig &inconfig, size_t fid, const std::vector<std::map<std::string, std::vector<eweight_type>>* >& thisfWeight,
            std::vector<SystStruct>& syst_vector);


    /* Function: given configuration, generate spectrum at central value. 
     * Note: assume the input config has SBNfit-style files, TODO: check if compatible with CAF-style
     */
    PROspec CreatePROspecCV(const PROconfig& configin);
};
#endif
