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


    struct SystStruct {

	//members
        std::string systname;
        int n_univ;
        std::string mode;
        std::string weight_formula;
        std::vector<float> knobval;
        int index;

	// pointer to cv spectrum and multi-universe spectrum from systematic variation
	std::unique_ptr<PROspec> p_cv;	
        std::vector<std::unique_ptr<PROspec>> p_multi_spec;

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
	    bool HasWeightFormula() const {return weight_formula == "1";}

        inline 
            const std::string& GetWeightFormula() const {return weight_formula;}

        std::vector<std::vector<eweight_type>> GetCovVec();
        std::vector<eweight_type> GetKnobs(int index, std::string variation);


        /* Function: check if num of universes of this systematics matches with its type 
 	 * Note: multisim mode can have many universes, while minmax mode can only have 2
 	 */
        void SanityCheck() const;

	/* Function: clean up all the member spectra (but ONLY spectra) */
        void CleanSpecs();

	/* Function: create EMPTY spectra with given length 
 	 */ 
        void CreateSpecs(long int num_bins);


	/* Function: given global bin index, and event weight, fill the central value spectrum */
	void FillCV(long int global_bin, double event_weight);

	/* Function: given global bin index, and event weight, fill the spectrum of given universe */
	void FillUniverse(int universe, long int global_bin, double event_weight);

    };


    /* Function: given config, read files in the xml, and grab all systematic variations 
     * TODO: not finished yet
     */
    int PROcess_SBNfit(const PROconfig &inconfig);
    int PROcess_CAFana(const PROconfig &inconfig);

    int PROcess_CAFana_Event(const PROconfig &inconfig, std::vector<std::unique_ptr<TTreeFormula>> & formulas, std::vector<SystStruct> &syst_vector, double reco_val, double add_weight, int global_bin);

    /* Function: given configuration, generate spectrum at central value. 
     * Note: assume the input config has SBNfit-style files, TODO: check if compatible with CAF-style
     */
    PROspec CreatePROspecCV(const PROconfig& configin);



    /* Function: assume currently reading one entry of a file, update systematic variation spectrum 
     * Note: designed to be called internally by PROcess_SBNfit() function
     *
     * Arguments: 
     * 		branch: pointer to branch variable, each corresponding to one subchannel 
     * 		eventweight_map: a map between systematic string to list of variation weights
     * 		subchannel_index: index associated with current branch/subchannel
     *		syst_vector: list of SystStruct TO BE UPDATED, each stores all variation spectra of one systematic
     *		syst_additional_weight: additional weight applied to systematic variation
     */
    void process_sbnfit_event(const PROconfig &inconfig, const std::shared_ptr<BranchVariable>& branch, const std::map<std::string, std::vector<eweight_type>>& eventweight_map, int subchannel_index, std::vector<SystStruct>& syst_vector, const std::vector<double>& syst_additional_weight);


    
};
#endif
