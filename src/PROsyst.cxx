#include "PROsyst.h"

namespace PROfit {

void PROsyst::AddMatrix(const SystStruct& syst){

    std::string sysname = syst.GetSysName();

    //generate matrix only if it's not already in the map 
    if(covmat_map.find(sysname) == covmat_map.end()){

	//generate fractional covariance matrix 	
        Eigen::MatrixXd frac_matrix = syst.GenerateCovarMatrix();
 	if(!SystStruct::isPositiveSemiDefinite_WithTolerance(frac_matrix)){
	    LOG<LOG_ERROR>(L"%1% || Fractional Covariance Matrix is not positive semi-definite!") % __func__;
	    log<LOG_ERROR>(L"Terminating.");
            exit(EXIT_FAILURE);
	}
	covmat_map[sysname] = frac_matrix;

	//genrate correlation matrix 
    }

    return;
}
};

