#include "PROspec.h"
using namespace PROfit;


PROspec::PROspec(PROconfig const & inConfig){

    spec = Eigen::VectorXd::Zero(inConfig.m_num_bins_total);
    error = Eigen::VectorXd::Zero(inConfig.m_num_bins_total);
    bins = Eigen::VectorXd::Zero(inConfig.m_num_bins_total+1);
    
    std::cout<<spec[0]<<std::endl;

}


TH1D PROspec::toTH1D(PROconfig const & inConfig){


    TH1D hSpec("","",10,0,1); 

    return hSpec;

}
