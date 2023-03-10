#include "PROspec.h"
using namespace PROfit;


PROspec::PROspec(PROconfig const & inConfig){
    configName = inConfig; 

    spec = Eigen::VectorXd::Zero(configName.num_modes);
    std::cout<<spec[0]<<std::endl;

}


TH1D PROspec::toTH1D(){



    TH1D hSpec("","",10,0,1); 


    return hSpec;

}
