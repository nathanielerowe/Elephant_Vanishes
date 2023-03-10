#include "PROspec.h"
using namespace PROfit;


PROspec::PROspec(PROconfig const & inconfig){
    config = inconfig; 

    spec = Eigen::VectorXd::Zero(config.num_modes,config.num_modes);
    std::cout<<spec[0]<<std::endl;

}
