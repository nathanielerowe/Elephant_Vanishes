#ifndef PROSPEC_H_
#define PROSPEC_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>

// ROOT
#include "TH1D.h"

// EIGEN
#include <Eigen/Dense>
#include <Eigen/SVD>

// PROfit
#include "PROconfig.h"

namespace PROfit{

    class PROspec {

        public:

            //Constructors
            PROspec() {};
            PROspec(PROconfig const & configin); //Load in config file EMPTY hists
            //PROspec(std::string &xmlname); //Load directly from XML 

            PROconfig configName;

            //Base
            Eigen::VectorXd spec;
            Eigen::VectorXd error;


            //Useful Plotting routines
            TH1D toTH1D();

    };

}


#endif
