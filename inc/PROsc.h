#ifndef PROSC_H_
#define PROSC_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include <cstdlib>

// TINYXML2
#include "tinyxml2.h"

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"

namespace PROfit{


    class PROsc {
        private:

        public:

            PROsc(){};

            double Pmue(double dmsq, double sinsq2thmue, double enu, double baseline)               

                {
                double sinterm = sin(1.27*dmsq*(baseline/enu));
	        double prob    = sinsq2thmue*sinterm*sinterm;

                if(prob<0.0 || prob >1.0){
                    log<LOG_ERROR>(L"Your probability is outside the bounds of math.");
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                    }
                return prob;
                }

            double Pmumu(double dmsq, double sinsq2thmumu, double enu, double baseline)

                {
                double sinterm = sin(1.27*dmsq*(baseline/enu));
                double prob    = 1.0 - (sinsq2thmumu*sinterm*sinterm);

                if(prob<0.0 || prob >1.0){
                    log<LOG_ERROR>(L"Your probability is outside the bounds of math.");
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                    }

	        return prob;
                }
        };

}
#endif
