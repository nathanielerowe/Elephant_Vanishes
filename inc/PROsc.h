#ifndef PROSC_H_
#define PROSC_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <cmath>

// TINYXML2
#include "tinyxml2.h"

//PROfit
#include "PROlog.h"

namespace PROfit{


    /* 
     * Class: Oscillation Class for oscillating event-by-event MC and prodice PROspec
     * Note:
     *  Currently 3+1 SBL approximation is only model and hardcoded. So a very simple class.
     * Todo: 
     *  Add Nuquids interface if desired
     *  Add interfacte for defining own model
     */

    class PROsc {
        private:
        public:

            PROsc(){};

            /* Function: 3+1 numu->nue apperance prob in SBL approx */
            float Pmue(float dmsq, float sinsq2thmue, float enu, float baseline) const{
                if(sinsq2thmue > 1) sinsq2thmue = 1;
                if(sinsq2thmue < 0) sinsq2thmue = 0;

                float sinterm = std::sin(1.27*dmsq*(baseline/enu));
                float prob    = sinsq2thmue*sinterm*sinterm;

                if(prob<0.0 || prob >1.0){
                    //Do we need this when its hardcoded simple cmath above?
                    log<LOG_ERROR>(L"%1% || Your probability %2% is outside the bounds of math.") % __func__ % prob;
                    log<LOG_ERROR>(L"%1% || Terminating.") % __func__;
                    exit(EXIT_FAILURE);
                }
                return prob;
            }

            /* Function: 3+1 numu->numue disapperance prob in SBL approx */
            float Pmumu(float dmsq, float sinsq2thmumu, float enu, float baseline) const{
                if(sinsq2thmumu > 1) sinsq2thmumu = 1;
                if(sinsq2thmumu < 0) sinsq2thmumu = 0;

                float sinterm = std::sin(1.27*dmsq*(baseline/enu));
                float prob    = 1.0 - (sinsq2thmumu*sinterm*sinterm);

                if(prob<0.0 || prob >1.0){
                    log<LOG_ERROR>(L"%1% || Your probability %2% is outside the bounds of math. dmsq = %3%, sinsq2thmumu = %4%, enu = %5%, baseline = %6%") % __func__ % prob % dmsq % sinsq2thmumu % enu % baseline;
                    log<LOG_ERROR>(L"%1% || Terminating.") % __func__;
                    exit(EXIT_FAILURE);
                }

                return prob;
            }
    };

}
#endif
