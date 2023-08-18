#ifndef PROCHI_H_
#define PROCHI_H_

// STANDARD

// OUR INCLUDES
#include "PROconfig.h"
#include "PROsyst.h"
#include "PROpeller.h"
#include "PROsc.h"

namespace PROfit{



    class PROchi
    {
        private:
            std::string model_tag;
            const PROconfig *config;
            const PROpeller *peller;
            const PROsyst *syst; 
            const PROsc *osc;

            Eigen::VectorXd last_param;
            float last_value;
        public:
            PROchi(const std::string tag, const PROconfig *conin, const PROpeller *pin, const PROsyst *systin, const PROsc *oscin) : model_tag(tag), config(conin), peller(pin), syst(systin), osc(oscin) {last_value = 0.0; last_param = Eigen::VectorXd::Zero(config->m_num_bins_total); }
            float operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient)
            {
              
                // Get Spectra from FillRecoSpectra

                // Calcuate Full Covariance matrix
                
                // Collapse Covariance and Spectra 
                
                // Invert Collaped Matrix Matrix 

                // Calculate Chi^2  value
                //Put actuall chi^2 here :)
                float value = param[0]+gradient[0];

                // Simple gradient here
                Eigen::VectorXd diff = param-last_param;
                gradient = (value-last_value)/diff.array();

                //Update last param
                last_param = param;
                last_value = value;
                return value;
            }
    };


}
#endif
