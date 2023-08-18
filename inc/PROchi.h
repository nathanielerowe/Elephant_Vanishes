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
               
                Eigen::VectorXd diff = param-last_param;
                float value = param[0]+gradient[0];

                //gradient = (last_value-value)/diff;

                last_param = param;
                last_value = value;
                return value;
            }
    };


}
#endif
