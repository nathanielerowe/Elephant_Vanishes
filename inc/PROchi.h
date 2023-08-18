#ifndef PROCHI_H_
#define PROCHI_H_

// STANDARD

// OUR INCLUDES
#include "PROconfig.h"
#include "PROsyst.h"
#include "PROpeller.h"

namespace PROfit{



    class PROchi
    {
        private:
            std::string model_tag;
            const PROconfig *config;
            const PROpeller *peller;
            const PROsyst *syst; 
        public:
            PROchi(const std::string tag, const PROconfig *conin, const PROpeller *pin, const PROsyst *systin) : model_tag(tag), config(conin), peller(pin), syst(systin) {}
            double operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient)
            {
                

                return 0;
            }
    };


}
#endif
