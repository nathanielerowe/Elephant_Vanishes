#ifndef PROCNP_H_
#define PROCNP_H_

// STANDARD
#include <string>
#include <vector>

#include <Eigen/Eigen>

// OUR INCLUDES
#include "PROconfig.h"
#include "PROsyst.h"
#include "PROpeller.h"
#include "PROsc.h"
#include "PROcess.h"
#include "PROmetric.h"

namespace PROfit{

    class PROCNP : public PROmetric
    {
        // TODO: How much of this should be in PROmetric instead?
        private:
            std::string model_tag;

            const PROconfig *config;
            const PROpeller *peller;
            const PROsyst *syst; 
            const PROsc *osc;
            const PROspec data;
            int nparams;
            int nsyst;
            EvalStrategy strat;
            std::vector<float> physics_param_fixed;
                        //Do we want to fix any param?
            int fixed_index;
            float fixed_val;

            //Save last values for gradient calculation
            Eigen::VectorXd last_param;
            float last_value;

            bool correlated_systematics;
            Eigen::MatrixXd prior_covariance;
            Eigen::MatrixXd collapsed_stat_covariance;

        public:

            /*Function: Constructor bringing all objects together*/
            PROCNP(const std::string tag, const PROconfig *conin, const PROpeller *pin, const PROsyst *systin, const PROsc *oscin, const PROspec &datain, int nparams, int nsyst, EvalStrategy strat = EventByEvent, std::vector<float> physics_param_fixed = std::vector<float>());

            /*Function: operator() is what is passed to minimizer.*/
            virtual double operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient);
            virtual double operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient, bool nograd);

            PROmetric *Clone() const {
                return new PROCNP(*this);
            }

            virtual void reset() {
                physics_param_fixed.clear();
                last_value = 0;
                last_param = Eigen::VectorXd::Constant(last_param.size(), 0);
            }

            virtual void set_physics_param_fixed(const std::vector<float> &physics_param) {
                physics_param_fixed = physics_param;
            }

            virtual void override_systs(const PROsyst &new_syst) {
                nparams -= nsyst;
                syst = &new_syst;
                nsyst = syst->GetNSplines();
                nparams += nsyst;
            }

            float Pull(const Eigen::VectorXd &systs);

            void fixSpline(int fix, double valin);

    };


}
#endif
