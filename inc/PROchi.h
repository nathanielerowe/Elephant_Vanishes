#ifndef PROCHI_H_
#define PROCHI_H_

// STANDARD

// OUR INCLUDES
#include "PROconfig.h"
#include "PROsyst.h"
#include "PROpeller.h"
#include "PROsc.h"
#include "PROcess.h"

namespace PROfit{



    class PROchi
    {
        private:
            std::string model_tag;
            const PROconfig *config;
            const PROpeller *peller;
            const PROsyst *syst; 
            const PROsc *osc;
            const PROspec data;

            Eigen::VectorXd last_param;
            float last_value;
        public:
            PROchi(const std::string tag, const PROconfig *conin, const PROpeller *pin, const PROsyst *systin, const PROsc *oscin, const PROspec &datain) : model_tag(tag), config(conin), peller(pin), syst(systin), osc(oscin), data(datain) {last_value = 0.0; last_param = Eigen::VectorXd::Zero(config->m_num_bins_total); }
            float operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient)
            {

                //std::vector<float> shifts = param(Eigen::SeqN(2,1)).array();
                //std::vector<float> fitparams = param(Eigen::SeqN(0,2)).array();
 
                // Get Spectra from FillRecoSpectra
                Eigen::VectorXd subvector1 = param.segment(0, 2);
                std::vector<float> fitparams(subvector1.data(), subvector1.data() + subvector1.size());
                Eigen::VectorXd subvector2 = param.segment(2,1);
                std::vector<float> shifts(subvector2.data(), subvector2.data() + subvector2.size());

                PROspec result = FillRecoSpectra(*config, *peller, *syst, *osc, shifts, fitparams);

                // Calcuate Full Covariance matrix
                Eigen::MatrixXd diag = result.Spec().array().matrix().asDiagonal(); 
                Eigen::MatrixXd full_covariance =  diag*(syst->fractional_covariance)*diag;
                
                // Collapse Covariance and Spectra 
                Eigen::MatrixXd collapsed_full_covariance =  CollapseMatrix(*config,full_covariance);  

                // Invert Collaped Matrix Matrix 
                Eigen::MatrixXd inverted_collapsed_full_covariance = collapsed_full_covariance.inverse();

                // Calculate Chi^2  value
                Eigen::VectorXd delta  = result.Spec() - data.Spec(); 
                float value = (delta)*inverted_collapsed_full_covariance*(delta.transpose())  + subvector2.array().square().sum();

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
