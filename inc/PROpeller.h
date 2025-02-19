#ifndef PROPELLER_H_
#define PROPELLER_H_

#include "PROconfig.h"
#include <Eigen/Eigen>

// STANDARD
#include <vector>

namespace PROfit{


    /*Class: The PROpeller, which moves the analysis forward. A class to keep all MC events for oscllation event-by-event.
    */
    class PROpeller {

        private:
            int nevents;

        public:

            //Empty Constructor
            PROpeller(){
                nevents = -1;
                trueLE.clear();
                added_weights.clear();
                bin_indices.clear();
                model_rule.clear();
                true_bin_indices.clear();
            }

            /*Function: Primary Constructor from raw std::vectors of MC values */ 
            PROpeller(const PROconfig &config, std::vector<float> &intruth, std::vector<float> &inadded_weights, std::vector<int> &inbin_indices, std::vector<int> &inmodel_rule, std::vector<int> &intrue_bin_indices) : trueLE(intruth), added_weights(inadded_weights), bin_indices(inbin_indices), model_rule(inmodel_rule), true_bin_indices(intrue_bin_indices){
                nevents = trueLE.size();
                hist = Eigen::MatrixXf::Constant(config.m_num_truebins_total, config.m_num_bins_total, 0);
                for(size_t i = 0; i < bin_indices.size(); ++i)
                    hist(true_bin_indices[i], bin_indices[i]) += added_weights[i];
	    }

            /* the Core MC is saved in these vectors.*/
            std::vector<float> trueLE;
            std::vector<float> added_weights;
            std::vector<int>   bin_indices;        /*Precalculated Bin index*/
            std::vector<int>   model_rule;
            std::vector<int>   true_bin_indices;
            Eigen::MatrixXf    hist;
            Eigen::VectorXf    histLE;

    };

}
#endif
