#ifndef PROPELLER_H_
#define PROPELLER_H_

#include "PROconfig.h"
#include "PROserial.h"
#include <Eigen/Eigen>
// STANDARD
#include <vector>
#include <chrono>
namespace PROfit{


    /*Class: The PROpeller, which moves the analysis forward. A class to keep all MC events for oscllation event-by-event.
    */
    class PROpeller {

        private:
            friend class boost::serialization::access;
            int nevents;

            // Serialization function for boost that will allow for save state of propeller
            template <class Archive>
                void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
                    ar & nevents;
                    ar & pcosth;
                    ar & trueLE;
                    ar & added_weights;
                    ar & bin_indices;
                    ar & model_rule;
                    ar & true_bin_indices;
                    ar & hist;
                    ar & histLE;
                }

        public:

            //Empty Constructor
            PROpeller(){
                nevents = -1;
                pmom.clear();
                pcosth.clear();
                trueLE.clear();
                added_weights.clear();
                bin_indices.clear();
                model_rule.clear();
                true_bin_indices.clear();
            };

            /*Function: Primary Constructor from raw std::vectors of MC values */ 
            PROpeller(const PROconfig &config, std::vector<float> &intruth, std::vector<float> &inpmom, std::vector<float> &inpcosth, std::vector<float> &inadded_weights, std::vector<int> &inbin_indices, std::vector<int> &inmodel_rule, std::vector<int> &intrue_bin_indices) : trueLE(intruth), pmom(inpmom), pcosth(inpcosth), added_weights(inadded_weights), bin_indices(inbin_indices), model_rule(inmodel_rule), true_bin_indices(intrue_bin_indices){
                nevents = trueLE.size();
                hist = Eigen::MatrixXf::Constant(config.m_num_truebins_total, config.m_num_bins_total, 0);
                for(size_t i = 0; i < bin_indices.size(); ++i)
                    hist(true_bin_indices[i], bin_indices[i]) += added_weights[i];
            };

            /* the Core MC is saved in these vectors.*/

            std::vector<float> trueLE;
            std::vector<float> added_weights;
            std::vector<int>   bin_indices;        /*Precalculated Bin index*/
            std::vector<int>   model_rule;
            std::vector<int>   true_bin_indices;
            std::vector<float> pmom;
            std::vector<float> pcosth;
            Eigen::MatrixXf    hist;
            Eigen::VectorXf    histLE;


            // boost serialize save to file
            void save(const std::string& filename) const {
                auto start = std::chrono::high_resolution_clock::now();
                std::ofstream ofs(filename);
                boost::archive::text_oarchive oa(ofs);
                oa << *this;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                log<LOG_INFO>(L"%1% || Serialization save of PROpeller into file  %2% took %3% seconds") % __func__ % filename.c_str() % elapsed.count();
            }

            // Load from file
            void load(const std::string& filename) {
                auto start = std::chrono::high_resolution_clock::now();
                std::ifstream ifs(filename);
                boost::archive::text_iarchive ia(ifs);
                ia >> *this;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                log<LOG_INFO>(L"%1% || Serialization load of PROpeller from file  %2% took %3% seconds") % __func__ % filename.c_str() %elapsed.count();
            }


    };

}
#endif
