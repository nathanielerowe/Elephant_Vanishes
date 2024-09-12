#ifndef PROPELLER_H_
#define PROPELLER_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace PROfit{


    /*Class: The PROpeller, which moves the analysis forward. A class to keep all MC events for oscllation event-by-event.
    */
    class PROpeller {

        private:
            int nevents;

        public:

            //Emmpty Constructor
            PROpeller(){
                nevents = -1;
                truth.clear();
                reco.clear();
                baseline.clear();
                pdg.clear();
                added_weights.clear();
                bin_indices.clear();
                model_rule.clear();
            }

            /*Function: Primary Constructor from raw std::vectors of MC values */ 
            PROpeller(std::vector<float> &intruth, std::vector<float> &inreco, std::vector<float> &inbaseline, std::vector<int> &inpdg, std::vector<float> &inadded_weights, std::vector<int> &inbin_indices, std::vector<int> &inmodel_rule) : truth(intruth), reco(inreco), baseline(inbaseline), pdg(inpdg), added_weights(inadded_weights), bin_indices(inbin_indices), model_rule(inmodel_rule){
	 	nevents = truth.size();
	    }

            /* the Core MC is saved in these six vectors.*/
            std::vector<float> truth;
            std::vector<float> reco;
            std::vector<float> baseline;
            std::vector<int>   pdg;
            std::vector<float> added_weights;
            std::vector<int>   bin_indices; /*Precaulcated Bin index*/
            std::vector<int>   model_rule;

    };

}
#endif
