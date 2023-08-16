#ifndef PROSPEC_H_
#define PROSPEC_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

// ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

// EIGEN
#include <Eigen/Dense>
#include <Eigen/SVD>

// PROfit
#include "PROconfig.h"

namespace PROfit{

    class PROspec {

	private:
            //Base
            Eigen::VectorXd spec;
            Eigen::VectorXd error_square;
            //Eigen::VectorXd bins;

        public:

            //Constructors
            PROspec() {}
            PROspec(const Eigen::VectorXd &in_spec, const Eigen::VectorXd &in_error) : spec(in_spec), error_square(in_error){}

	    /* Function: create PROspec of given size */
	    PROspec(long int num_bins);

            //PROspec(PROconfig const & configin); //Load in config file EMPTY hists
            //PROspec(std::string &xmlname); //Load directly from XML 


	    /* Function: given subchannel name/index, generate TH1D histogram of corresponding spectrum */
            TH1D toTH1D(const PROconfig& inconfig, int subchannel_index);
            TH1D toTH1D(const PROconfig& inconfig, const std::string& subchannel_fullname);


	    /* Function: save TH1Ds of all subchannels into a root file */
	    void toROOT(const PROconfig& inconfig, const std::string& output_name);


	    /* Function: fill given bin with provided weight 
 	     * Note: it does NOT check whether the given bin is out of range or not 
 	     */
	    void Fill(long int bin_index, double weight);

	    /* Function: zero out the spectrum and error, but keep the dimension */
	    void Zero();

            /* Function: Print out spec*/
	    void Print() const;


    };

}


#endif
