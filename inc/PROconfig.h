#ifndef PROCONFIG_H_
#define PROCONFIG_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

// TINYXML2
#include "tinyxml2.h"

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"

namespace PROfit{

    class PROconfig {
        protected:
        public:


            PROconfig() {}; //always have an empty?
            PROconfig(const std::string &xml);

            int LoadFromXML(const std::string & filename);

            std::vector<std::string> fullnames;

            int m_num_detectors;
            int m_num_channels;
            int m_num_modes;

            std::vector<int> m_num_bins;

            bool m_has_oscillation_patterns;


            //vectors of length num_channels
            std::vector<int> m_num_subchannels; 

            //the xml names are the way we track which channels and subchannels we want to use later
            std::vector<std::string> m_mode_names; 			
            std::vector<std::string> m_mode_plotnames; 			
           
            std::vector<std::string> m_detector_names; 		
            std::vector<std::string> m_detector_plotnames; 		
            
            std::vector<std::string> m_channel_names; 		
            std::vector<std::string> m_channel_plotnames; 		
            std::vector<std::string> m_channel_units; 		
           	std::vector<int> m_channel_num_bins;
            std::vector<std::vector<double> > m_channel_bin_edges;
        	std::vector<std::vector<double> > m_channel_bin_widths;


            std::vector<std::vector<std::string >> m_subchannel_names; 
            std::vector<std::vector<std::string >> m_subchannel_plotnames; 
            std::vector<std::vector<int >> m_subchannel_datas; 
            std::vector<std::vector<int> > m_subchannel_osc_patterns; 




            int m_num_bins_detector_block;
            int m_num_bins_mode_block;
            int m_num_bins_total;

            int m_num_bins_detector_block_collapsed;
            int m_num_bins_mode_block_collapsed;
            int m_num_bins_total_collapsed;

            std::string m_xmlname;	

            Eigen::MatrixXd collapsingVector;



    };




}
#endif
