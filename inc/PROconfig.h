#ifndef PROCONFIG_H_
#define PROCONFIG_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>

// TINYXML2
#include "tinyxml2.h"

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

//PROfit
#include "PROlog.h"

//ROOT
#include "TTreeFormula.h"

namespace PROfit{

    struct BranchVariable{
        std::string name;
        std::string type;
        std::string associated_hist;
        std::string associated_systematic;
        bool central_value;

        std::shared_ptr<TTreeFormula> branch_formula=NULL;
        std::shared_ptr<TTreeFormula> branch_true_value_formula=NULL;
        std::shared_ptr<TTreeFormula> branch_true_L_formula=NULL;

        bool oscillate;
        std::string true_param_name;
        std::string true_L_name;

        float value_f;
        float true_value_f;
        float true_L_f;

        double value_d;
        double true_value_d;
        double true_L_d;

        BranchVariable(std::string n, std::string t, std::string a) : name(n), type(t), associated_hist(a) {oscillate=false; associated_systematic=""; central_value =false;};
        BranchVariable(std::string n, std::string t, std::string a_hist, std::string a_syst, bool cv) : name(n), type(t), associated_hist(a_hist), associated_systematic(a_syst) { oscillate=false; central_value=cv;};
        virtual void* GetValue(){};
        virtual void* GetTrueValue(){};
        virtual void* GetTrueL(){};

        std::shared_ptr<TTreeFormula> GetFormula(){
            return branch_formula;
        }

        int SetOscillate(bool inbool){ oscillate = inbool;};
        bool GetOscillate(){ return oscillate;};

    };


    struct BranchVariable_d: public BranchVariable{
        BranchVariable_d(std::string n, std::string t, std::string a) : BranchVariable(n,t,a) {value_d=0;true_value_d=0; true_L_d = 0;};
        BranchVariable_d(std::string n, std::string t, std::string a_hist, std::string a_syst, bool cv) : BranchVariable(n,t,a_hist, a_syst, cv) {value_d=0;true_value_d=0; true_L_d = 0;};
        void* GetValue(){ 
            if(branch_formula == NULL) return &value_d;
            else{   
                branch_formula->GetNdata();
                value_d = (double)branch_formula->EvalInstance();
                return &value_d;
            }
        }

        void* GetTrueL(){
            if(branch_true_L_formula == NULL) return &true_L_d;
            else{
                branch_true_L_formula->GetNdata();
                true_L_d = (double)branch_true_L_formula->EvalInstance();
                return &true_L_d;
            }
        }

        void* GetTrueValue(){ 
            if(branch_true_value_formula == NULL) return &true_value_d;
            else{
                branch_true_value_formula->GetNdata();
                true_value_d = (double)branch_true_value_formula->EvalInstance();
                return &true_value_d;
            }
        }
    };

    struct BranchVariable_f: public BranchVariable{
        BranchVariable_f(std::string n, std::string t, std::string a) : BranchVariable(n,t,a) {value_f=0;true_value_f=0; true_L_f = 0;};
        void* GetValue(){ 
            if(branch_formula == NULL) return &value_f;
            else{
                branch_formula->GetNdata();
                value_f = (float)branch_formula->EvalInstance();
                return &value_f;
            }
        }

        void* GetTrueValue(){ 
            if(branch_true_value_formula == NULL) return &true_value_f;
            else{
                branch_true_value_formula->GetNdata();
                true_value_f = (float)branch_true_value_formula->EvalInstance();
                return &true_value_f;
            }
        }

        void* GetTrueL(){
            if(branch_true_L_formula == NULL) return &true_L_f;
            else{
                branch_true_L_formula->GetNdata();
                true_L_f = (float)branch_true_L_formula->EvalInstance();
                return &true_L_f;
            }
        }
    };


    class PROconfig {
        protected:
        public:


            PROconfig() {}; //always have an empty?
            PROconfig(const std::string &xml);

            int LoadFromXML(const std::string & filename);

            std::vector<std::string> m_fullnames;

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

            bool m_write_out_variation;
            bool m_form_covariance;
            std::string m_write_out_tag;


            int m_num_bins_detector_block;
            int m_num_bins_mode_block;
            int m_num_bins_total;

            int m_num_bins_detector_block_collapsed;
            int m_num_bins_mode_block_collapsed;
            int m_num_bins_total_collapsed;

            std::string m_xmlname;	

            Eigen::MatrixXd collapsingVector;


            //This section entirely for montecarlo generation of a covariance matrix or PROspec 
            //For generating a covariance matrix from scratch, this contains the number of montecarlos (weights in weight vector) and their names.
            int m_num_mcgen_files;
            std::vector<std::string> m_mcgen_tree_name;	
            std::vector<std::string> m_mcgen_file_name;	
            std::vector<std::string> m_mcgen_additional_weight_name;
            std::vector<bool> m_mcgen_additional_weight_bool;
            std::vector<double> m_mcgen_additional_weight;
            std::vector<int> m_mcgen_maxevents;	
            std::vector<double> m_mcgen_pot;	
            std::vector<double> m_mcgen_scale;	
            std::vector<bool> m_mcgen_fake;
            std::map<std::string,std::vector<std::string>> m_mcgen_file_friend_map;
            std::map<std::string,std::vector<std::string>> m_mcgen_file_friend_treename_map;


            //specific bits for covariancegeneration
            std::vector<std::string> m_mcgen_eventweight_branch_names;
            std::vector<std::string> m_mcgen_weightmaps_formulas;
            std::vector<std::string> m_mcgen_weightmaps_uses;
            std::vector<std::string> m_mcgen_weightmaps_patterns;
            std::vector<std::string> m_mcgen_weightmaps_mode;
            std::map<std::string,bool> m_mcgen_variation_allowlist;
            std::map<std::string,bool> m_mcgen_variation_denylist;
            std::map<std::string, std::vector<std::string>> m_mcgen_shapeonly_listmap; //a map of shape-only systematic and corresponding subchannels

            std::vector<std::vector<BranchVariable*>> m_branch_variables;



            //FIX skepic
            std::vector<std::string> systematic_name;



    };




}
#endif
