#include "PROcess.h"
#include "PROsc.h"
#include "PROspec.h"
#include "PROcreate.h"
#include "PROtocall.h"
#include "PROsyst.h"
#include "PROchi.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/SRWeightPSet.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"


#include <memory>
#include <string>

using namespace PROfit;
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[])
{

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.INPUT.xml"; 

    app.add_option("-x,--xml", xmlname, "XML file to read configuration from.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);

    //Initilize configuration from the XML;
    PROconfig config(xmlname);


    //PROfile development
    if(false){

        //Inititilize PROpeller to keep MC
        PROpeller prop;

        //Initilize objects for systematics storage
        std::vector<SystStruct> systsstructs;

        //Process the CAF files to grab and fill all SystStructs and PROpeller
        PROcess_CAFAna(config, systsstructs, prop);

        //Build a PROsyst to sort and analyze all systematics
        PROsyst systs(systsstructs);

        //Define the model (currently 3+1 SBL)
        PROsc osc(prop);

        //
        PROspec data = systsstructs.back().CV();

        std::vector<std::string> spline_names;
        int cnt=0;
        for(auto&s:systsstructs){
            //index,systname,mode
            log<LOG_INFO>(L"%1% || Starting Fixed fit %2% %3% %4% %5%") % __func__ % s.systname.c_str() % s.index % s.mode.c_str() % cnt ;
            if(s.mode=="multisigma"){
                spline_names.push_back(s.systname); 
                cnt++;
            }
   
        }


        LBFGSpp::LBFGSBParam<float> param;
        param.epsilon = 1e-6;
        param.max_iterations = 100;
        param.max_linesearch = 50;
        param.delta = 1e-6;

        LBFGSpp::LBFGSBSolver<float> solver(param);
        int nparams = systs.GetNSplines();
        std::vector<float> physics_params; 


        TCanvas *c =  new TCanvas("Profile", "Profile" , 400*4, 400*7);
        c->Divide(4,7);


        std::vector<std::unique_ptr<TGraph>> graphs; 

        //hack
        std::vector<float> priorX;
        std::vector<float> priorY;

       for(int i=0; i<=20;i++){
           float which_value = -2.0+0.2*i;
           priorX.push_back(which_value);
           priorY.push_back(which_value*which_value);

       }
       std::unique_ptr<TGraph> gprior = std::make_unique<TGraph>(priorX.size(), priorX.data(), priorY.data());



        for(int w=0; w<nparams;w++){
            int which_spline = w;


            std::vector<float> knob_vals;
            std::vector<float> knob_chis;

            for(int i=0; i<=20;i++){
                
                Eigen::VectorXf lb = Eigen::VectorXf::Constant(nparams, -3.0);
                Eigen::VectorXf ub = Eigen::VectorXf::Constant(nparams, 3.0);
                Eigen::VectorXf x = Eigen::VectorXf::Constant(nparams, 0.0);
                Eigen::VectorXf grad = Eigen::VectorXf::Constant(nparams, 0.0);
                Eigen::VectorXf bestx = Eigen::VectorXf::Constant(nparams, 0.0);


                float which_value = -2.0+0.2*i;
                float fx;
                knob_vals.push_back(which_value);

                lb[which_spline] = which_value;
                ub[which_spline] = which_value;
                x[which_spline] = which_value;


                PROchi chi("3plus1", config, prop, &systs, osc, data, PROchi::BinnedChi2, physics_params);
                chi.fixSpline(which_spline,which_value);

                log<LOG_INFO>(L"%1% || Starting Fixed fit ") % __func__  ;
                try {
                    x = Eigen::VectorXf::Constant(nparams, 0.012);
                    solver.minimize(chi, x, fx, lb, ub);
                } catch(std::runtime_error &except) {
                    log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
                }

                std::string spec_string = "";
                for(auto &f : x) spec_string+=" "+std::to_string(f); 
                log<LOG_INFO>(L"%1% || Fixed value of %2% for spline %3% was post  : %4% ") % __func__ % which_spline % which_value % fx;
                log<LOG_INFO>(L"%1% || BF splines @ %2%") % __func__ %  spec_string.c_str();

                knob_chis.push_back(fx);
            }            

            log<LOG_INFO>(L"%1% || Knob Values: %2%") % __func__ %  knob_vals;
            log<LOG_INFO>(L"%1% || Knob Chis: %2%") % __func__ %  knob_chis;

            c->cd(w+1);
            std::unique_ptr<TGraph> g = std::make_unique<TGraph>(knob_vals.size(), knob_vals.data(), knob_chis.data());
            std::string tit = spline_names[which_spline]+ ";#sigma Shift; #Chi^{2}";
            g->SetTitle(tit.c_str());
            graphs.push_back(std::move(g));
            graphs.back()->Draw("AL");
            graphs.back()->SetLineWidth(2);
            
            gprior->Draw("L same");
            gprior->SetLineStyle(2);
            gprior->SetLineWidth(1);
        
        }

        c->SaveAs("PROfile.pdf","pdf");

        delete c;


    }
    if(true){


        //Inititilize PROpeller to keep MC
        PROpeller prop;

        //Initilize objects for systematics storage
        std::vector<SystStruct> systsstructs;

        //Process the CAF files to grab and fill all SystStructs and PROpeller
        PROcess_CAFAna(config, systsstructs, prop);

        //Build a PROsyst to sort and analyze all systematics
        PROsyst systs(systsstructs);

        //Define the model (currently 3+1 SBL)
        PROsc osc(prop);

        PROspec data = systsstructs.back().CV();

        //CovarCheck
    //std::cout<<systs.fractional_covariance<<std::endl;

    //SplineCheck
    //PROspec test = systs.GetSplineShiftedSpectrum(config, prop, "GENIEReWeight_SBN_v1_multisigma_MaCCRES" , -2.0);
    //data.Print();
    //test.Print();


        std::vector<float> shift = {0.569111, 0.121685, 0.687138, 0.341982, 0.553843, 0.697239, 0.192107, 0.456814, 0.116422, 0.0806362, 0.501815, 0.781406, 0.321165, 0.525954, 0.0396827, 0.373735, 0.54393, 0.106545, 0.43473, 0.479401, 0.255056, 0.152222, 0.0812247, 0.153592, 0.256307, 0.0998781, 0.403397, 0.864244};
        std::vector<float> param = {0,0};
        std::vector<float> nulpram;

        //PROspec CVe = FillRecoSpectra(config, prop, systs, &osc, nulpram, nulpram, false);
        //CVe.plotSpectrum(config,"CVe");
        //CVe.Print();
        //PROspec CVb = FillRecoSpectra(config, prop, systs, &osc, nulpram, nulpram, true);
        //CVb.plotSpectrum(config,"CVb");
        //CVb.Print();

        //PROspec SYe = FillRecoSpectra(config, prop, systs, &osc, shift, nulpram, false);
        //SYe.plotSpectrum(config,"SYe");
        //SYe.Print();
        //PROspec SYb = FillRecoSpectra(config, prop, systs, &osc, shift, nulpram, true);
        //SYb.plotSpectrum(config,"SYb");
        //SYb.Print();

        //PROspec OSe = FillRecoSpectra(config, prop, systs, &osc, nulpram, param, false);
        //OSe.plotSpectrum(config,"OSe");
        //OSe.Print();
        //PROspec OSb = FillRecoSpectra(config, prop, systs, &osc, nulpram, param, true);
        //OSb.plotSpectrum(config,"OSb");
        //OSb.Print();

        //PROspec FUe = FillRecoSpectra(config, prop, systs, &osc, shift, param, false);
        //FUe.plotSpectrum(config,"FUe");
        //FUe.Print();
        //PROspec FUb = FillRecoSpectra(config, prop, systs, &osc, shift, param, true);
        //FUb.plotSpectrum(config,"FUb");
        //FUb.Print();
    }

    //test matrix generation 
    if(false){
        //fill in sysetamtic variations 
        std::vector<SystStruct> syst_vector;
        PROcess_SBNfit(config, syst_vector);
        std::cout << syst_vector.size() << std::endl;

        //generate covariance matrix 
        Eigen::MatrixXf fractional_matrix = Eigen::MatrixXf::Zero(15, 15);
        for(auto& s : syst_vector)
            fractional_matrix += PROsyst::GenerateFracCovarMatrix(s);
        std::cout << " Formed fractional covariance matrix by PROfit: " << std::endl; 
        std::cout << fractional_matrix << std::endl;

        //check is matrix is psd
        bool res1 =  PROsyst::isPositiveSemiDefinite(fractional_matrix), res2 = PROsyst::isPositiveSemiDefinite_WithTolerance(fractional_matrix);
        std::cout << "Matrix is positive semidefinite? " << res1 << " " << res2 << std::endl;


        Eigen::MatrixXf matrix1(3,3), matrix2(2,2), matrix3(2,2);
        matrix1 << 3,2,1,2,3,1,1,2,3;
        matrix2 << 1, 2, 2, 3;
        matrix3 << 1, 2, 2, 4;
        res1 =  PROsyst::isPositiveSemiDefinite(matrix1), res2 = PROsyst::isPositiveSemiDefinite_WithTolerance(matrix1);
        std::cout << "Matrix1 is positive semidefinite? " << res1 << " " << res2 <<  std::endl;
        std::cout << matrix1 << std::endl;
        res1 =  PROsyst::isPositiveSemiDefinite(matrix2), res2 = PROsyst::isPositiveSemiDefinite_WithTolerance(matrix2);
        std::cout << "Matrix2 is positive semidefinite? " << res1 << " " << res2 <<  std::endl;
        std::cout << matrix2 << std::endl;
        res1 =  PROsyst::isPositiveSemiDefinite(matrix3), res2 = PROsyst::isPositiveSemiDefinite_WithTolerance(matrix3);
        std::cout << "Matrix3 is positive semidefinite? " << res1 << " " << res2 <<  std::endl;
        std::cout << matrix3 << std::endl;


        //matrix collapsing 
        Eigen::MatrixXf collapsed_fraction = CollapseMatrix(config, fractional_matrix);
        std::cout << "Collapsed Matrix by PROfit: " << std::endl;
        std::cout << collapsed_fraction << std::endl;


        Eigen::MatrixXf correlation_matrix = PROsyst::GenerateCorrMatrix(fractional_matrix);
        std::cout << "Correlation Matrix" << std::endl;
        std::cout << correlation_matrix << std::endl;


    }
    return 0;
}
