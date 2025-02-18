#include "PROconfig.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
#include "PROsc.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROcess.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <Eigen/Eigen>
#include <random>

#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"


using namespace PROfit;

log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[])
{
    gStyle->SetOptStat(0);

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 100;
    bool oscillate = false;

    //floats
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_flag("-O,--Oscillate", oscillate, "Fit with oscillations.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);

    //Initilize configuration from the XML;
    PROconfig myConf(xmlname);

    //Inititilize PROpeller to keep MC
    PROpeller myprop;
    
    //Initilize objects for systematics storage
    std::vector<SystStruct> systsstructs;

    //Process the CAF files to grab and fill all SystStructs and PROpeller
    PROcess_CAFAna(myConf, systsstructs, myprop);

    //Build a PROsyst to sort and analyze all systematics
    PROsyst systs(systsstructs);

    //Define the model (currently 3+1 SBL)
    PROsc osc(myprop);

    //Setup minimization parameetrsM
    LBFGSpp::LBFGSBParam<float> param;  
    //LBFGSpp::LBFGSParam<float> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    param.max_linesearch = 50;
    LBFGSpp::LBFGSBSolver<float> solver(param); 
    //LBFGSpp::LBFGSSolver<float> solver(param); 

    Eigen::VectorXf data = systsstructs.back().CV().Spec();

    //int nparams = 2 + systs.GetNSplines();
    int nparams = 2*oscillate + systs.GetNSplines();

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<float> d;
    std::vector<float> throws;
    for(size_t i = 0; i < systs.GetNSplines(); i++)
        throws.push_back(d(gen));
    PROspec newSpec = PROspec::PoissonVariation(systs.GetSplineShiftedSpectrum(myConf, myprop, throws));

    std::cout << "CV: " << std::endl;
    systsstructs.back().CV().Print();
    std::cout << "Throw: " << std::endl;
    newSpec.Print();

    //Build chi^2 object
    PROchi chi("3plus1",&myConf,&myprop,&systs,oscillate ? &osc : NULL, newSpec, nparams, systs.GetNSplines());

    // Bounds
    Eigen::VectorXf lb = Eigen::VectorXf::Constant(nparams, -3.0);
    if(oscillate){
        lb(0) = 0.01; lb(1) = 0;
    }
    Eigen::VectorXf ub = Eigen::VectorXf::Constant(nparams, 3.0);
    if(oscillate){
        ub(0) = 100; ub(1) = 1;
    }
    // Initial guess
    Eigen::VectorXf x = Eigen::VectorXf::Constant(nparams, 0.2);

    // x will be overwritten to be the best point found
    float fx;
    int niter = -1;
    try {
        niter = solver.minimize(chi, x, fx, lb, ub);
    } catch(std::runtime_error &except) {
        log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
    }


    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    std::cout << "Throw: \n" << Eigen::VectorXf::Map(throws.data(), throws.size()).transpose() << std::endl;
    return 0;
}

