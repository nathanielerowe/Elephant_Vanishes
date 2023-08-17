#include "PROconfig.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
#include "PROpeller.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/NumericalDiff>

using namespace PROfit;

log_level_t GLOBAL_LEVEL = LOG_DEBUG;

class ChiTest
{
    private:
        int n;
    public:
        ChiTest(int n_) : n(n_) {}
        double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad)
        {
            double fx = 0.0;
            for(int i = 0; i < n; i += 2)
            {
                double t1 = 1.0 - x[i];
                double t2 = 10 * (x[i + 1] - x[i] * x[i]);
                grad[i + 1] = 20 * t2;
                grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
                fx += t1 * t1 + t2 * t2;
            }
            return fx;
        }
};


int main(int argc, char* argv[])
{

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 100;

    //doubles
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);



    PROconfig myConf(xmlname);
   /* Ibrahims Test Space
   PROpeller myprop;
    std::vector<SystStruct> systs;
    PROcess_CAFs(myConf, systs, myprop);
    std::cout<< "starting loop!!" << std::endl;
    
    for(int i=0; i<=myprop.truth.size(); ++i){
        std::cout<< myprop.truth[i] << " "<< myprop.reco[i] <<" "<< myprop.pdg[i] << std::endl;

    }
 */



    std::vector<SystStruct> systsstructs;
    PROcess_CAFana(myConf, systsstructs);
    PROsyst systs(systsstructs);
    PROspec p05 = systs.GetSplineShiftedSpectrum(systsstructs[0].CV(), "GENIEReWeight_ICARUS_v1_multisigma_MaCCRES", 0.5);
    p05.Print();
    std::cout << systs.GrabMatrix("piplus_Flux") << std::endl;

    //PROspec mySpec(myConf);
    //TH1D hmm = mySpec.toTH1D(myConf);


    return 0;
    LBFGSpp::LBFGSBParam<double> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    LBFGSpp::LBFGSBSolver<double> solver(param); 

    int n=78;
    ChiTest fun(n);

    // Bounds
    Eigen::VectorXd lb = Eigen::VectorXd::Constant(n, 0.0);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(n, std::numeric_limits<double>::infinity());

    // Initial guess
    Eigen::VectorXd x = Eigen::VectorXd::Constant(n, 2.0);


    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);


    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}

