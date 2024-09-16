#include "PROconfig.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROcess.h"

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

    //Initilize configuration from the XML;
    PROconfig myConf(xmlname);

    //Inititilize PROpeller to keep MC
    PROpeller myprop;
    
    //Initilize objects for systematics storage
    std::vector<SystStruct> systsstructs;

    //Process the CAF files to grab and fill all SystStructs and PROpeller
    PROcess_CAFs(myConf, systsstructs, myprop);

    //Build a PROsyst to sort and analyze all systematics
    PROsyst systs(systsstructs);

    //Define the model (currently 3+1 SBL)
    PROsc osc(myprop);

    //Setup minimization parameetrs
    LBFGSpp::LBFGSBParam<double> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    LBFGSpp::LBFGSBSolver<double> solver(param); 

    Eigen::VectorXd data = systsstructs.back().CV().Spec();

    int nparams = 3;

    //Build chi^2 object
    PROchi chi("3plus1",&myConf,&myprop,&systs,&osc, systsstructs.back().CV(), nparams, 1);

    // Bounds
    Eigen::VectorXd lb(3);
    lb << 0.01 , 0  , -3;
    Eigen::VectorXd ub(3);
    ub <<  100, 1, 3 ;

    // Initial guess
    Eigen::VectorXd x(3);
    x << 3, 0.5, 0.0;

    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(chi, x, fx, lb, ub);
    
    log<LOG_DEBUG>(L"%1% || FINISHED MINIMIZING: NITERATIONS %2%  and MINIMUM PARAMS  %3% %4% %5%" ) % __func__ % niter % x[0] % x[1] % x[2];

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}

