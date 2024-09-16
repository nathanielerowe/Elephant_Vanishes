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

int main(int argc, char* argv[])
{

    CLI::App app{"PROplot, plotting spectra"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    std::string output_filetag = "test"; 

    //doubles
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");
    app.add_option("-o,--output", output_filetag, "Output Filename.");

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
    PROsc osc;

    //GetSpectrum
    PROspec spec = systsstructs.back().CV();
   
    //Plot
    spec.plotSpectrum(myConf,"testplot");
    
    return 0;
}

