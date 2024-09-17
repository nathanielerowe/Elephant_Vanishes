#include "PROcess.h"
#include "PROsc.h"
#include "PROspec.h"
#include "PROcreate.h"
#include "PROtocall.h"
#include "PROsyst.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/SRWeightPSet.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include "TFile.h"

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

        std::vector<float> shift = {0.569111, 0.121685, 0.687138, 0.341982, 0.553843, 0.697239, 0.192107, 0.456814, 0.116422, 0.0806362, 0.501815, 0.781406, 0.321165, 0.525954, 0.0396827, 0.373735, 0.54393, 0.106545, 0.43473, 0.479401, 0.255056, 0.152222, 0.0812247, 0.153592, 0.256307, 0.0998781, 0.403397, 0.864244};
        std::vector<float> param = {0,0};
        std::vector<float> nulpram;

        PROspec CVe = FillRecoSpectra(config, prop, systs, &osc, nulpram, nulpram, false);
        CVe.plotSpectrum(config,"CVe");
        CVe.Print();
        PROspec CVb = FillRecoSpectra(config, prop, systs, &osc, nulpram, nulpram, true);
        CVb.plotSpectrum(config,"CVb");
        CVb.Print();

        PROspec SYe = FillRecoSpectra(config, prop, systs, &osc, shift, nulpram, false);
        SYe.plotSpectrum(config,"SYe");
        SYe.Print();
        PROspec SYb = FillRecoSpectra(config, prop, systs, &osc, shift, nulpram, true);
        SYb.plotSpectrum(config,"SYb");
        SYb.Print();

        PROspec OSe = FillRecoSpectra(config, prop, systs, &osc, nulpram, param, false);
        OSe.plotSpectrum(config,"OSe");
        OSe.Print();
        PROspec OSb = FillRecoSpectra(config, prop, systs, &osc, nulpram, param, true);
        OSb.plotSpectrum(config,"OSb");
        OSb.Print();

        PROspec FUe = FillRecoSpectra(config, prop, systs, &osc, shift, param, false);
        FUe.plotSpectrum(config,"FUe");
        FUe.Print();
        PROspec FUb = FillRecoSpectra(config, prop, systs, &osc, shift, param, true);
        FUb.plotSpectrum(config,"FUb");
        FUb.Print();
    }

    //test matrix generation 
    if(false){
        //fill in sysetamtic variations 
        std::vector<SystStruct> syst_vector;
        PROcess_SBNfit(config, syst_vector);
        std::cout << syst_vector.size() << std::endl;

        //generate covariance matrix 
        Eigen::MatrixXd fractional_matrix = Eigen::MatrixXd::Zero(15, 15);
        for(auto& s : syst_vector)
            fractional_matrix += PROsyst::GenerateFracCovarMatrix(s);
        std::cout << " Formed fractional covariance matrix by PROfit: " << std::endl; 
        std::cout << fractional_matrix << std::endl;

        //check is matrix is psd
        bool res1 =  PROsyst::isPositiveSemiDefinite(fractional_matrix), res2 = PROsyst::isPositiveSemiDefinite_WithTolerance(fractional_matrix);
        std::cout << "Matrix is positive semidefinite? " << res1 << " " << res2 << std::endl;


        Eigen::MatrixXd matrix1(3,3), matrix2(2,2), matrix3(2,2);
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
        Eigen::MatrixXd collapsed_fraction = CollapseMatrix(config, fractional_matrix);
        std::cout << "Collapsed Matrix by PROfit: " << std::endl;
        std::cout << collapsed_fraction << std::endl;


        Eigen::MatrixXd correlation_matrix = PROsyst::GenerateCorrMatrix(fractional_matrix);
        std::cout << "Correlation Matrix" << std::endl;
        std::cout << correlation_matrix << std::endl;


    }
    return 0;
}
