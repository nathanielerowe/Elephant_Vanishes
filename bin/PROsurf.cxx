#include "PROconfig.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROtocall.h"
#include "PROcreate.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROcess.h"
#include "PROsurf.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>
#include <algorithm>
#include <numeric>
#include <random>
#include <thread>

#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace PROfit;

//log_level_t GLOBAL_LEVEL = LOG_DEBUG;
log_level_t GLOBAL_LEVEL = LOG_ERROR;

int main(int argc, char* argv[])
{

    gStyle->SetOptStat(0);
    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 50000;
    size_t nfit = 1, nthread = 1;
    //Define a filename to save chisq values in
    std::string filename;
    int gridDim = 40;

    //doubles
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");
    app.add_option("-g,--grid",gridDim, "Grid Dim, NxN (default 40)");
    app.add_option("-n,--nfit",nfit, "Number of fits.");
    app.add_option("-t,--nthread",nthread, "Number of fits.");
    app.add_option("-o, --outfile", filename, "text file name (preferably absolute path otherwise behavior undefined)");

    CLI11_PARSE(app, argc, argv);

    if(nthread > nfit) nthread = nfit;

    //Initilize configuration from the XML;
    PROconfig config(xmlname);

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

    size_t nbinsx = gridDim, nbinsy = gridDim;

    PROspec data = systsstructs.back().CV();

    PROsurf surface(nbinsx, PROsurf::LogAxis, 1e-4, 1.0, nbinsy, PROsurf::LogAxis, 1e-2, 1e2);
    //surface.FillSurfaceFast(config, prop, systs, osc, data, filename);
    surface.FillSurface(config, prop, systs, osc, data, filename);


    std::vector<double> binedges_x, binedges_y;
    for(size_t i = 0; i < surface.nbinsx+1; i++)
        binedges_x.push_back(std::pow(10, surface.edges_x(i)));
    for(size_t i = 0; i < surface.nbinsy+1; i++)
        binedges_y.push_back(std::pow(10, surface.edges_y(i)));

    TH2D surf("surf", ";sin^{2}2#theta_{#mu#mu};#Deltam^{2}", surface.nbinsx, binedges_x.data(), surface.nbinsy, binedges_y.data());

    for(size_t i = 0; i < surface.nbinsx; i++) {
        for(size_t j = 0; j < surface.nbinsy; j++) {
            surf.SetBinContent(i+1, j+1, surface.surface(i, j));
        }
    }

    TCanvas c;
    c.SetLogy();
    c.SetLogx();
    c.SetLogz();
    surf.Draw("colz");
    c.Print("PROfit_surface.pdf");

    return 0;
}
