#include "PROconfig.h"
#include "PROlog.h"
#include "PROmetric.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROtocall.h"
#include "PROcreate.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROCNP.h"
#include "PROcess.h"
#include "PROsurf.h"

#include "CLI11.h"

#include <Eigen/Eigen>

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
    std::string xmlname, filename, xlabel, ylabel;
    int maxevents;
    size_t  nthread;
    float xlo, xhi, ylo, yhi;
    std::array<float, 2> xlims, ylims, injected_pt{0, 0};
    std::vector<int> grid_size;
    std::map<std::string, float> injected_systs;
    std::vector<std::string> syst_list, systs_excluded;
    bool eventbyevent=false, statonly = false, logx=true, logy=true;
    std::string chi2;

    app.add_option("-x, --xml",       xmlname, "Input PROfit XML config.")->required();
    app.add_option("-m, --max",       maxevents, "Max number of events to run over.")->default_val(50000);
    app.add_option("-v, --verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].")->default_val(LOG_ERROR);
    app.add_option("-t, --nthread",   nthread, "Number of fits.")->default_val(1);
    app.add_option("-o, --outfile",   filename, "If you want chisq to be dumped to text file, provide name")->default_val("");
    app.add_option("-g, --grid", grid_size, "Set grid size. If one dimension passed, grid assumed to be square, else rectangular")->expected(0, 2)->default_val(40);
    app.add_option("-c, --chi2", chi2, "Which chi2 function to use. Options are PROchi or PROCNP")->default_str("PROchi");
    CLI::Option *xlim_opt = app.add_option("--xlims", xlims, "Limits for x-axis");
    CLI::Option *ylim_opt = app.add_option("--ylims", ylims, "Limits for y-axis");
    app.add_option("--xlo", xlo, "Lower limit for x-axis")->excludes(xlim_opt)->default_val(1e-4);
    app.add_option("--xhi", xhi, "Upper limit for x-axis")->excludes(xlim_opt)->default_val(1);
    app.add_option("--ylo", ylo, "Lower limit for y-axis")->excludes(ylim_opt)->default_val(1e-2);
    app.add_option("--yhi", yhi, "Upper limit for y-axis")->excludes(ylim_opt)->default_val(1e2);
    app.add_option("--xlabel", xlabel, "X-axis label")->default_val("sin^{2}2#theta_{#mu#mu}");
    app.add_option("--ylabel", ylabel, "Y-axis label")->default_val("#Deltam^{2}_{41}");
    app.add_option("--inject", injected_pt, "Physics parameters to inject as true signal.")->default_str("0 0");
    app.add_option("--inject-systs", injected_systs, "Systematic shifts to inject. Map of name and shift value in sigmas. Only spline systs are supported right now.");
    app.add_option("--syst-list", syst_list, "Override list of systematics to use (note: all systs must be in the xml).");
    app.add_option("--exclude-systs", systs_excluded, "List of systematics to exclude.")->excludes("--syst-list"); 

    app.add_flag("--event-by-event",    eventbyevent, "Do you want to weight event-by-event?");
    app.add_flag("-s, --statonly", statonly, "Run a stats only surface instead of fitting systematics");
    app.add_flag("--logx,!--linx", logx, "Specify if x-axis is logarithmic or linear (default log)");
    app.add_flag("--logy,!--liny", logy, "Specify if y-axis is logarithmic or linear (default log)");

    CLI11_PARSE(app, argc, argv);


    if (grid_size.empty()) {
        grid_size = {40, 40};
    }
    if (grid_size.size() == 1) {
        grid_size.push_back(grid_size[0]); //make it square
    }

    if(*xlim_opt) {
        xlo = xlims[0];
        xhi = xlims[1];
    }
    if(*ylim_opt) {
        ylo = ylims[0];
        yhi = ylims[1];
    }

    bool savetoroot = filename.size() > 5 && filename.substr(filename.size() - 5) == ".root";

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
    
    std::vector<float> pparams = {std::log10(injected_pt[0]), std::log10(injected_pt[1])};
    log<LOG_INFO>(L"%1% || Injected point: sinsq2t = %2%, dmsq = %3%") % __func__ % injected_pt[0] % injected_pt[1];
    for(const auto& [name, shift]: injected_systs)
      log<LOG_INFO>(L"%1% || Injected syst: %2% shifted by %3%") % __func__ % name.c_str() % shift;
    //Grab Asimov Data
    PROspec data = injected_pt[0] != 0 && injected_pt[1] != 0 ? FillRecoSpectra(config, prop, systs, &osc, injected_systs, pparams, !eventbyevent) :
                   injected_systs.size() ? FillRecoSpectra(config, prop, systs, injected_systs, !eventbyevent) :
                   FillCVSpectrum(config, prop, !eventbyevent);
    Eigen::VectorXf data_vec = CollapseMatrix(config, data.Spec());
    Eigen::VectorXf err_vec_sq = data.Error().array().square();
    Eigen::VectorXf err_vec = CollapseMatrix(config, err_vec_sq).array().sqrt();
    data = PROspec(data_vec, err_vec);

    if(syst_list.size()) {
      systs = systs.subset(syst_list);
    } else if(systs_excluded.size()) {
      systs = systs.excluding(systs_excluded);
    }

    PROmetric *metric;
    if(chi2 == "PROchi") {
        metric = new PROchi("", &config, &prop, &systs, &osc, data, systs.GetNSplines(), systs.GetNSplines(), PROmetric::BinnedChi2);
    } else if(chi2 == "PROCNP") {
        metric = new PROCNP("", &config, &prop, &systs, &osc, data, systs.GetNSplines(), systs.GetNSplines(), PROmetric::BinnedChi2);
    } else {
        log<LOG_ERROR>(L"%1% || Unrecognized chi2 function %2%") % __func__ % chi2.c_str();
        abort();
    }

    //Define grid and Surface
    size_t nbinsx = grid_size[0], nbinsy = grid_size[1];
    PROsurf surface(*metric, nbinsx, logx ? PROsurf::LogAxis : PROsurf::LinAxis, xlo, xhi,
                    nbinsy, logy ? PROsurf::LogAxis : PROsurf::LinAxis, ylo, yhi);
    
    if(statonly)
        surface.FillSurfaceStat(config, !savetoroot ? filename : "");
    else
        surface.FillSurface(systs, !savetoroot ? filename : "", nthread);

    //And do a PROfile of pulls at the data also
    //PROfile(config,prop,systs,osc,data,filename+"_PROfile");

    std::vector<float> binedges_x, binedges_y;
    for(size_t i = 0; i < surface.nbinsx+1; i++)
        binedges_x.push_back(logx ? std::pow(10, surface.edges_x(i)) : surface.edges_x(i));
    for(size_t i = 0; i < surface.nbinsy+1; i++)
        binedges_y.push_back(logy ? std::pow(10, surface.edges_y(i)) : surface.edges_y(i));

    TH2D surf("surf", (";"+xlabel+";"+ylabel).c_str(), surface.nbinsx, binedges_x.data(), surface.nbinsy, binedges_y.data());

    for(size_t i = 0; i < surface.nbinsx; i++) {
        for(size_t j = 0; j < surface.nbinsy; j++) {
            surf.SetBinContent(i+1, j+1, surface.surface(i, j));
        }
    }

    if(savetoroot) {
      log<LOG_INFO>(L"%1% || Saving surface to %2% as TH2D named \"surf.\"") % __func__ % filename.c_str();
      TFile fout(filename.c_str(), "RECREATE");
      surf.Write();
    }

    TCanvas c;
    if(logy)
        c.SetLogy();
    if(logx)
        c.SetLogx();
    c.SetLogz();
    surf.Draw("colz");
    c.Print("PROfit_surface.pdf");

    delete metric;

    return 0;
}
