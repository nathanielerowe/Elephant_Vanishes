#include "PROconfig.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROcess.h"
#include "PROsurf.h"
#include "PROfitter.h"
#include "PROmodel.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <Eigen/Eigen>

#include "PROtocall.h"
#include "TH2D.h"
#include "TStyle.h"
#include <filesystem>

using namespace PROfit;

log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[])
{
    gStyle->SetOptStat(0);
    CLI::App app{"PROfit, a PROfessional, PROductive fitting and oscillation framework. Together let's minimize PROfit!"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    std::vector<float> osc_params;
    std::map<std::string, float> injected_systs;

    std::string analysis_tag = "PROfit"
    size_t nthread = 1;


    //Global Arguments for all PROfit enables subcommands.
    app.add_option("-x,--xml", xmlname, "Input PROfit XML configuration file.")->required();
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4]->[Error,Warning,Info,Debug].")->default_val(GLOBAL_LEVEL);
    app.add_option("-t,--tag", analysis_tag, "Analysis Tag used for output identification.")->default_str("PROfit");
    app.add_option("-n, --nthread",   nthread, "Number of threads to parallelize over.")->default_val(1);
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-c, --chi2", chi2, "Which chi2 function to use. Options are PROchi or PROCNP")->default_str("PROchi");
    app.add_option("--inject", osc_params, "Physics parameters to inject as true signal.")->expected(-1);// HOW TO
    app.add_option("--syst-list", syst_list, "Override list of systematics to use (note: all systs must be in the xml).");
    app.add_option("--exclude-systs", systs_excluded, "List of systematics to exclude.")->excludes("--syst-list"); 
    app.add_flag("--scale-by-width", cv_scale, "Scale histgrams by 1/(bin width).");
    app.add_flag("--event-by-event", eventbyevent, "Do you want to weight event-by-event?");
    app.add_flag("--statonly", statonly, "Run a stats only surface instead of fitting systematics");
    app.add_flag("--shapeonly", shapeonly, "Run a shape only analysis");
    app.add_flag("--rateonly", rateonly, "Run a rate only analysis");
    app.add_flag("--force",force,"Force loading binary data even if hash is incorrect (Be Careful!)");


    //PROcess, into binary data [Do this once first!]
    CLI::App *process_command = app.add_subcommand("process", "PROcess the MC and systematics in root files into binary data for future rapid loading.");

    //PROsurf, make a 2D surface scan of physics parameters
    CLI::App *surface_command = app.add_subcommand("surface", "Make a 2D surface scan of two physics parameters, profiling over all others.");
    surface_command->add_option("-g, --grid", grid_size, "Set grid size. If one dimension passed, grid assumed to be square, else rectangular")->expected(0, 2)->default_val(40);
    CLI::Option *xlim_opt = surface_command->add_option("--xlims", xlims, "Limits for x-axis");
    CLI::Option *ylim_opt = surface_command->add_option("--ylims", ylims, "Limits for y-axis");
    surface_command->add_option("--xlo", xlo, "Lower limit for x-axis")->excludes(xlim_opt)->default_val(1e-4);
    surface_command->add_option("--xhi", xhi, "Upper limit for x-axis")->excludes(xlim_opt)->default_val(1);
    surface_command->add_option("--ylo", ylo, "Lower limit for y-axis")->excludes(ylim_opt)->default_val(1e-2);
    surface_command->add_option("--yhi", yhi, "Upper limit for y-axis")->excludes(ylim_opt)->default_val(1e2);
    surface_command->add_option("--xlabel", xlabel, "X-axis label")->default_val("sin^{2}2#theta_{#mu#mu}");
    surface_command->add_option("--ylabel", ylabel, "Y-axis label")->default_val("#Deltam^{2}_{41}");
    surface_command->add_flag("--logx,!--linx", logx, "Specify if x-axis is logarithmic or linear (default log)");
    surface_command->add_flag("--logy,!--liny", logy, "Specify if y-axis is logarithmic or linear (default log)");

    //PROfile, make N profile'd chi^2 for each physics and nuisence parameters
    CLI::App *profile_command = app.add_subcommand("profile", "Make a 1D profiled chi2 for each physics and nuisence parameter.");

    //PROplot, plot things
    CLI::App *plot_command = app.add_subcommand("prot", "Make plots of CV, or injected point with error bars and covariance.");
    plot_command->add_flag("--with-splines", with_splines, "Include graphs of splines in output.");



    CLI11_PARSE(app, argc, argv);

    log<LOG_INFO>(L"%1% || PROfit commandline input arguments. xml: %2%, tag: %3%, nthread: %4% ") % __func__ % xmlname.c_str() % analysis_tag.c_str() % nthread ;

    //Initilize configuration from the XML;
    PROconfig config(xmlname);

    //Inititilize PROpeller to keep MC
    PROpeller prop;

    //Initilize objects for systematics storage
    std::vector<SystStruct> systsstructs;

    //input/output logic
    std::string propBinName = analysis_tag+"_prop.bin";
    std::string systBinName = analysis_tag+"_prop.bin";

    if((*process_command) || (!std::filesystem::exists(systBinName) || !std::filesystem::exists(propBinName))  ){
        log<LOG_INFO>(L"%1% || Processing PROpeller and PROsysts from XML defined root files, and saving to binary output also: %2%") % __func__ % propBinName.c_str();
        //Process the CAF files to grab and fill all SystStructs and PROpeller
        PROcess_CAFAna(config, systsstructs, prop);
        prop.save(propBinName);    
        saveSystStructVector(systsstructs,systBinName);
        log<LOG_INFO>(L"%1% || Done processing PROpeller and PROsysts from XML defined root files, and saving to binary output also: %2%") % __func__ % propBinName.c_str();
        return 0;
    }else{
        log<LOG_INFO>(L"%1% || Loading PROpeller and PROsysts from precalc binary input: %2%") % __func__ % propBinName.c_str();
        prop.load(propBinName);
        loadSystStructVector(systsstructs, systBinName);

        log<LOG_INFO>(L"%1% || Done loading. Config hash (%2%) and binary loaded PROpeller (%3%) or PROsyst hash(%4%) are here. ") % __func__ %  config.hash % prop.hash % systsstructs[0].hash;
        if(config.hash!=prop.hash && config.hash!=systsstructs.front().hash){
            if(force){
                log<LOG_WARNING>(L"%1% || WARNING config hash (%2%) and binary loaded PROpeller (%3%) or PROsyst hash(%4%) not compatable! ") % __func__ %  config.hash % prop.hash % systsstructs.front().hash;
                log<LOG_WARNING>(L"%1% || WARNING But we are forcing ahead, be SUPER clear and happy you understand what your doing.  ") % __func__;
            }else{
                log<LOG_ERROR>(L"%1% || ERROR config hash (%2%) and binary loaded PROpeller (%3%) or PROsyst hash(%4%) not compatable! ") % __func__ %  config.hash % prop.hash % systsstructs.front().hash;
                return 1;
            }
        }
    }

    //Build a PROsyst to sort and analyze all systematics
    PROsyst systs(systsstructs);
    std::unique_ptr<PROmodel> model = get_model_from_string(config.m_model_tag, prop);

    if(syst_list.size()) {
        systs = systs.subset(syst_list);
    } else if(systs_excluded.size()) {
        systs = systs.excluding(systs_excluded);
    }

    if(!osc_params.size()) {
        log<LOG_ERROR>(L"%1% || Expected %2% physics parameters to be provided for oscillation plot.")
            % __func__ % model->nparams;
        exit(EXIT_FAILURE);
    }
    Eigen::VectorXf pparams = Eigen::VectorXf::Constant(model->nparams + systs.GetNSplines(), 0);
    if(osc_params.size()) {
        if(osc_params.size() != model->nparams) {
            log<LOG_ERROR>(L"%1% || Incorrect number of physics parameters provided. Expected %2%, found %3%.")
                % __func__ % model->nparams % osc_params.size();
            exit(EXIT_FAILURE);
        }
        for(size_t i = 0; i < osc_params.size(); ++i) {
            pparams(i) = std::log10(osc_params[i]);
        }
    }

    Eigen::VectorXf allparams = Eigen::VectorXf::Constant(model->nparams + systs.GetNSplines(), 0);
    for(int i = 0; i < pparams.size(); ++i) allparams(i) = pparams(i);
    for(const auto& [name, shift]: injected_systs) {
        log<LOG_INFO>(L"%1% || Injected syst: %2% shifted by %3%") % __func__ % name.c_str() % shift;
        auto it = std::find(systs.spline_names.begin(), systs.spline_names.end(), name);
        if(it == systs.spline_names.end()) {
            log<LOG_ERROR>(L"%1% || Error: Unrecognized spline %2%. Ignoring this injected shift.") % __func__ % name.c_str();
            continue;
        }
        int idx = std::distance(systs.spline_names.begin(), it);
        allparams(idx) = shift;
    }

    //Create CV or injected data spectrum for all subsequent steps
    PROspec data = osc_params.size() || injected_systs.size() ? FillRecoSpectra(config, prop, systs, *model, allparams, !eventbyevent) :  FillCVSpectrum(config, prop, !eventbyevent);
    Eigen::VectorXf data_vec = CollapseMatrix(config, data.Spec());
    Eigen::VectorXf err_vec_sq = data.Error().array().square();
    Eigen::VectorXf err_vec = CollapseMatrix(config, err_vec_sq).array().sqrt();
    data = PROspec(data_vec, err_vec);


    //Some global minimizer params
    LBFGSpp::LBFGSBParam<float> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    param.max_linesearch = 250;
    param.delta = 1e-6;

    //Metric Time
    PROmetric *metric;
    if(chi2 == "PROchi") {
        metric = new PROchi("", config, prop, &systs, *model, data, eventbyevent ? PROmetric::EventByEvent : PROmetric::BinnedChi2);
    } else if(chi2 == "PROCNP") {
        metric = new PROCNP("", config, prop, &systs, *model, data, eventbyevent ? PROmetric::EventByEvent : PROmetric::BinnedChi2);
    } else {
        log<LOG_ERROR>(L"%1% || Unrecognized chi2 function %2%") % __func__ % chi2.c_str();
        abort();
    }


    //***********************************************************************
    //***********************************************************************
    //******************** PROfile PROfile PROfile **************************
    //***********************************************************************
    //***********************************************************************

    if(*profile_command){


        size_t nparams = 2 + systs.GetNSplines();
        Eigen::VectorXf lb = Eigen::VectorXf::Constant(nparams, -3.0);
        lb(0) = -2; lb(1) = -std::numeric_limits<float>::infinity();
        Eigen::VectorXf ub = Eigen::VectorXf::Constant(nparams, 3.0);
        ub(0) = 2; ub(1) = 0;
        for(size_t i = 2; i < nparams; ++i) {
            lb(i) = systs.spline_lo[i-2];
            ub(i) = systs.spline_hi[i-2];
        }
        PROfitter fitter(ub, lb, param);

        float chi2 = fitter.Fit(*metric); 
        Eigen::VectorXf best_fit = fitter.best_fit;
        Eigen::MatrixXf post_covar = fitter.Covariance();

        //std::string hname = "#chi^{2}/ndf = " + to_string(chi2) + "/" + to_string(config.m_num_bins_total_collapsed);
        std::string hname = "";

        Eigen::VectorXf subvector1 = best_fit.segment(0, 2);
        std::vector<float> fitparams(subvector1.data(), subvector1.data() + subvector1.size());

        Eigen::VectorXf subvector2 = best_fit.segment(2, systs.GetNSplines());
        std::vector<float> shifts(subvector2.data(), subvector2.data() + subvector2.size());

        Eigen::VectorXf post_fit = CollapseMatrix(config, FillRecoSpectra(config, prop, systs, *model, best_fit, true).Spec());
        TH1D post_hist("ph", hname.c_str(), config.m_num_bins_total_collapsed, config.m_channel_bin_edges[0].data());
        for(size_t i = 0; i < config.m_num_bins_total_collapsed; ++i) {
            post_hist.SetBinContent(i+1, post_fit(i));
        }

        PROfile(config, prop, systs, *model, data, chi, analysis_tag+".pdf", true, nthread, best_fit);

        return 0;

        //***********************************************************************
        //***********************************************************************
        //******************** PROsurf PROsurf PROsurf **************************
        //***********************************************************************
        //***********************************************************************
    }else if(*surface_command){

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

        //Define grid and Surface
        size_t nbinsx = grid_size[0], nbinsy = grid_size[1];
        PROsurf surface(*metric, 1, 0, nbinsx, logx ? PROsurf::LogAxis : PROsurf::LinAxis, xlo, xhi,
                nbinsy, logy ? PROsurf::LogAxis : PROsurf::LinAxis, ylo, yhi);

        if(statonly)
            surface.FillSurfaceStat(config, analysis_tag+"_statonly_surface.txt");
        else
            surface.FillSurface(analysis_tag+"_surface.txt");

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

        log<LOG_INFO>(L"%1% || Saving surface to %2% as TH2D named \"surf.\"") % __func__ % analysis_tag.c_str();
        TFile fout((analysis_tag+"_surf.root").c_str(), "RECREATE");
        surf.Write();

        TCanvas c;
        if(logy)
            c.SetLogy();
        if(logx)
            c.SetLogx();
        c.SetLogz();
        surf.Draw("colz");
        c.Print((analysis_tag+"_surface.pdf").c_str());
        delete metric;
        fout.close();
        return 0;

        //***********************************************************************
        //***********************************************************************
        //******************** PROplot PROplot PROplot **************************
        //***********************************************************************
        //***********************************************************************
    }else if(*proplot_command){





        //***************************** END *********************************
    }else{
        log<LOG_WARNING>(L"%1% || Please pass a subcommand to tell PROfit to do something! see --help for ideas.") % __func__);
        return 1;
    }


    return 0;
}

