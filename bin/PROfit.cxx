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

using namespace PROfit;

log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[])
{
    gStyle->SetOptStat(0);
    CLI::App app{"PROfit, a PROfessional, PROductive fitting and oscillation framework. Together let's minimize PROfit!"}; 

    // Define options
    std::string xmlname = "NULL.xml", filename = "profit"; 
    std::array<float, 2> injected_pt{0, 0};
    std::map<std::string, float> injected_systs;
   
    std::string analysis_tag;
    size_t nthread = 1;
   

    //Global Arguments for all PROfit enables subcommands.
    app.add_option("-x,--xml", xmlname, "Input PROfit XML configuration file.")->required();
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4]->[Error,Warning,Info,Debug].")->default_val(GLOBAL_LEVEL);
    app.add_option("-t,--tag", analysis_tag, "Analysis Tag used for output identification.")->default_str("PROfit");
    app.add_option("-n, --nthread",   nthread, "Number of threads to parallelize over.")->default_val(1);
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-c, --chi2", chi2, "Which chi2 function to use. Options are PROchi or PROCNP")->default_str("PROchi");
    app.add_option("--inject", injected_pt, "Physics parameters to inject as true signal.")->default_str("0 0");// HOW TO
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

    log<LOG_INFO>(L"%1% || PROfit commandline input arguments. xml: %2%, outfile: %3%, nthread: %4% ") % __func__ % xmlname.c_str() % filename.c_str() % nthread ;

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
    //PROsc osc(myprop);
    std::unique_ptr<PROmodel> model = get_model_from_string(myConf.m_model_tag, myprop);

    Eigen::VectorXf pparams{{std::log10(injected_pt[0]), std::log10(injected_pt[1])}};
    log<LOG_INFO>(L"%1% || PROfit Injected point: sinsq2t  %2% dmsq %3%") % __func__ % injected_pt[0] % injected_pt[1] ;
    PROspec data = injected_pt[0] != 0 && injected_pt[1] != 0 ? 
        FillRecoSpectra(myConf, myprop, systs, *model, pparams, true) :
        FillCVSpectrum(myConf, myprop, true);
    Eigen::VectorXf data_vec = CollapseMatrix(myConf, data.Spec());
    Eigen::VectorXf err_vec_sq = data.Error().array().square();
    Eigen::VectorXf err_vec = CollapseMatrix(myConf, err_vec_sq).array().sqrt();
    data = PROspec(data_vec, err_vec);

    Eigen::VectorXf cv_spec = CollapseMatrix(myConf, FillCVSpectrum(myConf, myprop, true).Spec());

    TH1D cv_hist("cv", "CV", myConf.m_num_bins_total_collapsed, myConf.m_channel_bin_edges[0].data());
    TH1D data_hist("dh", "Data", myConf.m_num_bins_total_collapsed, myConf.m_channel_bin_edges[0].data());
    for(size_t i = 0; i < myConf.m_num_bins_total_collapsed; ++i) {
        cv_hist.SetBinContent(i+1, cv_spec(i));
        data_hist.SetBinContent(i+1, data_vec(i));
        data_hist.SetBinError(i+1, std::sqrt(data_vec(i)));
    }

    //PROfile(myConf, myprop, systs, osc, data, "profit_test", true);

    PROchi chi("", myConf, myprop, &systs, *model, data, PROfit::PROchi::BinnedChi2);

    LBFGSpp::LBFGSBParam<float> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    param.max_linesearch = 250;
    param.delta = 1e-6;

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

    float chi2 = fitter.Fit(chi); 
    Eigen::VectorXf best_fit = fitter.best_fit;
    //Eigen::MatrixXd post_covar = fitter.ScaledCovariance(chi2, myConf.m_num_bins_total_collapsed);
    Eigen::MatrixXf post_covar = fitter.Covariance();

    //std::string hname = "#chi^{2}/ndf = " + to_string(chi2) + "/" + to_string(myConf.m_num_bins_total_collapsed);
    std::string hname = "";

    Eigen::VectorXf subvector1 = best_fit.segment(0, 2);
    std::vector<float> fitparams(subvector1.data(), subvector1.data() + subvector1.size());

    Eigen::VectorXf subvector2 = best_fit.segment(2, systs.GetNSplines());
    std::vector<float> shifts(subvector2.data(), subvector2.data() + subvector2.size());

    Eigen::VectorXf post_fit = CollapseMatrix(myConf, FillRecoSpectra(myConf, myprop, systs, *model, best_fit, true).Spec());
    TH1D post_hist("ph", hname.c_str(), myConf.m_num_bins_total_collapsed, myConf.m_channel_bin_edges[0].data());
    for(size_t i = 0; i < myConf.m_num_bins_total_collapsed; ++i) {
        post_hist.SetBinContent(i+1, post_fit(i));
    }

    if(false){//currently turn off these as not working as wanted
        TCanvas ch;
        cv_hist.SetTitle(hname.c_str());
        cv_hist.SetLineColor(kBlack);
        cv_hist.Draw("hist");
        post_hist.SetLineColor(kRed);
        post_hist.Draw("hist same");
        data_hist.Draw("E same");
        ch.Print((filename+"_hists.pdf").c_str(), "pdf");

        TH2D cov("cov", hname.c_str(), post_covar.rows(), 0, post_covar.rows(), post_covar.cols(), 0, post_covar.cols());
        for(size_t i = 0; i < post_covar.rows(); ++i) {
            for(size_t j = 0; j < post_covar.cols(); ++j) {
                cov.SetBinContent(i+1, j+1, post_covar(i,j));
            }
        }
        TCanvas c1;
        cov.Draw("colz");
        c1.Print((filename+"_cov.pdf").c_str(), "pdf");

        std::vector<std::string> names;
        for(size_t i = 0; i < model->nparams; ++i) names.push_back(model->param_names[i]);
        for(size_t i = 0; i < systs.GetNSplines(); ++i) names.push_back(systs.spline_names[i]);

        TH1D *hsyst_pre  = new TH1D("hp", hname.c_str(), nparams, 0, nparams);
        TH1D *hsyst_post = new TH1D("ho", hname.c_str(), nparams, 0, nparams);
        for(size_t i = 0; i < nparams; i++) {
            hsyst_pre->GetXaxis()->SetBinLabel(i+1, names[i].c_str());
            hsyst_pre->SetBinContent(i+1, 0.0);
            hsyst_pre->SetBinError(i+1, 1.0);
            hsyst_post->GetXaxis()->SetBinLabel(i+1, names[i].c_str());
            hsyst_post->SetBinContent(i+1, best_fit(i));
            hsyst_post->SetBinError(i+1, std::sqrt(post_covar(i,i)));
        }
        //hsyst_pre->Write(("pulls_pre_"+det).c_str());
        //hsyst_post->Write(("pulls_post_"+det).c_str());
        TCanvas c("c","", 1200, 600);
        TPad p("pp", "pp", 0, 0, 1, 1);
        p.SetBottomMargin(0.4);
        p.SetLeftMargin(0.1);
        p.SetRightMargin(0.05);
        p.cd();
        hsyst_pre->SetLineColor(kRed);
        hsyst_pre->SetMarkerColor(kRed);
        hsyst_pre->GetYaxis()->SetTitle("Parameter Pull [#sigma]");
        hsyst_pre->Draw("E");
        hsyst_post->SetLineColor(kBlack);
        hsyst_post->SetMarkerColor(kBlack);
        hsyst_post->Draw("E same");
        c.cd();
        p.Draw();
        c.Print((filename+"_pulls.pdf").c_str(), "pdf");

    }

    PROfile(myConf, myprop, systs, *model, data, chi,filename, true, nthread, best_fit);

    return 0;
}

