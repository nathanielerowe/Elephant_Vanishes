#include "PROconfig.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROtocall.h"
#include "PROcreate.h"
#include "PROpeller.h"
#include "PROchi.h"
#include "PROcess.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <Eigen/Eigen>
#include <algorithm>
#include <numeric>
#include <random>
#include <thread>

using namespace PROfit;

//log_level_t GLOBAL_LEVEL = LOG_DEBUG;
log_level_t GLOBAL_LEVEL = LOG_ERROR;

struct fc_out{
    int niter_syst, niter_osc;
    double chi2_syst, chi2_osc, dmsq, sinsq2tmm;
    Eigen::VectorXd best_fit_syst, best_fit_osc, syst_throw;
};

struct fc_args {
    size_t todo;
    std::vector<double>* dchi2s;
    std::vector<fc_out>* out;
    size_t *count;
    const PROconfig config;
    const PROpeller prop;
    PROsyst systs;
    Eigen::VectorXd cv, err;
    Eigen::MatrixXd L;
    int thread;
};

void fc_worker(fc_args args) {
    auto dchi2s = args.dchi2s;
    auto out = args.out;
    auto count = args.count;
    auto &config = args.config;
    auto &prop = args.prop;
    auto &systs = args.systs;
    auto &cv = args.cv;
    auto &err = args.err;
    auto &L = args.L;
    std::random_device rd{};
    std::mt19937 rng{rd()};
    PROsc osc(prop);
    while(dchi2s->size() < args.todo) {
        (*count)++;
        std::cout << "Thread #" << args.thread << " Size " << dchi2s->size() << " count " << *count << "\n";
        std::normal_distribution<float> d;
        std::vector<float> throws;
        Eigen::VectorXd throwC = Eigen::VectorXd::Constant(cv.size(), 0);
        for(size_t i = 0; i < systs.GetNSplines(); i++)
            throws.push_back(d(rng));
        for(size_t i = 0; i < cv.size(); i++)
            throwC(i) = d(rng);
        PROspec shifted = systs.GetSplineShiftedSpectrum(config, prop, throws);
        PROspec newSpec = PROspec::PoissonVariation(PROspec(shifted.Spec() + L * throwC, shifted.Error()));

        // No oscillations
        LBFGSpp::LBFGSBParam<double> param;  
        param.epsilon = 1e-6;
        param.max_iterations = 100;
        param.max_linesearch = 50;
        param.delta = 1e-6;
        LBFGSpp::LBFGSBSolver<double> solver(param); 
        int nparams = systs.GetNSplines();
        PROchi chi("3plus1",&config,&prop,&systs,NULL, newSpec, nparams, systs.GetNSplines(), PROchi::BinnedChi2);
        Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
        Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
        Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);

        double fx;
        int niter;
        std::vector<double> chi2s;
        int nfit = 0;
        do {
            nfit++;
            for(size_t i = 0; i < nparams; ++i)
                x(i) = 0.3*d(rng);
            // x will be overwritten to be the best point found
            log<LOG_INFO>(L"%1% || Fit without oscillations") % __func__;
            try {
                niter = solver.minimize(chi, x, fx, lb, ub);
            } catch(std::runtime_error &except) {
                log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
                continue;
            }
            chi2s.push_back(fx);
        } while(chi2s.size() < 10 && nfit < 100);

        if(chi2s.size() < 10) continue;
        fx = *std::min_element(chi2s.begin(), chi2s.end());
        //std::cout<<"Chis min: "<<fx<<" || ";
        //for(auto &c: chi2s)std::cout<<c<<" ";
        //std::cout<<std::endl;
        log<LOG_INFO>(L"%1% || Chi2 min %2% || %3%") % __func__ % fx % chi2s;

        // With oscillations
        LBFGSpp::LBFGSBParam<double> param_osc;  
        param_osc.epsilon = 1e-6;
        param_osc.max_iterations = 100;
        param_osc.max_linesearch = 250;
        param_osc.delta = 1e-6;
        LBFGSpp::LBFGSBSolver<double> solver_osc(param_osc); 
        nparams = 2 + systs.GetNSplines();
        PROchi chi_osc("3plus1",&config,&prop,&systs,&osc, newSpec, nparams, systs.GetNSplines(), PROchi::BinnedChi2);
        Eigen::VectorXd lb_osc = Eigen::VectorXd::Constant(nparams, -3.0);
        //lb_osc(0) = 0.01; lb_osc(1) = 0;
        lb_osc(0) = -2; lb_osc(1) = -std::numeric_limits<double>::infinity();
        Eigen::VectorXd ub_osc = Eigen::VectorXd::Constant(nparams, 3.0);
        //ub_osc(0) = 100; ub_osc(1) = 1;
        ub_osc(0) = 2; ub_osc(1) = 0;
        Eigen::VectorXd x_osc = Eigen::VectorXd::Constant(nparams, 0.0);
        x_osc(0) = 1.0;

        // x will be overwritten to be the best point found
        log<LOG_INFO>(L"%1% || Fit with oscillations") % __func__;
        double fx_osc;
        int niter_osc;
        std::vector<std::tuple<double, double, double>> res;
        nfit = 0;
        do {
            nfit++;
            for(size_t i = 0; i < nparams; ++i)
                x_osc(i) = d(rng);
            x_osc(0) = x_osc(0) < -2 ? -2 : x_osc(0) > 2 ? 2 : x_osc(0);
            x_osc(1) -= 4;
            try {
                niter_osc = solver_osc.minimize(chi_osc, x_osc, fx_osc, lb_osc, ub_osc);
            } catch(std::runtime_error &except) {
                log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
                //for(auto &t: throws) log<LOG_ERROR>(L"%1% || Throws, %2%") % __func__ % t;
                continue;
            }
            res.emplace_back(fx_osc, x(0), x(1));
        } while(res.size() < 10 && nfit < 100);
        if(res.size() < 10) continue;
        const auto [chi2_osc, ldmsq, lss2t] = *std::min_element(res.begin(), res.end(),
            [](const auto& l, const auto& r) {
              return std::get<0>(l) < std::get<0>(r);
            });
        //log<LOG_INFO>(L"%1% || Chi2 min %2% || %3%") % __func__ % fx_osc % chi2s;
        

        Eigen::VectorXd t = Eigen::VectorXd::Constant(throws.size(), 0);
        for(size_t i = 0; i < throws.size(); i++) t(i) = throws[i];

        out->push_back({niter, niter_osc, fx, chi2_osc, std::pow(10, ldmsq), std::pow(10, lss2t), x, x_osc.segment(2, nparams-2), t});

        dchi2s->push_back(std::abs(fx - chi2_osc ));
    }
}

int main(int argc, char* argv[])
{
    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 100;
    size_t nfit = 1, nthread = 1;

    //doubles
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");
    app.add_option("-n,--nfit",nfit, "Number of fits.");
    app.add_option("-t,--nthread",nthread, "Number of threads.");

    CLI11_PARSE(app, argc, argv);

    if(nthread > nfit) nthread = nfit;

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

    std::random_device rd{};
    std::mt19937 gen{rd()};

    Eigen::VectorXd data = systsstructs.back().CV().Spec();
    Eigen::VectorXd err = systsstructs.back().CV().Error();
    Eigen::MatrixXd diag = data.array().matrix().asDiagonal();
    Eigen::MatrixXd full_cov = diag * systs.fractional_covariance * diag;
    //Eigen::MatrixXd stat_cov = CollapseMatrix(myConf, diag);
    Eigen::LLT<Eigen::MatrixXd> llt(CollapseMatrix(myConf, full_cov));//+ stat_cov);
    
    std::vector<std::vector<double>> dchi2s;
    dchi2s.reserve(nthread);
    std::vector<std::vector<fc_out>> outs;
    outs.reserve(nthread);
    std::vector<size_t> counts;
    counts.reserve(nthread);
    std::vector<std::thread> threads;
    for(size_t i = 0; i < nthread; i++) {
        dchi2s.emplace_back();
        outs.emplace_back();
        counts.push_back(0);
        // TODO: nfit/nthread is not correct if nthread does not evenly divide nfit
        fc_args args{nfit/nthread, &dchi2s.back(), &outs.back(), &counts.back(), myConf, myprop, systs, data, err, llt.matrixL(), (int)i};
        threads.emplace_back(fc_worker, args);
    }
    for(auto&& t: threads) {
        t.join();
    }

    TH1D hdchi2s("hdchi2s", ";#Delta#chi^{2}", 50, 0, 50);
    std::vector<double> dmsqedges;
    for(size_t i = 0; i < 41; ++i)
        dmsqedges.push_back(std::pow(10.0, -2.0 + i * 4.0 / 40));
    TH1D hdmsqs("hdmsqs", ";#Deltam^{2}", 40, dmsqedges.data());
    TH1D hss2t("hss2t", ";sin^{2}2#theta_{#mu#mu}", 50, 0, 1);

    size_t count = std::accumulate(counts.begin(), counts.end(), 0);
    std::vector<double> flattened_dchi2s;
    for(const auto& v: dchi2s) {
        for(const auto& dchi2: v) {
            flattened_dchi2s.push_back(dchi2);
            hdchi2s.Fill(dchi2);
        }
    }
    for(const auto& out: outs) {
      for(const auto& fco: out) {
        hdmsqs.Fill(fco.dmsq);
        hss2t.Fill(fco.sinsq2tmm);
      }
    }

    TCanvas c;
    c.SetLogy();
    hdchi2s.Draw("hist");
    c.Print("dchi2s.pdf");
    hss2t.Draw("hist");
    c.Print("sinsq2tmm.pdf");
    c.SetLogx();
    hdmsqs.Draw("hist");
    c.Print("dmsqs.pdf");

    std::sort(flattened_dchi2s.begin(), flattened_dchi2s.end());

    std::cout << "90% delta chi2: " << flattened_dchi2s[0.9*flattened_dchi2s.size()] << std::endl;
    std::cout << "Total attempted throws: " << count << std::endl;

    return 0;
}

