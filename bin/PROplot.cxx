#include "PROconfig.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
#include "PROmodel.h"
#include "PROpeller.h"
#include "PROcess.h"
#include "PROtocall.h"

#include "CLI11.h"

#include <Eigen/Eigen>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

#include <cstdlib>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace PROfit;

log_level_t GLOBAL_LEVEL;

std::map<std::string, std::unique_ptr<TH1D>> getCVHists(const PROspec & spec, const PROconfig& inconfig, bool scale = false);
std::map<std::string, std::unique_ptr<TH2D>> covarianceTH2D(const PROsyst &syst, const PROconfig &config, const PROspec &cv);
std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> getSplineGraphs(const PROsyst &systs, const PROconfig &config);
std::unique_ptr<TGraphAsymmErrors> getErrorBand(const PROconfig &config, const PROpeller &prop, const PROsyst &syst, bool scale = false);

int main(int argc, char* argv[])
{
    gStyle->SetOptStat(0);
    CLI::App app{"PROplot, plotting spectra using ROOT"}; 

    // Define options
    std::string xmlname, filename; 
    std::vector<float> osc_params;
    std::map<std::string, float> apply_shift;
    std::vector<std::string> syst_list, systs_excluded;
    bool eventbyevent = false, with_splines = false, all_scale = false, cv_scale = false, err_scale = false, osc_scale = false;

    CLI::App *cv_command = app.add_subcommand("cv", "Make CV histograms");
    CLI::App *cov_command = app.add_subcommand("covariance", "Make 2D histogram of covariance matrix");
    CLI::App *spline_command = app.add_subcommand("splines", "Make graphs of the splines");
    CLI::App *errband_command = app.add_subcommand("error-band", "Make error band");
    CLI::App *osc_command = app.add_subcommand("osc", "Plot spectrum with an injected signal.");
    CLI::App *all = app.add_subcommand("all", "Make all plots (splines optional) using default options.");

    // Options for whole application
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.")->required();
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].")->default_val(LOG_WARNING);
    app.add_option("-o,--output", filename, "Output Filename.")->default_val("PROplot.pdf");
    app.add_option("--syst-list", syst_list, "Override list of systematics to use (note: all systs must be in the xml).");
    app.add_option("--exclude-systs", systs_excluded, "List of systematics to exclude.")->excludes("--syst-list"); 

    app.add_flag("--event-by-event",    eventbyevent, "Do you want to weight event-by-event?");

    cv_command->add_flag("--scale-by-width", cv_scale, "Scale histgrams by 1/(bin width).");
    errband_command->add_flag("--scale-by-width", err_scale, "Scale histgrams by 1/(bin width).");
    osc_command->add_flag("--scale-by-width", osc_scale, "Scale histograms by 1/(bin width).");
    osc_command->add_option("--signal", osc_params, "Oscillation parameters to use.")->expected(-1);

    all->add_flag("--with-splines", with_splines, "Include graphs of splines in output.");
    all->add_flag("--scale-by-width", all_scale, "Scale histgrams by 1/(bin width).");
    all->add_option("--signal", osc_params, "Oscillation parameters to use.")->expected(-1);


    //app.add_option("--apply-osc", apply_osc, "Apply oscillations with these paramters to the CV spectrum.")->default_str("0 0");
    //app.add_option("--apply-shift", apply_shift, "Apply these shifts to the CV spectrum.");

    CLI11_PARSE(app, argc, argv);

    bool savetoroot = filename.size() > 5 && filename.substr(filename.size() - 5) == ".root";
    bool savetopdf = filename.size() > 4 && filename.substr(filename.size() - 4) == ".pdf";
    if(!savetoroot && !savetopdf) {
        log<LOG_ERROR>(L"%1% || Output file specified, %2%, is not a pdf or root file.") % __func__ % filename.c_str();
        return 1;
    }

    //Initilize configuration from the XML;
    PROconfig config(xmlname);

    PROpeller prop1;
    std::vector<SystStruct> systsstructs1;
    PROcess_CAFAna(config, systsstructs1, prop1);

    

    prop1.save("mytest.txt");
    saveSystStructVector(systsstructs1,"structs.txt");

    PROpeller prop;
    std::vector<SystStruct> systsstructs;
    
    prop.load("mytest.txt");
    loadSystStructVector(systsstructs,"structs.txt");

    PROsyst systs(systsstructs);
    std::unique_ptr<PROmodel> model = get_model_from_string(config.m_model_tag, prop);

    if(syst_list.size()) {
      systs = systs.subset(syst_list);
    } else if(systs_excluded.size()) {
      systs = systs.excluding(systs_excluded);
    }

    if(*osc_command && !osc_params.size()) {
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

    PROspec spec = FillCVSpectrum(config, prop, !eventbyevent);


    if(savetopdf) {
        TCanvas c;
        
        c.Print((filename + "[").c_str(), "pdf");
        if(*all || *cv_command) {
            std::map<std::string, std::unique_ptr<TH1D>> cv_hists = getCVHists(spec, config, cv_scale || all_scale);
            size_t global_subchannel_index = 0;
            size_t global_channel_index = 0;
            for(size_t im = 0; im < config.m_num_modes; im++){
                for(size_t id =0; id < config.m_num_detectors; id++){
                    for(size_t ic = 0; ic < config.m_num_channels; ic++){
                        std::unique_ptr<THStack> s = std::make_unique<THStack>((std::to_string(global_channel_index)).c_str(),(std::to_string(global_channel_index)).c_str());
                        std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.59,0.89,0.59,0.89);
                        leg->SetFillStyle(0);
                        leg->SetLineWidth(0);
                        for(size_t sc = 0; sc < config.m_num_subchannels[ic]; sc++){
                            const std::string& subchannel_name  = config.m_fullnames[global_subchannel_index];
                            s->Add(cv_hists[subchannel_name].get());
                            leg->AddEntry(cv_hists[subchannel_name].get(), config.m_subchannel_plotnames[ic][sc].c_str() ,"f");
                            ++global_subchannel_index;
                        }

                        s->Draw("hist");
                        leg->Draw();

                        s->SetTitle((config.m_mode_names[im]  +" "+ config.m_detector_names[id]+" "+ config.m_channel_names[ic]).c_str());
                        s->GetXaxis()->SetTitle(config.m_channel_units[ic].c_str());
                        if(cv_scale || all_scale)
                            s->GetYaxis()->SetTitle("Events/GeV");
                        else
                            s->GetYaxis()->SetTitle("Events");

                        c.Print(filename.c_str(), "pdf");
                    }
                }
            }
        }

        if((*all && osc_params.size()) || *osc_command) {
            PROspec osc_spec = FillRecoSpectra(config, prop, systs, *model, pparams, !eventbyevent);
            std::map<std::string, std::unique_ptr<TH1D>> cv_hists = getCVHists(spec, config, osc_scale || all_scale);
            std::map<std::string, std::unique_ptr<TH1D>> osc_hists = getCVHists(osc_spec, config, osc_scale || all_scale);
            size_t global_subchannel_index = 0;
            size_t global_channel_index = 0;
            for(size_t im = 0; im < config.m_num_modes; im++){
                for(size_t id =0; id < config.m_num_detectors; id++){
                    for(size_t ic = 0; ic < config.m_num_channels; ic++){
                        TH1D* osc_hist = NULL;
                        TH1D* cv_hist = NULL;
                        for(size_t sc = 0; sc < config.m_num_subchannels[ic]; sc++){
                            const std::string& subchannel_name  = config.m_fullnames[global_subchannel_index];
                            const auto &h = cv_hists[subchannel_name];
                            const auto &o = osc_hists[subchannel_name];
                            if(sc == 0) {
                                cv_hist = (TH1D*)h->Clone();
                                osc_hist = (TH1D*)o->Clone();
                            } else {
                                cv_hist->Add(&*h);
                                osc_hist->Add(&*o);
                            }
                            ++global_subchannel_index;
                        }
                        if(osc_scale || all_scale)
                            cv_hist->GetYaxis()->SetTitle("Events/GeV");
                        else
                            cv_hist->GetYaxis()->SetTitle("Events");
                        cv_hist->SetTitle((config.m_mode_names[im]  +" "+ config.m_detector_names[id]+" "+ config.m_channel_names[ic]).c_str());
                        cv_hist->GetXaxis()->SetTitle("");
                        cv_hist->SetLineColor(kBlack);
                        cv_hist->SetFillColor(kWhite);
                        cv_hist->SetFillStyle(0);
                        osc_hist->SetLineColor(kBlue);
                        osc_hist->SetFillColor(kWhite);
                        osc_hist->SetFillStyle(0);
                        cv_hist->SetLineWidth(3);
                        osc_hist->SetLineWidth(3);
                        TH1D *rat = (TH1D*)osc_hist->Clone();
                        rat->Divide(cv_hist);
                        rat->SetTitle("");
                        rat->GetYaxis()->SetTitle("Ratio");
                        TH1D *one = (TH1D*)rat->Clone();
                        one->Divide(one);
                        one->SetLineColor(kBlack);
                        one->GetYaxis()->SetTitle("Ratio");

                        std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.59,0.89,0.59,0.89);
                        leg->SetFillStyle(0);
                        leg->SetLineWidth(0);
                        leg->AddEntry(cv_hist, "No Oscillations", "l");
                        std::string oscstr = "";//"#splitline{Oscilations:}{";
                        for(int j=0;j<model->nparams;j++){
                            oscstr+=model->pretty_param_names[j]+ " : "+ to_string_prec(osc_params[j],2) + (j==0 ? ", " : "" );
                        }
                        //oscstr+="}";

                        leg->AddEntry(osc_hist, oscstr.c_str(), "l");

                        TPad p1("p1", "p1", 0, 0.25, 1, 1);
                        p1.SetBottomMargin(0);
                        p1.cd();
                        cv_hist->Draw("hist");
                        osc_hist->Draw("hist same");
                        leg->Draw("same");

                        TPad p2("p2", "p2", 0, 0, 1, 0.25);
                        p2.SetTopMargin(0);
                        p2.SetBottomMargin(0.3);
                        p2.cd();
                        one->GetYaxis()->SetTitleSize(0.1);
                        one->GetYaxis()->SetLabelSize(0.1);
                        one->GetXaxis()->SetTitleSize(0.1);
                        one->GetXaxis()->SetLabelSize(0.1);
                        one->GetYaxis()->SetTitleOffset(0.5);
                        one->Draw("hist");
                        one->SetMaximum(rat->GetMaximum()*1.2);
                        one->SetMinimum(rat->GetMinimum()*0.8);
                        rat->Draw("hist same");
                        
                        c.cd();
                        p1.Draw();
                        p2.Draw();

                        c.Print(filename.c_str(), "pdf");

                        delete cv_hist;
                        delete osc_hist;
                    }
                }
            }
        }

        if(*all || *cov_command) {
            std::map<std::string, std::unique_ptr<TH2D>> matrices = covarianceTH2D(systs, config, spec);
            for(const auto &[name, mat]: matrices) {
                mat->Draw("colz");
                c.Print(filename.c_str(), "pdf");
            }
        }

        if(*all || *errband_command) {
            std::map<std::string, std::unique_ptr<TH1D>> cv_hists = getCVHists(spec, config, cv_scale || all_scale);
            size_t global_subchannel_index = 0;
            size_t global_channel_index = 0;
            for(size_t im = 0; im < config.m_num_modes; im++){
                for(size_t id =0; id < config.m_num_detectors; id++){
                    for(size_t ic = 0; ic < config.m_num_channels; ic++){
                        std::unique_ptr<THStack> s = std::make_unique<THStack>((std::to_string(global_channel_index)).c_str(),(std::to_string(global_channel_index)).c_str());
                        std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.59,0.89,0.59,0.89);
                        leg->SetFillStyle(0);
                        leg->SetLineWidth(0);
                        for(size_t sc = 0; sc < config.m_num_subchannels[ic]; sc++){
                            const std::string& subchannel_name  = config.m_fullnames[global_subchannel_index];
                            s->Add(cv_hists[subchannel_name].get());
                            leg->AddEntry(cv_hists[subchannel_name].get(), config.m_subchannel_plotnames[ic][sc].c_str() ,"f");
                            ++global_subchannel_index;
                        }
                        std::unique_ptr<TGraphAsymmErrors> err_band = getErrorBand(config, prop, systs, err_scale || all_scale);
                        err_band->SetLineColor(kRed+1);                        
                        err_band->GetYaxis()->SetTitle("Events/GeV");
                        err_band->Draw("AP");
                        err_band->GetXaxis()->SetRangeUser(0.3, 3.0);


                        s->Draw("hist SAME");
                        leg->Draw("SAME");
                        err_band->SetTitle((config.m_mode_names[im]  +" "+ config.m_detector_names[id]+" "+ config.m_channel_names[ic]).c_str());
                        err_band->GetXaxis()->SetTitle(config.m_channel_units[ic].c_str());
                         TH1* dummy = new TH1F("", "", 1, 0, 1);
                         dummy->SetLineColor(kRed+1);
                        leg->AddEntry(dummy->Clone(), "Syst", "l");
                        err_band->Draw("SAME P");
                        c.Print(filename.c_str(), "pdf");
                    }
                }
            }
        }

        if((*all && with_splines) || *spline_command) {
            std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> spline_graphs = getSplineGraphs(systs, config);
            c.Clear();
            c.Divide(4,4);
            for(const auto &[syst_name, syst_bins]: spline_graphs) {
                int bin = 0;
                bool unprinted = true;
                for(const auto &[fixed_pts, curve]: syst_bins) {
                    unprinted = true;
                    c.cd(bin%16+1);
                    fixed_pts->SetMarkerColor(kBlack);
                    fixed_pts->SetMarkerStyle(kFullCircle);
                    fixed_pts->GetXaxis()->SetTitle("#sigma");
                    fixed_pts->GetYaxis()->SetTitle("Weight");
                    fixed_pts->SetTitle((syst_name+" - True Bin "+std::to_string(bin)).c_str());
                    fixed_pts->Draw("PA");
                    curve->Draw("C same");
                    ++bin;
                    if(bin % 16 == 0) {
                        c.Print(filename.c_str(), "pdf");
                        unprinted = false;
                    }
                }
                if(unprinted)
                    c.Print(filename.c_str(), "pdf");
            }
        }

        c.Print((filename + "]").c_str(), "pdf");
    } else {
        TFile fout(filename.c_str(), "RECREATE");

        if(*all || *cv_command) {
            std::map<std::string, std::unique_ptr<TH1D>> cv_hists = getCVHists(spec, config, cv_scale || all_scale);
            fout.mkdir("CV_hists");
            fout.cd("CV_hists");
            for(const auto &[name, hist]: cv_hists) {
                hist->Write(name.c_str());
            }
        }

        if((*all && osc_params.size()) || *osc_command) {
            PROspec osc_spec = FillRecoSpectra(config, prop, systs, *model, pparams, !eventbyevent);
            std::map<std::string, std::unique_ptr<TH1D>> osc_hists = getCVHists(osc_spec, config, osc_scale || all_scale);
            fout.mkdir("Osc_hists");
            fout.cd("Osc_hists");
            for(const auto &[name, hist]: osc_hists) {
                hist->Write(name.c_str());
            }
        }

        if(*all || *cov_command) {
            std::map<std::string, std::unique_ptr<TH2D>> matrices = covarianceTH2D(systs, config, spec);
            fout.mkdir("Covariance");
            fout.cd("Covariance");
            for(const auto &[name, mat]: matrices)
                mat->Write(name.c_str());
        }

        if(*all || *errband_command) {
            std::unique_ptr<TGraphAsymmErrors> err_band = getErrorBand(config, prop, systs, err_scale || all_scale);
            fout.mkdir("ErrorBand");
            fout.cd("ErrorBand");
            err_band->Write("err_band");

        }
        if((*all && with_splines) || *spline_command) {
            std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> spline_graphs = getSplineGraphs(systs, config);
            fout.mkdir("Splines");
            fout.cd("Splines");
            for(const auto &[name, syst_splines]: spline_graphs) {
                size_t bin = 0;
                for(const auto &[fixed_pts, curve]: syst_splines) {
                    fixed_pts->Write((name+"_fixedpts_"+std::to_string(bin)).c_str());
                    curve->Write((name+"_curve_"+std::to_string(bin)).c_str());
                    bin++;
                }
            }
        }
    }
    
    return 0;
}

std::map<std::string, std::unique_ptr<TH1D>> getCVHists(const PROspec &spec, const PROconfig& inconfig, bool scale) {
    std::map<std::string, std::unique_ptr<TH1D>> hists;  

    size_t global_subchannel_index = 0;
    size_t global_channel_index = 0;
    for(size_t im = 0; im < inconfig.m_num_modes; im++){
        for(size_t id =0; id < inconfig.m_num_detectors; id++){
            for(size_t ic = 0; ic < inconfig.m_num_channels; ic++){
                for(size_t sc = 0; sc < inconfig.m_num_subchannels[ic]; sc++){
                    const std::string& subchannel_name  = inconfig.m_fullnames[global_subchannel_index];
                    const std::string& color = inconfig.m_subchannel_colors[ic][sc];
                    int rcolor = color == "NONE" ? kRed - 7 : inconfig.HexToROOTColor(color);
                    std::unique_ptr<TH1D> htmp = std::make_unique<TH1D>(spec.toTH1D(inconfig, global_subchannel_index));
                    htmp->SetLineWidth(1);
                    htmp->SetLineColor(kBlack);
                    htmp->SetFillColor(rcolor);
                    if(scale) htmp->Scale(1,"width");
                    hists[subchannel_name] = std::move(htmp);

                    log<LOG_DEBUG>(L"%1% || Printot %2% %3% %4% %5% %6% : Integral %7% ") % __func__ % global_channel_index % global_subchannel_index % subchannel_name.c_str() % sc % ic % hists[subchannel_name]->Integral();
                    ++global_subchannel_index;
                }//end subchan
                ++global_channel_index;
            }//end chan
        }//end det
    }//end mode
    return hists;
}

std::map<std::string, std::unique_ptr<TH2D>> covarianceTH2D(const PROsyst &syst, const PROconfig &config, const PROspec &cv) {
    std::map<std::string, std::unique_ptr<TH2D>> ret;
    Eigen::MatrixXf fractional_cov = syst.fractional_covariance;
    Eigen::MatrixXf diag = cv.Spec().array().matrix().asDiagonal(); 
    Eigen::MatrixXf full_covariance =  diag*fractional_cov*diag;
    Eigen::MatrixXf collapsed_full_covariance =  CollapseMatrix(config,full_covariance);  
    Eigen::VectorXf collapsed_cv = CollapseMatrix(config, cv.Spec());
    Eigen::MatrixXf collapsed_cv_inv_diag = collapsed_cv.asDiagonal().inverse();
    Eigen::MatrixXf collapsed_frac_cov = collapsed_cv_inv_diag * collapsed_full_covariance * collapsed_cv_inv_diag;

    std::unique_ptr<TH2D> cov_hist = std::make_unique<TH2D>("cov", "Fractional Covariance Matrix;Bin # ;Bin #", config.m_num_bins_total, 0, config.m_num_bins_total, config.m_num_bins_total, 0, config.m_num_bins_total);
    std::unique_ptr<TH2D> collapsed_cov_hist = std::make_unique<TH2D>("ccov", "Collapsed Fractional Covariance Matrix;Bin # ;Bin #", config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed, config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed);

    std::unique_ptr<TH2D> cor_hist = std::make_unique<TH2D>("cor", "Correlation Matrix;Bin # ;Bin #", config.m_num_bins_total, 0, config.m_num_bins_total, config.m_num_bins_total, 0, config.m_num_bins_total);
    std::unique_ptr<TH2D> collapsed_cor_hist = std::make_unique<TH2D>("ccor", "Collapsed Correlation Matrix;Bin # ;Bin #", config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed, config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed);

    for(size_t i = 0; i < config.m_num_bins_total; ++i)
        for(size_t j = 0; j < config.m_num_bins_total; ++j){
            cov_hist->SetBinContent(i+1,j+1,fractional_cov(i,j));
            cor_hist->SetBinContent(i+1,j+1,fractional_cov(i,j)/(sqrt(fractional_cov(i,i))*sqrt(fractional_cov(j,j))));
        }

    for(size_t i = 0; i < config.m_num_bins_total_collapsed; ++i)
        for(size_t j = 0; j < config.m_num_bins_total_collapsed; ++j){
            collapsed_cov_hist->SetBinContent(i+1,j+1,collapsed_frac_cov(i,j));
            collapsed_cor_hist->SetBinContent(i+1,j+1,collapsed_frac_cov(i,j)/(sqrt(collapsed_frac_cov(i,i))*sqrt(collapsed_frac_cov(j,j))));
        }

    ret["total_frac_cov"] = std::move(cov_hist);
    ret["collapsed_total_frac_cov"] = std::move(collapsed_cov_hist);
    ret["total_cor"] = std::move(cor_hist);
    ret["collapsed_total_cor"] = std::move(collapsed_cor_hist);

    for(const auto &name: syst.covar_names) {
        const Eigen::MatrixXf &covar = syst.GrabMatrix(name);
        const Eigen::MatrixXf &corr = syst.GrabCorrMatrix(name);

        std::unique_ptr<TH2D> cov_h = std::make_unique<TH2D>(("cov"+name).c_str(), (name+" Fractional Covariance;Bin # ;Bin #").c_str(), config.m_num_bins_total, 0, config.m_num_bins_total, config.m_num_bins_total, 0, config.m_num_bins_total);
        std::unique_ptr<TH2D> corr_h = std::make_unique<TH2D>(("cor"+name).c_str(), (name+" Correlation;Bin # ;Bin #").c_str(), config.m_num_bins_total, 0, config.m_num_bins_total, config.m_num_bins_total, 0, config.m_num_bins_total);
        for(size_t i = 0; i < config.m_num_bins_total; ++i){
            for(size_t j = 0; j < config.m_num_bins_total; ++j){
                cov_h->SetBinContent(i+1,j+1,covar(i,j));
                corr_h->SetBinContent(i+1,j+1,corr(i,j));
            }
        }

        ret[name+"_cov"] = std::move(cov_h);
        ret[name+"_corr"] = std::move(corr_h);
    }

    return ret;
}

std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> 
getSplineGraphs(const PROsyst &systs, const PROconfig &config) {
    std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> spline_graphs;

    for(size_t i = 0; i < systs.GetNSplines(); ++i) {
        const std::string &name = systs.spline_names[i];
        const PROsyst::Spline &spline = systs.GrabSpline(name);
        //using Spline = std::vector<std::vector<std::pair<float, std::array<float, 4>>>>;
        std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>> bin_graphs;
        bin_graphs.reserve(config.m_num_truebins_total);
        
        for(size_t j = 0; j < config.m_num_truebins_total; ++j) {
            const std::vector<std::pair<float, std::array<float, 4>>> &spline_for_bin = spline[j];
            std::unique_ptr<TGraph> curve = std::make_unique<TGraph>();
            std::unique_ptr<TGraph> fixed_pts = std::make_unique<TGraph>();
            for(size_t k = 0; k < spline_for_bin.size(); ++k) {
                //const auto &[lo, coeffs] = spline_for_bin[k];
                float lo = spline_for_bin[k].first;
                std::array<float, 4> coeffs = spline_for_bin[k].second;
                float hi = k < spline_for_bin.size() - 1 ? spline_for_bin[k+1].first : systs.spline_hi[i];
                auto fn = [coeffs](float shift){
                    return coeffs[0] + coeffs[1]*shift + coeffs[2]*shift*shift + coeffs[3]*shift*shift*shift;
                };
                fixed_pts->SetPoint(fixed_pts->GetN(), lo, fn(0)); 
                if(k == spline_for_bin.size() - 1)
                    fixed_pts->SetPoint(fixed_pts->GetN(), hi, fn(hi - lo));
                float width = (hi - lo) / 20;
                for(size_t l = 0; l < 20; ++l)
                    curve->SetPoint(curve->GetN(), lo + l * width, fn(l * width));
            }
            bin_graphs.push_back(std::make_pair(std::move(fixed_pts), std::move(curve)));
        }
        spline_graphs[name] = std::move(bin_graphs);
    }

    return spline_graphs;
}

std::unique_ptr<TGraphAsymmErrors> getErrorBand(const PROconfig &config, const PROpeller &prop, const PROsyst &syst, bool scale) {
    //TODO: Only works with 1 mode/detector/channel
    Eigen::VectorXf cv = CollapseMatrix(config, FillCVSpectrum(config, prop, true).Spec());
    std::vector<float> edges = config.GetChannelBinEdges(0);
    std::vector<float> centers;
    for(size_t i = 0; i < edges.size() - 1; ++i)
        centers.push_back((edges[i+1] + edges[i])/2);
    std::vector<Eigen::VectorXf> specs;
    for(size_t i = 0; i < 1000; ++i)
        specs.push_back(FillSystRandomThrow(config, prop, syst).Spec());
        //specs.push_back(CollapseMatrix(config, FillSystRandomThrow(config, prop, syst).Spec()));
    TH1D tmphist("th", "", cv.size(), edges.data());
    for(size_t i = 0; i < cv.size(); ++i)
        tmphist.SetBinContent(i+1, cv(i));
    if(scale) tmphist.Scale(1, "width");
    //std::unique_ptr<TGraphAsymmErrors> ret = std::make_unique<TGraphAsymmErrors>(cv.size(), centers.data(), cv.data());
    std::unique_ptr<TGraphAsymmErrors> ret = std::make_unique<TGraphAsymmErrors>(&tmphist);
    for(size_t i = 0; i < cv.size(); ++i) {
        std::array<float, 1000> binconts;
        for(size_t j = 0; j < 1000; ++j) {
            binconts[j] = specs[j](i);
        }
        float scale_factor = tmphist.GetBinContent(i+1)/cv(i);
        std::sort(binconts.begin(), binconts.end());
        float ehi = std::abs((binconts[840] - cv(i))*scale_factor);
        float elo = std::abs((cv(i) - binconts[160])*scale_factor);
        ret->SetPointEYhigh(i, ehi);
        ret->SetPointEYlow(i, elo);
    }
    return ret;
}

