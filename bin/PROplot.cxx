#include "PROconfig.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROcreate.h"
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

#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace PROfit;

log_level_t GLOBAL_LEVEL;

std::map<std::string, std::unique_ptr<TH1D>> getCVHists(const PROspec & spec, const PROconfig& inconfig);
std::unique_ptr<TH2D> covarianceTH2D(const PROsyst &syst, const PROconfig &config);
std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> getSplineGraphs(const PROsyst &systs, const PROconfig &config);

int main(int argc, char* argv[])
{
    gStyle->SetOptStat(0);
    CLI::App app{"PROplot, plotting spectra using ROOT"}; 

    // Define options
    std::string xmlname, filename; 
    std::array<float, 2> apply_osc{0, 0};
    std::map<std::string, float> apply_shift;
    std::vector<std::string> syst_list, systs_excluded;
    bool eventbyevent = false;

    bool covariance_plot = false, spline_plots = false, err_band_plot = false;

    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.")->required();
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].")->default_val(LOG_WARNING);
    app.add_option("-o,--output", filename, "Output Filename.")->default_val("PROplot.pdf");
    app.add_option("--apply-osc", apply_osc, "Apply oscillations with these paramters to the CV spectrum.")->default_str("0 0");
    app.add_option("--apply-shift", apply_shift, "Apply these shifts to the CV spectrum.");
    app.add_option("--syst-list", syst_list, "Override list of systematics to use (note: all systs must be in the xml).");
    app.add_option("--exclude-systs", systs_excluded, "List of systematics to exclude.")->excludes("--syst-list"); 

    app.add_flag("--event-by-event",    eventbyevent, "Do you want to weight event-by-event?");

    CLI::Option *all_flag = app.add_flag("-a,--all", "Make all plots");
    app.add_flag("--covariance", covariance_plot, "Make TH2 representation of total fractional covariance matrix.")->excludes(all_flag);
    app.add_flag("--splines", spline_plots, "Make TGraphs of the splines in each true bin")->excludes(all_flag);
    app.add_flag("--error-band", err_band_plot, "Plot an error band on the CV. If saving to root file, will save as a TGraphAsymmErrors")->excludes(all_flag);

    CLI11_PARSE(app, argc, argv);

    bool savetoroot = filename.size() > 5 && filename.substr(filename.size() - 5) == ".root";
    bool savetopdf = filename.size() > 4 && filename.substr(filename.size() - 4) == ".pdf";
    if(!savetoroot && !savetopdf) {
        log<LOG_ERROR>(L"%1% || Output file specified, %2%, is not a pdf or root file.") % __func__ % filename.c_str();
        return 1;
    }

    //Initilize configuration from the XML;
    PROconfig config(xmlname);

    PROpeller prop;
    std::vector<SystStruct> systsstructs;
    PROcess_CAFAna(config, systsstructs, prop);

    PROsyst systs(systsstructs);
    PROsc osc(prop);

    std::vector<float> pparams = {std::log10(apply_osc[0]), std::log10(apply_osc[1])};
    std::cout << "Injected point: sinsq2t = " << apply_osc[0] << " dmsq = " << apply_osc[1] << std::endl;
    for(const auto& [name, shift]: apply_shift)
      std::cout << "Injected syst: " << name << " shifted by " << shift << " sigma\n";

    PROspec spec = 
        apply_osc[0] != 0 && apply_osc[1] != 0 ? 
            FillRecoSpectra(config, prop, systs, &osc, apply_shift, pparams, !eventbyevent) :
        apply_shift.size() ? 
            FillRecoSpectra(config, prop, systs, apply_shift, !eventbyevent) :
            FillCVSpectrum(config, prop, !eventbyevent);

    if(syst_list.size()) {
      systs = systs.subset(syst_list);
    } else if(systs_excluded.size()) {
      systs = systs.excluding(systs_excluded);
    }
   
    std::map<std::string, std::unique_ptr<TH1D>> cv_hists = getCVHists(spec, config);
    std::unique_ptr<TH2D> covariance = covarianceTH2D(systs, config);
    std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> spline_graphs = getSplineGraphs(systs, config);

    if(savetopdf) {
        TCanvas c;
        
        c.Print((filename + "[").c_str(), "pdf");
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
                    s->GetYaxis()->SetTitle("Events/GeV");

                    c.Print(filename.c_str(), "pdf");
                }
            }
        }

        if(*all_flag || covariance_plot) {
            covariance->Draw("colz");
            c.Print(filename.c_str(), "pdf");
        }

        if(*all_flag || spline_plots) {
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

        fout.mkdir("CV_hists");
        fout.cd("CV_hists");
        for(const auto &[name, hist]: cv_hists) {
            hist->Write(name.c_str());
        }

        if(*all_flag || covariance_plot) {
            fout.mkdir("Covariance");
            fout.cd("Covariance");
            covariance->Write("collapsed_frac_cov");
        }

        if(*all_flag || spline_plots) {
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

std::map<std::string, std::unique_ptr<TH1D>> getCVHists(const PROspec &spec, const PROconfig& inconfig) {
    bool div_bin = true;

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
                    if(div_bin) htmp->Scale(1,"width");
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

std::unique_ptr<TH2D> covarianceTH2D(const PROsyst &syst, const PROconfig &config) {
    Eigen::MatrixXd fractional_cov = CollapseMatrix(config, syst.fractional_covariance);

    std::unique_ptr<TH2D> cov_hist = std::make_unique<TH2D>("cov", "Fractional Covariance Matrix;Bin #;Bin #", config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed, config.m_num_bins_total_collapsed, 0, config.m_num_bins_total_collapsed);

    for(size_t i = 0; i < config.m_num_bins_total_collapsed; ++i)
        for(size_t j = 0; j < config.m_num_bins_total_collapsed; ++j)
            cov_hist->SetBinContent(i+1,j+1,fractional_cov(i,j));

    return cov_hist;
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
                const auto &[lo, coeffs] = spline_for_bin[k];
                float hi = k < spline_for_bin.size() - 1 ? spline_for_bin[k].first : systs.spline_hi[i];
                auto fn = [coeffs](float shift){
                    return coeffs[0] + coeffs[1]*shift + coeffs[2]*shift*shift + coeffs[3]*shift*shift*shift;
                };
                fixed_pts->SetPoint(fixed_pts->GetN(), lo, fn(0)); 
                if(k == spline_for_bin.size() - 1)
                    fixed_pts->SetPoint(fixed_pts->GetN(), hi, fn(hi - lo));
                float width = (hi - lo) / 10;
                for(size_t l = 0; l < 10; ++l)
                    curve->SetPoint(curve->GetN(), lo + l * width, fn(l * width));
            }
            bin_graphs.push_back(std::make_pair(std::move(fixed_pts), std::move(curve)));
        }
        spline_graphs[name] = std::move(bin_graphs);
    }

    return spline_graphs;
}

