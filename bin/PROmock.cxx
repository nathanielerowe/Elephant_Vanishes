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
#include "PROfitter.h"
#include "PROmodel.h"

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
#include "TGraph.h"
#include "TStyle.h"
#include "TFile.h"
#include "TRatioPlot.h"
#include "TGraphAsymmErrors.h"

using namespace PROfit;

//Run a single fit at a given point in parameter space (nominally no osc) with mock data
log_level_t GLOBAL_LEVEL = LOG_ERROR;

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);
  CLI::App app{"Test for PROfit"}; 
  
  // Define options
  std::string xmlname = "NULL.xml"; 
  int maxevents = 50000;
  size_t  nthread = 1;
  //Define a filename to save chisq values in
  std::string filename;
  bool binned=false;
  bool plotonly=false;
  bool floatosc=false;
  std::vector<int> grid_size;
  std::vector<std::string> mockparams;
  std::vector<float> mockshifts;
  std::string reweights_file;
  std::vector<std::string> mockreweights;
  std::vector<TH2D*> weighthists;
  std::vector<float> physics_params_in;
  
  app.add_option("-x, --xml",       xmlname, "Input PROfit XML config.");
  app.add_option("-m, --max",       maxevents, "Max number of events to run over.");
  app.add_option("-v, --verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].")->default_val(2);
  app.add_option("-t, --nthread",   nthread, "Number of fits.")->default_val(1);
  app.add_option("-o, --outfile",   filename, "If you want chisq to be dumped to text file, provide name")->default_str("profit");
  app.add_option("-s, --mocksys",   mockparams, "Vector of systematics parameter names to vary for mock data");
  app.add_option("-u, --mockvals",  mockshifts, "Vector of size of shifts.");
  app.add_option("-f, --rwfile", reweights_file, "File containing histograms for reweighting");
  app.add_option("-r, --mockrw",   mockreweights, "Vector of reweights to use for mock data");
  app.add_option("-p, --pparams",   physics_params_in, "deltam^2, sin^22thetamumu, default no osc");
  app.add_flag("-q, --plotonly", plotonly, "Skip the fit and just produce the plots");
  app.add_flag("-c, --floatosc", floatosc, "Let oscillation parameters float in the fit");

  app.add_flag(  "-b, --binned",    binned, "Do you want to weight event-by-event?");

  CLI11_PARSE(app, argc, argv);


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
  log<LOG_DEBUG>(L"%1% || model name %2% ") % __func__ % config.m_model_tag.c_str();
  std::unique_ptr<PROmodel> model = get_model_from_string(config.m_model_tag, prop);
  log<LOG_DEBUG>(L"%1% || model size %2% by %3% ") % __func__ % model->hists[0].rows() % model->hists[0].cols();

  //Convert to eigen and to log
  if (physics_params_in.size() < model->nparams) {
    for (size_t i=physics_params_in.size(); i<model->nparams; ++i) {
      physics_params_in.push_back(0.0);
      log<LOG_WARNING>(L"%1% || Model expects %2% params but only %3% params are given") % __func__ % model->nparams % physics_params_in.size();
    }
  }
  std::vector<float> physics_params_log;
  for (size_t i=0; i<physics_params_in.size(); i++) {
    physics_params_log.push_back(std::log10(physics_params_in[i]));
  }

  Eigen::VectorXf physics_params = Eigen::VectorXf::Map(physics_params_log.data(), physics_params_log.size());

  for (size_t i = 0; i < mockparams.size(); ++i) {
    log<LOG_INFO>(L"%1% || Mock data parameters %2%") % __func__  % mockparams[i].c_str();
  }
  for (size_t i = 0; i < mockreweights.size(); ++i) {
    log<LOG_INFO>(L"%1% || Mock reweight %2%") % __func__  % mockreweights[i].c_str();
  }
  
  PROspec data, cv;

  cv = FillCVSpectrum(config, prop, binned);
  if (mockparams.empty() && mockreweights.empty()) {
    log<LOG_INFO>(L"%1% || Will use CV MC (with any requested oscillations) as data for this study") % __func__  ;
    data = physics_params_in[0] != 0 && physics_params_in[1] != 0 ? 
      FillRecoSpectra(config, prop, systs, *model, physics_params, binned) :
      FillCVSpectrum(config, prop, binned);
  }
  else if (!mockreweights.empty()) {
    log<LOG_INFO>(L"%1% || Will use reweighted MC (with any requested oscillations) as data for this study") % __func__  ;
        log<LOG_INFO>(L"%1% || Any parameter shifts requested will be ignored (fix later?)") % __func__  ;
    auto file = std::make_unique<TFile>(reweights_file.c_str());
    log<LOG_DEBUG>(L"%1% || Set file to : %2% ") % __func__ % reweights_file.c_str();
    log<LOG_DEBUG>(L"%1% || Size of reweights vector : %2% ") % __func__ % mockreweights.size() ;
    for (size_t i=0; i < mockreweights.size(); ++i) {
      log<LOG_DEBUG>(L"%1% || Mock reweight i : %2% ") % __func__ % mockreweights[i].c_str() ;
      TH2D* rwhist = (TH2D*)file->Get(mockreweights[i].c_str());
      weighthists.push_back(rwhist);
      log<LOG_DEBUG>(L"%1% || Read in weight hist ") % __func__ ;      
    }
    data = FillWeightedSpectrumFromHist(config,prop,weighthists,*model,physics_params,binned);
  }
  else{
    if (mockshifts.size() != mockparams.size()) {
      log<LOG_INFO>(L"%1% || No shifts provided, setting all parameter shifts to 1 sigma") % __func__  ;
      for (size_t i = 0; i < mockparams.size(); ++i) {
	mockshifts.push_back(1.0);
      }
    }
    log<LOG_INFO>(L"%1% || Will use systematic shifted MC (no osc) as data for this study") % __func__  ;    
    data = systs.GetSplineShiftedSpectrum(config,prop,mockparams,mockshifts);
  }

  //stupid hack, must be a better way to do this
  //Set up options:
  std::vector<const char*> xlabel(4);
  xlabel[0] = "Reconstructed Neutrino Energy";
  xlabel[1] = "True Leading Proton Momentum";
  xlabel[2] = "True Leading Proton Cos(Theta)";
  xlabel[3] = "Check what variable you are plotting!";
  int xi;
  if (xmlname.find("standard") != std::string::npos) {
    xi = 0;
  }
  else if (xmlname.find("pmom") != std::string::npos) {
    xi = 1;
  }
  else if (xmlname.find("costh") != std::string::npos) {  
    xi = 2;
  }
  else {
    xi = 3;
  }

  TH1D hcv = cv.toTH1D(config,0);
  TH1D hmock = data.toTH1D(config,0);
  hcv.Scale(1, "width");
  hmock.Scale(1, "width");
  hcv.GetYaxis()->SetTitle("Events/GeV");
  hmock.GetYaxis()->SetTitle("Events/GeV");
  hcv.GetXaxis()->SetTitle(xlabel[xi]);
  hmock.GetXaxis()->SetTitle(xlabel[xi]);
  hcv.SetTitle("");
  hmock.SetTitle("");

  TCanvas *c = new TCanvas((filename+"_spec_cv").c_str(), (filename+"_spec_cv").c_str(), 800, 800);
  hcv.SetLineColor(kBlack);
  hmock.SetLineColor(5);
  hmock.SetFillColor(5);
  TRatioPlot * rp = new TRatioPlot(&hmock,&hcv);
  rp->Draw();
  rp->GetLowerRefGraph()->SetMarkerStyle(21);
  TGraphAsymmErrors *lowerGraph = dynamic_cast<TGraphAsymmErrors*>(rp->GetLowerRefGraph());
  if (lowerGraph) {
    int nPoints = lowerGraph->GetN();
    for (int i = 0; i < nPoints; i++) {
      lowerGraph->SetPointError(i, 0, 0, 0, 0); // Set both x and y errors to zero
    }
  }
  std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.35,0.7,0.89,0.89);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->AddEntry(&hcv,"CV","l");
  leg->AddEntry(&hmock,"Mock data: ", "f");
  TObject *null = new TObject(); 
  int i=0;
  for (const auto& m : mockparams) {
    char ns[6];
    snprintf(ns, sizeof(ns),"%.2f", mockshifts[i]);
    leg->AddEntry(null, (m+": "+ns+ " sigma").c_str(),"");
    i++;
  }
  for (const auto& m : mockreweights) {
      leg->AddEntry(null, m.c_str(),"");
    i++;
  }
  for (const auto& m : physics_params_in) {
    leg->AddEntry(null, ("param: "+std::to_string(m)).c_str(),"");
    i++;
  }
  
  leg->Draw();
  c->SaveAs((filename+"_spec.pdf").c_str());

  if (!plotonly) {
    Eigen::VectorXf data_vec = CollapseMatrix(config, data.Spec());
    Eigen::VectorXf err_vec_sq = data.Error().array().square();
    Eigen::VectorXf err_vec = CollapseMatrix(config, err_vec_sq).array().sqrt();
    data = PROspec(data_vec, err_vec);

    
    //Define a metric
    PROchi chi("", config, prop, &systs, *model, data, PROfit::PROchi::BinnedChi2);
    log<LOG_DEBUG>(L"%1% || Run fit to get inputs for PROfile") % __func__ ;
    
    //Run a fit to get best-fit input for PROfile
    LBFGSpp::LBFGSBParam<float> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    param.max_linesearch = 250;
    param.delta = 1e-6;

    size_t nparams = model->nparams + systs.GetNSplines();
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

    log<LOG_DEBUG>(L"%1% || Fit finished, now running PROfile") % __func__ ;    

    PROfile(config, prop, systs, *model, data, chi, param, filename, floatosc, nthread, best_fit, physics_params);
  }
  
  return 0;
}
