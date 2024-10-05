#include "PROsurf.h"

using namespace PROfit;

std::vector<std::vector<double>> latin_hypercube_sampling(size_t num_samples, size_t dimensions, std::uniform_real_distribution<float>&dis, std::mt19937 &gen) {
    std::vector<std::vector<double>> samples(num_samples, std::vector<double>(dimensions));

    for (size_t d = 0; d < dimensions; ++d) {

        std::vector<double> perm(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            perm[i] = (i + dis(gen)) / num_samples;  
        }
        std::shuffle(perm.begin(), perm.end(), gen);  
        for (size_t i = 0; i < num_samples; ++i) {
            samples[i][d] = perm[i]; 
        }
    }

    return samples;
}

std::vector<int> sorted_indices(const std::vector<double>& vec) {
    std::vector<int> indices(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&vec](int i1, int i2) { return vec[i1] < vec[i2]; });
    return indices;
}



PROsurf::PROsurf(size_t nbinsx, LogLin llx, double x_lo, double x_hi, size_t nbinsy, LogLin lly, double y_lo, double y_hi) : nbinsx(nbinsx), nbinsy(nbinsy), edges_x(Eigen::VectorXd::Constant(nbinsx + 1, 0)), edges_y(Eigen::VectorXd::Constant(nbinsy + 1, 0)), surface(nbinsx, nbinsy) {
    if(llx == LogAxis) {
        x_lo = std::log10(x_lo);
        x_hi = std::log10(x_hi);
    }
    if(lly == LogAxis) {
        y_lo = std::log10(y_lo);
        y_hi = std::log10(y_hi);
    }
    for(size_t i = 0; i < nbinsx + 1; i++)
        edges_x(i) = x_lo + i * (x_hi - x_lo) / nbinsx;
    for(size_t i = 0; i < nbinsy + 1; i++)
        edges_y(i) = y_lo + i * (y_hi - y_lo) / nbinsy;
}

void PROsurf::FillSurfaceSimple(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, bool binned_weighting) {
    std::random_device rd{};
    std::mt19937 rng{rd()};
    std::normal_distribution<float> d;

    std::ofstream chi_file;
    PROchi::EvalStrategy strat = binned_weighting ? PROchi::BinnedChi2 : PROchi::EventByEvent;

    if(!filename.empty()){
        chi_file.open(filename);
    }

    std::cout << "In FillSurfaceFast\n";

    for(size_t i = 0; i < nbinsx; i++) {
        for(size_t j = 0; j < nbinsy; j++) {
            std::cout << "Filling point " << i << " " << j << std::endl;
            LBFGSpp::LBFGSBParam<double> param;  
            param.epsilon = 1e-6;
            param.max_iterations = 100;
            param.max_linesearch = 50;
            param.delta = 1e-6;
            LBFGSpp::LBFGSBSolver<double> solver(param); 
            size_t nparams = systs.GetNSplines();
            std::vector<float> physics_params = {(float)edges_y(j), (float)edges_x(i)};//deltam^2, sin^22thetamumu
            PROchi chi("3plus1",&config,&prop,&systs,&osc, data, nparams, systs.GetNSplines(), strat, physics_params);
            Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
            Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
            Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);

            double fx;
            int niter=-1;
            std::vector<double> chi2s;
            int nfit = 0;
            do {
                nfit++;
                for(size_t i = 0; i < nparams; ++i)
                    x(i) = 0.3*d(rng);
                // x will be overwritten to be the best point found
                try {
                    niter = solver.minimize(chi, x, fx, lb, ub);
                } catch(std::runtime_error &except) {
                    log<LOG_ERROR>(L"%1% || Fit failed on iter %2%,  %3%") % __func__ % niter % except.what();
                    continue;
                }
                chi2s.push_back(fx);
            } while(chi2s.size() < 10 && nfit < 100);
            fx = *std::min_element(chi2s.begin(), chi2s.end());
            surface(i, j) = fx;
            if(!filename.empty()){
                chi_file<<"\n"<<edges_x(i)<<" "<<edges_y(j)<<" "<<fx<<std::flush;
            }
        }
    }
}

std::vector<surfOut> PROsurf::PointHelper(const PROconfig *config, const PROpeller *prop, const PROsyst *systs, const PROsc *osc, const PROspec *data, std::vector<surfOut> multi_physics_params, PROchi::EvalStrategy strat, bool binned_weighting, int start, int end){

    std::random_device rd{};
    std::mt19937 rng{rd()};
    std::normal_distribution<float> d;
    std::uniform_real_distribution<float> d_uni(-2.0, 2.0);
    strat = binned_weighting ? PROchi::BinnedChi2 : PROchi::EventByEvent;

    std::vector<surfOut> outs;

    for(int i=start; i<end;i++){

        surfOut output;
        std::vector<float> physics_params = multi_physics_params[i].grid_val;
        output.grid_val = physics_params;
        output.grid_index = multi_physics_params[i].grid_index;

        LBFGSpp::LBFGSBParam<double> param;
        param.epsilon = 1e-6;
        param.max_iterations = 100;
        param.max_linesearch = 50;
        param.delta = 1e-6;

        LBFGSpp::LBFGSBSolver<double> solver(param);
        int nparams = systs->GetNSplines();
        PROchi chi("3plus1",config,prop,systs,osc,*data, nparams, systs->GetNSplines(), strat, physics_params);

        Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
        Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
        Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);
        Eigen::VectorXd grad = Eigen::VectorXd::Constant(nparams, 0.0);
        Eigen::VectorXd bestx = Eigen::VectorXd::Constant(nparams, 0.0);

        //First do 100 simple function calls suing LATIN hypercube setup
        double fx;
        int niter =-1;
        int N_multistart = 100;
        std::vector<double> chi2s_multistart;
        std::vector<std::vector<double>> latin_samples = latin_hypercube_sampling(N_multistart, nparams,d_uni,rng);


        log<LOG_INFO>(L"%1% || Starting MultiGlobal runs : %2%") % __func__ % N_multistart ;
        for(int s=0; s<N_multistart; s++){


            x = Eigen::Map<Eigen::VectorXd>(latin_samples[s].data(), latin_samples[s].size());
            fx =  chi(x,grad,false);
            chi2s_multistart.push_back(fx);

        }
        //Sort so we can take the best N_localfits for further zoning
        std::vector<int> best_multistart = sorted_indices(chi2s_multistart);    

        log<LOG_INFO>(L"%1% || Ending MultiGlobal Best two are : %2% and %3%") % __func__ % chi2s_multistart[best_multistart[0]] %   chi2s_multistart[best_multistart[1]];
        log<LOG_INFO>(L"%1% || Best Points is  : %2% ") % __func__ % latin_samples[best_multistart[0]];

        int N_localfits = 5;
        std::vector<double> chi2s_localfits;
        double chimin = 9999999;

        log<LOG_INFO>(L"%1% || Starting Local Gradients runs : %2%") % __func__ % N_localfits ;
        for(int s=0; s<N_localfits; s++){
            //Get the nth
            x = Eigen::Map<Eigen::VectorXd>( latin_samples[best_multistart[s]].data(), latin_samples[best_multistart[s]].size());   
            try {
                niter = solver.minimize(chi, x, fx, lb, ub);
            } catch(std::runtime_error &except) {
                log<LOG_ERROR>(L"%1% || Fit failed on niter,%2% : %3%") % __func__ % niter % except.what();
            }
            chi2s_localfits.push_back(fx);
            if(fx<chimin){
                bestx = x;
                chimin=fx;
            }
            log<LOG_INFO>(L"%1% ||  LocalGrad Run : %2% has a chi %3%") % __func__ % s % fx;
            std::string spec_string = "";
            for(auto &f : x) spec_string+=" "+std::to_string(f); 
            log<LOG_INFO>(L"%1% || Best Point is  : %2% ") % __func__ % spec_string.c_str();


        }


        // and do CV
        log<LOG_INFO>(L"%1% || Starting CV fit ") % __func__  ;
        try {
            x = Eigen::VectorXd::Constant(nparams, 0.012);
            niter = solver.minimize(chi, x, fx, lb, ub);
        } catch(std::runtime_error &except) {
            log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
        }
        chi2s_localfits.push_back(fx);
        if(fx<chimin){
            bestx = x;
            chimin=fx;
        }

        log<LOG_INFO>(L"%1% ||  CV Run has a chi %2%") % __func__ %  fx;
        std::string spec_string = "";
        for(auto &f : x) spec_string+=" "+std::to_string(f); 
        log<LOG_INFO>(L"%1% || Best Point post CV is  : %2% ") % __func__ % spec_string.c_str();



        log<LOG_INFO>(L"%1% || FINAL has a chi %2%") % __func__ %  chimin;
        spec_string = "";
        for(auto &f : bestx) spec_string+=" "+std::to_string(f); 
        log<LOG_INFO>(L"%1% || FINAL is  : %2% ") % __func__ % spec_string.c_str();

        output.chi = chimin;
        outs.push_back(output);

    }

    return outs;

}


void PROsurf::FillSurface(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename, bool binned_weighting, int nThreads) {
    std::random_device rd{};
    std::mt19937 rng{rd()};
    std::normal_distribution<float> d;
    std::uniform_real_distribution<float> d_uni(-2.0, 2.0);
    PROchi::EvalStrategy strat = binned_weighting ? PROchi::BinnedChi2 : PROchi::EventByEvent;

    int nparams = systs.GetNSplines();
    Eigen::VectorXd lastx = Eigen::VectorXd::Constant(nparams, 0.01);

    std::ofstream chi_file;

    if(!filename.empty()){
        chi_file.open(filename);
    }

    data.plotSpectrum(config,"TTCV");

    std::vector<surfOut> grid;
    for(size_t i = 0; i < nbinsx; i++) {
        for(size_t j = 0; j < nbinsy; j++) {
            std::vector<int> grid_pts = {(int)i,(int)j};
            std::vector<float> physics_params = {(float)edges_y(j), (float)edges_x(i)};  //deltam^2, sin^22thetamumu
            surfOut pt; pt.grid_val = physics_params; pt.grid_index = grid_pts;
            grid.push_back(pt);
        }
    }

    int loopSize = grid.size();
    int chunkSize = loopSize / nThreads;

    //std::vector<surfOut> PointHelper(const PROconfig *config, const PROpeller *prop, const PROsyst *systs, const PROosc *osc, const PROspec *data, std::vector<surfOut>> multi_physics_params, PROchi::EvalStrategy strat, bool binned_weighting, int start, int end);

    std::vector<std::future<std::vector<surfOut>>> futures; 

    log<LOG_INFO>(L"%1% || Starting THREADS  : %2% , Loops %3%, Chunks %4%") % __func__ % nThreads % loopSize % chunkSize;

    for (int t = 0; t < nThreads; ++t) {
        int start = t * chunkSize;
        int end = (t == nThreads - 1) ? loopSize : start + chunkSize;
        futures.emplace_back(std::async(std::launch::async, [&, start, end]() {
                    return this->PointHelper(&config, &prop, &systs, &osc, &data, grid, strat, binned_weighting, start, end);
                    }));

    }

    std::vector<surfOut> combinedResults;
    for (auto& fut : futures) {
        std::vector<surfOut> result = fut.get();
        combinedResults.insert(combinedResults.end(), result.begin(), result.end());
    }

    for (const auto& item : combinedResults) {
        log<LOG_INFO>(L"%1% || Finished  : %2% %3% %4%") % __func__ % item.grid_val[1] % item.grid_val[0] %item.chi ;
        surface(item.grid_index[0], item.grid_index[1]) = item.chi;
        chi_file<<"\n"<<item.grid_val[1]<<" "<<item.grid_val[0]<<" "<<item.chi<<std::flush;
    }

}

std::vector<double> findMinAndBounds(TGraph* g, double val,double range) {
    int n = g->GetN();
    double minY = 1e9, minX = 0;
    for (int i = 0; i < n; ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        if (y < minY) {
            minY = y;
            minX = x;
        }
    }

    double leftX = minX, rightX = minX;
    
    // Search to the left of the minimum
    for (double x = minX; x >= -range; x -= 0.01) {
        double y = g->Eval(x);
        if (y >= val) {
            leftX = x;
            break;
        }
    }

    // Search to the right of the minimum
    for (double x = minX; x <= range; x += 0.01) {
        double y = g->Eval(x);
        if (y >= val) {
            rightX = x;
            break;
        }
    }
    
    return {minY,leftX,rightX};
}


int PROfit::PROfile(const PROconfig &config, const PROpeller &prop, const PROsyst &systs, const PROsc &osc, const PROspec &data, std::string filename) {


    LBFGSpp::LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    param.max_linesearch = 50;
    param.delta = 1e-6;

    LBFGSpp::LBFGSBSolver<double> solver(param);
    int nparams = systs.GetNSplines();
    std::vector<float> physics_params; 



    TCanvas *c =  new TCanvas(filename.c_str(), filename.c_str() , 400*4, 400*7);
    c->Divide(4,7);


    std::vector<std::unique_ptr<TGraph>> graphs; 

    //hack
    std::vector<double> priorX;
    std::vector<double> priorY;

    for(int i=0; i<=20;i++){
        double which_value = -2.0+0.2*i;
        priorX.push_back(which_value);
        priorY.push_back(which_value*which_value);

    }
    std::unique_ptr<TGraph> gprior = std::make_unique<TGraph>(priorX.size(), priorX.data(), priorY.data());



    for(int w=0; w<nparams;w++){
        int which_spline = w;


        std::vector<double> knob_vals;
        std::vector<double> knob_chis;

        for(int i=0; i<=20;i++){

            Eigen::VectorXd lb = Eigen::VectorXd::Constant(nparams, -3.0);
            Eigen::VectorXd ub = Eigen::VectorXd::Constant(nparams, 3.0);
            Eigen::VectorXd x = Eigen::VectorXd::Constant(nparams, 0.0);
            Eigen::VectorXd grad = Eigen::VectorXd::Constant(nparams, 0.0);
            Eigen::VectorXd bestx = Eigen::VectorXd::Constant(nparams, 0.0);


            double which_value = -2.0+0.2*i;
            double fx;
            knob_vals.push_back(which_value);

            lb[which_spline] = which_value;
            ub[which_spline] = which_value;
            x[which_spline] = which_value;


            PROchi chi("3plus1", &config, &prop, &systs, &osc, data, nparams, systs.GetNSplines(), PROchi::BinnedChi2, physics_params);
            chi.fixSpline(which_spline,which_value);

            log<LOG_INFO>(L"%1% || Starting Fixed fit ") % __func__  ;
            try {
                x = Eigen::VectorXd::Constant(nparams, 0.012);
                solver.minimize(chi, x, fx, lb, ub);
            } catch(std::runtime_error &except) {
                log<LOG_ERROR>(L"%1% || Fit failed, %2%") % __func__ % except.what();
            }

            std::string spec_string = "";
            for(auto &f : x) spec_string+=" "+std::to_string(f); 
            log<LOG_INFO>(L"%1% || Fixed value of %2% for spline %3% was post  : %4% ") % __func__ % which_spline % which_value % fx;
            log<LOG_INFO>(L"%1% || BF splines @ %2%") % __func__ %  spec_string.c_str();

            knob_chis.push_back(fx);
        }            

        log<LOG_INFO>(L"%1% || Knob Values: %2%") % __func__ %  knob_vals;
        log<LOG_INFO>(L"%1% || Knob Chis: %2%") % __func__ %  knob_chis;

        c->cd(w+1);
        std::unique_ptr<TGraph> g = std::make_unique<TGraph>(knob_vals.size(), knob_vals.data(), knob_chis.data());
        std::string tit = systs.spline_names[which_spline]+ ";#sigma Shift; #Chi^{2}";
        g->SetTitle(tit.c_str());
        graphs.push_back(std::move(g));
        graphs.back()->Draw("AL");
        graphs.back()->SetLineWidth(2);

        gprior->Draw("L same");
        gprior->SetLineStyle(2);
        gprior->SetLineWidth(1);

    }

    c->SaveAs((filename+".pdf").c_str(),"pdf");

    delete c;

    //Next version
    TCanvas *c2 =  new TCanvas((filename+"1sigma").c_str(), (filename+"1sigma").c_str() , 40*nparams, 400);
    c2->cd();
    
    //plot 2sigma also? default no, as its messier
    bool twosig = false;

    std::vector<double> bfvalues;
    std::vector<double> values1_up;
    std::vector<double> values1_down;

    std::vector<double> values2_up;
    std::vector<double> values2_down;

    for(auto &g:graphs){
        std::vector<double> tmp1 = findMinAndBounds(g,1.0,1.0);
        bfvalues.push_back(tmp[0]);
        values1_down.push_back(tmp[1]);
        values1_up.push_back(tmp[2]);

        if(twosig){
            std::vector<double> tmp2 = findMinAndBounds(g,4.0,2.0);
            values2_down.push_back(tmp2[1]);
            values2_up.push_back(tmp2[2]);
        }
    }

    
    int nBins = systs.spline_names.size();


    double range = twosig? 2 : 1;
    double wid = twosig? 0.4  : 0.8;
    double off1 = twosig?0.1   : 0.0;
    double off2 = twosig?0.5  : 0.0;

    TH1D *h1up = new TH1D("hup", "Bar Chart", nBins, 0, nBins);
    TH1D *h1down = new TH1D("hdown", "Bar Chart", nBins, 0, nBins);

    TH1D *h2up = new TH1D("h2up", "Bar Chart", nBins, 0, nBins);
    TH1D *h2down = new TH1D("h2down", "Bar Chart", nBins, 0, nBins);


    h1up->SetFillColor(kBlue-7);
    h1up->SetBarWidth(wid);
    h1up->SetBarOffset(off1);
    h1up->SetStats(0);
    h1up->SetMinimum(-range);
    h1up->SetMaximum(range);

    h1down->SetFillColor(kBlue-7);
    h1down->SetBarWidth(wid);
    h1down->SetBarOffset(off1);
    h1down->SetStats(0);

    h2up->SetFillColor(38);
    h2up->SetBarWidth(wid);
    h2up->SetBarOffset(off2);
    h2up->SetStats(0);

    h2down->SetFillColor(38);
    h2down->SetBarWidth(wid);
    h2down->SetBarOffset(off2);
    h2down->SetStats(0);


    // Fill the histogram with values from the vector
    for (int i = 0; i < nBins; ++i) {
        h1up->SetBinContent(i+1, values1_up[i]); 
        h1down->SetBinContent(i+1, values1_down[i]); 

        h2up->SetBinContent(i+1, values2_up[i]); 
        h2down->SetBinContent(i+1, values2_down[i]); 

        h1up->GetXaxis()->SetBinLabel(i+1,systs.spline_names[i].c_str());

    }
    h1up->SetTitle("");
    h1up->Draw("b");
    h1up->GetYaxis()->SetTitle("#sigma Shift");

    h1down->Draw("b same");
    if(twosig){
        h2up->Draw("b same");
        h2down->Draw("b same");
    }

    TLine l(0,0,nBins,0);
    l.SetLineStyle(2);
    l.SetLineColor(kBlack);
    l.SetLineWidth(1);
    l.Draw();

    for (int i = 0; i < nBins; ++i) {
        std::unique_ptr<TMarker> star = std::make_unique<TMarker>(i + wid/2.0, bfvalues[i], 29);
        star->SetMarkerSize(2); 
        star->SetMarkerColor(kBlack); 
        star->Draw();
    }


    c2->SaveAs((filename+"_1sigma.pdf").c_str(),"pdf");
    delete c2;

    return 0;
}

