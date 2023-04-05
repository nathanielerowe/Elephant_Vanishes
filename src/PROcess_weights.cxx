#include "PROcess_weights.h"
#include "TTree.h"
#include "TFile.h"

using namespace PROfit;

void SystStruct::CleanSpecs(){
    multi_vecspec.clear();
    return;
}

void void SetSpecDimension(int row, int col){
}
void SystStruct::SanityCheck() const{
    if(mode == "minmax" && n_univ != 2){
        log<LOG_ERROR>(L"%1% || Systematic variation %2% is tagged as minmax mode, but has %3% universes (can only be 2)") % __func__ % systname.c_str() % n_univ;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }

    log<LOG_INFO>(L"%1% || Systematic variation %2% passed sanity check!!") % __func__ % systname.c_str();
    log<LOG_INFO>(L"%1% || Systematic variation %2% has %3% universes, and is in %4% mode with weight formula: %5%") % __func__ % systname.c_str() % n_univ % mode.c_str() % weight_formula.c_str();
    return;
}



int PROcess(const PROconfig &inconfig){
 
    log<LOG_DEBUG>(L"%1% || Starting to construct CovarianceMatrixGeneration in EventWeight Mode  ") % __func__ ;

    int universes_used = 0;

    std::vector<std::string> variations;
    std::vector<std::string> variations_tmp;

    int num_files = inconfig.m_num_mcgen_files;

    log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;

    std::vector<long int> nentries(num_files,0);
    std::vector<int> used_montecarlos(num_files,0);

    std::vector<std::unique_ptr<TFile>> files(num_files);
    std::vector<TTree*> trees(num_files,nullptr);//keep as bare pointers because of ROOT :(
    std::vector<std::vector<std::map<std::string, std::vector<eweight_type>>* >> f_weights(num_files);
    std::map<std::string, int> map_systematic_num_universe;

    //inconfig.m_mcgen_additional_weight.resize(num_files,1.0); its const, not allowed

    int good_event = 0;

    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = inconfig.m_mcgen_file_name.at(fid);

        files[fid] = std::make_unique<TFile>(fn.c_str(),"read");
        trees[fid] = (TTree*)(files[fid]->Get(inconfig.m_mcgen_tree_name.at(fid).c_str()));
        nentries[fid]= (long int)trees.at(fid)->GetEntries();

	if(files[fid]->IsOpen()){
    	    log<LOG_INFO>(L"%1% || Root file succesfully opened: %2%") % __func__  % fn.c_str();
	}else{
    	    log<LOG_ERROR>(L"%1% || Fail to open root file: %2%") % __func__  % fn.c_str();
	    exit(EXIT_FAILURE);
	}
    	log<LOG_INFO>(L"%1% || Total Entries: %2%") % __func__ %  nentries[fid];

        //Some POT counting (FIX)
        //Guanqun: POT counting not needed for covariance matrix generation 
        //double pot_scale = 1.0;
        //if(inconfig.m_mcgen_pot[fid]!=-1){
        //    pot_scale = FIX_plot_pot/inconfig.m_mcgen_pot[fid];
        //}
        //mcgen_scale[fid] = inconfig.m_mcgen_scale[fid]*pot_scale;

    	//first, grab friend trees
        auto mcgen_file_friend_treename_iter = inconfig.m_mcgen_file_friend_treename_map.find(fn);
        if (mcgen_file_friend_treename_iter != inconfig.m_mcgen_file_friend_treename_map.end()) {

            auto mcgen_file_friend_iter = inconfig.m_mcgen_file_friend_map.find(fn);
            if (mcgen_file_friend_iter == inconfig.m_mcgen_file_friend_map.end()) {
    	    	log<LOG_ERROR>(L"%1% || Friend TTree provided but no friend file??") % __func__;
		log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }

            for(int k=0; k < mcgen_file_friend_treename_iter->second.size(); k++){

                std::string treefriendname = mcgen_file_friend_treename_iter->second.at(k);
                std::string treefriendfile = mcgen_file_friend_iter->second.at(k);
                trees[fid]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());
            }
        }

        // grab branches 
        int num_branch = inconfig.m_branch_variables[fid].size();
	f_weights[fid].resize(num_branch);
        for(int ib = 0; ib != num_branch; ++ib) {

            const auto& branch_variable = inconfig.m_branch_variables[fid][ib];

            //quick check that this branch associated subchannel is in the known chanels;
            int is_valid_subchannel = 0;
            for(const auto &name: inconfig.m_fullnames){
                if(branch_variable->associated_hist==name){
                    log<LOG_DEBUG>(L"%1% || Found a valid subchannel for this branch %2%") % __func__  % name.c_str();
                    ++is_valid_subchannel;
                }
            }
            if(is_valid_subchannel==0){
    	    	log<LOG_ERROR>(L"%1% || This branch did not match one defined in the .xml : %2%") % __func__ % inconfig.m_xmlname.c_str();
    	    	log<LOG_ERROR>(L"%1% || There is probably a typo somehwhere in xml!") % __func__;
		log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);

            }else if(is_valid_subchannel>1){
    	    	log<LOG_ERROR>(L"%1% || This branch matched more than 1 subchannel!: %2%") % __func__ %  branch_variable->associated_hist.c_str();
		log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }

            branch_variable->branch_formula = std::make_shared<TTreeFormula>(("branch_form_"+std::to_string(fid) +"_" + std::to_string(ib)).c_str(), branch_variable->name.c_str(), trees[fid]);
    	    log<LOG_INFO>(L"%1% || Setting up reco variable for this branch: %2%") % __func__ %  branch_variable->name.c_str();


	    //grab monte carlo weight
	    if(m_mcgen_additional_weight_bool[fid][ib]){
                branch_variable->branch_monte_carlo_weight_formula  =  std::make_shared<TTreeFormula>(("branch_add_weight_"+std::to_string(fid)+"_" + std::to_string(ib)).c_str(),inconfig.m_mcgen_additional_weight_name[fid][id].c_str(),trees[fid]);
    	    	log<LOG_INFO>(L"%1% || Setting up additional monte carlo weight for this branch: %2%") % __func__ %  inconfig.m_mcgen_additional_weight_name[fid][id].c_str();
            }


	    //grab eventweight branch
    	    log<LOG_INFO>(L"%1% || Setting up eventweight map for this branch: %2%") % __func__ %  inconfig.m_mcgen_eventweight_branch_names[fid][ib].c_str();
	    trees[fid]->SetBranchAddress(inconfig.m_mcgen_eventweight_branch_names[fid][ib].c_str(), &(f_weights[fid][ib]));

	    if(!f_weights[fid][ib]){
            	log<LOG_ERROR>(L"%1% || Could not read eventweight branch for file=%2%") % __func__ % fid ;
		log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
	    }
        } //end of branch loop


	//calculate how many "universes" each systematoc has.
    	log<LOG_INFO>(L"%1% || Start calculating number of universes for systematics") % __func__;
        trees.at(fid)->GetEntry(good_event);
	for(int ib = 0; ib != num_branch; ++ib) {
            const auto& branch_variable = inconfig.m_branch_variables[fid][ib];
	    auto& f_weight = f_weights[fid][ib];

            for(const auto& it : *f_weight){
    	    	log<LOG_INFO>(L"%1% || On systematic: %2%") % __func__ % it.first.c_str();

            	if(inconfig.m_mcgen_variation_allowlist.count(it.first)==0){
    	    	    log<LOG_INFO>(L"%1% || Skip systematic: %2% as its not in the AllowList!!") % __func__ % it.first.c_str();
                    continue;
                }

            	if(inconfig.m_mcgen_variation_denylist.count(it.first)>0){
    	    	    log<LOG_INFO>(L"%1% || Skip systematic: %2% as it is in the DenyList!!") % __func__ % it.first.c_str();
                    continue;
            	}

            	log<LOG_INFO>(L"%1% || %2% has %3% montecarlo variations in branch %4%") % __func__ % it.first.c_str() % it.second.size() % branch_variable->associated_hist.c_str();

		map_systematic_num_universe[it.first] = std::max(map_systematic_num_universe[it.first], it.second.size());
	    }
        }
    } // end fid

    log<LOG_INFO>(L"%1% || Found %2% unique variations") % __func__ % map_systematic_num_universe.size();
    log<LOG_INFO>(L"%1% || Now start to grab related weightmaps") % __func__;

    int total_num_systematics = map_systematic_num_universe.size();
    std::vector<SystStruct> syst_vector;
    for(auto& sys_pair : map_systematic_num_universe){

	const std::string& sys_name = sys_pair.first;
        syst_vector.emplace_back(sys_name, sys_pair.second);

            
	// Check to see if pattern is in this variation
	std::string sys_weight_formula = "1", sys_mode ="";

	for(int i = 0 ; i != inconfig.m_mcgen_weightmaps_patterns.size(); ++i){
            if (inconfig.m_mcgen_weightmaps_uses[i] && sys_name.find(inconfig.m_mcgen_weightmaps_patterns[i]) != std::string::npos) {
                sys_weight_formula = sys_weight_formula + "*(" + inconfig.m_mcgen_weightmaps_formulas[i]+")";
                sys_mode=inconfig.m_mcgen_weightmaps_mode[i];

		log<LOG_INFO>(L"%1% || Systematic variation %2% is a match for patten %3%") % __func__ % sys_name.c_str() % inconfig.m_mcgen_weightmaps_patterns[i].c_str();
		log<LOG_INFO>(L"%1% || Corresponding weight is : %2%") % __func__ % inconfig.m_mcgen_weightmaps_formulas[i].c_str();
		log<LOG_INFO>(L"%1% || Corresponding mode is : %2%") % __func__ % inconfig.m_mcgen_weightmaps_mode[i].c_str();
            }
	}

	if(sys_weight_formula != "1" || sys_mode !=""){
	    syst_vector.back().SetWeightFormula(sys_weight_formula);
	    syst_vector.back().SetMode(sys_mode);
	}
    }


    //sanity check 
    for(const auto& s : syst_vector)
	s.SanityCheck();

    std::vector<std::vector<std::unique_ptr<TTreeFormula>>> additional_weight_formulas(num_files);//(num_files, std::vector<std::unique_ptr<TTreeFormula>>(variations.size()));FIX why does this fail?
    for (auto& innerVector : additional_weight_formulas) {
        innerVector.resize(variations.size());
    }

    for(int fid=0; fid < num_files; ++fid) {
        files[fid]->cd();
        for(int vid = 0; vid < variations.size(); vid++){ 
            additional_weight_formulas[fid][vid] =  std::make_unique<TTreeFormula>(("weightMapsFormulas_"+std::to_string(fid)+"_"+std::to_string(vid)).c_str(), s_formulas[vid].c_str(),trees[fid]);
        }
    }


    //CHeck This FIX
    std::vector<int> num_universes_per_variation;
    std::map<int, std::string> map_universe_to_var;
    std::vector<int> vec_universe_to_var;
    std::map<std::string, int> map_var_to_num_universe;

    map_universe_to_var.clear();
    vec_universe_to_var.clear();
    num_universes_per_variation.clear();


    //But in reality we want the max universes to be the sum of all max variaitons across all files, NOT the sum over all files max variations.
    universes_used = num_universes_per_variation.size();


    std::cout << " -------------------------------------------------------------" << std::endl;
    std::cout << " Initilizing " << universes_used << " universes." << std::endl;
    std::cout << " -------------------------------------------------------------" << std::endl;

    std::vector<double> base_vec (inconfig.m_num_bins_total,0.0);

    std::cout << " Full concatanated vector has : " << inconfig.m_num_bins_total << std::endl;

    std::vector<std::vector<double>> multi_vecspec; ///FIX REplace with Eigen
    multi_vecspec.clear();
    multi_vecspec.resize(universes_used,base_vec);

    std::cout << " multi_vecspec now initilized of size :" << multi_vecspec.size() << std::endl;
    std::cout << " Reading the data files" << std::endl;
    //watch.Reset();
    //watch.Start();

    for(int j=0; j < num_files; j++){
        int nevents = std::min(inconfig.m_mcgen_maxevents[j], nentries[j]);
        std::cout << " Starting @ data file=" << files[j]->GetName() <<" which has "<<nevents<<" Events. "<<std::endl;
        size_t nbytes = 0;
        for(int i=0; i < nevents; i++) {
            if(i%100==0)std::cout<<" -- uni :"<<i<<" / "<<nevents<<std::endl;
            nbytes+= trees[j]->GetEntry(i);
            ProcessEvent(inconfig, *(f_weights[j]),j,i);
            //INPUT PROCESS FIX

        } //end of entry loop
        std::cout << " nbytes read=" << nbytes << std::endl;

    } //end of file loop



    //watch.Stop();
    //std::cout << " done CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

    /***************************************************************
     *		Now some clean-up and Writing
     * ************************************************************/

    //FIX check?
    //NOPE In fact, calling Close() on a TFile object managed by a std::unique_ptr can lead to undefined behavior because the std::unique_ptr may delete the object before the Close() method finishes executing. So, it is best to rely on the std::unique_ptr to manage the lifetime of the TFile object and not call Close() explicitly.
    //for(auto f: files){
    //    std::cout << " TFile::Close() file=" << f->GetName() <<  std::endl;
    //    f->Close();
    //}
    std::cout << " End" << std::endl;


    return 0;
}


/*
void ProcessEvent(const PROconfig &inconfig,
        const std::map<std::string, 
        std::vector<eweight_type> >& thisfWeight,
        size_t fileid,
        int entryid) {

    double abnormally_large_weight = 1e3;//1e20;//20.0;
    double global_weight = 1.0;

    if( montecarlo_additional_weight_bool[fileid]){
        montecarlo_additional_weight_formulas[fileid]->GetNdata();
        global_weight = montecarlo_additional_weight_formulas[fileid]->EvalInstance();
    };//this will be 1.0 unless specifi
    global_weight *= montecarlo_scale[fileid];

    double additional_CV_weight = 1.0;

    const auto bnbcorr_iter = thisfWeight.find(bnbcorrection_str);
    if (bnbcorr_iter != thisfWeight.end())    additional_CV_weight *= (*bnbcorr_iter).second.front();

    if(std::isinf(global_weight) or (global_weight != global_weight)){
        std::stringstream ss;
        ss << "ProcessEvent\t||\tERROR  error @ " << entryid
            << " in File " << montecarlo_file.at(fileid) 
            << " as its either inf/nan: " << global_weight << std::endl;
        throw std::runtime_error(ss.str());
    }

    if(!EventSelection(fileid)) return;

    // precompute the weight size
    std::vector<double> weights(universes_used,global_weight);

    //Loop over all variations
    std::map<std::string, std::vector<eweight_type> >::const_iterator var_iter;
    int woffset = 0;
    int vid = 0;
    for(const auto& var : variations){


        //check if variation is in this file, if it isn't: then just push back 1's of appropiate number to keep universes consistent
        //this is of length of whatever the maximum length that was found in ANY file
        auto expected_num_universe_sz = map_var_to_num_universe.at(var); 

        //is  
        var_iter = thisfWeight.find(var);
        int quick_fix = 0;
	double indiv_variation_weight = 1.0;

        if (var_iter == thisfWeight.end()) {
            //This we need to drop this for new version (where we add 1's)
            //std::cout<<var<<" is not in this universe, adding "<<expected_num_universe_sz<<" to woffset "<<woffset<<std::endl;
            //woffset += expected_num_universe_sz;
            //continue;
        }else {
            //first one is what we expect, second is whats directly inside the map.

            if (expected_num_universe_sz != (*var_iter).second.size()) {
                std::stringstream ss;
                //std::cout<< "Number of universes is not the max in this file" << std::endl;
                //throw std::runtime_error(ss.str());

                if( (*var_iter).second.size() > expected_num_universe_sz && var_iter != thisfWeight.end()){
                    ss << ". REALLY BAD!!  iter.size() " <<(*var_iter).second.size()<<" expected "<<expected_num_universe_sz<<" on var "<<var<<std::endl;
                    throw std::runtime_error(ss.str());
                }
            }
            //so if this file contains smaller number of variations
            quick_fix = (*var_iter).second.size();
            //std::cout<< "So setting quick fix to the true number of universes in this varaiion : "<<quick_fix<< std::endl;

	    if(!montecarlo_fake[fileid]){
            	//Grab newwer variation specfic weights;
             	m_variation_weight_formulas[fileid][vid]->GetNdata();
            	indiv_variation_weight = m_variation_weight_formulas[fileid][vid]->EvalInstance();
            	if((indiv_variation_weight!= indiv_variation_weight || indiv_variation_weight <0) && !montecarlo_fake[fileid]){
		    	std::cout << "varaition: " << var << std::endl;
                    	std::cout<<"ERROR! the additional wight is nan or negative "<<indiv_variation_weight<<" Breakign!"<<std::endl;
                    	exit(EXIT_FAILURE);
            	}
  	    }
            //std::cout<<var<<" "<<indiv_variation_weight<<" "<<fileid<<" "<<vid<<std::endl;

        }

        if (woffset >= weights.size()) {
            std::stringstream ss;
            ss << "woffset=" << woffset << " weights sz=" << weights.size() << " !" << std::endl;
            throw std::runtime_error(ss.str());
        }

        for(size_t wid0 = 0, wid1 = woffset; wid1 < (woffset + expected_num_universe_sz); ++wid0, ++wid1) {
            double wei = 1.0;
            //            std::cout<<"wid0 "<<wid0<<"/ "<<expected_num_universe_sz<<"  wid1 "<<wid1<<" / "<<weights.size()<<" woffset "<<woffset<<" quick_fix "<<quick_fix<<std::endl;

            if(wid0<quick_fix){
                wei = (*var_iter).second[wid0];
            }
       
            bool is_inf = std::isinf(wei);
            bool is_nan = (wei != wei);

            if(is_inf or is_nan){
                std::stringstream ss;
                ss << "ProcessEvent\t||\t ERROR! Killing :: event # " << entryid
                    << " in File " << montecarlo_file.at(fileid) << " weight: " << wei << " global bnb: " << global_weight << " in " << var << std::endl;
                //                throw std::runtime_error(ss.str());
                wei = 1.0;                
            }

            if(wei > abnormally_large_weight){
                std::cout<<"ProcessEvent\t||\tATTENTION!! HUGE weight: "<<wei<<" at "<<var<<" event "<<entryid<<" file "<<fileid<<std::endl;
                wei=1.0;
            }
		

	    weights[wid1] *= wei*indiv_variation_weight;
        }

        woffset += expected_num_universe_sz;
        vid++;
    }//end of all variations

    if (woffset != weights.size()) {
        std::stringstream ss;
        ss << "woffset=" << woffset << " weights sz=" << weights.size() << " !" << std::endl;
        throw std::runtime_error(ss.str());
    }

    //So the size of weights must equal global universes ya?
    if(universes_used != num_universes_per_variation.size()){
        std::stringstream ss;
        ss <<otag<<" ERROR "<<std::endl;
        ss <<"weights.size() "<<weights.size()<<std::endl;
        ss <<"universes_used "<<universes_used<<std::endl;
        ss <<"multi_vecspec.size() "<<multi_vecspec.size()<<std::endl;
        ss <<"num_universes_per_variation.size() "<<num_universes_per_variation.size()<<std::endl;
        throw std::runtime_error(ss.str());
    }
    for(int t=0; t < branch_variables[fileid].size(); t++) {


        const auto branch_var_jt = branch_variables[fileid][t];
        int ih = spec_central_value.map_hist.at(branch_var_jt->associated_hist);

        branch_var_jt->GetFormula()->GetNdata();
        double reco_var = branch_var_jt->GetFormula()->EvalInstance();
        // double reco_varold = *(static_cast<double*>(branch_var_jt->GetValue()));
        int reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);
        spec_central_value.hist[ih].Fill(reco_var, global_weight*additional_CV_weight);
        //std::cout<<reco_var<<" "<<reco_bin<<" "<<ih<<std::endl;

        for(int m=0; m<weights.size(); m++){
            if(reco_bin<0) continue;
            multi_vecspec[m][reco_bin] += weights[m];
        }
    }

    return;*/
//}

