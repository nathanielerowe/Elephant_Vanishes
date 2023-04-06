#include "PROcess_weights.h"
#include "TTree.h"
#include "TFile.h"

namespace PROfit {

void SystStruct::CleanSpecs(){
    if(p_multi_spec)
        p_multi_spec.reset(nullptr);
    return;
}

void SystStruct::CreateSpecs(int row, int col){
    this->CleanSpecs();
    log<LOG_INFO>(L"%1% || Creating multi-universe spectrum with dimension (%2% x %3%)") % __func__ % row % col;
    p_multi_spec = std::make_unique<multi_spec>(row, col);
    p_multi_spec->setZero(row, col);
    return;
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



int PROcess_SBNfit(const PROconfig &inconfig){
 
    log<LOG_DEBUG>(L"%1% || Starting to construct CovarianceMatrixGeneration in EventWeight Mode  ") % __func__ ;

    int universes_used = 0;

    std::vector<std::string> variations;
    std::vector<std::string> variations_tmp;

    int num_files = inconfig.m_num_mcgen_files;

    log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;

    std::vector<long int> nentries(num_files,0);
    std::vector<std::unique_ptr<TFile>> files(num_files);
    std::vector<TTree*> trees(num_files,nullptr);//keep as bare pointers because of ROOT :(
    std::vector<std::vector<std::map<std::string, std::vector<eweight_type>>* >> f_event_weights(num_files);
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
	f_event_weights[fid].resize(num_branch);
        for(int ib = 0; ib != num_branch; ++ib) {

            std::shared_ptr<BranchVariable> branch_variable = inconfig.m_branch_variables[fid][ib];

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
	    if(inconfig.m_mcgen_additional_weight_bool[fid][ib]){
                branch_variable->branch_monte_carlo_weight_formula  =  std::make_shared<TTreeFormula>(("branch_add_weight_"+std::to_string(fid)+"_" + std::to_string(ib)).c_str(),inconfig.m_mcgen_additional_weight_name[fid][ib].c_str(),trees[fid]);
    	    	log<LOG_INFO>(L"%1% || Setting up additional monte carlo weight for this branch: %2%") % __func__ %  inconfig.m_mcgen_additional_weight_name[fid][ib].c_str();
            }


	    //grab eventweight branch
    	    log<LOG_INFO>(L"%1% || Setting up eventweight map for this branch: %2%") % __func__ %  inconfig.m_mcgen_eventweight_branch_names[fid][ib].c_str();
	    trees[fid]->SetBranchAddress(inconfig.m_mcgen_eventweight_branch_names[fid][ib].c_str(), &(f_event_weights[fid][ib]));

	    if(!f_event_weights[fid][ib]){
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
	    auto& f_weight = f_event_weights[fid][ib];

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

		map_systematic_num_universe[it.first] = std::max((int)map_systematic_num_universe[it.first], (int)it.second.size());
	    }
        }
    } // end fid

    log<LOG_INFO>(L"%1% || Found %2% unique variations") % __func__ % map_systematic_num_universe.size();
    for(auto& sys_pair : map_systematic_num_universe){
    	log<LOG_DEBUG>(L"%1% || Variation: %2% --> %3% universes") % __func__ % sys_pair.first.c_str() % sys_pair.second;
    }

    //constuct object for each systematic variation, and grab weight maps
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


    //create 2D multi-universe spec.
    for(auto& s : syst_vector){
	int nrow = s.GetNUniverse(), ncol = inconfig.m_num_bins_total;
 	s.CreateSpecs(nrow, ncol);	
    }


    time_t start_time = time(nullptr);
    log<LOG_INFO>(L"%1% || Start reading the files..") % __func__;
    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = inconfig.m_mcgen_file_name.at(fid);
        long int nevents = std::min(inconfig.m_mcgen_maxevents[fid], nentries[fid]);
	log<LOG_DEBUG>(L"%1% || Start @files: %2% which has %3% events") % __func__ % fn.c_str() % nevents;


    	std::vector<std::unique_ptr<TTreeFormula>> sys_weight_formula;
	for(const auto& s : syst_vector){
	    sys_weight_formula.push_back(std::make_unique<TTreeFormula>(("weightMapsFormulas_"+std::to_string(fid)+"_"+ s.GetSysName()).c_str(), s.GetWeightFormula().c_str(),trees[fid]));
  	}
	log<LOG_DEBUG>(L"%1% || Finished setting up systematic weight formula") % __func__;


        for(long int i=0; i < nevents; ++i) {

	    trees[fid]->GetEntry(i);
            if(i%100==0)
		log<LOG_DEBUG>(L"%1% || -- uni : %2% / %3%") % __func__ % i % nevents;
            
            ProcessEvent(inconfig, fid, f_event_weights[fid], syst_vector);
            //INPUT PROCESS FIX

        } //end of entry loop

    } //end of file loop

    time_t time_took = time(nullptr) - start_time;
    log<LOG_INFO>(L"%1% || Finish reading files, it took %2% seconds..") % __func__ % time_took;


    //FIX check?
    //NOPE In fact, calling Close() on a TFile object managed by a std::unique_ptr can lead to undefined behavior because the std::unique_ptr may delete the object before the Close() method finishes executing. So, it is best to rely on the std::unique_ptr to manage the lifetime of the TFile object and not call Close() explicitly.
    //for(auto f: files){
    //  log<LOG_DEBUG>(L"%1% || Closing file: %2%") %__func__ % f->GetName();
    //    f->Close();
    //}
    log<LOG_INFO>(L"%1% || DONE") %__func__ ;

    return 0;
}

int PROcess_CAFana(const PROconfig &inconfig){

    log<LOG_DEBUG>(L"%1% || Starting to construct CovarianceMatrixGeneration in EventWeight Mode  ") % __func__ ;

    int universes_used = 0;
    int num_files = inconfig.m_num_mcgen_files;
    std::vector<SystStruct> syst_vector;
    
    log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;

    std::vector<int> nentries(num_files,0);
    std::vector<int> used_montecarlos(num_files,0);


    std::vector<std::unique_ptr<TFile>> files(num_files);
    std::vector<TTree*> trees(num_files,nullptr);//keep as bare pointers because of ROOT :(

    //std::vector<std::map<std::string, std::vector<eweight_type> >* > f_weights(num_files,nullptr);
    std::vector<int> pset_indices_tmp;
    std::vector<int> pset_indices;
    std::map<std::string, int> map_systematic_num_universe;
    std::vector<std::vector<float>> knobvals;
    std::vector<std::string> cafana_pset_names;

    int good_event = 0;


    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = inconfig.m_mcgen_file_name.at(fid);

    
        log<LOG_DEBUG>(L"%1% || On file %2% - %3% Getting SRglobal ") % __func__ % fid % inconfig.m_mcgen_file_name.at(fid).c_str()  ;
        caf::SRGlobal* global = NULL;
        files[fid] = std::make_unique<TFile>(fn.c_str(),"read");
        trees[fid] = (TTree*)(files[fid]->Get(inconfig.m_mcgen_tree_name.at(fid).c_str()));
        TTree * globalTree = (TTree*)(files[fid]->Get("globalTree"));
        nentries[fid]= (int)trees.at(fid)->GetEntries();

        globalTree->SetBranchAddress("global", &global);
        log<LOG_DEBUG>(L"%1% || On file %2% - %3% Getting first entry from global tree ") % __func__ % fid % fn.c_str()  ;
        globalTree->GetEntry(0);

        log<LOG_DEBUG>(L"%1% || On file %2% - %3% Grabbing Weights ") % __func__ % fid % fn.c_str()  ;
        for(unsigned int i = 0; i < global->wgts.size(); ++i) {
            const caf::SRWeightPSet& pset = global->wgts[i];
            pset_indices.push_back(i);
            cafana_pset_names.push_back(pset.name);
            map_systematic_num_universe[pset.name] = std::max(map_systematic_num_universe[pset.name], pset.nuniv);
            knobvals.push_back(pset.map.at(0).vals);
        }

        auto mcgen_file_friend_treename_iter = inconfig.m_mcgen_file_friend_treename_map.find(fn);
        if (mcgen_file_friend_treename_iter != inconfig.m_mcgen_file_friend_treename_map.end()) {

            auto mcgen_file_friend_iter = inconfig.m_mcgen_file_friend_map.find(fn);
            if (mcgen_file_friend_iter == inconfig.m_mcgen_file_friend_map.end()) {
                std::stringstream ss;
                throw std::runtime_error(ss.str());
            }

            for(int k=0; k < (*mcgen_file_friend_iter).second.size(); k++){

                std::string treefriendname = (*mcgen_file_friend_treename_iter).second.at(k);
                std::string treefriendfile = (*mcgen_file_friend_iter).second.at(k);
                trees[fid]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());

                log<LOG_DEBUG>(L"%1% || On file %2% - %3% Adding FriendTrees : %3% ") % __func__ % fid % fn.c_str()  % treefriendname.c_str();
            }
        }

        //if(inconfig.m_mcgen_additional_weight_bool[fid]){
        //we have an additional weight we want to apply at run time, otherwise its just set at 1.
        //    std::cout<<"Setting Additional weight of : "<< inconfig.m_mcgen_additional_weight_name[fid].c_str()<<std::endl; 
        //FIX FIX
        //additional_weight_formulas[fid] =  std::make_shared<TTreeFormula>(("a_w"+std::to_string(fid)).c_str(),inconfig.m_mcgen_additional_weight_name[fid].c_str(),trees[fid]);
        //}

        std::cout<<"Total Entries: "<<trees.at(fid)->GetEntries()<<" good event "<<good_event<<std::endl;
        trees.at(fid)->GetEntry(good_event);

        //This bit will calculate how many "universes" the file has. if ALL default is the inputted xml value



        log<LOG_DEBUG>(L"%1% || On file %2% - %3% Starting on cafana pset loop to build SysVec ") % __func__ % fid % fn.c_str()  ;
        for(int v =0; v< cafana_pset_names.size(); v++){
            std::string varname = cafana_pset_names[v];

            log<LOG_DEBUG>(L"%1% || starting %2% ") % __func__ % varname.c_str()  ;

            if(inconfig.m_mcgen_variation_allowlist.size()> 0 ){
                if(inconfig.m_mcgen_variation_allowlist.count(varname)==0){
                    log<LOG_INFO>(L"%1% || Skipping %2% as its not in allowlist ") % __func__ % varname.c_str();
                    continue;
                }
            }

            if(inconfig.m_mcgen_variation_denylist.size()> 0 ){
                if(inconfig.m_mcgen_variation_denylist.count(varname)>0){
                    log<LOG_INFO>(L"%1% || Skipping %2% as it in denylist ") % __func__ % varname.c_str();
                    continue;
                }
            }

            //Variation Weight Maps Area
            int n_wei = inconfig.m_mcgen_weightmaps_patterns.size();
            std::string variation_mode = "multisim"; 
            std::string s_formula = "1";

            for(int i=0; i< n_wei; i++){
                // Check to see if pattern is in this variation
                if (varname.find(inconfig.m_mcgen_weightmaps_patterns[i]) != std::string::npos) {
                    log<LOG_DEBUG>(L"%1% || Variation %2% is a match for pattern %3%") % __func__ % varname.c_str() % inconfig.m_mcgen_weightmaps_patterns[i].c_str() ;

                    s_formula = s_formula + "*(" + inconfig.m_mcgen_weightmaps_formulas[i]+")";
                    log<LOG_DEBUG>(L"%1% || Weight is thus %2% and mode is %3% ") % __func__ % s_formula.c_str() % inconfig.m_mcgen_weightmaps_mode[i].c_str();
                }
            }


            log<LOG_DEBUG>(L"%1% || %2% has %3% montecarlos in fie %4% ") % __func__ % varname.c_str() % map_systematic_num_universe[varname] % fid  ;
            used_montecarlos[fid] += map_systematic_num_universe[varname];

            //Some code to check if this varname is already in the syst_vector. If it is, check if things are the same, otherwise PANIC!
            log<LOG_DEBUG>(L"%1% || emplace syst_vector") % __func__  ;
            log<LOG_DEBUG>(L"%1% || varname: %2% , map: %3% , sform: %5% , varmode %4%, knobvals[v]size %6%, psetindex %7%  ") % __func__  % varname.c_str() % map_systematic_num_universe[varname] % variation_mode.c_str() % s_formula.c_str() % knobvals[v].size() % pset_indices[v] ;
             
            syst_vector.emplace_back(varname, map_systematic_num_universe[varname], variation_mode, s_formula, knobvals[v], pset_indices[v]);
        }
    } // end fid
    
    std::vector<std::vector<std::unique_ptr<TTreeFormula>>> additional_weight_formulas(num_files);
    for (auto& innerVector : additional_weight_formulas) {
        innerVector.resize(syst_vector.size());
    }

    for(int fid=0; fid < num_files; ++fid) {
        files[fid]->cd();
        for(int vid = 0; vid < syst_vector.size(); vid++){ 
            additional_weight_formulas[fid][vid] =  std::make_unique<TTreeFormula>(("weightMapsFormulas_"+std::to_string(fid)+"_"+std::to_string(vid)).c_str(), syst_vector[vid].formula.c_str(),trees[fid]);
        }
    }

    //But in reality we want the max universes to be the sum of all max variaitons across all files, NOT the sum over all files max variations.
    //universes_used = num_universes_per_variation.size();
    //FIX Replace With Eigen
    

    for(int j=0; j < num_files; j++){
        int nevents = nentries[j];//std::min(inconfig.m_mcgen_maxevents[j], nentries[j]);
        std::cout << " Starting @ data file=" << files[j]->GetName() <<" which has "<<nevents<<" Events. "<<std::endl;
        for(int i=0; i < nevents; i++) {
            if(i%100==0)std::cout<<" -- uni :"<<i<<" / "<<nevents<<std::endl;
            trees[j]->GetEntry(i);
            //ProcessEvent(inconfig, *(f_weights[j]),j,i);

        } //end of entry loop
    } //end of file loop
    std::cout << " End" << std::endl;

    return 0;
}




}
