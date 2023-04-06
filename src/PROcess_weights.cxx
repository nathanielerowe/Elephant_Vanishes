#include "PROcess_weights.h"
#include "TTree.h"
#include "TFile.h"

namespace PROfit {

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

};

