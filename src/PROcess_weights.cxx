#include "PROcess_weights.h"
#include "TTree.h"
#include "TFile.h"

using namespace PROfit;

int PROcess_CAFana(const PROconfig &inconfig){

    log<LOG_DEBUG>(L"%1% || Starting to construct CovarianceMatrixGeneration in EventWeight Mode  ") % __func__ ;

    int universes_used = 0;
    //Is there a Global weight to be applied to ALL weights, CV and otherwise, inside the eventweight class? 
    std::string bnbcorrection_str = "NAN";//this has mainly beeen superseeded by XML creation

    std::vector<std::string> variations;
    std::vector<std::string> variations_tmp;

    int num_files = inconfig.m_num_mcgen_files;

    log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;

    std::vector<int> nentries(num_files,0);
    std::vector<int> used_montecarlos(num_files,0);

    std::vector<std::unique_ptr<TFile>> files(num_files);
    std::vector<TTree*> trees(num_files,nullptr);//keep as bare pointers because of ROOT :(

    //std::vector<std::map<std::string, std::vector<eweight_type> >* > f_weights(num_files,nullptr);
    std::vector<int> pset_indices_tmp;
    std::vector<int> pset_indices;
    std::map<std::string, int> map_systematic_num_universe;
    std::vector<std::vector<float>> knobvals_tmp;
    std::vector<std::vector<float>> knobvals;
    std::vector<std::string> cafana_pset_names;

    int good_event = 0;

    std::vector<SystStruct> syst_vector;
    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = inconfig.m_mcgen_file_name.at(fid);

        caf::SRGlobal* global = NULL;
        files[fid] = std::make_unique<TFile>(fn.c_str(),"read");
        trees[fid] = (TTree*)(files[fid]->Get(inconfig.m_mcgen_tree_name.at(fid).c_str()));
        TTree * globalTree = (TTree*)(files[fid]->Get("globalTree"));
        nentries[fid]= (int)trees.at(fid)->GetEntries();

        globalTree->SetBranchAddress("global", &global);
        globalTree->GetEntry(0);

        for(unsigned int i = 0; i < global->wgts.size(); ++i) {
            const caf::SRWeightPSet& pset = global->wgts[i];
            pset_indices.push_back(i);
            cafana_pset_names.push_back(pset.name);
            map_systematic_num_universe[pset.name] = std::max(map_systematic_num_universe[pset.name], pset.nuniv);
            knobvals_tmp[i] = pset.map.at(0).vals;
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



        for(int v =0; v< cafana_pset_names.size(); v++){
            std::string varname = cafana_pset_names[v];

            std::cout<<"On : "<<varname<<std::endl;

            if(inconfig.m_mcgen_variation_allowlist.size()> 0 ){
                if(inconfig.m_mcgen_variation_allowlist.count(varname)==0){
                    std::cout<<"Skipping "<<varname<<" as its not in the AllowList!!"<<std::endl;
                    continue;
                }
            }

            if(inconfig.m_mcgen_variation_denylist.size()> 0 ){
                if(inconfig.m_mcgen_variation_denylist.count(varname)>0){
                    std::cout<<"Skipping "<<varname<<" as it is the DenyList!!"<<std::endl;
                    continue;
                }
            }

            //Variation Weight Maps Area
            int n_wei = inconfig.m_mcgen_weightmaps_patterns.size();
            std::string variation_mode;
            std::string s_formula = "1";

            for(int i=0; i< n_wei; i++){
                // Check to see if pattern is in this variation
                if (varname.find(inconfig.m_mcgen_weightmaps_patterns[i]) != std::string::npos) {
                    std::cout << "Variation "<<varname<<" is a match for pattern "<<inconfig.m_mcgen_weightmaps_patterns[i]<<std::endl;
                    s_formula = s_formula + "*(" + inconfig.m_mcgen_weightmaps_formulas[i]+")";
                    std::cout<<" -- weight is thus "<<s_formula<<std::endl;
                    std::cout<<" -- mode is "<<inconfig.m_mcgen_weightmaps_mode[i]<<std::endl;
                    variation_mode=inconfig.m_mcgen_weightmaps_mode[i];
                }
            }

            std::cout << varname << " has " << map_systematic_num_universe[varname] << " montecarlos in file " << fid << std::endl;
            used_montecarlos[fid] += map_systematic_num_universe[varname];
            syst_vector.emplace_back(varname, map_systematic_num_universe[varname], variation_mode, s_formula, knobvals[v], pset_indices[v]);
        }
    } // end fid
    
    std::vector<std::vector<std::unique_ptr<TTreeFormula>>> additional_weight_formulas(num_files);
    for (auto& innerVector : additional_weight_formulas) {
        innerVector.resize(variations.size());
    }

    for(int fid=0; fid < num_files; ++fid) {
        files[fid]->cd();
        for(int vid = 0; vid < syst_vector.size(); vid++){ 
            additional_weight_formulas[fid][vid] =  std::make_unique<TTreeFormula>(("weightMapsFormulas_"+std::to_string(fid)+"_"+std::to_string(vid)).c_str(), syst_vector[vid].formula.c_str(),trees[fid]);
        }
    }

    std::cout << " File: 0 | " << inconfig.m_mcgen_file_name.at(0) << " has " << used_montecarlos.at(0) << " montecarlos" << std::endl;
    for(int i=1; i<num_files; i++){
        std::cout << " File: " << i <<" |  "<<inconfig.m_mcgen_file_name[i]<< " has " << used_montecarlos.at(i) << " montecarlos" << std::endl;
        if(used_montecarlos.at(i)!= used_montecarlos.at(i-1)){
            std::cerr << " Warning, number of universes for are different between files" << std::endl;
            std::cerr << " The missing universes are Set to weights of 1. Make sure this is what you want!" << std::endl;
            for(int j=0; j<num_files; j++){
                if(universes_used < used_montecarlos.at(j)) 
                    universes_used = used_montecarlos.at(j);
                std::cerr <<"File " << j << " montecarlos: " << used_montecarlos.at(j) << std::endl;
            }
        }
    }
    //But in reality we want the max universes to be the sum of all max variaitons across all files, NOT the sum over all files max variations.
    //WHAT
    //universes_used = num_universes_per_variation.size();

    std::cout << " -------------------------------------------------------------" << std::endl;
    std::cout << " Initilizing " << universes_used << " universes." << std::endl;
    std::cout << " -------------------------------------------------------------" << std::endl;

    
    std::cout << " Full concatanated vector has : " << inconfig.m_num_bins_total << std::endl;
    //FIX Replace With Eigen
    std::vector<std::vector<double>> multi_vecspec;
    std::vector<double> base_vec (inconfig.m_num_bins_total,0.0);
    multi_vecspec.clear();
    multi_vecspec.resize(universes_used,base_vec);

    std::cout << " multi_vecspec now initilized of size :" << multi_vecspec.size() << std::endl;
    std::cout << " Reading the data files" << std::endl;

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

