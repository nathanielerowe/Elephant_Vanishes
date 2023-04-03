#include "PROconfig.h"
using namespace PROfit;

PROconfig::PROconfig(const std::string &xml){
    m_xmlname = xml;
    PROconfig::LoadFromXML(m_xmlname);

    m_num_modes = 1;
    m_num_detectors=1;
    m_num_channels = 2;
    m_num_bins_total = 100;
    m_num_bins_total_collapsed = 10;

    //A matrix for collapsing the full-vector
    //left multiply this matrix by the full-vector to get collapsed vector
    //TODO: check corr=init for multi detector


    /*collapsingVector = Eigen::MatrixXd::Zero(num_bins_total,num_bins_total_collapsed);

      for(int im = 0; im < num_modes; im++){
      for(int id =0; id < num_detectors; id++){
      int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
      int corr = edge;
      for(int ic = 0; ic < num_channels; ic++){
      int corner=edge;
      for(int j=0; j< num_bins.at(ic); j++){
      for(int sc = 0; sc < num_subchannels.at(ic); sc++){
      int place = j+sc*num_bins[ic]+corner;
      collapsingVector(place,corr)=1;
      }
      }
      corr++;
      }
      }
      }*/

}


int PROconfig::LoadFromXML(const std::string &filename){


    //Setup TiXml documents
    tinyxml2::XMLDocument doc;
    doc.LoadFile(filename.c_str());
    bool loadOkay = !doc.ErrorID();

    bool use_universe = 1; //FIX
    try{
        if(loadOkay) log<LOG_INFO>(L"%1% || Correctly loaded and parsed XML, continuing") % __func__;
        else throw 404;    
    }
    catch (int ernum) {
        log<LOG_ERROR>(L"%1% || ERROR: Failed to load XML configuration file names %4%. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % filename.c_str();
        log<LOG_ERROR>(L"This generally means broken brackets or attribute syntax in xml itself.");
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }

    tinyxml2::XMLHandle hDoc(&doc);

    tinyxml2::XMLElement *pMode, *pDet, *pChan;

    //max subchannels 100? Can we avoid this
    m_subchannel_plotnames.resize(100);
    m_subchannel_datas.resize(100);
    m_subchannel_names.resize(100);
    m_subchannel_osc_patterns.resize(100);
    char *end;

    //Grab the first element. Note very little error checking here! make sure they exist.
    pMode = doc.FirstChildElement("mode");
    pDet =  doc.FirstChildElement("detector");
    pChan = doc.FirstChildElement("channel");

    std::cout<<pMode<<std::endl;
    if(!pMode){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 mode defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }else{
        while(pMode){
            // What modes are we running in (e.g nu, nu bar, horn current=XXvolts....) Can have as many as we want
            const char* mode_name= pMode->Attribute("name");
            if(mode_name==NULL){
                log<LOG_ERROR>(L"%1% || Modes need a name! Please define a name attribute for all modes. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                m_mode_names.push_back(mode_name);
            }

            const char* mode_plotname= pMode->Attribute("plotname");
            if(mode_plotname==NULL){
                m_mode_plotnames.push_back(m_mode_names.back());
            }else{
                m_mode_plotnames.push_back(mode_plotname);
            }

            pMode = pMode->NextSiblingElement("mode");
            log<LOG_DEBUG>(L"%1% || Loading Mode %2%  ") % __func__ % m_mode_names.back().c_str() ;

        }
    }

    // How many detectors do we want!
    if(!pDet){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 detector defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);

    }else{

        while(pDet){

            const char* detector_name= pDet->Attribute("name");
            if(detector_name==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: Need all detectors to have a name attribute @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                m_detector_names.push_back(detector_name);
            }

            const char* detector_plotname = pDet->Attribute("plotname");
            if(detector_plotname==NULL){
                m_detector_plotnames.push_back(m_detector_names.back());
            }else{
                m_detector_plotnames.push_back(detector_plotname);
            }

            pDet = pDet->NextSiblingElement("detector");
            log<LOG_DEBUG>(L"%1% || Loading Det %2%  ") % __func__ % m_detector_names.back().c_str();

        }
    }

    //How many channels do we want! At the moment each detector must have all channels
    int nchan = 0;
    if(!pChan){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 channel defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }else{


        while(pChan){
            // Read in how many bins this channel uses

            const char* channel_name= pChan->Attribute("name");
            if(channel_name==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: Need all channels to have names in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                m_channel_names.push_back(channel_name);
            }

            const char* channel_plotname= pChan->Attribute("plotname");
            if(channel_plotname==NULL){
                m_channel_plotnames.push_back(m_channel_names.back());
            }else{
                m_channel_plotnames.push_back(channel_plotname);
            }


            const char* channel_unit= pChan->Attribute("unit");
            if(channel_unit==NULL){
                m_channel_units.push_back("");
            }else{
                m_channel_units.push_back(channel_unit);
            }

            log<LOG_DEBUG>(L"%1% || Loading Channel %2% with   ") % __func__ % m_channel_names.back().c_str() ;


            // What are the bin edges and bin widths (bin widths just calculated from edges now)
            tinyxml2::XMLElement *pBin = pChan->FirstChildElement("bins");
            std::stringstream iss(pBin->Attribute("edges"));

            double number;
            std::vector<double> binedge;
            std::vector<double> binwidth;
            std::string binstring = "";
            while ( iss >> number ){
                binedge.push_back( number );
                binstring+=" "+std::to_string(number);
            }

            log<LOG_DEBUG>(L"%1% || Loading Bins with edges %2%  ") % __func__ % binstring.c_str();

            for(int b = 0; b<binedge.size()-1; b++){
                binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
            }

            m_channel_num_bins.push_back(binedge.size()-1);

            m_channel_bin_edges.push_back(binedge);
            m_channel_bin_widths.push_back(binwidth);


            // Now loop over all this channels subchanels. Not the names must be UNIQUE!!
            tinyxml2::XMLElement *pSubChan;

            pSubChan = pChan->FirstChildElement("subchannel");
            int nsubchan=0;
            while(pSubChan){

                const char* subchannel_name= pSubChan->Attribute("name");
                if(subchannel_name==NULL){
                    log<LOG_ERROR>(L"%1% || ERROR: Subchannels need a name in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);

                }else{
                    m_subchannel_names[nchan].push_back(subchannel_name);
                    log<LOG_DEBUG>(L"%1% || Subchannel Starting:  %2%") % __func__ % m_subchannel_names.at(nchan).back().c_str() ;

                }


                const char* subchannel_plotname= pSubChan->Attribute("plotname");
                if(subchannel_plotname==NULL){
                    m_subchannel_plotnames[nchan].push_back(m_subchannel_names[nchan].back());
                }else{
                    m_subchannel_plotnames[nchan].push_back(subchannel_plotname);
                }

                const char* subchannel_data= pSubChan->Attribute("data");
                if(subchannel_data==NULL){
                    m_subchannel_datas[nchan].push_back(0);
                }else{
                    m_subchannel_datas[nchan].push_back(1);
                }

                //0 means dont oscillate, 11 means electron disapearance, -11 means antielectron dis..etc..
                if(pSubChan->Attribute("osc"))
                {
                    m_has_oscillation_patterns = true;
                    m_subchannel_osc_patterns.at(nchan).push_back(strtod(pSubChan->Attribute("osc"), &end));
                }else{
                    m_has_oscillation_patterns = false;
                    m_subchannel_osc_patterns.at(nchan).push_back(0);
                }

                log<LOG_DEBUG>(L"%1% || Subchannel %2% with and osc pattern %3% and isdata %4%") % __func__ % m_subchannel_names.at(nchan).back().c_str() % m_subchannel_osc_patterns.at(nchan).back() % m_subchannel_datas.at(nchan).back();

                nsubchan++;
                pSubChan = pSubChan->NextSiblingElement("subchannel");
            }
            m_num_subchannels.push_back(nsubchan);

            nchan++;
            pChan = pChan->NextSiblingElement("channel");
        }
    }//end channel loop



    //Now onto mcgen, for CV specs or for covariance generation
    tinyxml2::XMLElement *pMC, *pWeiMaps, *pList, *pSpec, *pShapeOnlyMap;
    pMC   = doc.FirstChildElement("MCFile");
    pWeiMaps = doc.FirstChildElement("WeightMaps");
    pList = doc.FirstChildElement("variation_list");
    pSpec = doc.FirstChildElement("varied_spectrum");
    pShapeOnlyMap = doc.FirstChildElement("ShapeOnlyUncertainty");

    if(pMC){//Skip if not in XML
        while(pMC)
        {
            const char* tree = pMC->Attribute("treename");
            if(tree==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: You must have an associated root TTree name for all MonteCarloFile tags.. eg. treename='events' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                m_mcgen_tree_name.push_back(tree);
            }

            const char* file = pMC->Attribute("filename");
            if(file==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: You must have an associated root TFile name for all MonteCarloFile tags.. eg. filename='my.root' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                m_mcgen_file_name.push_back(file);
            }


            const char* maxevents = pMC->Attribute("maxevents");
            if(maxevents==NULL){
                m_mcgen_maxevents.push_back(1e16);
            }else{
                m_mcgen_maxevents.push_back(strtod(maxevents,&end) );
            }

            //arbitray scaling you can have
            const char* scale = pMC->Attribute("scale");
            if(scale==NULL){
                m_mcgen_scale.push_back(1.0);
            }else{
                m_mcgen_scale.push_back(strtod(scale,&end) );
            }

            const char* inpot = pMC->Attribute("pot");
            if(inpot==NULL){
                m_mcgen_pot.push_back(-1.0);
            }else{
                m_mcgen_pot.push_back(strtod(inpot,&end) );
            }

            //Is this useful? 
            const char* isfake = pMC->Attribute("fake");
            if(isfake==NULL){
                m_mcgen_fake.push_back(false);
            }else{
                m_mcgen_fake.push_back(true);
            }

            log<LOG_DEBUG>(L"%1% || MultisimFile %2%, treename: %3%  ") % __func__ % m_mcgen_file_name.back().c_str() % m_mcgen_tree_name.back().c_str();


            //Here we can grab some friend tree information
            tinyxml2::XMLElement *pFriend;
            pFriend = pMC->FirstChildElement("friend");
            while(pFriend){


                std::string ffname;
                const char* friend_filename = pFriend->Attribute("filename");
                if(friend_filename==NULL){
                    ffname = m_mcgen_file_name.back();
                }else{
                    ffname = friend_filename;
                }

                if(m_mcgen_file_friend_treename_map.count(m_mcgen_file_name.back())>0){

                    (m_mcgen_file_friend_treename_map[m_mcgen_file_name.back()]).push_back( pFriend->Attribute("treename") );
                    (m_mcgen_file_friend_map[m_mcgen_file_name.back()]).push_back(ffname);

                }else{
                    std::vector<std::string> temp_treename = {pFriend->Attribute("treename")};
                    std::vector<std::string> temp_filename = {ffname};

                    m_mcgen_file_friend_treename_map[m_mcgen_file_name.back()] = temp_treename;
                    m_mcgen_file_friend_map[m_mcgen_file_name.back()] = temp_filename;
                }
                pFriend = pFriend->NextSiblingElement("friend");
            }//END of friend loop


            tinyxml2::XMLElement *pBranch;
            pBranch = pMC->FirstChildElement("branch");

            std::vector<BranchVariable*> TEMP_branch_variables;
            while(pBranch){

                const char* bnam = pBranch->Attribute("name");
                const char* btype = pBranch->Attribute("type");
                const char* bhist = pBranch->Attribute("associated_subchannel");
                const char* bsyst = pBranch->Attribute("associated_systematic");
                const char* bcentral = pBranch->Attribute("central_value");
                const char* bwname = pBranch->Attribute("eventweight_branch_name");
                const char* badditional_weight = pBranch->Attribute("additional_weight");

                if(bwname== NULL){
                    log<LOG_WARNING>(L"%1% || WARNING: No eventweight branch name passed, defaulting to 'weights' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    m_mcgen_eventweight_branch_names.push_back("weights");
                }else{
                    log<LOG_DEBUG>(L"%1% || Setting eventweight branch name %2%") %__func__ % bnam;
                    m_mcgen_eventweight_branch_names.push_back(std::string(bwname));
                }

                if(bnam == NULL){
                    log<LOG_ERROR>(L"%1% || ERROR!: Each branch must include the name of the branch to use. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"%1% || ERROR!: e.g name = 'ereco' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }

                if(btype == NULL){
                    log<LOG_WARNING>(L"%1% || WARNING: No branch type has been specified, assuming double. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    btype= "double";
                }

                if(bhist == NULL){
                    log<LOG_ERROR>(L"%1% || Each branch must have an associated_subchannel to fill! On branch %4% : @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % bnam;
                    log<LOG_ERROR>(L"%1% || e.g associated_subchannel='mode_det_chan_subchannel ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }


                if(bsyst == NULL){
                    if(use_universe == false){
                        log<LOG_WARNING>(L"%1% || WARNING: No root file with unique systematic variation is provided ") % __func__;
                        log<LOG_ERROR>(L"%1% || ERROR! please provide what systematic variation this file correpsonds to!") % __func__;
                        log<LOG_ERROR>(L"Terminating.");
                        exit(EXIT_FAILURE);
                    }
                    systematic_name.push_back("");
                }else{
                    systematic_name.push_back(bsyst);	

                }


                if(badditional_weight == NULL){
                    m_mcgen_additional_weight_bool.push_back(0);
                    m_mcgen_additional_weight_name.push_back("");
                }else{
                    m_mcgen_additional_weight_name.push_back(badditional_weight);
                    m_mcgen_additional_weight_bool.push_back(1);
                    log<LOG_DEBUG>(L"%1% || Setting an additional weight for branch %2% using the branch %3% as a reweighting.") % __func__ % bnam %badditional_weight;

                }



                if((std::string)btype == "double"){
                    if(use_universe){
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist ) );
                    } else  if((std::string)bcentral == "true"){
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist,bsyst, true) );
                        log<LOG_DEBUG>(L"%1% || Setting as  CV for det sys.") % __func__ ;
                    } else {
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist,bsyst, false) );
                        log<LOG_DEBUG>(L"%1% || Setting as individual (not CV) for det sys.") % __func__ ;
                    }
                }else{
                    log<LOG_ERROR>(L"%1% || ERROR: currently only double, allowed for input branch variables (sorry!) i @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % bnam;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }

                std::string oscillate = "false";
                if(pBranch->Attribute("oscillate")!=NULL){
                    oscillate =pBranch->Attribute("oscillate");
                }	

                if(oscillate == "false"){
                    log<LOG_DEBUG>(L"%1% || Oscillations are OFF ") % __func__ ;
                    TEMP_branch_variables.back()->SetOscillate(false);
                }else if(oscillate == "true"){
                    log<LOG_DEBUG>(L"%1% || Oscillations are Set to  ON ") % __func__;
                    TEMP_branch_variables.back()->SetOscillate(true);
                    TEMP_branch_variables.back()->true_param_name = pBranch->Attribute("true_param_name");
                    if(pBranch->Attribute("true_L_name") != NULL){
                        //for oscillation that needs both E and L
                        TEMP_branch_variables.back()->true_L_name = pBranch->Attribute("true_L_name");
                        log<LOG_DEBUG>(L"%1% || Oscillations using true param name:   %2% and baseline %3% ") % __func__ % pBranch->Attribute("true_param_name") % pBranch->Attribute("true_L_name") ;
                    }else{
                        //for oscillations that only needs E, such as an energy-dependent scaling for single photon NCpi0!
                        log<LOG_DEBUG>(L"%1% || Oscillations using  Energy only dependent oscillation ( or shift/normalization)  %2% ") % __func__ % pBranch->Attribute("true_param_name") ;
                    }
                }else{
                    log<LOG_DEBUG>(L"%1% || Do Not Oscillate  ") % __func__  ;
                    TEMP_branch_variables.back()->SetOscillate(false);
                }

                log<LOG_DEBUG>(L"%1% || Associated subchannel: %2% ") % __func__ % bhist;

                pBranch = pBranch->NextSiblingElement("branch");
            }
            m_branch_variables.push_back(TEMP_branch_variables);
            //next file
            pMC=pMC->NextSiblingElement("MCFile");
        }
    }

    if(!pList){
        log<LOG_DEBUG>(L"%1% || No Allowlist or Denylist set, including ALL variations by default.") % __func__  ;
    }else{
        while(pList){

            tinyxml2::XMLElement *pAllowList = pList->FirstChildElement("allowlist");
            while(pAllowList){
                std::string wt = std::string(pAllowList->GetText());
                m_mcgen_variation_allowlist[wt] = true; 
                log<LOG_DEBUG>(L"%1% || Allowlisting variations: %2%") % __func__ % wt.c_str() ;
                pAllowList = pAllowList->NextSiblingElement("allowlist");
            }

            tinyxml2::XMLElement *pDenyList = pList->FirstChildElement("denylist");
            while(pDenyList){
                std::string bt = std::string(pDenyList->GetText());
                m_mcgen_variation_denylist[bt] = true; 
                log<LOG_DEBUG>(L"%1% || Denylisting variations: %2%") % __func__ % bt.c_str() ;
                pDenyList = pDenyList->NextSiblingElement("denylist");
            }
            pList = pList->NextSiblingElement("variation_list");
        }
    }
    //weightMaps
    if(!pWeiMaps){
        log<LOG_DEBUG>(L"%1% || WeightMaps not set, all weights for all variations are 1 (individual branch weights still apply)") % __func__  ;
    }else{
        while(pWeiMaps){


            tinyxml2::XMLElement *pVariation;
            pVariation = pWeiMaps->FirstChildElement("variation");
    
            while(pVariation){


                const char* w_pattern = pVariation->Attribute("pattern");
                const char* w_formula = pVariation->Attribute("weight_formula");
                const char* w_use = pVariation->Attribute("use");
                const char* w_mode = pVariation->Attribute("mode");

                if(w_pattern== NULL){
                    log<LOG_ERROR>(L"%1% || ERROR! No pattern passed for this variation in WeightMaps. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }else{
                    log<LOG_DEBUG>(L"%1% || Loading WeightMaps Variation Pattern: %2%") %__func__ % w_pattern;
                    m_mcgen_weightmaps_patterns.push_back(std::string(w_pattern));
                }


                if(w_formula== NULL){
                    log<LOG_WARNING>(L"%1% || Warning, No formula passed for this variation in WeightMaps. Setting to 1. Make sure this is wanted behaviour.") %__func__ ;
                    m_mcgen_weightmaps_formulas.push_back("1");
                }else{
                    log<LOG_DEBUG>(L"%1% || Loading WeightMaps Variation Formula: %2%") %__func__ % w_formula;
                    m_mcgen_weightmaps_formulas.push_back(std::string(w_formula));
                }

                if(w_use== NULL){
                    m_mcgen_weightmaps_uses.push_back("true");
                }else{
                    m_mcgen_weightmaps_uses.push_back(std::string(w_use));
                }

                if(w_mode== NULL){
                    log<LOG_WARNING>(L"%1% || Warning, No mode passed for this variaiton in  WeightMaps. Assuming default multisim.  Make sure this is wanted behaviour.") %__func__ ;
                    m_mcgen_weightmaps_mode.push_back("multisim");
                }else{

                    std::string mode = std::string(w_mode);
                    if(mode=="multisim" || mode=="minmax"){
                        m_mcgen_weightmaps_mode.push_back(mode);
                    }else{
                        log<LOG_ERROR>(L"%1% || ERROR! The mode passed in is %4% but only allowed is multisim or minmax. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % w_mode;
                        log<LOG_ERROR>(L"Terminating.");
                        exit(EXIT_FAILURE);
                    }
                }

                pVariation = pVariation->NextSiblingElement("variation");
            }

            pWeiMaps=pWeiMaps->NextSiblingElement("WeightMaps");
        }
    }


    if(!pShapeOnlyMap){

    }else{
        while(pShapeOnlyMap){

            log<LOG_WARNING>(L"%1% || Warning!  Not setting up for shape-only covariance matrix generation. MAKE SURE this is what you want if you're generating covariance matrix!!!") % __func__;
            
            std::string pshapeonly_systematic_name = std::string(pShapeOnlyMap->Attribute("name"));
            const char* pshapeonly_systematic_use = pShapeOnlyMap->Attribute("use");
            bool pshapeonly_systematic_use_bool = true;

            if(pshapeonly_systematic_use == NULL || std::string(pshapeonly_systematic_use) == "true"){
                std::cout << "" << pshapeonly_systematic_name << std::endl;
                log<LOG_DEBUG>(L"%1% || Setting up shape-only covariance matrix for systematic: %2% ") % __func__ % pshapeonly_systematic_name.c_str();

            }else if(std::string(pshapeonly_systematic_use) == "false"){
                log<LOG_DEBUG>(L"%1% || Setting up shape-only covariance matrix for systematic: %2% ? False ") % __func__ % pshapeonly_systematic_name.c_str();
                pshapeonly_systematic_use_bool = false;
            }else{
                log<LOG_WARNING>(L"%1% || INVALID argument received for Attribute use of ShapeOnlyUncertainty element for systematic: %2% . Default it to true ") % __func__ % pshapeonly_systematic_name.c_str();
            }

            tinyxml2::XMLElement *pSubchannel;
            pSubchannel = pShapeOnlyMap->FirstChildElement("subchannel");	

            while(pshapeonly_systematic_use_bool && pSubchannel){

                const char* pshapeonly_subchannel_name = pSubchannel->Attribute("name");
                const char* pshapeonly_subchannel_use = pSubchannel->Attribute("use");

                if(pshapeonly_subchannel_use && std::string(pshapeonly_subchannel_use) == "false" ){
                    std::cout << " Subchannel " << std::string(pshapeonly_subchannel_name) << " is not included " << std::endl;
                }else{
                    std::cout <<  " Subchannel " << std::string(pshapeonly_subchannel_name) << " is included " << std::endl;
                    m_mcgen_shapeonly_listmap[pshapeonly_systematic_name].emplace_back(pshapeonly_subchannel_name);
                }

                pSubchannel = pSubchannel->NextSiblingElement("subchannel");
            }

            pShapeOnlyMap = pShapeOnlyMap->NextSiblingElement("ShapeOnlyUncertainty");
        }

        std::cout<< " Finish setting up systematics and subchannels for shape-only covariance matrix generation " << std::endl;
    }


    while(pSpec){
        const char* swrite_out = pSpec->Attribute("writeout");
        const char* swrite_out_tag = pSpec->Attribute("writeout_tag");
        const char* sform_matrix = pSpec->Attribute("form_matrix");	

        if( std::string(swrite_out) == "true") m_write_out_variation = true;
        if(m_write_out_variation){
            if(swrite_out_tag == NULL) m_write_out_tag="UNSET";
            else m_write_out_tag = std::string(swrite_out_tag);
        }

        if( std::string(sform_matrix) == "false") m_form_covariance = false;

        pSpec = pSpec->NextSiblingElement("varied_spectrum");
    }












    return 0;
}
