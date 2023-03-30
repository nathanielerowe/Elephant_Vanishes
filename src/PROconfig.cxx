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













    return 0;
}
