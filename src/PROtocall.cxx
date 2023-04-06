#include "PROtocall.h"

using namespace PROfit;

IndexMapper::IndexMapper(const PROconfig &in_config): mp_config(&in_config){
    this->generate_index_map();
}

void IndexMapper::ChangePROconfig(const PROconfig &in_config){
    this->clean();
    mp_config = &in_config;
    this->generate_index_map();
    return;
}

int IndexMapper::FindLocalBin(double reco_value, int channel_index) const{
    this->sanity_check();

    //find local bin 
    const std::vector<double>& bin_edges = mp_config->m_channel_bin_edges.at(channel_index);
    auto pos_iter = std::upper_bound(bin_edges.begin(), bin_edges.end(), reco_value);

    //over/under-flow, don't care for now
    if(pos_iter == bin_edges.end() || pos_iter == bin_edges.begin())
	return -1; 
    return pos_iter - bin_edges.begin() - 1; 
}

int IndexMapper::FindLocalBin(const PROconfig &inconfig, double reco_value, int channel_index){
	
    log<LOG_WARNING>(L"%1% || Updating PROconfig object associated with IndexMapper") % __func__;
    log<LOG_WARNING>(L"%1% || Provided channel index %2% might be invalid.. ") % __func__ % channel_index;

    this->ChangePROconfig(inconfig);
    return this->FindLocalBin(reco_value, channel_index);
}

long int IndexMapper::FindGlobalBin(double reco_value, int sub_index) const {
    long int global_bin_start = m_map_subchannel_index_to_global_index_start.at(sub_index);
    int channel_index = m_map_subchannel_index_to_channel_index.at(sub_index);

    long int local_bin = FindLocalBin(reco_value, channel_index);

    return local_bin == -1 ? -1 : global_bin_start + local_bin;
}

long int IndexMapper::FindGlobalBin(const PROconfig &inconfig, double reco_value, int sub_index){
    this->ChangePROconfig(inconfig);
    return this->FindGlobalBin(reco_value, sub_index);
}

int IndexMapper::GetSubchannelIndex(const std::string& fullname) const{
    return m_map_fullname_subchannel_index.at(fullname);
}

int IndexMapper::GetSubchannelIndex(const PROconfig &inconfig, const std::string& fullname){
    this->ChangePROconfig(inconfig);
    return m_map_fullname_subchannel_index.at(fullname);
}

void IndexMapper::generate_index_map(){
    this->sanity_check();
    this->clean();

    int global_subchannel_index = 0;
    for(int im = 0; im < mp_config->m_num_modes; im++){

        long int mode_bin_start = im*mp_config->m_num_bins_mode_block;

        for(int id =0; id < mp_config->m_num_detectors; id++){

            long int detector_bin_start = id*mp_config->m_num_bins_detector_block;
            long int channel_bin_start = 0;

            for(int ic = 0; ic < mp_config->m_num_channels; ic++){
                for(int sc = 0; sc < mp_config->m_num_subchannels[ic]; sc++){

                    std::string temp_name  = mp_config->m_mode_names[im] +"_" +mp_config->m_detector_names[id]+"_"+mp_config->m_channel_names[ic]+"_"+mp_config->m_subchannel_names[ic][sc];
                    long int global_bin_index = mode_bin_start + detector_bin_start + channel_bin_start + sc*mp_config->m_channel_num_bins[ic];

		    m_map_fullname_subchannel_index[temp_name] = global_subchannel_index;
                    m_map_subchannel_index_to_global_index_start[global_subchannel_index] = global_bin_index;
                    m_map_subchannel_index_to_channel_index[global_subchannel_index] = ic;

                    ++global_subchannel_index;
                }
                channel_bin_start += mp_config->m_channel_num_bins[ic]*mp_config->m_num_subchannels[ic];
            }
        }
    }
    return;
}


void IndexMapper::clean(){
    //clean all maps 
    if(mp_config){
	m_map_fullname_subchannel_index.clear();
	m_map_subchannel_index_to_global_index_start.clear();
	m_map_subchannel_index_to_channel_index.clear();
    }

    return;
}

void IndexMapper::sanity_check(){
    if(!mp_config){
	log<LOG_ERROR>(L"%1% || IndexMapper not associated with PROconfig object, check if it exists or got destroyed ") % __func__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }
    log<LOG_DEBUG>(L"%1% || IndexMapper is associated with a object!!") % __func__;
    return;
}
