#include "PROtocall.h"

namespace PROfit{

int FindGlobalBin(const PROconfig &inconfig, double reco_value, const std::string& subchannel_fullname){
    int subchannel_index = inconfig.GetSubchannelIndex(subchannel_fullname);
    return FindGlobalBin(inconfig, reco_value, subchannel_index);
}

int FindGlobalBin(const PROconfig &inconfig, double reco_value, int subchannel_index){
    int global_bin_start = inconfig.GetGlobalBinStart(subchannel_index);
    int channel_index = inconfig.GetChannelIndex(subchannel_index);
    int local_bin = FindLocalBin(inconfig, reco_value, channel_index);
    return local_bin == -1 ? -1 : global_bin_start + local_bin;
}


int FindLocalBin(const PROconfig &inconfig, double reco_value, int channel_index){
    
    //find local bin 
    const std::vector<double>& bin_edges = inconfig.GetChannelBinEdges(channel_index);
    auto pos_iter = std::upper_bound(bin_edges.begin(), bin_edges.end(), reco_value);

    //over/under-flow, don't care for now
    if(pos_iter == bin_edges.end() || pos_iter == bin_edges.begin()){
	log<LOG_DEBUG>(L"%1% || Reco value: %2% is in underflow or overflow bins, return bin of -1") % __func__ % reco_value;
	log<LOG_DEBUG>(L"%1% || Channel %2% has bin lower edge: %3% and bin upper edge: %4%") % __func__ % channel_index % *bin_edges.begin() % bin_edges.back();
	return -1; 
    }
    return pos_iter - bin_edges.begin() - 1; 
}

int FindGlobalTrueBin(const PROconfig &inconfig, double true_value, const std::string& subchannel_fullname){
    int subchannel_index = inconfig.GetSubchannelIndex(subchannel_fullname);
    return FindGlobalTrueBin(inconfig, true_value, subchannel_index);
}

int FindGlobalTrueBin(const PROconfig &inconfig, double true_value, int subchannel_index){
    int global_bin_start = inconfig.GetGlobalTrueBinStart(subchannel_index);
    int channel_index = inconfig.GetChannelIndex(subchannel_index);
    int local_bin = FindLocalTrueBin(inconfig, true_value, channel_index);
    return local_bin == -1 ? -1 : global_bin_start + local_bin;
}


int FindLocalTrueBin(const PROconfig &inconfig, double true_value, int channel_index){
    
    //find local bin 
    const std::vector<double>& bin_edges = inconfig.GetChannelTrueBinEdges(channel_index);
    auto pos_iter = std::upper_bound(bin_edges.begin(), bin_edges.end(), true_value);

    //over/under-flow, don't care for now
    if(pos_iter == bin_edges.end() || pos_iter == bin_edges.begin()){
	log<LOG_DEBUG>(L"%1% || True value: %2% is in underflow or overflow bins, return bin of -1") % __func__ % true_value;
	log<LOG_DEBUG>(L"%1% || Channel %2% has bin lower edge: %3% and bin upper edge: %4%") % __func__ % channel_index % *bin_edges.begin() % bin_edges.back();
	return -1; 
    }
    return pos_iter - bin_edges.begin() - 1; 
}

int FindSubchannelIndexFromGlobalBin(const PROconfig &inconfig, int global_bin, bool reco_bin ){
   if(reco_bin)
   	return inconfig.GetSubchannelIndexFromGlobalBin(global_bin);
   else
	return inconfig.GetSubchannelIndexFromGlobalTrueBin(global_bin);
}

};
