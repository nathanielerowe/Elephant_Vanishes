#ifndef PRO_TO_CALL_H
#define PRO_TO_CALL_H

// C++ include 
#include <algorithm>
#include <unordered_map>
#include <string>

// PROfit include 
#include "PROlog.h"
#include "PROconfig.h"

namespace PROfit{

    class IndexMapper{
	private:
	    const PROconfig* mp_config;
            std::unordered_map<std::string, int> m_map_fullname_subchannel_index;
            std::unordered_map<int, long int> m_map_subchannel_index_to_global_index_start;
            std::unordered_map<int, int> m_map_subchannel_index_to_channel_index;


	    /* Function: fill in mapping between subchannel name/index to global indices */
	    void generate_index_map();

	    /* Function: reset everything in the class */
	    void clean();

	    /* Function: simply check if mp_config points to anything */
	    void sanity_check();
	public:

	    IndexMapper():mp_config(nullptr){}
	    IndexMapper(const PROconfig &inconfig);


	    /* Function: associate IndexMapper with another PROconfig, and update mapping
 	     * Warning: every binning mapping for previous config will be lost in this case 
 	     */
	    void ChangePROconfig(const PROconfig &in_config);


	    /* Function: given a value for reconstructed variable, figure out which local bin in the histogram it belongs to
 	     * Note: bin index start from 0, not 1
 	     * Note: return value of -1 means the reco value is out of range
	     */
	    int FindLocalBin(double reco_value, int channel_index) const;
	    int FindLocalBin(const PROconfig &inconfig, double reco_value, int channel_index);

	    /* Function: given a value for reconstructed variable, figure out which global bin it belongs to
 	     * Note: bin index start from 0, not 1
 	     * Note: return value of -1 means the reco value is out of range
	     */
	    long int FindGlobalBin(double reco_value, int subchannel_index) const;
	    long int FindGlobalBin(const PROconfig &inconfig, double reco_value, int subchannel_index);


	    /* Function: given subchannel full name, return global subchannel index 
 	     * Note: index start from 0, not 1
	     */
	    int GetSubchannelIndex(const std::string& fullname) const;
	    int GetSubchannelIndex(const PROconfig &inconfig, const std::string& fullname);
    };

};

#endif
