#ifndef PRO_TO_CALL_H
#define PRO_TO_CALL_H

namespace PROfit{

    /* Function: given a value for reconstructed variable, figure out which local bin in the histogram it belongs to
     */
    int FindLocalBin(const PROconfig &inconfig, double value, int channel_index);

    /* Function: given a value for reconstructed variable, figure out which global bin it belongs to
     */
    int FindLocalBin(const PROconfig &inconfig, double value, int channel_index);

};

#endif
