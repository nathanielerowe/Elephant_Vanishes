#ifndef PROCONFIG_H_
#define PROCONFIG_H_

// STANDARD
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

// TINYXML2
#include "tinyxml2.h"


namespace PROfit{

    class PROconfig {
        protected:
        public:

        
            PROconfig() {}; //always have an empty?
            PROconfig(const std::string &xml);
            
            int LoadFromXML(const char * filedata);

            std::vector<std::string> fullnames;
            std::vector<int> vec_is_data;

            int num_detectors;
            int num_channels;
            int num_modes;

            std::string xmlname;	


    };


    

}
#endif
