#include "PROconfig.h"
using namespace PROfit;

PROconfig::PROconfig(const std::string &xml){
    xmlname = xml;

    std::string line, text;
    std::ifstream in(xmlname);
    while(std::getline(in, line))  text += line + "\n";
    const char* xmldata = text.c_str();
    PROconfig::LoadFromXML(xmldata);
    num_modes = 10;
}


int PROconfig::LoadFromXML(const char * filedata){

    //Setup TiXml documents
    tinyxml2::XMLDocument doc;
    bool loadOkay = doc.Parse(filedata, 0);

    try{
        if(loadOkay) log<LOG_INFO>(L"%1% || Correctly loaded and parsed XML, continuing") % __func__;
        else throw 404;    
    }
    catch (int ernum) {
        log<LOG_ERROR>(L"%1% || ERROR: Failed to load XML configuration file. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"This generally means broken brackets or attribute syntax in xml itself.");
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }

    tinyxml2::XMLHandle hDoc(&doc);

    return 0;
}
