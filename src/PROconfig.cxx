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

    tinyxml2::XMLHandle hDoc(&doc);


    return 0;
}
