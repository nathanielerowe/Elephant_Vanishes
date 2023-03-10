#include "PROconfig.h"
#include "PROspec.h"

#include "CLI11.h"

using namespace PROfit;

int main(int argc, char* argv[])
{

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 100;

    //doubles
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");

    CLI11_PARSE(app, argc, argv);



    PROconfig myConf(xmlname);
    PROspec mySpec(myConf);

    TH1D hmm = mySpec.toTH1D();


    return 0;
}

