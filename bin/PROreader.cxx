#include "PROconfig.h"
#include "PROspec.h"
#include "PROcovariancegen.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/SRWeightPSet.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include "TFile.h"

#include <memory>
#include <string>

using namespace PROfit;
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[])
{

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string filename = "NULL.flat.caf.root"; 
    
    app.add_option("-f,--file", filename, "File to read systematics from.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);

    auto file = std::make_unique<TFile>(filename.c_str());
    TTree* globalTree = (TTree*)file->Get("globalTree");

    //auto global = std::make_unique<SRGlobal>();
    caf::SRGlobal* global = NULL;
    globalTree->SetBranchAddress("global", &global);

    globalTree->GetEntry(0);

    for(const auto& pset: global->wgts) {
    //for(unsigned int i = 0; i < global->wgts.size(); ++i) {
        //const caf::SRWeightPSet& pset = global->wgts[i];
        //std::cout << "i is: " << i << std::endl;
        if(pset.map.size() != 1) continue;
        std::cout << pset.name << " (type " << pset.type << "): with " << pset.nuniv << " universes\n";
        std::cout << pset.map.at(0).param.name << std::endl;
        for(const auto& val: pset.map.at(0).vals) {
            std::cout << val << ' ';
        }
        std::cout << std::endl;
    }
    

}


