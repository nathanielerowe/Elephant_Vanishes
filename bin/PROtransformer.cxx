#include "PROconfig.h"
#include "PROspec.h"
#include "PROcovariancegen.h"
#include "PROcreate.h"
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
    std::string flat_filename = "NULL.flat.caf.root"; 
    std::string selection_filename = "NULL.flat.caf.root"; 
    std::string destination_file = "NULLDEFAULT";

    app.add_option("-f,--flatfile", flat_filename, "File to read systematics from.");
    app.add_option("-s,--selectionfile", selection_filename, "File after neutrino selection.");
    app.add_option("-o,--output", destination_file, "Output file");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);
 
    PROcess_WeightMap(flat_filename, selection_filename, destination_file);

    return 0;
}
