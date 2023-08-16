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
#include <cstdlib>
#include <highfive/H5File.hpp>

using namespace PROfit;
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[]){
 
  CLI::App app{"Test for PROfit"};
  // Define options
  std::string infilename = "NULL.flat.caf.root";
  std::string outfilename = "null.h5";
  app.add_option("-o,--outfile", outfilename, "H5 file to output to.");
  app.add_option("-i, --infile", infilename, "root file input");
  app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");
 
  CLI11_PARSE(app, argc, argv);

  //open flat caf file
  auto flat_caf = std::make_unique<TFile>(infilename.c_str(),"read");
  std::cout<< "opened file!" << std::endl;
  TTree *mytree = (TTree*)flat_caf->Get("nutree");
  std::cout<< "defined tree!"<<std::endl;
  double recoE = 0;
  double trueE =0;
  double  baseline=0;
  
  mytree->SetBranchAddress("trueE", &trueE);
  mytree->SetBranchAddress("recoE", &recoE);
  mytree->SetBranchAddress("baseline", &baseline);

  std::cout << "set branch addresses" << std::endl;
  //std::vector<float> v_trueE;
  std::vector<float> v_recoE;
  //std::vector<float> v_baseline;
  std::vector<float> L_over_E;

  for(int i = 0 ; i != mytree->GetEntries(); ++i){
      mytree->GetEntry(i);
      //v_trueE.push_back(trueE);
      v_recoE.push_back(recoE);
      //v_baseline.push_back(baseline);
      L_over_E.push_back(baseline/recoE);
  }
  //}

  HighFive::File myfile(outfilename, HighFive::File::Truncate);
  //myfile.createDataSet("nugroup/trueE", trueE);
  myfile.createDataSet("nugroup/recoE", v_recoE);
  //myfile.createDataSet("nugroup/baseline", baseline);
  myfile.createDataSet("nugroup/LE", L_over_E);
  /*auto dset1 = myfile.getDataSet("/nugroup/trueE/")
  auto trueE = dset1.read<std::vector<float>>

  //}
  //PROcess_WeightMap(flat_filename, selection_filename, destination_file);
  
    return 0;
    
}
