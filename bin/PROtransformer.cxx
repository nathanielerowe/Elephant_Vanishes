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
#include <hdf5.h>

using namespace PROfit;
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

int main(int argc, char* argv[]){
 
  CLI::App app{"Test for PROfit"};
  // Define options
  std::string filename = "NULL.flat.caf.root";
  std::string loadtype = "NONE";
  app.add_option("-f,--file", filename, "File to read systematics from.");
  app.add_option("-t,--type", loadtype, "h5 or root?");
  app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");
 
  CLI11_PARSE(app, argc, argv);


  if(loadtype=="root"){
  //open flat caf file
  auto flat_caf = std::make_unique<TFile>(filename.c_str(),"read");
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

  for(int i = 0 ; i != mytree->GetEntries(); ++i){
      mytree->GetEntry(i);
  }
  }

  else if(loadtype=="h5"){

    const char* groupName = "nugroup"; // Group containing the dataset
    const char* datasetName = "trueE";

    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_id = H5Gopen2(file_id, groupName, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen2(group_id, datasetName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    int num_dimensions = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t dimensions[num_dimensions];
    H5Sget_simple_extent_dims(dataspace_id, dimensions, NULL);
    size_t data_size = H5Tget_size(H5Dget_type(dataset_id));
    void* data_buffer = malloc(data_size);

    for (hsize_t i = 0; i < dimensions[0]; ++i) {
        hsize_t start[] = {i, 0, 0};
        //H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, dimensions, NULL);
        H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, dataspace_id, H5P_DEFAULT, data_buffer);
        printf("Entry %llu: %s\n", (unsigned long long)i, (char*)data_buffer);
    }

    free(data_buffer);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

  }
  //PROcess_WeightMap(flat_filename, selection_filename, destination_file);
  
    return 0;
    
}
