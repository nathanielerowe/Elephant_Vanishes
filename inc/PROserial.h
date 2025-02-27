#ifndef PROSERIEL_H
#define PROSERIEL_H


#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include <Eigen/Dense>

#include "PROcreate.h"

//extrnd boost to handle read/write serialization of Eigen::MatrixXf and Eigen::VectorXF
namespace boost {
    namespace serialization {

        template <class Archive>
            void serialize(Archive& ar, Eigen::MatrixXf& mat, [[maybe_unused]] const unsigned int version) {
                int rows = mat.rows(), cols = mat.cols();
                ar & rows & cols;
                if (Archive::is_loading::value) mat.resize(rows, cols);
                for (int i = 0; i < rows; ++i)
                    for (int j = 0; j < cols; ++j)
                        ar & mat(i, j);
            }

        template <class Archive>
            void serialize(Archive& ar, Eigen::VectorXf& vec, [[maybe_unused]] const unsigned int version) {
                int size = vec.size();
                ar & size;
                if (Archive::is_loading::value) vec.resize(size);
                for (int i = 0; i < size; ++i)
                    ar & vec(i);
            }

    } 
} 


namespace PROfit {


    void saveSystStructVector(const std::vector<SystStruct> &structs, const std::string &filename) {
            std::ofstream ofs(filename);
            boost::archive::text_oarchive oa(ofs);
            oa & structs;  
    }

    void loadSystStructVector(std::vector<SystStruct> &structs, const std::string &filename) {
            std::ifstream ifs(filename);
            boost::archive::text_iarchive ia(ifs);
            ia & structs;  
    }






}
#endif

