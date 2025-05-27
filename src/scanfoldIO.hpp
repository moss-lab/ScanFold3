#ifndef scanfoldIO
#define scanfoldIO

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "basepair.hpp"
#include "shared.hpp"
namespace io 
{
        namespace py = pybind11;
        void list_to_ct(std::vector<basepair::BasePair>& base_pair_list, std::string out_file, double filter,
            char strand, std::string name, size_t start_coordinate, size_t end_coordinate);
        void py_list_to_ct(py::list& base_pair_list, std::string out_file, double filter,
            char strand, std::string name, size_t start_coordinate, size_t end_coordinate);
        void py_makedbn(py::list& pairs, std::string sequence, std::string header, std::string file_name);
}

#endif