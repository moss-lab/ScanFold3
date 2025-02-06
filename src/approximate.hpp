#ifndef APPROXIMATE
#define APPROXIMATE
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "shared.hpp"
#include "basepair.hpp"
#include "matrix.hpp"
/*
Greedy algorithm to find approximate lowest z-score global structure in O(nlogn)+O(1.5n) time
*/
namespace approximate
{
    namespace py = pybind11;
    typedef std::shared_ptr<basepair::BasePair> base_pair_pointer;
    std::vector<base_pair_pointer> greedy_approximation(matrix::BasePairMatrix& mat);
    //python wrapper
    py::list py_greedy_approximation(matrix::BasePairMatrix& mat);
}
#endif