//compiled w/ g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
//requires boost in include path, use 'sudo apt-get install libboost-all-dev'
/*
library header for ScanFold-Fold, contains python bindings
*/

#ifndef FOLD
#define FOLD
#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "basepair.hpp"
#include "flow.hpp"
#include "window.hpp"
#include "matrix.hpp"
#include "approximate.hpp"
#include "scanfoldIO.hpp"
#include "shared.hpp"
namespace py = pybind11;

//python bindings
PYBIND11_MODULE(fold, var)
{
    var.doc() = "Classes and functions needed for ScanFold-Fold";
    var.def("readScanTSV", &window::readScanTSV, "read the .tsv file produced by scanfold-scan into a C++ data structure");
    var.def("greedy_approximation", &approximate::py_greedy_approximation);
    var.def("findPairsInDotBracket", &shared::py_findPairsInDotBracket);
    var.def("getDBFromPairs", &shared::py_getDBFromPairs);
    var.def("filter_base_pairs", &basepair::py_filterBasePairs);
    py::class_<basepair::BasePair, std::shared_ptr<basepair::BasePair>>(var, "BasePair")
        .def(py::init<size_t, size_t>())
        .def(py::init<int, int, char, char, size_t>())
        .def(py::init<int, int, char, char, double, double, double, double, size_t, size_t, size_t>())
        .def(py::init<int, int, char, char, double, double, double, double, size_t, size_t>())
        .def("getZNorm", py::overload_cast<>(&basepair::BasePair::getZNorm))
        .def("getZNorm", py::overload_cast<>(&basepair::BasePair::getZNorm, py::const_))
        .def("getAvgZScore", &basepair::BasePair::getAvgZScore)
        .def("getMFENorm", &basepair::BasePair::getMFENorm)
        .def("getAvgMFE", &basepair::BasePair::getAvgMFE)  
        .def("getPValNorm", &basepair::BasePair::getPValNorm)
        .def("getAvgPVal", &basepair::BasePair::getAvgPVal)
        .def("getEDNorm", &basepair::BasePair::getEDNorm)
        .def("getAvgED", &basepair::BasePair::getAvgED)
        .def("getptr", &basepair::BasePair::getptr)
        //.def("print", py::overload_cast<>(&basepair::BasePair::print))
        .def("print", &basepair::BasePair::py_print)
        .def_readwrite("i_coord", &basepair::BasePair::icoord)
        .def_readwrite("j_coord", &basepair::BasePair::jcoord)
        .def_readwrite("i_nucleotide", &basepair::BasePair::inuc)
        .def_readwrite("j_nucleotide", &basepair::BasePair::jnuc)
    ;
    py::class_<window::ScanFoldWindow>(var, "ScanFoldWindow")
        .def(py::init<std::string&, size_t>())
        .def("getPairs", &window::ScanFoldWindow::getPairs)
        .def("print", &window::ScanFoldWindow::print)
        .def_readonly("Start", &window::ScanFoldWindow::Start)
        .def_readonly("End", &window::ScanFoldWindow::End)
        .def_readonly("Temperature", &window::ScanFoldWindow::Temperature)
        .def_readonly("NativeMFE", &window::ScanFoldWindow::NativeMFE)
        .def_readonly("Zscore", &window::ScanFoldWindow::Zscore)
        .def_readonly("pvalue", &window::ScanFoldWindow::pvalue)
        .def_readonly("Sequence", &window::ScanFoldWindow::Sequence)
        .def_readonly("Structure", &window::ScanFoldWindow::Structure)
        .def_readonly("centroid", &window::ScanFoldWindow::centroid)
    ; 
    py::class_<matrix::BasePairMatrix>(var, "BasePairMatrix")
        .def(py::init<std::string>())
        .def("getZNorm", &matrix::BasePairMatrix::getZNorm)
        .def("getAvgZScore", &matrix::BasePairMatrix::getAvgZScore)
        .def("getBestPairing", py::overload_cast<py::list &>(&matrix::BasePairMatrix::getBestPairing))
        .def("toCSV", &matrix::BasePairMatrix::toCSV)
        .def("getWinSize", &matrix::BasePairMatrix::getWinSize)
        .def("getStepSize", &matrix::BasePairMatrix::getStepSize)
        .def("get", &matrix::BasePairMatrix::get)
        .def("getSequenceLength", &matrix::BasePairMatrix::getSequenceLength)
        .def("get_all_pairs", &matrix::BasePairMatrix::py_getAllPairs)
        .def("getSequence", &matrix::BasePairMatrix::getSequence)
    ;
    py::class_<flow::FlowGraph>(var, "FlowGraph")
        .def(py::init<matrix::BasePairMatrix&>())
        .def("minimum_weighted_matching", py::overload_cast<py::list &, bool>(&flow::FlowGraph::minimum_weighted_matching))
    ;
}
#endif
