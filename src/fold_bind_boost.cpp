#include "fold.h"
//python bindings
BOOST_PYTHON_MODULE(fold)
{
    using namespace boost::python;
    Py_Initialize();
    python::class_<BasePair>("BasePair")
        .def(init<size_t win, size_t len>())
        .def(init<int i, int j, char ival, char jval, size_t len>())
        .def(init<int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len, size_t num>())
        .def(init<int i, int j, char ival, char jval, double z, double m, double e, double p, size_t win, size_t len>())
        .def("update", &BasePair::update)
        .def("getZNorm", &BasePair::getZNorm)
        .def("getAvgZScore", &BasePair::getAvgZScore)
        .def("getMFENorm", &BasePair::getMFENorm)
        .def("getAvgMFE", &BasePair::getAvgMFE)
        .def("getPValNorm", &BasePair::getPValNorm)
        .def("getAvgPVal", &BasePair::getAvgPVal)
        .def("getEDNorm", &BasePair::getEDNorm)
        .def("getAvgED", &BasePair::getAvgED)
        .def_readonly("i_coord", &BasePair::icoord)
        .def_readonly("j_coord", &BasePair::jcoord)
        .def_readonly("i_nucleotide", &BasePair::inuc)
        .def_readonly("j_nucleotide", &BasePair::jnuc)
    ;
    python::class_<BasePairMatrix>("BasePairMatrix")
        .def(init<std::string tsv_name>()) 
        .ded("update", &BasePairMatrix::update)
        .def("getZNorm", &BasePairMatrix::getZNorm)
        .def("getAvgZScore", &BasePairMatrix::getAvgZScore)
        .def("getMaxZNorm", &BasePairMatrix::getMaxZNorm)
        .def("get", &BasePairMatrix::get)
        .def<void (BasePairMatrix::*)(python::list)>("getBestPairing", &BasePairMatrix::getBestPairing)
        .def("toCSV", &BasePairMatrix::toCSV)
    ;
}