#include "fold.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;
//python bindings
PYBIND11_MODULE(fold, var)
{
    var.doc() = "Classes and functions needed for ScanFold-Fold";

}