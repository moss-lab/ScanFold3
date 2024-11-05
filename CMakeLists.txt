cmake_minimum_required(VERSION 3.15...3.30)

# project version
project(fold VERSION 1.0)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# get required packages
find_package(pybind11 REQUIRED)
find_package(PythonInterp 3)
find_package(PythonLibs   3)
include_directories(${PYTHON_INCLUDE_DIRS})

# create pybind11 module
pybind11_add_module(fold src/fold.cpp)

# set output directory
set_target_properties(fold PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib/fold"
)
