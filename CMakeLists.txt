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

# set source files
set(SOURCES ${SOURCES}
    src/test.cpp
    src/shared.cpp 
    src/window.cpp 
    src/basepair.cpp 
    src/matrix.cpp 
    src/flow.cpp 
)

# set binding file
set(BINDINGS
    src/fold.cpp
)
# create pybind11 module
#pybind11_add_module(fold src/fold.cpp)
add_library(fold MODULE "${SOURCES}" "${BINDINGS}")
target_link_libraries(fold PRIVATE pybind11::module)
# set output directory
set_target_properties(fold PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib/fold"
)
