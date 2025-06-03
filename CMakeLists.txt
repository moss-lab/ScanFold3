cmake_minimum_required(VERSION 3.15...3.30)

# project version
project(fold VERSION 1.0)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(PYBIND11_PYTHON_VERSION 3)
#set(PYBIND11_FINDPYTHON ON)
#find_package(pybind11 CONFIG)
# get required packages
#add_library(pybind11::pybind11 INTERFACE IMPORTED)
#target_include_directories(pybind11::pybind11 INTERFACE ${PATH_TO_PYBIND11})
#find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
#include_directories(${PYTHON_INCLUDE_DIRS})
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs   REQUIRED)
#find_package(pybind11 REQUIRED)

# get required packages
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
#find_package(PythonInterp 3)
#find_package(PythonLibs   3)
#include_directories(${PYTHON_INCLUDE_DIRS})

# set source files
set(SOURCES ${SOURCES} 
    src/shared.cpp 
    src/window.cpp 
    src/basepair.cpp 
    src/matrix.cpp 
    src/flow.cpp 
    src/approximate.cpp 
    src/scanfoldIO.cpp
)

# set binding file
set(BINDINGS
    src/fold.cpp
)
# create pybind11 module
#pybind11_add_module(fold src/fold.cpp)
add_library(fold MODULE "${SOURCES}" "${BINDINGS}")
target_link_libraries(fold pybind11::module)
#target_link_libraries(fold PRIVATE pybind11::pybind11)
# set output directory
set_target_properties(fold PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib/fold"
)
