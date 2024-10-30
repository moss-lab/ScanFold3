cmake_minimum_required(VERSION 3.15...3.30)

# project version
project(fold VERSION 3.0
	     LANGUAGES CXX
       )

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set output directory
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin )

# packages 
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

#find_package(pybind11 CONFIG REQUIRED)
#message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

# project files
add_executable(fold lib/fold.cpp lib/fold.hpp)
#add_library(fold MODULE lib/fold.cpp lib/fold.hpp)
#PYTHON_ADD_MODULE(fold lib/fold.cpp lib/fold.h)
target_link_libraries(fold ${PYTHON_LIBRARIES})
#target_link_libraries(fold pybind11::module)

