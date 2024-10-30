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

# project files
add_executable(fold lib/fold.cpp lib/fold.hpp)
target_link_libraries(fold ${PYTHON_LIBRARIES})

