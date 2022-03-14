cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project( FRINGE
         VERSION 0.1
         LANGUAGES CXX)

####!!!!-- Checking prerequisites --!!!!

###Add modules
list(APPEND CMAKE_MODULE_PATH ${FRINGE_SOURCE_DIR}/cmake)

#Set C++ Standard to use
set(CMAKE_CXX_STANDARD 11)

###Find gdal
FIND_PACKAGE(GDAL 3.0 REQUIRED)

###Find Armadillo
FIND_PACKAGE(Armadillo 8.100 REQUIRED)

###Find openmp
find_package(OpenMP COMPONENTS CXX)

###Find linalg
FIND_PACKAGE(LAPACK REQUIRED)

###Find python3
FIND_PACKAGE(Python 3.7 COMPONENTS Interpreter Development NumPy)
include(useCython)

find_package(pybind11 REQUIRED)

###Include FRINGE specific paths
include (ConfigProject)

####!!!!-- End of checking prerequisites --!!!!


####!!!!--Start installation--!!!!
add_subdirectory (src)
add_subdirectory (python)
