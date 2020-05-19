###CMake for building FRInGE


###Make sure that symlinks are not being used.
###From GMT's cmake related files
get_filename_component (srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (bindir "${CMAKE_BINARY_DIR}" REALPATH)

###Minimum cmake version
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)


###Project
project( FRINGE
         VERSION 0.1
         LANGUAGES C CXX)

####!!!!-- Checking prerequisites --!!!!

###Add modules
list(APPEND CMAKE_MODULE_PATH ${FRINGE_SOURCE_DIR}/cmake)

#Set C++ Standard to use
set(CMAKE_CXX_STANDARD 11)

###Find gdal
FIND_PACKAGE(GDAL 3.0 REQUIRED)

###Find Armadillo
FIND_PACKAGE(Armadillo REQUIRED)

###Find openmp
find_package(OpenMP COMPONENTS CXX)

###Find linalg
FIND_PACKAGE(LAPACK REQUIRED)

###Find python3
FIND_PACKAGE(Python 3.7 COMPONENTS Interpreter Development NumPy)
include(useCython)

###Include FRINGE specific paths
include (ConfigProject)

####!!!!-- End of checking prerequisites --!!!!


####!!!!--Start installation--!!!!
add_subdirectory (src)
add_subdirectory (python)
