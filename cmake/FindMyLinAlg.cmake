# Find the Lapack and Blas libraries
#
# The following variables are optionally searched for defaults
#  MYLINALG_ROOT_DIR:            Base directory where all components are found
#
# The following are set after configuration is done:
#  MYLINALG_FOUND
#  MYLINALG_INCLUDE_DIRS
#  MYLINALG_LIBRARIES

if (APPLE)
    set (MYLINALG_INCLUDE_SEARCH_PATHS
        /opt/local/include
        /opt/local/atlas/include
        $ENV{MYLINALG_ROOT_DIR}
        $ENV{MYLINALG_ROOT_DIR}/include
        )

    set (MYLINALG_LIB_SEARCH_PATHS
        /opt/local/lib
        /opt/local/lib/atlas
        $ENV{MYLINALG_ROOT_DIR}
        $ENV{MYLINALG_ROOT_DIR}/lib
        )
else()
    set(MYLINALG_INCLUDE_SEARCH_PATHS
        /usr/include/atlas
        /usr/include/atlas-base
        $ENV{MYLINALG_ROOT_DIR}
        $ENV{MYLINALG_ROOT_DIR}/include
        )

    set (MYLINALG_LIB_SEARCH_PATHS
        /usr/lib64
        /usr/lib
        /usr/lib/atlas
        /usr/lib/atlas-base
        $ENV{MYLINALG_ROOT_DIR}
        $ENV{MYLINALG_ROOT_DIR}/lib
        )
        
endif()

message (STATUS "Searching: ${MYLINALG_LIB_SEARCH_PATHS}")

find_path(MYLINALG_CBLAS_INCLUDE_DIR   NAMES cblas.h   PATHS ${MYLINALG_INCLUDE_SEARCH_PATHS})
find_path(MYLINALG_CLAPACK_INCLUDE_DIR NAMES clapack.h PATHS ${MYLINALG_INCLUDE_SEARCH_PATHS})

find_library(MYLINALG_CBLAS_LIBRARY NAMES  cblas_r.a cblas.a cblas_r cblas PATHS ${MYLINALG_LIB_SEARCH_PATHS})
find_library(MYLINALG_PTCBLAS_LIBRARY NAMES ptcblas_r.a ptcblas.a ptcblas_r ptcblas PATHS ${MYLINALG_LIB_SEARCH_PATH})

find_library(MYLINALG_BLAS_LIBRARY NAMES   f77blas_r.a f77blas.a f77blas_r f77blas  PATHS ${MYLINALG_LIB_SEARCH_PATHS})
find_library(MYLINALG_PTBLAS_LIBRARY NAMES ptf77blas_r.a ptf77blas.a ptf77blas_r ptf77blas PATHS ${MYLINALG_LIB_SEARCH_PATHS})

find_library(MYLINALG_LAPACK_LIBRARY NAMES lapack.a alapack_r.a alapack.a lapack_atlas.a atllapack.a lapack alapack_r alapack lapack_atlas atllapack PATHS ${MYLINALG_LIB_SEARCH_PATHS})

find_library(MYLINALG_ATLAS_LIBRARY NAMES atlas_r.a atlas.a atlas_r atlas PATHS ${MYLINALG_LIB_SEARCH_PATHS})

set(LOOKED_FOR
  MYLINALG_CBLAS_INCLUDE_DIR
  MYLINALG_CLAPACK_INCLUDE_DIR

  MYLINALG_CBLAS_LIBRARY
  MYLINALG_PTCBLAS_LIBRARY
  MYLINALG_BLAS_LIBRARY
  MYLINALG_PTBLAS_LIBRARY
  MYLINALG_LAPACK_LIBRARY
  MYLINALG_ATLAS_LIBRARY
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MYLINALG DEFAULT_MSG ${LOOKED_FOR})

if(MYLINALG_FOUND)
  set(MYLINALG_INCLUDE_DIR ${MYLINALG_CBLAS_INCLUDE_DIR} ${MYLINALG_CLAPACK_INCLUDE_DIR})
  set(MYLINALG_LIBRARIES ${MYLINALG_LAPACK_LIBRARY} ${MYLINALG_CBLAS_LIBRARY} ${MYLINALG_BLAS_LIBRARY} ${MYLINALG_ATLAS_LIBRARY})
  mark_as_advanced(${LOOKED_FOR})

  message(STATUS "Found MYLINALG (include: ${MYLINALG_CBLAS_INCLUDE_DIR} library: ${MYLINALG_BLAS_LIBRARY} lapack: ${MYLINALG_LAPACK_LIBRARY}")
endif(MYLINALG_FOUND)

