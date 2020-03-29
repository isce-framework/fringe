###List all the sources to build nmap
###None in this case as we only have one .cpp file

###Main file - with main() function
set( NMAP_MAIN nmap.cpp)

###List of all files that get included
set ( NMAP_INCFILES ../../include/fringe/fringe_common.hpp
                    ../../include/fringe/args.hxx
                    ../../include/fringe/ulongmask.hpp
                    KS2sample.hpp
                    AD2unique.hpp
                    nmap_cuda.h)

###List all includes
set( NMAP_INCS  ../../include
                ${ARMADILLO_INCLUDE_DIRS}
                ${GDAL_INCLUDE_DIRS} )

###List all libs
set( NMAP_LIBS  ${ARMADILLO_LIBRARIES}
                ${GDAL_LIBRARIES}
                m)

###List all CUDA stuff
#set( NMAP_CUDA_INCFILES ../../include/fringe/cuda/cudaUtils.h
#                        ../../include/fringe/cuda/ulongmask.h
#                        KS2sample_cuda.h)

#set( NMAP_CUDA_SRCS nmap_cuda.cu )


##Section for install C++ executable
#if (WITH_CUDA)
#    add_executable(nmap ${NMAP_INCFILES} ${NMAP_CUDA_INCFILES} ${NMAP_CUDA_SRCS} ${NMAP_MAIN})
#    set_source_files_properties(${NMAP_CUDA_SRCS} PROPERTIES COMPILE_FLAGS "-Xcompiler -fopenmp")
#    set_source_files_properties(${NMAP_MAIN}  PROPERTIES COMPILE_FLAGS "-fopenmp -O3")
#    message(STATUS "Building nmap with CUDA support")
#    target_compile_definitions(nmap PRIVATE BUILD_WITH_CUDA)
#else()
    add_executable(nmap ${NMAP_INCFILES}  ${NMAP_MAIN})
    message(STATUS "Building nmap without CUDA support")
#endif()

target_include_directories(nmap PUBLIC  ${NMAP_INCS})
target_link_libraries(nmap PUBLIC ${NMAP_LIBS})
target_link_libraries(nmap PRIVATE OpenMP::OpenMP_CXX)

install (TARGETS 
        nmap
        DESTINATION ${FRINGE_BINDIR}
        COMPONENT exec)


###Section for cython extension
set_source_files_properties(
    nmaplib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)

#if (WITH_CUDA)
#    message(STATUS "Building nmaplib with CUDA support")

    ##First build a library with CUDA
    #    add_library(nmapcudalib SHARED
    #    ${NMAP_CUDA_INCFILES} ${NMAP_CUDA_SRCS})
    #target_include_directories(nmapcudalib PUBLIC ${NMAP_INCS})

    ###Then build cython module using library
    #cython_add_module(nmaplib nmaplib.pyx)
    #target_compile_definitions(nmaplib PRIVATE BUILD_WITH_CUDA)
#else()
    cython_add_module(nmaplib nmaplib.pyx) 
#endif()

target_include_directories(nmaplib PUBLIC ${NMAP_INCS} . )
target_link_libraries(nmaplib PUBLIC ${NMAP_LIBS}) 
target_link_libraries(nmaplib PRIVATE OpenMP::OpenMP_CXX)

###Final installation
install (TARGETS
        nmaplib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    nmap.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)

if (WITH_CUDA)
    install (TARGETS
            nmapcudalib
            DESTINATION ${FRINGE_LIBDIR}
            COMPONENT pycudalib)
endif()
