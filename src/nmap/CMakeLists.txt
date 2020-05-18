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

###Section for cython extension
set_source_files_properties(
    nmaplib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(nmaplib nmaplib.pyx) 

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
