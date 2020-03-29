###List of srcs
set(DESPECK_SRCS    ../../include/fringe/ulongmask.hpp
                    ../../include/fringe/args.hxx
                    despeck.hpp
                    despeck.cpp )

###List of includes
set(DESPECK_INCS    ../../include
                    ${GDAL_INCLUDE_DIRS}
                    ${ARMADILLO_INCLUDE_DIRS} )

###List of libs
set(DESPECK_LIBS    ${ARMADILLO_LIBRARIES}
                    ${GDAL_LIBRARIES}
                    m)

###Section for C++ executable
add_executable(despeck ${DESPECK_SRCS})
target_include_directories(despeck PUBLIC ${DESPECK_INCS})
target_link_libraries(despeck PUBLIC ${DESPECK_LIBS})
target_link_libraries(despeck PRIVATE OpenMP::OpenMP_CXX)


install (TARGETS 
        despeck
        DESTINATION ${FRINGE_BINDIR}
        COMPONENT exec)

###Section for installing cython extension
set_source_files_properties(
    despecklib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(despecklib despecklib.pyx)
target_include_directories(despecklib PUBLIC ${DESPECK_INCS} .)
target_link_libraries(despecklib PUBLIC ${DESPECK_LIBS})
target_link_libraries(despecklib PRIVATE OpenMP::OpenMP_CXX)
install (TARGETS
        despecklib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    despeck.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
