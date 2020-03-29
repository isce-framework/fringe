#List of srcs
set(AMPDISPERSION_SRCS  ../../include/fringe/args.hxx
                        ampdispersion.hpp
                        ampdispersion.cpp)

#List of includes
set(AMPDISPERSION_INCS  ../../include
                        ${GDAL_INCLUDE_DIRS})

#List of libs
set(AMPDISPERSION_LIBS  ${GDAL_LIBRARIES}
                        m)

##Section for C++ executable
add_executable(ampdispersion ${AMPDISPERSION_SRCS})
target_include_directories(ampdispersion PUBLIC  ${AMPDISPERSION_INCS})
target_link_libraries(ampdispersion PUBLIC ${AMPDISPERSION_LIBS})
target_link_libraries(ampdispersion PRIVATE OpenMP::OpenMP_CXX)

install (TARGETS 
        ampdispersion
        DESTINATION ${FRINGE_BINDIR}
        COMPONENT exec)


###Section for cython extension
set_source_files_properties(
    ampdispersionlib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(ampdispersionlib ampdispersionlib.pyx)
target_include_directories(ampdispersionlib PUBLIC ${AMPDISPERSION_INCS} .)
target_link_libraries(ampdispersionlib PUBLIC ${AMPDISPERSION_LIBS})
target_link_libraries(ampdispersionlib PRIVATE OpenMP::OpenMP_CXX)


###Final installation
install (TARGETS
        ampdispersionlib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    ampdispersion.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
