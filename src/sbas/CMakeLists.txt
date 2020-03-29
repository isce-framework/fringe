#List of srcs
set(SBAS_SRCS           sbas.hpp
                        sbas.cpp)

#List of includes
set(SBAS_INCS  ../../include
                ${GDAL_INCLUDE_DIRS}
                ${ARMADILLO_INCLUDE_DIRS})

#List of libs
set(SBAS_LIBS  ${GDAL_LIBRARIES}
               ${ARMADILLO_LIBRARIES}
               m)

###Section for cython extension
set_source_files_properties(
    sbaslib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(sbaslib sbaslib.pyx)
target_include_directories(sbaslib PUBLIC ${SBAS_INCS} .)
target_link_libraries(sbaslib PUBLIC ${SBAS_LIBS})
target_link_libraries(sbaslib PRIVATE OpenMP::OpenMP_CXX)


###Final installation
install (TARGETS
        sbaslib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    sbas.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
