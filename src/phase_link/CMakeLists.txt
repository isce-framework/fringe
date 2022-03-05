######Start of actual cmake part
set(PHASELINK_SRCS    ../../include/fringe/ulongmask.hpp
                ../../include/fringe/args.hxx
                phase_link.hpp
                phase_link.cpp)

	set(PHASELINK_INCS    ../../include
                ${ARMADILLO_INCLUDE_DIRS}
                ${GDAL_INCLUDE_DIRS})

	set(PHASELINK_LIBS    ${ARMADILLO_LIBRARIES}
                ${GDAL_LIBRARIES}
                ${LAPACK_LIBRARIES}
                m)

###Section for cython extension
set_source_files_properties(
    phase_linklib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(phase_linklib phase_linklib.pyx)
target_include_directories(phase_linklib PUBLIC ${PHASELINK_INCS} .)
target_link_libraries(phase_linklib PUBLIC ${PHASELINK_LIBS}) 
target_link_libraries(phase_linklib PRIVATE OpenMP::OpenMP_CXX)

###Final installation
install (TARGETS
        phase_linklib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    phase_link.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
