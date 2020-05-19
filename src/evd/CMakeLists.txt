######Start of actual cmake part
set(EVD_SRCS    ../../include/fringe/ulongmask.hpp
                ../../include/fringe/args.hxx
                evd.hpp
                evd.cpp)

set(EVD_INCS    ../../include
                ${ARMADILLO_INCLUDE_DIRS}
                ${GDAL_INCLUDE_DIRS})

set(EVD_LIBS    ${ARMADILLO_LIBRARIES}
                ${GDAL_LIBRARIES}
                ${LAPACK_LIBRARIES}
                m)

###Section for cython extension
set_source_files_properties(
    evdlib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(evdlib evdlib.pyx)
target_include_directories(evdlib PUBLIC ${EVD_INCS} .)
target_link_libraries(evdlib PUBLIC ${EVD_LIBS}) 
target_link_libraries(evdlib PRIVATE OpenMP::OpenMP_CXX)

###Final installation
install (TARGETS
        evdlib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    evd.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
