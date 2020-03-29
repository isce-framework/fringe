#List of srcs
set(FIT_SRCS    Fit.h
                Fit.cc)

#List of includes
set(FIT_INCS  ../../include
                ${GDAL_INCLUDE_DIRS}
                ${ARMADILLO_INCLUDE_DIRS})

#List of libs
set(FIT_LIBS  ${GDAL_LIBRARIES}
              ${ARMADILLO_LIBRARIES}
              m )

###Section for cython extension
set(library_sources SHARED Fit.cc)
add_library(fit SHARED ${library_sources})
#set_target_properties(fit PROPERTIES COMPILE_FLAGS "-DARMA_EXTRA_DEBUG")
target_link_libraries(fit ${FIT_LIBS})
target_include_directories(fit PUBLIC ${GDAL_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

###Final installation
install (TARGETS fit
        DESTINATION ${FRINGE_LIBDIR}
        COMPONENT lib)


###Section for cython extension
set_source_files_properties(
    fitlib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(fitlib fitlib.pyx)
target_include_directories(fitlib PUBLIC ${FIT_INCS} .)
target_link_libraries(fitlib PUBLIC fit)


###Final installation
install (TARGETS
        fitlib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install (PROGRAMS
        TimeFunction.py 
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    fit_timeseries.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)


