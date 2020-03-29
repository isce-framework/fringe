#List of srcs
set(UNWRAP_SRCS           Unwrap.h
                       Unwrap.cc)

#List of includes
set(UNWRAP_INCS  ../../include
                ${GDAL_INCLUDE_DIRS}
                ${ARMADILLO_INCLUDE_DIRS})

#List of libs
set(UNWRAP_LIBS  ${GDAL_LIBRARIES}
               ${ARMADILLO_LIBRARIES}
               m )

###Section for cython extension
set(library_sources SHARED Unwrap.cc)
add_library(unwrap SHARED ${library_sources})
#set_target_properties(unwrap PROPERTIES COMPILE_FLAGS "-DARMA_EXTRA_DEBUG")
target_link_libraries(unwrap PUBLIC ${UNWRAP_LIBS})
target_link_libraries(unwrap PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(unwrap PUBLIC ${GDAL_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

###Final installation
install (TARGETS unwrap
        DESTINATION ${FRINGE_LIBDIR}
        COMPONENT lib)


###Section for cython extension
set_source_files_properties(
    Unwrap.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(unwraplib unwraplib.pyx)
target_include_directories(unwraplib PUBLIC ${UNWRAP_INCS} .)
target_link_libraries(unwraplib PUBLIC unwrap)


###Final installation
install (TARGETS
        unwraplib
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

#install(PROGRAMS
#    unwrap_error.py
#    DESTINATION ${FRINGE_BINDIR}
#    COMPONENT pyexec)


