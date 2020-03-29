###List of source files
set(CALAMP_SRCS ../../include/fringe/args.hxx
                calamp.hpp
                calamp.cpp )

###List of includes
set(CALAMP_INCS ../../include
                ${GDAL_INCLUDE_DIRS})

###List of libraries
set(CALAMP_LIBS ${GDAL_LIBRARIES})

##Section for C++ executable
add_executable(calamp ${CALAMP_SRCS})
target_include_directories(calamp PUBLIC  ${CALAMP_INCS})
target_link_libraries(calamp PUBLIC ${CALAMP_LIBS})
target_link_libraries(calamp PRIVATE OpenMP::OpenMP_CXX)

install (TARGETS 
        calamp
        DESTINATION ${FRINGE_BINDIR}
        COMPONENT exec)


###Section for cython extension
set_source_files_properties(
    calamplib.pyx
    PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(calamplib calamplib.pyx)
target_include_directories(calamplib PUBLIC ${CALAMP_INCS} .)
target_link_libraries(calamplib PUBLIC ${CALAMP_LIBS})
target_link_libraries(calamplib PRIVATE OpenMP::OpenMP_CXX)

##Final installation
install (TARGETS
         calamplib
         DESTINATION ${FRINGE_PYDIR}
         COMPONENT pylib)

install (PROGRAMS 
         calamp.py
         DESTINATION ${FRINGE_BINDIR}
         COMPONENT pyexec)
