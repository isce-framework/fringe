#List of srcs
set(TOPORESIDUAL_SRCS  ../../include/fringe/args.hxx
                        toporesidual.hpp
                        toporesidual.cpp)

#List of includes
set(TOPORESIDUAL_INCS  ../../include
                       ../../include/fringe
                        ${GDAL_INCLUDE_DIRS})

#List of libs
set(TOPORESIDUAL_LIBS  ${GDAL_LIBRARIES}
                        gomp m)

##Section for C++ executable
add_executable(toporesidual ${TOPORESIDUAL_SRCS})
set_target_properties(toporesidual PROPERTIES COMPILE_FLAGS "-fopenmp" LINK_FLAGS "-fopenmp")
target_include_directories(toporesidual PUBLIC  ${TOPORESIDUAL_INCS})
target_link_libraries(toporesidual ${TOPORESIDUAL_LIBS})

install (TARGETS 
        toporesidual
        DESTINATION ${fringe_BINDIR}
        COMPONENT exec)


###Section for cython extension
#set_source_files_properties(
#    toporesiduallib.pyx
#    PROPERTIES CYTHON_IS_CXX TRUE)
#cython_add_module(toporesiduallib toporesiduallib.pyx)
#target_include_directories(toporesiduallib PUBLIC ${TOPORESIDUAL_INCS} .)
#set_target_properties(toporesiduallib PROPERTIES COMPILE_FLAGS "-fopenmp")
#target_link_libraries(toporesiduallib ${TOPORESIDUAL_LIBS})


###Final installation
#install (TARGETS
#        toporesiduallib
#        DESTINATION ${fringe_PYDIR}
#        COMPONENT pylib)

#install(PROGRAMS
#    toporesidual.py
#    DESTINATION ${fringe_BINDIR}
#    COMPONENT pyexec)
