###Final installation
file(GLOB files "*.py")
install (PROGRAMS
        ${files} 
        DESTINATION ${FRINGE_PYDIR}/timefn
        COMPONENT pylib)
