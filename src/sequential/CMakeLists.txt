###Final installation
install (PROGRAMS
        Stack.py
        DESTINATION ${FRINGE_PYDIR}
        COMPONENT pylib)

install(PROGRAMS
    sequential.py
    DESTINATION ${FRINGE_BINDIR}
    COMPONENT pyexec)
