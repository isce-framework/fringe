install (PROGRAMS   uavsar2vrt.py
                    stripmap2vrt.py
                    tops2vrt.py
                    adjustMiniStacks.py
                    unwrap_fringe.py 
                    unwrapStack.py 
                    Network.py
                    StackVRT.py
                    integratePS.py
        DESTINATION ${FRINGE_BINDIR}
        COMPONENT python)
