option(ENABLE_LAPACK_QDIAG "Enable LAPACK diagonalization of real matrices" OFF)
if(ENABLE_LAPACK_QDIAG)
    add_definitions(-DLAPACK_QDIAG)
endif()
