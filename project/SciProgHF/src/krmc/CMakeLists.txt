include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/krmc
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libkrmc.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_KRMC_FORTRAN_SOURCES
    krmc_gradient.F
    krmccan.F
    krmcgascip.F
    krmcgp.F
    krmcinp.F
    krmcmlf.F
    krmcopt.F
    krmcout.F
    krmcsav.F
    krmctra.F
    krmcvec.F
)

set(FIXED_KRMC_FORTRAN_SOURCES_noruntimecheck
    krmcluci.F
    krci_ctl.F # Miro: DKRCI_THRGCI not defined
)

set(FREE_KRMC_FORTRAN_SOURCES
   mc_energies.F90
   mcscf_routines.F90
)

set(KRMC_C_SOURCES
   gendet.c
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_KRMC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_KRMC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_KRMC_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'krmc' module with some exceptions")
    set_source_files_properties(${FREE_KRMC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_KRMC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${KRMC_C_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS_runtimecheck})
endif()

add_library(
    krmc
    OBJECT
    ${FIXED_KRMC_FORTRAN_SOURCES}
    ${FIXED_KRMC_FORTRAN_SOURCES_noruntimecheck}
    ${FREE_KRMC_FORTRAN_SOURCES}
    ${KRMC_C_SOURCES}
)

#
# add intermodule dependencies to satisfy compilation and linking
#
add_dependencies(krmc gp)
add_dependencies(krmc dirac)
add_dependencies(krmc main)
add_dependencies(krmc interface_mpi)
if(ENABLE_QCORR)
#   it's module qcorr_cfg
    add_dependencies(krmc qcorr)
endif()

