include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/lucita
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For liblucita.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_LUCITA_FORTRAN_SOURCES
    doubledir.F
    fock.F
    io.F
    lucia_cc.F
    lucia_icci.F
    property.F
    transform.F
)

#miro: specify source files where one can not apply runtime checks due to old-fashion non-rigorous programming style
# this is due to WORK(:) array variable which can not undergo runtime checks
set(FIXED_LUCITA_FORTRAN_SOURCES_noruntimecheck
    vecmat.F
    lucita_diagonalize.F
    density.F
    sigma.F
    diagonal.F
    strings.F
    ciinfo.F
    interface.F
    memory.F
    program.F
    utils.F
    par_lucia.F #Miro: Subscript #1 of the array WORK has value 5905696052053 which is greater than the upper bound of 1
    par_solver.F
)

set(FREE_LUCITA_FORTRAN_SOURCES
)

set(LUCITA_C_SOURCES
 flush_stdout.c
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_LUCITA_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_LUCITA_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_LUCITA_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

#
if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'lucita' module with many exceptions")
    set_source_files_properties(${FREE_LUCITA_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_LUCITA_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${LUCITA_C_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS_runtimecheck})
endif()

add_library(
    lucita
    OBJECT
    ${FIXED_LUCITA_FORTRAN_SOURCES}
    ${FIXED_LUCITA_FORTRAN_SOURCES_noruntimecheck}
    ${FREE_LUCITA_FORTRAN_SOURCES}
    ${LUCITA_C_SOURCES}
)

add_dependencies(lucita qcorr)
add_dependencies(luciarel interface_mpi)
