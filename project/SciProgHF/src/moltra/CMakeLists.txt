include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/moltra
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libmoltra.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_MOLTRA_FORTRAN_SOURCES
   tradr1.F
   tradr2.F
   tradr4.F
   tradr6.F
   trainp.F
   traone.F
   traout.F
   trapam.F
   traprp.F
   # traone.F depends on spfao
   ${PROJECT_SOURCE_DIR}/src/dirac/dirspf.F
)

set(FREE_MOLTRA_FORTRAN_SOURCES
    dcbtra_h.F90
    tradr3.F90
    tradr5.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_MOLTRA_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_MOLTRA_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_MOLTRA_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'moltra' module")
    set_source_files_properties(${FREE_MOLTRA_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_MOLTRA_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    moltra
    OBJECT
    ${FREE_MOLTRA_FORTRAN_SOURCES}
    ${FIXED_MOLTRA_FORTRAN_SOURCES}
)

# deal with intermodule dependencies
# moltra depends on:
add_dependencies(moltra main)
add_dependencies(moltra dirac)
add_dependencies(moltra dirrci)
if(ENABLE_EXATENSOR)
add_dependencies(moltra exacorr)
endif()
add_dependencies(moltra interface_mpi)
