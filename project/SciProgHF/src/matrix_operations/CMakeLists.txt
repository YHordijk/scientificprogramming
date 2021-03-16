include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libmatrix_operations.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_MATRIX_OPERATIONS_FORTRAN_SOURCES
)

set(FREE_MATRIX_OPERATIONS_FORTRAN_SOURCES
    dirac_interface.F90
    matrix_defop.F90
    matrix_genop.F90
    matrix_module.F90
    matrix_operations.F90
    matrix_stat.F90
    memory_bouncer.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_MATRIX_OPERATIONS_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_MATRIX_OPERATIONS_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'matrix_operations' module without exceptions")
    set_source_files_properties(${FREE_MATRIX_OPERATIONS_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_MATRIX_OPERATIONS_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    matrix_operations
    OBJECT
    ${FIXED_MATRIX_OPERATIONS_FORTRAN_SOURCES}
    ${FREE_MATRIX_OPERATIONS_FORTRAN_SOURCES}
)
#
# deal with intermodule dependencies
#
add_dependencies(matrix_operations grid)
add_dependencies(matrix_operations gp)
