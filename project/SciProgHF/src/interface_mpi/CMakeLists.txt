include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libinterface_mpi.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_INTERFACE_MPI_FORTRAN_SOURCES
    interface_to_mpi_f77.F
)

set(FREE_INTERFACE_MPI_FORTRAN_SOURCES
    integer_kind_mpilib.F90
    integer_model.F90
    interface_to_mpi.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_INTERFACE_MPI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_INTERFACE_MPI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'interface_mpi' module without exceptions")
    set_source_files_properties(${FREE_INTERFACE_MPI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_INTERFACE_MPI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()


add_library(
    interface_mpi
    OBJECT
    ${FIXED_INTERFACE_MPI_FORTRAN_SOURCES}
    ${FREE_INTERFACE_MPI_FORTRAN_SOURCES}
)
