include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

add_definitions(-DBUILD_PELIB)

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libpelib_interface.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_PELIB_FORTRAN_SOURCES
    pe_dirac_interfaces.F
)
set(FREE_PELIB_FORTRAN_SOURCES
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_PELIB_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_PELIB_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'pelib_interface' module without exceptions")
    set_source_files_properties(${FREE_PELIB_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_PELIB_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
          pelib_interface
          ${FIXED_PELIB_FORTRAN_SOURCES}
          ${FREE_PELIB_FORTRAN_SOURCES}
)

# deal with intedependencies
# needs memory allocator etc..
add_dependencies(gen1int_interface gen1int)
add_dependencies(pelib gen1int_interface)
add_dependencies(pelib_interface pelib)
add_dependencies(pelib_interface interface_mpi)
