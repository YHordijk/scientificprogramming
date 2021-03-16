include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libinterest.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_INTEREST_FORTRAN_SOURCES
)

set(FREE_INTEREST_FORTRAN_SOURCES
    module_interest_interface.F90
    src/git_revision_info.f90
    src/module_interest_eri.f90
    src/module_interest_hrr.f90
    src/module_interest_one.f90
    src/module_interest_osr.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_INTEREST_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_INTEREST_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'interest' module without exceptions")
    set_source_files_properties(${FREE_INTEREST_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_INTEREST_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    interest
    OBJECT
    ${FREE_INTEREST_FORTRAN_SOURCES}
)
