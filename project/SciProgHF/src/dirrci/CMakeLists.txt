include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/dirrci
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdirrci.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_DIRRCI_FORTRAN_SOURCES
   cidirr.F
   cigenp.F
   cigosc.F
   cimain.F
)

set(FREE_DIRRCI_FORTRAN_SOURCES
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_DIRRCI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_DIRRCI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'dirrci' module without exceptions")
    set_source_files_properties(${FREE_DIRRCI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_DIRRCI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

# Create the module library
add_library(
    dirrci 
    OBJECT
    ${FIXED_DIRRCI_FORTRAN_SOURCES}
    ${FREE_DIRRCI_FORTRAN_SOURCES}
)

add_dependencies(dirrci gp)
