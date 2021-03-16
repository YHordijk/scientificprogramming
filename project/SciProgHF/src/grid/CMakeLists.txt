include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/grid
    ${PROJECT_SOURCE_DIR}/src/cfun # yes, depends 
)

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libgrid.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_GRID_FORTRAN_SOURCES
    Lebedev-Laikov.F
    dftgrd.F
    grtest.F
)

set(FREE_GRID_FORTRAN_SOURCES
  num_grid_cfg.F90
  num_grid_gen.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_GRID_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_GRID_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'grid' module without exceptions")
    set_source_files_properties(${FREE_GRID_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_GRID_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    grid
    OBJECT
    ${FIXED_GRID_FORTRAN_SOURCES}
    ${FREE_GRID_FORTRAN_SOURCES}
)

#
# deal with module intedependencies
#
add_dependencies(grid xcint)
add_dependencies(grid gp)
