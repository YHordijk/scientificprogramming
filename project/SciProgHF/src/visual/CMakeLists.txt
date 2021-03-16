include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libvisual.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_VISUAL_FORTRAN_SOURCES
)

set(FREE_VISUAL_FORTRAN_SOURCES
    integration.F90
    visual_cfg.F90
    visual.F90
    visual_in_point.F90
    visual_london.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_VISUAL_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_VISUAL_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'visual' module without exceptions")
    set_source_files_properties(${FREE_VISUAL_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_VISUAL_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    visual 
    OBJECT
    ${FIXED_VISUAL_FORTRAN_SOURCES}
    ${FREE_VISUAL_FORTRAN_SOURCES}
)

# add intermodule dependencies
add_dependencies(visual gp)
add_dependencies(visual matrix_operations)
