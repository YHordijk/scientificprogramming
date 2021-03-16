include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libinput.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_INPUT_FORTRAN_SOURCES
)

set(FREE_INPUT_FORTRAN_SOURCES
   input_reader.F90
   input_reader_sections.F90
   parse_input.F90
)
if(ENABLE_EXATENSOR)
set(FREE_INPUT_FORTRAN_SOURCES
   ${FREE_INPUT_FORTRAN_SOURCES}
   ${PROJECT_SOURCE_DIR}/src/exacorr/exacorr_cfg.F90
)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_INPUT_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_INPUT_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'input' module without exceptions")
    set_source_files_properties(${FREE_INPUT_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_INPUT_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    input
    OBJECT
    ${FIXED_INPUT_FORTRAN_SOURCES}
    ${FREE_INPUT_FORTRAN_SOURCES}
)

# deal with intermodule dependencies
# input reader needs cfg-modules of many DIRAC functionalities
add_dependencies(input reladc)
if (ENABLE_QCORR)
    add_dependencies(input qcorr)
endif()
add_dependencies(input gp)
add_dependencies(input dirac)
add_dependencies(input main)
add_dependencies(input visual)
add_dependencies(input relccsd)
if(ENABLE_EXATENSOR)
add_dependencies(input exacorr)
endif()
add_dependencies(input grid)
