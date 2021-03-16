include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/dirac
)

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdirac.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_DIRAC_FORTRAN_SOURCES
   atom2molecule_utils.F
   dftaux.F
   dirana.F
   diraos.F
   dirbss.F
   dircnv.F
   dirden.F
   dirgeo.F
   dirgp.F
   dirgrd.F
   dirmp2no.F
   dirone.F
   diropen.F
   dirout.F
   dirrdn.F
   dirscf.F
   dirset.F
   dirsol.F
   dirspf.F
   dirtra.F
   dirtwo.F
   dirv1c.F
   huckel.F
)

set(FREE_DIRAC_FORTRAN_SOURCES
   dirimpmos.F90
#   dirmp2no.F needs 
    ${PROJECT_SOURCE_DIR}/src/relccsd/relcc_cfg.F90
# needs also input reader..
    ${PROJECT_SOURCE_DIR}/src/input/parse_input.F90
    ${PROJECT_SOURCE_DIR}/src/input/input_reader.F90
    ${PROJECT_SOURCE_DIR}/src/input/input_reader_sections.F90
)


if(ENABLE_INTEREST)
    set(FREE_DIRAC_FORTRAN_SOURCES ${PROJECT_SOURCE_DIR}/src/interest/module_interest_interface.F90 ${FREE_DIRAC_FORTRAN_SOURCES})
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()


if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'dirac' module without exceptions")
    set_source_files_properties(${FREE_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    dirac
    OBJECT
    ${FIXED_DIRAC_FORTRAN_SOURCES}
    ${FREE_DIRAC_FORTRAN_SOURCES}
)

#
# Oh, intermodule dependendencies, worth to reduce if possible
#
add_dependencies(dirac main)
add_dependencies(dirac abacus)
add_dependencies(dirac gp)
add_dependencies(dirac xcfun_fortran_bindings)
add_dependencies(dirac embedding)
add_dependencies(dirac localization)
add_dependencies(dirac visual)
add_dependencies(dirac x2c)
add_dependencies(dirac interface_mpi)
if (ENABLE_INTEREST)
    add_dependencies(dirac interest)
endif()
if (ENABLE_PCMSOLVER)
    add_dependencies(dirac pcm_dirac)
endif()
if (ENABLE_AOOSOC)
    add_dependencies(dirac aoosoc)
endif()
if (ENABLE_SRDFT)
    add_dependencies(dirac srdft)
endif()
add_dependencies(dirac reladc)
if (ENABLE_QCORR)
    add_dependencies(dirac qcorr)
endif()
if (ENABLE_KRCC)
    add_dependencies(dirac krcc)
endif()
if(ENABLE_EXATENSOR)
    add_dependencies(dirac exatensor)
endif()


