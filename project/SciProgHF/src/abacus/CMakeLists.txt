include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/abacus
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For librabacus.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_ABACUS_FORTRAN_SOURCES
    abadip.F
    abalnr.F
    abamolden.F
    abanuc.F
    abaop2.F
    abaopt.F
    abaortvec.F
    abarint.F
    abatro.F
    abavib.F
    abavrml.F
    abawalk.F
    her1ave.F
    her1car.F
    her1cat.F
    her1drv.F
    her1int.F
    her1odc.F
    her1out.F
    her1pro.F
    her1sym.F
    her2ave.F
    her2dir.F
    her2drv.F
    her2el1.F
    her2el2.F
    her2fck.F
    her2gab.F
    her2her.F
    her2odc.F
    her2out.F
    her2sym.F
    herbas.F
    hercen.F
    herdirect.F
    herdrv.F
    hergam.F
    hergeo.F
    hergp.F
    hergroup.F
    hermem.F
    herrdn.F
    herrdn_dirac.F
    hersol.F
    hersym.F
    hertst.F
    vibnor.F
    herpar.F 
)

set(FREE_ABACUS_FORTRAN_SOURCES
   carsph.F90
   hermodules.F90
   filemodule.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_ABACUS_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_ABACUS_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()


if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the whole 'abacus' module without exceptions")
    set_source_files_properties(${FREE_ABACUS_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_ABACUS_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

# fix for Windows with 64bit integers on and build type release
# when using GNU (gfortran) compilers
# there is one source file for which -O3/-O2 is not good (I.Hrasko)
if(     (CMAKE_SYSTEM_NAME MATCHES "Windows") 
    AND (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    AND (ENABLE_64BIT_INTEGERS)
    AND (CMAKE_BUILD_TYPE MATCHES "release"))
    set_source_files_properties(
        ${PROJECT_SOURCE_DIR}/src/abacus/hermodules.F90
        PROPERTIES COMPILE_FLAGS "-O1" # added compiler flag overrides previous optimization flags
        )
    message(STATUS "src/abacus/hermodules.F90: assigned -O1 optimization flag to fix the compilation")
endif()

add_library(
    abacus
    OBJECT
    ${FIXED_ABACUS_FORTRAN_SOURCES}
    ${FREE_ABACUS_FORTRAN_SOURCES}
)

#
# deal with intermodule dependencies 
# 
add_dependencies(abacus main)
add_dependencies(abacus gp)
add_dependencies(abacus ecp)
add_dependencies(abacus xcint)  
add_dependencies(abacus grid)  
if (ENABLE_PELIB)
    add_dependencies(abacus pelib)  
endif()
