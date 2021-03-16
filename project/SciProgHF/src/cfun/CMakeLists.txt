include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libcfun.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_CFUN_FORTRAN_SOURCES
    dftsetcam.F
)

set(FREE_CFUN_FORTRAN_SOURCES
    fortran_io.F90
)

set(CFUN_C_SOURCES
   emb-functional.c
   fun-becke.c
   fun-cam.c
   fun-gga.c
   fun-gllb_hole.c
   fun-hjsx.c
   fun-kin_pw91k.c
   fun-kin_tf.c
   fun-kin_vw.c
   fun-kt.c
   fun-lb94.c
   fun-lbalpha.c
   fun-lyp.c
   fun-optx.c
   fun-p86c.c
   fun-pbec.c
   fun-pbex.c
   fun-pw86x.c
   fun-pw91c.c
   fun-pw91x.c
   fun-pz81.c
   fun-revpbex.c
   fun-rpbex.c
   fun-slater.c
   fun-srpbec.c
   fun-srvwn5.c
   fun-vwn.c
   functionals.c
   general.c
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_CFUN_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_CFUN_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'cfun' module")
    set_source_files_properties(${FREE_CFUN_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_CFUN_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${CFUN_C_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS_runtimecheck})
endif()

add_library(
    cfun
    OBJECT
    ${FIXED_CFUN_FORTRAN_SOURCES}
    ${FREE_CFUN_FORTRAN_SOURCES}
    ${CFUN_C_SOURCES}
)
