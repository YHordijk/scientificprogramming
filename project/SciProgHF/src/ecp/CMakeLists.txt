include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/ecp
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libecp.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_ECP_FORTRAN_SOURCES
    recp1.F
    recp_lnk.F
)

set(FREE_ECP_FORTRAN_SOURCES
   recp_chk.F90
   recp_cso.F90
   recp_cso_ang.F90
   recp_cso_ps.F90
   recp_cso_rad.F90
   recp_cso_sub.F90
   recp_cso_table.F90
   recp_function1.F90
   recp_function2.F90
   recp_inp.F90
   recp_inp_orb.F90
   recp_inp_print.F90
   recp_inp_read.F90
   recp_ipt.F90
   recp_ntr.F90
   recp_one.F90
   recp_one_stvcz.F90
   recp_one_sub.F90
   recp_output.F90
   recp_socfpd.F90
   recp_write.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_ECP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_ECP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'ecp' module without exceptions")
    set_source_files_properties(${FREE_ECP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_ECP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    ecp
    OBJECT
    ${FIXED_ECP_FORTRAN_SOURCES}
    ${FREE_ECP_FORTRAN_SOURCES}
)
