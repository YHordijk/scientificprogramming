include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
   )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libxcint.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_XCINT_FORTRAN_SOURCES
    )

set(FREE_XCINT_FORTRAN_SOURCES
   # miro: can not make dependency on the entire gp, source-file pick is necessary
   ${PROJECT_SOURCE_DIR}/src/gp/quaternion_algebra.F90
   density_eval.F90
   dft_cfg.F90
   dirac_ao_eval.F90
   dirac_ao_eval_sub.F90
   extra.F90
   file_units.F90
   interface_ao.F90
   interface_ao_specific.F90
   interface_file_io.F90
   interface_fortran.F90
   interface_functional_read.F90
   interface_grid.F90
   interface_mo.F90
   interface_mo_specific.F90
   overlap_diagnostic.F90
   sigma_prefactor_setting.F90
   xc_ac.F90
   xc_blas3.F90
   xc_derv.F90
   xc_geometric.F90
   xc_london_c1.F90
   xc_max_block_length.F90
   xc_mpi.F90
   xc_response.F90
   xcint_main.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_XCINT_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_XCINT_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'xcint' module without exceptions")
    set_source_files_properties(${FREE_XCINT_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_XCINT_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    xcint
    OBJECT
    ${FIXED_XCINT_FORTRAN_SOURCES}
    ${FREE_XCINT_FORTRAN_SOURCES}
)

# add dependencies 
if(ENABLE_XCFUN)
    add_dependencies(xcint xcfun_fortran_bindings)
    add_dependencies(xcint interface_mpi)
endif()
