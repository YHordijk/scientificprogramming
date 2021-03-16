if(ENABLE_EXATENSOR)
  include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

  set(FIXED_EXACORR_FORTRAN_SOURCES
    exacorr_dirac.F
  )

  set(FREE_EXACORR_FORTRAN_SOURCES
   exacorr_main.F90
   exacorr_datatypes.F90
   exacorr_global.F90
   exacorr_cfg.F90
   exacorr_utils.F90
   exacorr_tensor_methods.F90
   intermediates.F90
   exacorr_ao_to_mo.F90
   exacorr_gradient.F90
   exacorr_mo.F90
   exacorr_respect.F90
   exacorr_cholesky.F90
   exacorr_cc.F90
   talsh_ao_to_mo.F90
   talsh_common_routines.F90
   talsh_gradient.F90
   talsh_cc.F90
   mobasis_hartreefock.F90
  )

  if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_EXACORR_FORTRAN_SOURCES}  EXACORR COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_EXACORR_FORTRAN_SOURCES} EXACORR COMPILE_FLAGS "-qfixed")
  endif()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "20")
    if (CMAKE_BUILD_TYPE MATCHES "release") 
      set_source_files_properties(${PROJECT_SOURCE_DIR}/src/exacorr/exacorr_cc.F90 PROPERTY COMPILE_FLAGS "-O0")
      message(STATUS "src/exacorr/exacorr_cc.F90: disabling inteprocedural optimization for Intel compiler")
    endif()
  endif()


  add_library(
    exacorr
    OBJECT
    ${FIXED_EXACORR_FORTRAN_SOURCES}
    ${FREE_EXACORR_FORTRAN_SOURCES}
  )

# deal with dependencies
  add_dependencies(exacorr  main)
  add_dependencies(exacorr  exatensor)
endif()
