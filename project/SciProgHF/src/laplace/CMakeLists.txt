if (ENABLE_LAPLACE)
  include_directories(
      ${PROJECT_SOURCE_DIR}/src/include
      ${PROJECT_BINARY_DIR}/external/include/laplace-minimax
      )

  if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For liblaplace.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
  endif()

  # source files
  set(FREE_F95_LAPLACE_FORTRAN_SOURCES
      pamlaplce.F90
     )
 
   set(FREE_F03_LAPLACE_FORTRAN_SOURCES
      lap_driver.F90
     )

   # compile with 2003 standard
   if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
     set_source_files_properties(${FREE_F03_LAPLACE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-std=f2003")
   endif()
   
   if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
     set_source_files_properties(${FREE_F03_LAPLACE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-stand f03")
   endif()
   
   # check XL
   if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
     set_source_files_properties(${FREE_LAPLACE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
     set_source_files_properties(${FIXED_LAPLACE_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
   endif()
   
  if(ENABLE_RUNTIMECHECK)
     message(STATUS "runtime-check flags activated for the 'LAPLACE' module without exceptions")
     set_source_files_properties(${FREE_LAPLACE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_LAPLACE_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
  endif()

  add_library(
      laplace
      OBJECT
      ${FREE_F95_LAPLACE_FORTRAN_SOURCES}
      ${FREE_F03_LAPLACE_FORTRAN_SOURCES}
      )

  add_dependencies(laplace laplace-minimax)
  add_dependencies(laplace gp)
  add_dependencies(laplace dirac)
  add_dependencies(laplace main)

endif()
