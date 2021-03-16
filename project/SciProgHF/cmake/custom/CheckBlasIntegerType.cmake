#.rst:
# CheckBlasIntegerType
# --------------------
#
# Check integer setting of the blas library by compiling and running small Fortran program (*.F90).
# If program crashes, then we have integer*8 blas library.
#
#    CHECK_BLAS_INTEGER_TYPE(<var_compile> <var_run> <var_type>)
#
# ::
#
#  <var_compile>    - variable to store the result of the compilation (1 for success, empty for failure)
#                     Will be created as an internal cache variable.
#
#  <var_run>        - variable to store the result of the run (1 for integer*4 blas, 0 for integer*8 blas)
#                     <var_run> is valid only for <var_compile>=1 !
#                     Will be created as an internal cache variable.
#
#  <var_type>       - variable to store the blas integer type: 1 for integer*4 blas, 2 for integer*8 blas 
#                     and 0 for not-tested (or test failed).
#                     Will be created as an internal cache variable.
#
#
#
# Program is calling following macros:
#
#  ::
#
#   CHECK_Fortran_SOURCE_COMPILES1 (our macro, not CMake's)
#   CHECK_Fortran_SOURCE_RUNS1  (our macro, not CMake's)
#
# Program is compiling and running short Fortran source code of Hans Jorgen Aa Jensen, test-blas-integer-compatibility.F90.
#
# The following variables are to be set before calling this macro to modify
# the way the check is run:
#
# ::
#
#   CMAKE_REQUIRED_FLAGS = string of compile command line flags
#   CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#   CMAKE_REQUIRED_INCLUDES = list of include directories
#   CMAKE_REQUIRED_LIBRARIES = list of libraries to link - must contain blas !
#   CMAKE_REQUIRED_QUIET = execute quietly without messages
#

macro(CHECK_BLAS_INTEGER_TYPE VAR_COMPILE VAR_RUN VAR_TYPE)

  include(CheckFortranSourceCompiles1)
  include(CheckFortranSourceRuns1)

  if(NOT CMAKE_REQUIRED_QUIET)
    message(STATUS "Performing compile and run test of blas integer size (can deactivate with -D ENABLE_MATH_INT_TEST=OFF) ")
  endif()

   #  try if Fortran source compiles (Miro's adapted macro)
  file(READ "${PROJECT_SOURCE_DIR}/cmake/custom/test-blas-integer-compatibility.F90" blas_source)
  CHECK_Fortran_SOURCE_COMPILES1(${blas_source} ${VAR_COMPILE} )

  # try if Fortran source runs (Miro's adapted macro)
  CHECK_Fortran_SOURCE_RUNS1(${blas_source} ${VAR_RUN})

  if (${VAR_COMPILE} AND ${VAR_RUN})
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Test of blas integer size SUCCESS .. blas is integer*4 ")
    endif()
    set(${VAR_TYPE} 1)
  elseif (${VAR_COMPILE} AND NOT ${VAR_RUN})
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Test of blas integer size SUCCESS .. blas is integer*8 ")
    endif()
    set(${VAR_TYPE} 2)
  elseif (NOT ${VAR_COMPILE})
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "INFO: Test of blas integer size failed - could not compile or link source code...")
    endif()
    set(${VAR_TYPE} 0)
  endif()

endmacro()
