#
# Check the of integer size compatibility of math libraries (BLAS).
#
# These DIRAC CMake internal variables are utilized inside:
#
#  MATH_LIBS 
#  EXPLICIT_LIBS 
#  ENABLE_BUILTIN_BLAS
#  ENABLE_64BIT_INTEGERS 
#  MKL_FLAG
#  CMAKE_Fortran_FLAGS
#  CMAKE_Fortran_COMPILER
#
# This (own) macro is called: CHECK_BLAS_INTEGER_TYPE
#
# Hint: to check it manually, compile the short "test-blas-integer-compatibility.F90" source code
# and link it against desired BLAS library, for example:
#
#        ifort -i8 test-blas-integer-compatibility.F90 -lblas
#        ifort -mkl=parallel test-blas-integer-compatibility.F90
#        ifort -i8 -mkl=sequential test-blas-integer-compatibility.F90
#        gfortran test-blas-integer-compatibility.F90 -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
#        ifort  test-blas-integer-compatibility.F90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
#
# The executable, a.out,  always crashes when is linked against integer*8 (or ilp64) blas library
#

option(ENABLE_MATH_INT_TEST "Enable checking of math libs integer size compatibility" ON)

if (ENABLE_MATH_INT_TEST)

#
# Miro: we have multiple combinations of blas-library linking flags in MATH_LIBS,EXPLICIT_LIBS,MKL_FLAG
#
# The blas-integer-size test is always deactivated for ENABLE_BUILTIN_BLAS=ON, because on takes builtin blas source FIRST.
#
#

  if ( (MATH_LIBS OR EXPLICIT_LIBS OR MKL_FLAG) AND NOT ENABLE_BUILTIN_BLAS)

    message(STATUS "Configure-time Linked math libs integer compatibility test ENABLED")

    include(CheckBlasIntegerType) # own macro

    set(CMAKE_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})

    if (MKL_FLAG)
      # MKL_FLAG belong to CMAKE_REQUIRED_FLAGS !
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=${MKL_FLAG}")
    else()
      if (MATH_LIBS AND EXPLICIT_LIBS)
        #Miro: both MATH_LIBS and EXPLICIT_LIBS are active !
        # FYI: you can not combine, for example MKL and SYSTEM_NATIVE libraries
        set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBS};${EXPLICIT_LIBS}")
      elseif(MATH_LIBS AND NOT EXPLICIT_LIBS)
        set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBS}")
      elseif(NOT MATH_LIBS AND EXPLICIT_LIBS)
        set(CMAKE_REQUIRED_LIBRARIES "${EXPLICIT_LIBS}")
      else()
        message(STATUS "Have you any blas library linked to Dirac ? The linking step may crash ...")
      endif()
    endif()

#   perform BLAS-integer test only, LAPACK-integer test is done inside executable, dirac.x
    CHECK_BLAS_INTEGER_TYPE(blas_test_compiled blas_test_run blas_integer_type)

    if ("${blas_integer_type}" EQUAL 1)
      message(STATUS "We have BLAS-integer*4")
      if (ENABLE_64BIT_INTEGERS)
       # Miro: message does not show up with STATUS before printing WARNING !
       # all previous STATUS output are swallowed ! see issue #188
        message("Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
        message(WARNING "Incompatibitility between blas integer*4 and required integer*8 variables")
      endif()
    elseif ("${blas_integer_type}" EQUAL 2)
      message(STATUS "We have BLAS-integer*8")
      if (NOT ENABLE_64BIT_INTEGERS)
       # Miro: message does not show up with STATUS before printing WARNING !
       # all previous STATUS outputs are swallowed ! issue #188
        message("Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
        message(WARNING "Incompatibitility between blas integer*8 and the default integer*4 variables")
      endif()
    else()
      message(STATUS "Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "Blas integer*4/8 status not determined.")
      message(STATUS "WATCH THAT YOU HAVE COMPABILE INTEGER SIZES IN LINKED BLAS,LAPACK LIBRARIES!")
    endif()
  endif()
else()
  message(STATUS "Linked math libs integer size compatibility test DISABLED")
  message(STATUS "WATCH THAT YOU HAVE COMPABILE INTEGER SIZES IN LINKED BLAS,LAPACK LIBRARIES!")
endif()
