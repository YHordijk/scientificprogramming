
# Description
#     We use a Fortran test program
#     to find out whether we can compile code with ieee_arithmeticand run with -xHost.

#-------------------------------------------------------------------------------

if(CMAKE_Fortran_COMPILER_WORKS)

    set(_ieee_binary ${CMAKE_BINARY_DIR}/test_ieee_isnan_support.x)
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_SOURCE_DIR}/cmake/compilers/test_ieee_isnan_support.f90 -o ${_ieee_binary} ERROR_QUIET)
    execute_process(
        COMMAND ${_ieee_binary}
        OUTPUT_VARIABLE _ieee_binary_output
        ERROR_VARIABLE _ieee_binary_output
    )

    set(_intrinsic_binary ${CMAKE_BINARY_DIR}/test_isnan_support.x)
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_SOURCE_DIR}/cmake/compilers/test_isnan_support.f90 -o ${_intrinsic_binary} ERROR_QUIET)
    execute_process(
        COMMAND ${_intrinsic_binary}
        OUTPUT_VARIABLE _intrinsic_binary_output
        ERROR_VARIABLE _intrinsic_binary_output
    )

    # depending on result set XHOST_FLAG_AVAILABLE
    if(_ieee_binary_output MATCHES "ieee_arithmetic supported by compiler")
        message(STATUS "-- IEEE Exception Handling module supported by Fotran compiler")
        add_definitions(-DHAVE_IEEE_ISNAN)
    elseif(_intrinsic_binary_output  MATCHES "isnan intrinsic supported by compiler")
        message(STATUS "-- IEEE Exception Handling module no supported by Fotran compiler, but intrinsic ISNAN exists")
        add_definitions(-DHAVE_INTRINSIC_ISNAN)
    else()
        message(STATUS "-- IEEE Exception Handling AND intrisic ISNAN module not supported by Fotran compiler")
    endif()

endif()
