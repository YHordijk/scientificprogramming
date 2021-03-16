#.rst:
#
# Creates the build_info.h file in the build directory.
# This file can be included into sources to print out
# build information variables to the executable program output.

set(_user_name "unknown")
execute_process(
    COMMAND whoami
    TIMEOUT 1
    OUTPUT_VARIABLE _user_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
message(STATUS "User name: ${_user_name}")

set(_host_name "unknown")
execute_process(
    COMMAND hostname
    TIMEOUT 1
    OUTPUT_VARIABLE _host_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
message(STATUS "Host name: ${_host_name}")

set(_system "unknown")
if(CMAKE_SYSTEM)
    set(_system ${CMAKE_SYSTEM})
    message(STATUS "Operating system: ${CMAKE_SYSTEM}")
endif()

set(_cmake_version "unknown")
if(CMAKE_VERSION)
    set(_cmake_version ${CMAKE_VERSION})
    message(STATUS "CMake version: ${CMAKE_VERSION}")
endif()

set(_cmake_generator "unknown")
if(CMAKE_GENERATOR)
    set(_cmake_generator ${CMAKE_GENERATOR})
    message(STATUS "CMake generator: ${CMAKE_GENERATOR}")
endif()

set(_cmake_build_type "unknown")
if(CMAKE_BUILD_TYPE)
    set(_cmake_build_type ${CMAKE_BUILD_TYPE})
    message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
endif()

set(_configuration_time "unknown")
# get configuration time in UTC
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import datetime; print(datetime.datetime.utcnow())"
    OUTPUT_VARIABLE _configuration_time
    )
string(STRIP ${_configuration_time} _configuration_time) # delete newline
message(STATUS "Configuration time: ${_configuration_time}")

set(_python_version "unknown")
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -V
    OUTPUT_VARIABLE _python_version
    ERROR_VARIABLE _python_version
    )
string(REGEX MATCH "Python ([0-9].[0-9].[0-9])" temp "${_python_version}")
set(_python_version ${CMAKE_MATCH_1})
message(STATUS "Python version: ${_python_version}")

foreach(_lang Fortran C CXX)
    set(_${_lang}_compiler "unknown")
    if(CMAKE_${_lang}_COMPILER)
        set(_${_lang}_compiler ${CMAKE_${_lang}_COMPILER})
        message(STATUS "${_lang} compiler: ${CMAKE_${_lang}_COMPILER}")
    endif()

    set(_python_${_lang}_version "unknown")
    if(PYTHON_${_lang}_VERSION)
        set(_python_${_lang}_version ${PYTHON_${_lang}_VERSION})
        message(STATUS "${_lang} compiler version: ${CMAKE_${_lang}_COMPILER_ID} ${PYTHON_${_lang}_VERSION}")
    endif()

    set(_${_lang}_compiler_flags "unknown")
    if(CMAKE_${_lang}_FLAGS)
        set(_${_lang}_compiler_flags ${CMAKE_${_lang}_FLAGS})
        message(STATUS "${_lang} compiler flags: ${CMAKE_${_lang}_FLAGS}")
        if(ENABLE_RUNTIMECHECK)
            message(STATUS "${_lang} appended runtime-check compiler flags: ${CMAKE_${_lang}_FLAGS_runtimecheck}")
        endif()
    endif()
endforeach()

set(_static_linking ${ENABLE_STATIC_LINKING})
message(STATUS "Static linking: ${ENABLE_STATIC_LINKING}")

set(_enable_64bit_integers ${ENABLE_64BIT_INTEGERS})
message(STATUS "64-bit integers: ${ENABLE_64BIT_INTEGERS}")

set(_enable_mpi ${ENABLE_MPI})
message(STATUS "MPI parallelization: ${ENABLE_MPI}")

set(_mpi_launcher "unknown")
if(MPI_FOUND)
    set(_mpi_launcher ${MPIEXEC})
    message(STATUS "MPI launcher: ${MPIEXEC}")
endif()

set(_math_libs "unknown")
if(MATH_LIBS)
    set(_math_libs ${MATH_LIBS})
    message(STATUS "Math libraries: ${MATH_LIBS}")
endif()

set(_enable_builtin_blas ${ENABLE_BUILTIN_BLAS})
if(ENABLE_BUILTIN_BLAS)
    message(STATUS "Builtin BLAS: ${ENABLE_BUILTIN_BLAS}")
endif()

set(_enable_builtin_lapack ${ENABLE_BUILTIN_LAPACK})
if(ENABLE_BUILTIN_LAPACK)
    message(STATUS "Builtin LAPACK: ${ENABLE_BUILTIN_LAPACK}")
endif()

set(_explicit_libs "unknown")
if(EXPLICIT_LIBS)
    set(_explicit_libs ${EXPLICIT_LIBS})
    message(STATUS "Explicit libraries: ${EXPLICIT_LIBS}")
endif()

if(MKL_FLAG)
    message(STATUS "Intel MKL flag: ${MKL_FLAG}")
endif()

get_directory_property(_list_of_definitions DIRECTORY ${PROJECT_SOURCE_DIR}/src COMPILE_DEFINITIONS)
message(STATUS "Compile definitions: ${_list_of_definitions}")

# add definitions for external libraries used by dirac
# exatensor
set(_exatensor_git_repo_location "unknown")
set(_exatensor_configuration "unknown")
set(_exatensor_git_hash "unknown")

set(_enable_exatensor ${ENABLE_EXATENSOR})
message(STATUS "Exacorr module enabled : ${ENABLE_EXATENSOR}")

if(ENABLE_EXATENSOR)
   message(STATUS "The Exacorr module will be included in the Dirac executable and to the standalone exacorr.x")
   message(STATUS "Exacorr employs the ExaTensor library (https://github.com/ORNL-QCI/ExaTENSOR) for tensor operations")
   message(STATUS "Please read carefully the Dirac documentation guide for setting up the Exacorr runtime environment")

   set(_exatensor_git_repo_location ${EXATENSOR_GIT_REPO_LOCATION})
   message(STATUS "ExaTensor source code repository: ${EXATENSOR_GIT_REPO_LOCATION}")

   set(_exatensor_git_hash ${EXATENSOR_GIT_HASH})
   message(STATUS "Exatensor source code git hash: ${EXATENSOR_GIT_HASH}")

   set(_exatensor_configuration ${EXATENSOR_ENV})
   message(STATUS "ExaTensor build environment: ${EXATENSOR_ENV}")
endif()

# generate the build_info.h include file with defined set of variables
get_filename_component(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/build_info.h.in
    ${PROJECT_BINARY_DIR}/build_info.h
    @ONLY
    )

add_custom_target(
    build_info
    ALL DEPENDS ${PROJECT_BINARY_DIR}/build_info.h
    )
