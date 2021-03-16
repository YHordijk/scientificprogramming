cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(pelib Fortran)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
   )

option(ENABLE_MPI               "Enable MPI"                    OFF)
option(ENABLE_MPIF              "Use mpif.h instead of mpi.mod" OFF)
option(ENABLE_64BIT_INTEGERS    "Enable 64-bit integers"        OFF)
option(ENABLE_STATIC_LINKING    "Enable static linking"         OFF)

if(ENABLE_64BIT_INTEGERS)
    message("-- Using 64-bit integers")
    if(ENABLE_MPI)
        set(ENABLE_MPIF ON)
        message("-- WARNING: Using mpif.h instead of mpi.mod, which can cause conversion")
        message("            issues if MPI is not built with 64-bit integers.")
    endif()
endif()

if(ENABLE_MPI AND NOT DEFINED HOST_PROGRAM)
    find_package(MPI REQUIRED)
    if(MPI_Fortran_FOUND)
        if(NOT ${MPI_Fortran_COMPILER} MATCHES ${CMAKE_Fortran_COMPILER})
            message("-- WARNING: Detected MPI Fortran compiler does not match Fortran compiler")
#            set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
        endif()
        include_directories(${MPI_Fortran_INCLUDE_PATH})
        add_definitions(-DVAR_MPI)
        if(NOT ENABLE_MPIF)
            add_definitions(-DUSE_MPI_MOD_F90)
        endif()
    # the following use of MPI_ variables is deprecated
    elseif(MPI_FOUND)
        include_directories(${MPI_INCLUDE_PATH})
        add_definitions(-DVAR_MPI)
        if(NOT ENABLE_MPIF)
            add_definitions(-DUSE_MPI_MOD_F90)
        endif()
    else()
        message(FATAL_ERROR "-- MPI Fortran not found")
    endif()
elseif(ENABLE_MPI AND DEFINED HOST_PROGRAM)
    add_definitions(-DVAR_MPI)
    if(NOT ENABLE_MPIF)
        add_definitions(-DUSE_MPI_MOD_F90)
    endif()
endif()

if("${INTEGRAL_LIBRARY}" STREQUAL "GEN1INT")
    add_definitions(-DGEN1INT_INTEGRALS)
    message("-- Using Gen1Int integral library")
elseif("${INTEGRAL_LIBRARY}" STREQUAL "OFF")
    add_definitions(-DDUMMY_INTEGRALS)
    message("-- Using dummy integral interfaces")
else()
    message(FATAL_ERROR "-- No integral library enabled")
endif()

if("${HOST_PROGRAM}" STREQUAL "DIRAC")
    add_definitions(-DPRG_DIRAC)
    message("-- Host program set to DIRAC")
elseif("${HOST_PROGRAM}" STREQUAL "DALTON")
    message("-- Host program set to Dalton")
else()
    message("-- No host program set")
endif()

if(DEFINED PARENT_INCLUDE_DIR)
    include_directories(${PARENT_INCLUDE_DIR})
endif()

if(NOT DEFINED PARENT_MODULE_DIR)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PARENT_MODULE_DIR})
endif()
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

if(DEFINED PARENT_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${PARENT_INSTALL_PREFIX})
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT DEFINED PARENT_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE release CACHE STRING "-- Release build (options are: none release debug profile)" FORCE)
elseif(PARENT_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${PARENT_BUILD_TYPE})
endif()

include(ConfigCompilerFlags)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/pe_mpi.F90
    ${PROJECT_SOURCE_DIR}/src/pe_precision.F90
    ${PROJECT_SOURCE_DIR}/src/pe_constants.F90
    ${PROJECT_SOURCE_DIR}/src/pe_variables.F90
    ${PROJECT_SOURCE_DIR}/src/pe_utils.F90
    ${PROJECT_SOURCE_DIR}/src/pe_operators.F90
    ${PROJECT_SOURCE_DIR}/src/pe_multipole_tools.F90
    ${PROJECT_SOURCE_DIR}/src/pe_analysis_tools.F90
    ${PROJECT_SOURCE_DIR}/src/pe_potential_derivatives.F90
    ${PROJECT_SOURCE_DIR}/src/pe_induced_moments.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pe_blas_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pe_lapack_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pe_integral_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/polarizable_embedding.F90
    )

#add_library(pelib SHARED ${SOURCES})
#install(TARGETS pelib LIBRARY DESTINATION lib)

add_library(pelib STATIC ${SOURCES})
install(TARGETS pelib ARCHIVE DESTINATION lib)
