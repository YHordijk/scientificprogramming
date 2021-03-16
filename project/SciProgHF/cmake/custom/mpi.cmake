#
#
# depends on CMake's default module FindMPI (on Ubuntu, /usr/share/cmake-2.8/Modules/FindMPI.cmake)
#
# Also, check current, up-to-date FindMPI module, https://cmake.org/cmake/help/v3.0/module/FindMPI.html
#
#
if(MPI_FOUND)
    if (${MPI_Fortran_VERSION_MAJOR})
        message(STATUS "MPI version:  ${MPI_Fortran_VERSION_MAJOR}.${MPI_Fortran_VERSION_MINOR}")
    endif()

    include(CheckFortranSourceCompiles1)

#Miro: verify that Fortran "MPI"-compiler agrees with "CMake" compiler
    if (NOT ("${CMAKE_Fortran_COMPILER}" STREQUAL "${MPI_Fortran_COMPILER}"))
      message(STATUS "Info: MPI Fortran and CMake Fortran compilers differ, see ")
      message(STATUS "...CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
      message(STATUS ".....MPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}")
    endif()

    add_definitions(-DVAR_MPI)
    add_definitions(-DVAR_MPI2)

    # test whether MPI module is compatible with compiler

     message(STATUS "Compiling ${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-compiler-compatibility.F90 ...")
    file(READ "${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-compiler-compatibility.F90" _source)
    check_fortran_source_compiles1(
        ${_source}
        MPI_COMPILER_MATCHES
        )
    if(MPI_COMPILER_MATCHES)
        message(STATUS "mpi.mod matches the current MPI compiler, employing 'use mpi' and setting -DUSE_MPI_MOD_F90")
        add_definitions(-DUSE_MPI_MOD_F90)
    else()
        message(STATUS "mpi.mod is compiled with different compiler or compiler flags, using '#include \"mpif.h\"' instead")
    endif()

    # test whether MPI integer type matches
    message(STATUS "Compiling ${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-itype-compatibility.F90 ...")
    file(READ "${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-itype-compatibility.F90" _source)
    check_fortran_source_compiles1(
        ${_source}
        MPI_ITYPE_MATCHES
        )
    if(NOT MPI_ITYPE_MATCHES)
        if(ENABLE_64BIT_INTEGERS)
            message(STATUS "No 64-bit integer MPI interface found, will use 32-bit integer MPI interface")
            add_definitions(-DVAR_MPI_32BIT_INT)
            set(USE_32BIT_MPI_INTERFACE TRUE)
        else()
            message(STATUS "Cannot determine whether MPI is built for 32bit integers, check it yourself ")
        endif()
    else()
        if(ENABLE_64BIT_INTEGERS)
            message(STATUS "Well, program got compiled with 64-bit integers, but the internal MPI integer compatibility test will be performed inside dirac.x.")
        else()
            message(STATUS "Well, program got compiled with 32-bit integers, but the internal MPI integer compatibility test will be performed inside dirac.x.")
        endif()
    endif()

    if(ENABLE_64BIT_INTEGERS AND FORCE_32BIT_MPI_INTERFACE)
        message(STATUS "32-bit integer MPI interface activated by the user for 64-bit DIRAC !")
        add_definitions(-DVAR_MPI_32BIT_INT)
        set(USE_32BIT_MPI_INTERFACE TRUE)
    endif()

endif()
