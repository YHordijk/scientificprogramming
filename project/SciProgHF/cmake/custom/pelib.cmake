option(ENABLE_PELIB "Enable polarizabile embedding library" ON)

if(ENABLE_PELIB)
    include_directories(${PROJECT_BINARY_DIR}/external/pelib-build)

    add_definitions(-DHAS_PELIB)
    set(PE_HOST_PROGRAM "DIRAC")
    if(ENABLE_GEN1INT)
        set(PE_INTEGRAL_LIBRARY "GEN1INT")
    else()
        message(FATAL_ERROR "-- PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
    endif()

    set(PE_INCLUDE_DIR)

    # here we control whether PElib will use mpi.mod or mpif.h
    # MPI_COMPILER_MATCHES is defined in cmake/custom/mpi.cmake
    # TODO: we should use generator expressions instead but this probably requires
    # moving away from add_external() which we should also do at some point
    if(MPI_COMPILER_MATCHES)
      set(_enable_mpif "0")
    else()
      set(_enable_mpif "1")
    endif()

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DPARENT_INCLUDE_DIR=${PE_INCLUDE_DIR}
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DINTEGRAL_LIBRARY=${PE_INTEGRAL_LIBRARY}
	-DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_MPIF=${_enable_mpif}
        -DHOST_PROGRAM=${PE_HOST_PROGRAM}
    )

    unset(_enable_mpif)

    add_external(pelib)


    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libpelib.a
        ${EXTERNAL_LIBS})
endif()

message(STATUS "PElib module: ${ENABLE_PELIB}")
