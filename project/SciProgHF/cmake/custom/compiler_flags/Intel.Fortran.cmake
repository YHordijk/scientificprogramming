if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -assume byterecl -g -traceback -DVAR_IFORT")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif()
