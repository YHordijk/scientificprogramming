#.rst:
#
# Enables static linking by appending corresponding compiler flags.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS
#   CMAKE_C_FLAGS
#   CMAKE_CXX_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--static Enable static linking [default: False]."
#   define: "'-DENABLE_STATIC_LINKING={0}'.format(arguments['--static'])"

option(ENABLE_STATIC_LINKING "Enable static libraries linking" OFF)

if(ENABLE_STATIC_LINKING)
    if(DEFINED CMAKE_Fortran_COMPILER_ID)
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static -static-libgcc -static-intel")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Bstatic")
        endif()
    endif()

    if(DEFINED CMAKE_C_COMPILER_ID)
        if(CMAKE_C_COMPILER_ID MATCHES GNU)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static -fpic")
        endif()
    endif()

    if(DEFINED CMAKE_CXX_COMPILER_ID)
        if(CMAKE_C_COMPILER_ID MATCHES GNU)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static -fpic")
        endif()
    endif()
endif()
