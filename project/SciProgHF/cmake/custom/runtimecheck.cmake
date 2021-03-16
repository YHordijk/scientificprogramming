#.rst:
#
# Enables run-time checking by defining corresponding compiler flags.
#
# Variables defined (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS_runtimecheck
#   CMAKE_C_FLAGS_runtimecheck
#   CMAKE_CXX_FLAGS_runtimecheck
#
# autocmake.yml configuration::
#
#   docopt: "--check Enable runtime checking of compiled executables [default: False]"
#   define: "'-DENABLE_RUNTIMECHECK={0}'.format(arguments['--check'])"
#
# Important change, Miro, July 2018: do not append corresponding compiler flags here !
# Rather append the runtime-check flags separately to source files in corresponding DIRAC modules 
#
#
option(ENABLE_RUNTIMECHECK "Enable run-time checking" OFF)

message(STATUS "Enable run-time checking: ${ENABLE_RUNTIMECHECK}")

if(ENABLE_RUNTIMECHECK)
    if(DEFINED CMAKE_Fortran_COMPILER_ID)
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(CMAKE_Fortran_FLAGS_GNU_runtimecheck "-fbounds-check")
            set(CMAKE_Fortran_FLAGS_runtimecheck  ${CMAKE_Fortran_FLAGS_GNU_runtimecheck})
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            set(CMAKE_Fortran_FLAGS_Intel_runtimecheck "-check all -fpstkchk")
            set(CMAKE_Fortran_FLAGS_runtimecheck  ${CMAKE_Fortran_FLAGS_Intel_runtimecheck})
        endif()
    endif()

    if(DEFINED CMAKE_C_COMPILER_ID)
        if(CMAKE_C_COMPILER_ID MATCHES GNU)
            set(CMAKE_C_FLAGS_GNU_runtimecheck "-fbounds-check")
            set(CMAKE_C_FLAGS_runtimecheck  ${CMAKE_C_FLAGS_GNU_runtimecheck})
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES Intel)
            set(CMAKE_C_FLAGS_Intel_runtimecheck "-fp-stack-check -check-uninit")
            set(CMAKE_C_FLAGS_runtimecheck  ${CMAKE_C_FLAGS_Intel_runtimecheck})
        endif()
    endif()

    if(DEFINED CMAKE_CXX_COMPILER_ID)
        if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
            set(CMAKE_CXX_FLAGS_GNU_runtimecheck "-fbounds-check")
            set(CMAKE_CXX_FLAGS_runtimecheck  ${CMAKE_CXX_FLAGS_GNU_runtimecheck})
            #set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_GNU_runtimecheck}")
        endif()
        if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
            set(CMAKE_CXX_FLAGS_Intel_runtimecheck "-fp-stack-check -check-uninit -early-template-check")
            set(CMAKE_CXX_FLAGS_runtimecheck  ${CMAKE_CXX_FLAGS_Intel_runtimecheck})
            #set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_Intel_runtimecheck}")
        endif()
    endif()
    # control printout - is done in custom/build_info.cmake
endif()
