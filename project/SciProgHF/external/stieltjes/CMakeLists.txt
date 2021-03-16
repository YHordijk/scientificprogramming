cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(stieltjes)

include(CTest)

enable_language(Fortran C)
enable_testing()

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"       OFF)

# Compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -qintsize=8 -q64"
            )
    endif()
endif()

# you can add options like this
# option(ENABLE_RABOOF "Enable this and that" OFF)

# they can then be queried
# if(ENABLE_RABOOF)
#     add_definitions(-DRABOOF)
# endif()

# you can add definitions like this
# add_definitions(-DADAPTDSYEV)

add_library(
    stieltjes
    stieltjesmod.f90
    stieltjes_ql_diag.f90
    herm_spline.c
    )

add_executable(
    stieltjes.x
    driver.f90
    input_routines.f90
    )

target_link_libraries(
    stieltjes.x
    stieltjes
    )

add_test(
    main
    ${CMAKE_SOURCE_DIR}/test/test.py)

install(TARGETS stieltjes ARCHIVE DESTINATION lib)
