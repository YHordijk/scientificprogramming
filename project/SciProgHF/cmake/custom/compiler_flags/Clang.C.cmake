if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    if(CMAKE_C_COMPILER_ID MATCHES Clang)
        # hjaaj March 2018: cmake by mistake identifies xlc compiler as Clang 4.0.1
        #  ( probably fixed in cmake version > 3.11.x; see https://gitlab.kitware.com/cmake/cmake/issues/17784 )
        message(STATUS "INFO: C compiler ID is changed from Clang to XL because of cmake bug for XL C compiler")
        set(CMAKE_C_COMPILER_ID "XL")
        unset(CMAKE_C_COMPILER_VERSION) # will be wrong
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Clang)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
    set(CMAKE_C_FLAGS_DEBUG "-O0")
endif()
