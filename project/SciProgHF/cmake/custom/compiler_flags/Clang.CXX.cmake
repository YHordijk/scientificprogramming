if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
        # hjaaj March 2018: cmake by mistake identifies xlc++ compiler as Clang 4.0.1
        #  ( probably fixed in cmake version > 3.11.x; see https://gitlab.kitware.com/cmake/cmake/issues/17784 )
        message(STATUS "INFO: CXX compiler ID is changed from Clang to XL because of cmake bug for XL CXX compiler")
        set(CMAKE_CXX_COMPILER_ID "XL")
        unset(CMAKE_CXX_COMPILER_VERSION) # will be wrong
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG -Wno-unused")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -DDEBUG")
endif()
