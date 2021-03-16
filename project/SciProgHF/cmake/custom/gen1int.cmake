option(ENABLE_GEN1INT "Enable Gen1Int library" ON)

if(ENABLE_GEN1INT)

add_definitions(-DBUILD_GEN1INT)

include_directories(${PROJECT_BINARY_DIR}/external/gen1int-build)

    add_library(
        gen1int_interface
        src/gen1int/gen1int_cube.F90
        src/gen1int/gen1int_api.F90
        src/gen1int/gen1int_matrix.F90
        src/gen1int/gen1int_host.F90
        src/gen1int/gen1int_shell.F90
        )   

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_INCLUDE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    )

# added -DPARENT_INCLUDE_DIR line above;
# we need to add the module directory to the include directories explicitly,
# otherwise the IBM xlf compiler will not find the .mod files. /hjaaj Aug 2017

    add_external(gen1int)

    set(EXTERNAL_LIBS
        gen1int_interface
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${EXTERNAL_LIBS})
endif()

message(STATUS "Gen1Int module: ${ENABLE_GEN1INT}")
