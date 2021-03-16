#include <stdio.h>

#include "git_info.h"
#include "build_info.h"

/*
 * Print the configuration/build info generated at build time by the executable.
 * We do this from C and not from Fortran because Fortran has trouble with long lines.
 */

#if defined VAR_XLF
/* Miro: xlf90 linker requires underscore */
void print_info_(void)
#else
void print_info(void)
#endif
{
    printf("\n");
    printf("\n");
    printf("Version information\n");
    printf("-------------------\n");
    printf("\n");
    printf("Branch        | %s\n", GIT_BRANCH);
    printf("Commit hash   | %s\n", GIT_COMMIT_HASH);
    printf("Commit author | %s\n", GIT_COMMIT_AUTHOR);
    printf("Commit date   | %s\n", GIT_COMMIT_DATE);

    printf("\n");
    printf("\n");
    printf("Configuration and build information\n");
    printf("-----------------------------------\n");
    printf("\n");
    printf("Who compiled             | %s\n", USER_NAME);
    printf("Compiled on server       | %s\n", HOST_NAME);
    printf("Operating system         | %s\n", SYSTEM);
    printf("CMake version            | %s\n", CMAKE_VERSION);
    printf("CMake generator          | %s\n", CMAKE_GENERATOR);
    printf("CMake build type         | %s\n", CMAKE_BUILD_TYPE);
    printf("Configuration time       | %s\n", CONFIGURATION_TIME);
    printf("Python version           | %s\n", PYTHON_VERSION);
    printf("Fortran compiler         | %s\n", FORTRAN_COMPILER);
    printf("Fortran compiler version | %s\n", FORTRAN_COMPILER_VERSION);
    printf("Fortran compiler flags   | %s\n", FORTRAN_COMPILER_FLAGS);
    printf("C compiler               | %s\n", C_COMPILER);
    printf("C compiler version       | %s\n", C_COMPILER_VERSION);
    printf("C compiler flags         | %s\n", C_COMPILER_FLAGS);
    printf("C++ compiler             | %s\n", CXX_COMPILER);
    printf("C++ compiler version     | %s\n", CXX_COMPILER_VERSION);
    printf("C++ compiler flags       | %s\n", CXX_COMPILER_FLAGS);
    printf("Static linking           | %s\n", STATIC_LINKING);
    printf("64-bit integers          | %s\n", ENABLE_64BIT_INTEGERS);
    printf("MPI parallelization      | %s\n", ENABLE_MPI);
    printf("MPI launcher             | %s\n", MPI_LAUNCHER);
    printf("Math libraries           | %s\n", MATH_LIBS);
    printf("Builtin BLAS library     | %s\n", ENABLE_BUILTIN_BLAS);
    printf("Builtin LAPACK library   | %s\n", ENABLE_BUILTIN_LAPACK);
    printf("Explicit libraries       | %s\n", EXPLICIT_LIBS);
    printf("Compile definitions      | %s\n", DEFINITIONS);
#if defined ENABLE_EXATENSOR
    printf("\n");
    printf("EXACORR dependencies\n");
    printf("--------------------\n");
    printf("Exatensor source repo    | %s\n", EXATENSOR_REPO);
    printf("Exatensor git hash       | %s\n", EXATENSOR_HASH); 
    printf("Exatensor configuration  | %s\n", EXATENSOR_CONFIG);
#endif
    printf("\n");
    printf("\n");

    fflush(stdout);
}
