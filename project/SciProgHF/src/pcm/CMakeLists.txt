# Set up compilation of PCMSolver
set(CUSTOM_PCMSolver_LOCATION ${PROJECT_BINARY_DIR}/external/pcmsolver)
get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)

# PCMSolver does not know profile
if(CMAKE_BUILD_TYPE MATCHES "profile")
  set(PCM_BUILD_TYPE "release")
else()
  set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

list(APPEND PCMSolverCMakeArgs
  -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
  -DENABLE_CXX11_SUPPORT=OFF
  -DBOOST_INCLUDEDIR=${BOOST_INCLUDE_DIRS}
  -DBOOST_LIBRARYDIR=${BOOST_LIBRARIES}
  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
  -DENABLE_TESTS=OFF
  -DENABLE_LOGGER=OFF
  -DENABLE_TIMER=OFF
  -DBUILD_STANDALONE=OFF
  -DENABLE_FORTRAN_API=ON
  -DSTATIC_LIBRARY_ONLY=ON
  -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
  -DZLIB_ROOT=${ZLIB_ROOT}
  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
  -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  -DCMAKE_INSTALL_LIBDIR=lib
  )

ExternalProject_Add(pcmsolver
  PREFIX ${CUSTOM_PCMSolver_LOCATION}
  GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
  GIT_TAG v1.1.1
  CMAKE_ARGS "${PCMSolverCMakeArgs}"
  INSTALL_DIR "${CUSTOM_PCMSolver_LOCATION}/install"
  LOG_INSTALL 1
  LOG_DOWNLOAD 1
  LOG_UPDATE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_TEST 1
  LOG_INSTALL 1
  )

# Set also variables usually set by find_package
ExternalProject_Get_Property(pcmsolver INSTALL_DIR)
set(PCMSolver_LIBRARY "${INSTALL_DIR}/lib/libpcm.a")
file(MAKE_DIRECTORY ${INSTALL_DIR}/include/pcmsolver)  # note [1] below
set(PCMSolver_INCLUDE_DIRS "${INSTALL_DIR}/include" ${ZLIB_INCLUDE_DIRS})
set(PCMSolver_LIBRARIES ${PCMSolver_LIBRARY} ${ZLIB_LIBRARIES})

# Set target for DALTON to depend upon as set by find_package
add_library(PCMSolver::PCMSolver STATIC IMPORTED GLOBAL)
add_dependencies(PCMSolver::PCMSolver pcmsolver)
set_target_properties(PCMSolver::PCMSolver PROPERTIES
  IMPORTED_LOCATION "${PCMSolver_LIBRARY}"
  INTERFACE_LINK_LIBRARIES "${PCMSolver_LIBRARIES}"
  INTERFACE_INCLUDE_DIRECTORIES "${PCMSolver_INCLUDE_DIRS}"
  )

include_directories(SYSTEM "${PCMSolver_INCLUDE_DIRS}")

link_directories(${INSTALL_DIR}/lib)
set(PCMSOLVER_LIBS
  ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
  ${PCMSolver_LIBRARIES}
  ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  PARENT_SCOPE)
set(PCMSolver_PARSE_DIR ${INSTALL_DIR}/bin PARENT_SCOPE)

include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libpcm_dirac.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_PCM_DIRAC_FORTRAN_SOURCES
)

set(FREE_PCM_DIRAC_FORTRAN_SOURCES
    pcmmod_cfg.F90
    pcm_gp.F90
    pcm_integrals.F90
    pcm_linear_response.F90
    pcm_scf.F90
    pcm_write.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_PCM_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_PCM_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'pcm_dirac' module without exceptions")
    set_source_files_properties(${FREE_PCM_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_PCM_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    pcm_dirac
    OBJECT
    ${FIXED_PCM_DIRAC_FORTRAN_SOURCES}
    ${FREE_PCM_DIRAC_FORTRAN_SOURCES}
)

add_dependencies(pcm_dirac interface_mpi)
add_dependencies(pcm_dirac pcmsolver)

# [1] It's nice to have a full PCMSolver::PCMSolver target that has embedded
# the library, the linking library paths with dependencies, and the include
# paths with dependencies, just like FindPCMSolver supplies.
# Problem is that conventional target
# derived from ExternalProject_Add is a highly sought but not quite
# certified cmake pattern. Hence INTERFACE_INCLUDE_DIRECTORIES complains
# that the directories don't exist at configure time. Hence the hack to
# create an empty directory.
