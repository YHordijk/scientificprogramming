include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/luciarel # local include files
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libluciarel.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_LUCIAREL_FORTRAN_SOURCES
   io_r.F
   par_lushmem_r.F
   par_solver_r.F
   utils_r.F
   vecmat_r.F
)

#Miro: I had to hide sources from runtimechecking
set(FIXED_LUCIAREL_FORTRAN_SOURCES_noruntimecheck
   par_lucia_r.F
   property_r.F
   diagonalize_r.F
   diagonal_r.F
   gasci_main_driver_r.F
   memory_r.F
   strings_r.F
   ciinfo_r.F
   sigden_r.F
   interface_r.F
   resolution_utils_r.F
)

set(FREE_LUCIAREL_FORTRAN_SOURCES
   io_ttss_r.F90
   krci_cfg.F90
   krci_prop.F90
   mospinor_info.F90
   string_info.F90
   symmetry_setup.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_LUCIAREL_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_LUCIAREL_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_LUCIAREL_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'luciarel' module, with many exceptions")
    set_source_files_properties(${FREE_LUCIAREL_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_LUCIAREL_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()


# Create the module library
add_library(
    luciarel
    OBJECT
    ${FIXED_LUCIAREL_FORTRAN_SOURCES}
    ${FIXED_LUCIAREL_FORTRAN_SOURCES_noruntimecheck}
    ${FREE_LUCIAREL_FORTRAN_SOURCES}
)

# deal with intermodule dependencies
add_dependencies(luciarel gp)
add_dependencies(luciarel interface_mpi)
