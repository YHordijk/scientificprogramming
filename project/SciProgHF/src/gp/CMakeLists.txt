include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_BINARY_DIR}/src/gp
    ${PROJECT_BINARY_DIR}
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libgp.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_GP_FORTRAN_SOURCES
   finesort_scheme6.F
   global_file.F
   gphjj.F
   gpjeppe.F
   gpjth.F
   gpkrmc.F
   gplupar.F
   gpluuk.F
   gpsaue.F
   gptrygve.F
   gpulfek.F
   graphic.F
   io.F
   matexport.F
   math.F
   mempkg.F
   mpi2macro.F
   mxmpack.F
   order2.F
   qalgebra.F
   qpack.F
   spline.F
   time.F
)

set(FIXED_GP_FORTRAN_SOURCES_noruntimecheck
# Miro: I moved problematic source files here
   gplucia_io.F
   gplucia.F
)

set(FREE_GP_FORTRAN_SOURCES
    codata.F90
    debug_mem_alloc.F90
    dirac_cfg.F90
    dirac_mocoefficients.F90
    dynamic_string.F90
    electrostatic_potential.F90
    exception_handling.F90
    file_handler.F90
    f77_work_interface.F90
    gp_q_diag.F90
    gpsk.F90
    heapsort.F90
    include_h_modules.F90
    kinds.F90
    machine_parameters.F90
    memory_allocator.F90
    memory_control.F90
    memory_errorhandler.F90
    memory_internal_special_i8.F90
    memory_internal_ada_1d.F90
    memory_internal_ada_2d.F90
    memory_internal_ada_3d.F90
    memory_internal_ada_4d.F90
    memory_internal_init.F90
    memory_legacy_dirac.F90
    memory_parameters.F90
    memory_track_if.F90
    memory_tracker.F90
    moltra_labeling.F90
    numerical_operators.F90
    num_der.F90
    orbital_rotation_indices.F90
    osutils.F90
    picture_change_operations.F90
    quaternion_algebra.F90
    selftests.F90
    string_manipulations.f90
    type_laplace.F90
    unit_test_generator.f90
    unit_testing.F90
    vector_functions.F90
    xmlout.F90
    xml/xml_file.F90
    xml/xml_structure.F90
    xml/xml_parser.F90
)

set(GP_C_SOURCES
    execinfo_bsd.c
    gpc.c
    memory_sys_if.c
    print_info.c
    selftest_fortran_c_interoperability.c
    stacktraverse_bsd.c
    trace.c
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_GP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_GP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_GP_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'gp' module with some source files excluded")
    set_source_files_properties(${FREE_GP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_GP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${GP_C_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS_runtimecheck})
endif()

add_library(
    gp
    OBJECT
    ${FIXED_GP_FORTRAN_SOURCES}
    ${FIXED_GP_FORTRAN_SOURCES_noruntimecheck}
    ${FREE_GP_FORTRAN_SOURCES}
    ${GP_C_SOURCES}
)

# add intermodule dependencies
add_dependencies(gp interface_mpi)
add_dependencies(gp xcint)
add_dependencies(gp pdpack)
generate_git_info_header(${PROJECT_BINARY_DIR} git_info.h)
