include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libx2c.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_X2C_FORTRAN_SOURCES
)

set(FREE_X2C_FORTRAN_SOURCES
    fragment_x2c_cfg.F90
    x2c_2e_soc_interface.F90
    x2c_cfg.F90
    x2c_common_block_interface.F90
    x2c_construct_decoupling_matrices_onmo.F90
    x2c_construct_h1_inf_2c_matrix_onmo.F90
    x2c_driver.F90
    x2c_fragment_x2c_interface.F90
    x2c_import_export_C1_matrices_ao.F90
    x2c_io_interface.F90
    x2c_picture_change_transform_mo_coefficients_ao.F90
    x2c_picture_change_transform_operator_ao.F90
    x2c_setup_basic_h1_4c_matrix_ao2onmo.F90
    x2c_setup_defining_h1_4c_matrix_ao2onmo.F90
    x2c_transform_h1_inf_2c_matrix_onmo2ao.F90
    x2c_utility_routines.F90
)


if (ENABLE_AOOSOC)
# X2C depends on some aoosoc files and vice versa
    set(FREE_X2C_FORTRAN_SOURCES
        ${FREE_X2C_FORTRAN_SOURCES}
        ${PROJECT_SOURCE_DIR}/src/aoosoc/atomic_oo_order_so_correction_cfg.F90
        ${PROJECT_SOURCE_DIR}/src/aoosoc/atomic_oo_order_so_correction.F90
        ${PROJECT_SOURCE_DIR}/src/aoosoc/atomic_oo_order_so_correction_utils.F90
        ${PROJECT_SOURCE_DIR}/src/aoosoc/module_aoosoc_pfg.F90
        ${PROJECT_SOURCE_DIR}/src/aoosoc/module_aoosoc_init.F90
    )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_X2C_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_X2C_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'x2c' module without exceptions")
    set_source_files_properties(${FREE_X2C_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_X2C_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    x2c
    OBJECT
    ${FIXED_X2C_FORTRAN_SOURCES}
    ${FREE_X2C_FORTRAN_SOURCES}
)

# add intermodule dependencies
add_dependencies(x2c gp)
if (ENABLE_INTEREST) 
    add_dependencies(x2c interest)
endif()

