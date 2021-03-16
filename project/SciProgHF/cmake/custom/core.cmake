include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
    )

include(autocmake_CMakeDependentOption)

# FIXME re-introduce
# set_source_files_properties(${PROJECT_BINARY_DIR}/config_info.F90 PROPERTIES GENERATED 1)


#Miro: remove problematic linking library
list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES  "pthread")

# Ivan Hrasko's determination of compiler versions - very good tool !
include(FindCompilerVersion)

# Miro Ilias restored CMake printout, associated with build_info.h generation
include(build_info)

set(_list_of_executables)
macro(compile_standalone _executable _source _list_of_executables)
    add_executable(
        ${_executable}
        ${_source}
    )

    # each executable is of Fortran language ...
    set_property(TARGET ${_executable} PROPERTY LINKER_LANGUAGE Fortran)

    # extend list of executables
    set(_list_of_executables
        ${_list_of_executables}
        ${_executable}
    )
endmacro()

compile_standalone(dirac.x src/main/main.F90 "${_list_of_executables}")

#synchronize basis set directories of trunk/basis* into the installation directory
add_custom_command(TARGET dirac.x PRE_LINK COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis ${PROJECT_BINARY_DIR}/basis COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_ecp ${PROJECT_BINARY_DIR}/basis_ecp  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/basis_dalton ${PROJECT_BINARY_DIR}/basis_dalton COMMAND ${CMAKE_COMMAND} -E echo "INFO:basis set directories, basis*, synchronized into current installation directory")

# If exatensor is enabled we also make a stand-alone binary for convenience.
if (ENABLE_EXATENSOR)
   compile_standalone(exacorr.x src/exacorr/exacorr_main.F90   "${_list_of_executables}")
endif()

option(ENABLE_UTILITIES "Enable utilities" ON)
if(ENABLE_UTILITIES)
    compile_standalone(dirac_mointegral_export.x utils/dirac_mointegral_export.F90   "${_list_of_executables}")
    compile_standalone(rspread.x                 utils/rspread.F90                   "${_list_of_executables}")
    compile_standalone(labread.x                 utils/labread.F90                   "${_list_of_executables}")
    compile_standalone(twofit.x                  utils/twofit.F90                    "${_list_of_executables}")
    compile_standalone(cfread.x                  utils/cfread.F90                    "${_list_of_executables}")
    compile_standalone(cf_addlabels.x            utils/cf_addlabels.F90              "${_list_of_executables}")
    compile_standalone(pcmo_addlabels.x          utils/pcmo_addlabels.F90            "${_list_of_executables}")    
    compile_standalone(vibcal.x                  utils/vibcal.F90                    "${_list_of_executables}")
    compile_standalone(polfit.x                  utils/polfit.F90                    "${_list_of_executables}")
    compile_standalone(mx2fit.x                  utils/mx2fit.F90                    "${_list_of_executables}")
    compile_standalone(diag.x                    utils/diag.F90                      "${_list_of_executables}")
    compile_standalone(test_allocator.x          utils/test_allocator.F90            "${_list_of_executables}")
#   compile_standalone(test_davidson.x           utils/test_davidson.F90             "${_list_of_executables}")
endif()

if(ENABLE_GEN1INT)
    add_dependencies(dirac gen1int)
endif()

if(ENABLE_PELIB)
    add_dependencies(gen1int_interface gen1int)
    add_dependencies(pelib gen1int_interface)
    add_dependencies(dirac pelib)
endif()

if(ENABLE_RELCCSD_STANDALONE)
    if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
        set_source_files_properties(src/relccsd/ccmain.F PROPERTIES COMPILE_FLAGS "-qfixed")
    endif()
    compile_standalone(relccsd.x                 src/relccsd/ccmain.F                "${_list_of_executables}")
# miro: take care of object lib
    target_link_libraries(relccsd.x objlib.relccsd.x)
endif()

option(ENABLE_HSCC "Enable compilation of the standalone hsmrcc.x" OFF)
if(ENABLE_HSCC) # note: this functionality depends on the pam flag
    compile_standalone(hsmrcc.x src/cc_external/hsmrcc/hsmrmain.c "${_list_of_executables}")
    set_property(TARGET hsmrcc.x PROPERTY LINKER_LANGUAGE C)
endif()

option(ENABLE_HSFS "Enable compilation of the standalone hsfscc.x" OFF)
if(ENABLE_HSFS) # note: this functionality depends on the pam flag
    compile_standalone(hsfscc.x src/cc_external/hsfscc/hsfs_main.F90 "${_list_of_executables}")
endif()

foreach(
    _executable
    ${_list_of_executables}
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${_executable} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        SET_TARGET_PROPERTIES(${_executable} PROPERTIES LINK_FLAGS "-Wl,-E")
    endif()

#miro: libraries are OBJECT type
    add_library(objlib.${_executable} ${EXTERNAL_OBJECTS})
    SET_TARGET_PROPERTIES(objlib.${_executable}  PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries( ${_executable} objlib.${_executable} ${EXTERNAL_LIBS} )

    if(ENABLE_STATIC_LINKING)
        if(ENABLE_OLD_LINKER)
           set(_no_dyn_export "")
        else()
           # Miro: the "-Wl,--no-export-dynamic" flag is only for newer
           # ld-linkers, fix for older linkers
           set(_no_dyn_export "-Wl,--no-export-dynamic")
        endif()
        target_link_libraries(${_executable} ${_no_dyn_export})
        unset(_no_dyn_export)
        if (BLAS_TYPE MATCHES Atlas)
            # activate extra linking flags for dirac.x only if this bug occurs, please do not delete
            # http://gcc.gnu.org/bugzilla/show_bug.cgi?id=30471#c7
            # this bug is: signal 11 (SIGSEGV):Segmentation fault, function __restore_rt (0x1613DA0),from file sigaction.c
            target_link_libraries(${_executable} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
        endif()
        if ( CMAKE_Fortran_COMPILER_ID MATCHES GNU )
            # miro: GNU+NoneLibs/MKL/OPENBLAS-static needs this extra lib !
            target_link_libraries(${_executable} -ldl)
        endif()
        if ( BLAS_TYPE MATCHES OPENBLAS )
            # miro: OPENBLAS-static needs this extra lib !
            target_link_libraries(${_executable} -lgfortran)
        endif()
        if ( BLAS_TYPE MATCHES OPENBLAS AND CMAKE_Fortran_COMPILER_ID MATCHES Intel )
            # miro:  needs extra lib , plus remove -lpthread from CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS !
            target_link_libraries(${_executable} -lquadmath)
        endif()
      
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES PGI)
    # added -libstdc++ linking flag for each ${_executable} due to pgCC and xcfun
        target_link_libraries(${_executable} -lstdc++)
    endif()

    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    # add -lC linking flag for each ${_executable} to satisfy missing objects - sqrtf,logf,cosf,powf 
        target_link_libraries(${_executable} -lC)
    endif()

    if(ENABLE_EXATENSOR)
        # added -libstdc++ linking flag for each ${_executable} due to EXATENSOR
        target_link_libraries(objlib.${_executable} -lstdc++) 
	target_link_libraries(${_executable} -lstdc++)
    endif()
endforeach()

#printout of all linked libraries to dirac.x executable
get_target_property(OUT dirac.x  LINK_LIBRARIES)
message(STATUS "For checking, linked libraries to dirac.x: ${OUT} ")

# radovan: no! we cannot unset
#          these are used in the make install target
# unset(_list_of_executables)
