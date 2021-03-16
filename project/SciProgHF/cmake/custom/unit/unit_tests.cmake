option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)

if(ENABLE_UNIT_TESTS)
    message(STATUS "Unit control tests ENABLED")
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/unit)
    execute_process(
        COMMAND python ${PROJECT_SOURCE_DIR}/cmake/custom/unit/gather.py
                       ${PROJECT_SOURCE_DIR}/external/unit-tests
                       ${PROJECT_BINARY_DIR}/unit/list_of_unit_tests.cmake
        )
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR}/unit)

    #miro: OBJECT libraries throughout the DIRCA suite
    add_library(objlib.UNIT_TEST_LIBS ${EXTERNAL_OBJECTS})

    macro(add_unit_test _path _name)
        file(GLOB COPIED_FILES "${_path}/*")
        file(COPY ${COPIED_FILES} DESTINATION ${PROJECT_BINARY_DIR}/unit-tests/${_name})
        add_custom_target(
            ${_name}_files_copied WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} )
        add_executable(
            ${_name}.x
            ${_path}/test.F90
            )
        set_target_properties(
            ${_name}.x
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY
            "unit-tests/${_name}"
            )

        if (CMAKE_Fortran_COMPILER_ID MATCHES GNU AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        # add -lC linking flag for each unit test executable to satisfy missing objects sqrtf,logf...\
        # add -Wl,-bbigtoc flag to deal with TOC overflow
            target_link_libraries(
                ${_name}.x
                objlib.UNIT_TEST_LIB ${EXTERNAL_LIBS} -lC -Wl,-bbigtoc
                )
        else()
            target_link_libraries(
                ${_name}.x
                objlib.UNIT_TEST_LIBS ${EXTERNAL_LIBS}
                )
        endif()

        if(ENABLE_EXATENSOR)
          # added -libstdc++ linking flag for each ${_executable} due to EXATENSOR
          target_link_libraries(${_name}.x -lstdc++)
        endif()

        set_property(
            TARGET ${_name}.x
            PROPERTY
            LINKER_LANGUAGE Fortran
            )
        add_dependencies(
            ${_name}.x
            ${_name}_files_copied
            )
        add_test(
            ${_name}
            python ${PROJECT_SOURCE_DIR}/cmake/custom/unit/wrap.py --build-dir=${PROJECT_BINARY_DIR} --name=${_name}
            )
    endmacro()

    include(list_of_unit_tests)
else()
    message(STATUS "Unit control tests DISABLED")
endif()
