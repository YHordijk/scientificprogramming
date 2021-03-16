find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
    configure_file(
        ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile
    )
    add_custom_target(
        doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} -d validate -d filteroutput -d commentscan
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
endif()
