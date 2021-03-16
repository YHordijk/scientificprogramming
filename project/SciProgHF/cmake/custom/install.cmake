file(COPY ${PROJECT_SOURCE_DIR}/basis        DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/basis_dalton DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/basis_ecp    DESTINATION ${PROJECT_BINARY_DIR})
message(STATUS "Copied DIRAC basis set directories into the build directory")

foreach(
    EXECUTABLE
    ${_list_of_executables}
    )
    install(
        TARGETS ${EXECUTABLE}
        DESTINATION share/dirac
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

install(
    FILES ${PROJECT_BINARY_DIR}/pam
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

# workaround to install pam-dirac symlink:
# 1) copy pam to pam-dirac
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/pam ${PROJECT_BINARY_DIR}/pam-dirac)")
# 2) install real file pam-dirac
install(
    FILES ${PROJECT_BINARY_DIR}/pam-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )
# 3) remove real file
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")
# 4) create symlink
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/pam ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")

# write git hash to build dir
file(WRITE ${PROJECT_BINARY_DIR}/GIT_HASH "${_git_last_commit_hash}")
# copy version info to install dir
install(
    FILES ${PROJECT_BINARY_DIR}/GIT_HASH ${PROJECT_SOURCE_DIR}/VERSION
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_SOURCE_DIR}/basis_dalton
    ${PROJECT_SOURCE_DIR}/basis_ecp
    DESTINATION share/dirac
    PATTERN .git EXCLUDE
    )
