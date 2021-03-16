set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "R. Bast, H. J. Aa. Jensen, T. Saue, and L. Visscher, and contributors")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "R. Bast")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/doc/installation/general.rst")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_FILE_NAME "DIRAC")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "DIRAC-${PROGRAM_VERSION}-Source")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "DIRAC ${PROGRAM_VERSION}")

set(EXPORT_DIR ${PROJECT_BINARY_DIR}/export)
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "${EXPORT_DIR};/")

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include(CPack)
    message(STATUS "Adding target release")
    add_custom_target(release
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/fetch_external_sources ${PROJECT_SOURCE_DIR}
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/src                  ${EXPORT_DIR}
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_unreleased_code  ${EXPORT_DIR} ${CODE_REMOVAL_FLAGS}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/basis                ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/basis_dalton         ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/basis_ecp            ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/test                 ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CMakeLists.txt       ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/LICENSE              ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CHANGELOG.rst        ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/CTestConfig.cmake    ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/atomic-start-x2c-pam ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/cmake                ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/external             ${EXPORT_DIR}
        COMMAND cp -RL ${PROJECT_SOURCE_DIR}/doc                  ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/pam.in               ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/setup                ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/utils                ${EXPORT_DIR}
        COMMAND cp -R  ${PROJECT_SOURCE_DIR}/src/xcfun            ${EXPORT_DIR}/src
        COMMAND cp     ${PROJECT_SOURCE_DIR}/src/main/main.F90    ${EXPORT_DIR}/src/main
        COMMAND cp     ${PROJECT_SOURCE_DIR}/VERSION              ${EXPORT_DIR}
        COMMAND cp     ${PROJECT_SOURCE_DIR}/.gitignore           ${EXPORT_DIR}
        COMMAND echo "${GIT_COMMIT_HASH}" > ${EXPORT_DIR}/cmake/GIT_HASH
        COMMAND echo "set(IS_EXPORTED_TARBALL TRUE)" > ${EXPORT_DIR}/cmake/custom/exported.cmake
        COMMAND rm -f  ${EXPORT_DIR}/external/pcmsolver/.git
        COMMAND rm -f  ${EXPORT_DIR}/external/unit-tests/.git
        COMMAND rm -f  ${EXPORT_DIR}/external/stieltjes/.git
        COMMAND rm -f  ${EXPORT_DIR}/external/laplace-minimax/.git
        COMMAND rm -rf ${EXPORT_DIR}/src/openrsp/
        COMMAND rm -rf ${EXPORT_DIR}/src/aoosoc/
        COMMAND rm -rf ${EXPORT_DIR}/src/eri/
        COMMAND rm -rf ${EXPORT_DIR}/src/densfit/
        COMMAND rm -rf ${EXPORT_DIR}/src/prp/esr/
        COMMAND rm -rf ${EXPORT_DIR}/src/krcc/
        COMMAND rm -rf ${EXPORT_DIR}/src/srdft/
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_copyright_tags ${EXPORT_DIR}
        COMMAND ${PROJECT_SOURCE_DIR}/maintenance/remove_basis_sets ${EXPORT_DIR}
        COMMAND make package_source
        COMMENT "Packaging source files (this may take several minutes)"
        VERBATIM
        )
endif()
