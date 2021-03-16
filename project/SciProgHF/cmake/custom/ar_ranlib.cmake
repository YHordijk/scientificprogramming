# ar_ranlib.sh; hjaaj January 2017 
#

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")

    set(CMAKE_AR     ${PROJECT_SOURCE_DIR}/cmake/custom/ar.sh)
    set(CMAKE_RANLIB ${PROJECT_SOURCE_DIR}/cmake/custom/ranlib.sh)

    message(STATUS "INFO: ar and ranlib commands modified for Darwin")
    message(STATUS "-- new CMAKE_AR: ${CMAKE_AR}")
    message(STATUS "-- new CMAKE_RANLIB: ${CMAKE_RANLIB}")

endif()
