#.rst:
#
# placeholder for handling detection of laplace module. now will only be used to enable/disable module
#

option(ENABLE_LAPLACE "Enable Laplace transform code" ON)

if(ENABLE_LAPLACE)
  message(STATUS "Laplace transform code ENABLED")
  add_definitions(-DHAS_LAPLACE)
else()
  message(STATUS "Laplace transform code DISABLED")
endif()
