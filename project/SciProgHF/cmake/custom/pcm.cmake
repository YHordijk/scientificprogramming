#.rst:
#
# Detects PCMSolver module.
#
# Variables used::
#
#   PCMSOLVER_ROOT         - User-set path to the PCMSolver module
#   ENABLE_PCMSOLVER
#
# Variables defined::
#
#   PCMSOLVER_FOUND         - Was the PCMSolver module found
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--pcmsolver=<ENABLE_PCMSOLVER> Toggle use of PCMSolver <ON/OFF> [default: ON]."
#     - "--pcmsolver-dir=<PCMSOLVER_ROOT> The PCMSolver module location. [default: '']."
#   define:
#     - "'-DENABLE_PCMSOLVER={0}'.format(arguments['--pcmsolver'])"
#     - "'-DPCMSOLVER_ROOT={0}'.format(arguments['--pcmsolver-dir'])"

option(ENABLE_PCMSOLVER "Enable Polarizable Continuum Model through PCMSolver" ON)

if(ENABLE_PCMSOLVER)
  # User signal to try pre-built PCMSolver
  if(PCMSOLVER_ROOT)
    find_package(PCMSolver)
  endif()
  # Build PCMSolver as external package if pre-built failed or not signaled
  if(NOT PCMSolver_FOUND)
    message(STATUS "PCMSolver not found. The pre-packaged version will be built.")
    find_package(ZLIB QUIET)
    if(NOT ZLIB_FOUND)
      set(ENABLE_PCMSOLVER OFF)
      message(STATUS "Polarizable Continuum Model via PCMSolver DISABLED")
      message(STATUS " PCMSolver dependencies NOT satisfied:")
      if(NOT ZLIB_FOUND)
        IF (ENABLE_STATIC_LINKING)
          message(STATUS "  - install static Zlib development libraries")
        else()
          message(STATUS "  - install dynamic Zlib development libraries")
        endif()
      endif()
    else()
      message(STATUS "Polarizable Continuum Model via PCMSolver ENABLED")
      add_definitions(-DHAS_PCMSOLVER)
    endif()
  endif()
else()
  message(STATUS "Polarizable Continuum Model via PCMSolver DISABLED")
endif()
