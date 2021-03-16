#
# Miro: We must list here all external arguments passed from the principal DIRAC CMake module !
#
#  CMAKE_<Lang>_CIOMPILER_ID
#  CMAKE_<Lang>_Flags
#  ENABLE_64BIT_INTEGERS
#  .
#  .
#  .
#
#

cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        )
endif()

project(XCFun CXX C Fortran)

set (XCFun_VERSION_MAJOR 1)
set (XCFun_VERSION_MINOR 0)

option(XCFUN_NO_STDC++ "Disable libsrdc++" OFF)
option(XCFUN_ENABLE_TESTS "Enable tests" ON)

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
       if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL AND ENABLE_64BIT_INTEGERS)
          set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -fno-rtti -fno-exceptions -maix64 ")
          message(STATUS "Added g++ -maix64 flag for xcfun due to IBM AIX XL Fortran and integer*8")
        else()
	  set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -fno-rtti -fno-exceptions ")
        endif()
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-g -O3 -DNDEBUG -Wno-unused")

	if (XCFUN_NO_STDC++)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-threadsafe-statics")
	endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    # Intel has special flags on Windows
    # ihrasko: here I tried to translate all flags from Linux to Windows flags
    if (CMAKE_SYSTEM_NAME MATCHES "Windows")
        set (CMAKE_CXX_FLAGS "/Qwd981 /Qwd279 /Qwd383 /Qwd1572 /Qwd177")
        set (CMAKE_CXX_FLAGS_DEBUG "/Z7 /Od")
        set (CMAKE_CXX_FLAGS_RELEASE "/Z7 /O2 /DNDEBUG")
        set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS}")
    else ()
        set (CMAKE_CXX_FLAGS "-wd981 -wd279 -wd383 -wd1572 -wd177 -fno-rtti -fno-exceptions")
        set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
        # radovan: do not use -fast with Intel 11 (possibly 12, to be checked)
        # -fast causes a miscompile with intel 11
        set (CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG")
        set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
    endif ()
endif ()

if (CMAKE_C_COMPILER_ID MATCHES GNU)
      if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL AND ENABLE_64BIT_INTEGERS)
          set (CMAKE_C_FLAGS "-Wall -maix64")
          message(STATUS "Added gcc -maix64 flag for xcfun due to IBM AIX XL Fortran and integer*8")
       else()
      set (CMAKE_C_FLAGS "-Wall")
        endif()
      set (CMAKE_C_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
      set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG -Wno-unused")
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
    # Intel has special flags on Windows
    # ihrasko: here I tried to translate all flags from Linux to Windows flags
    if (CMAKE_SYSTEM_NAME MATCHES "Windows")
        set (CMAKE_C_FLAGS "/Qwd981 /Qwd279 /Qwd383 /Qwd1572 /Qwd177")
        set (CMAKE_C_FLAGS_DEBUG "/Z7 /Od")
        set (CMAKE_C_FLAGS_RELEASE "/Z7 /O2 /DNDEBUG")
        set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS}")
    else ()
        set (CMAKE_C_FLAGS "-wd981 -wd279 -wd383 -wd1572 -wd177")
        set (CMAKE_C_FLAGS_DEBUG "-g -O0")
        set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG")
        set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
    endif ()
endif ()

if (CMAKE_COMPILER_IS_GNUFORTRAN)
	set (CMAKE_Fortran_FLAGS "-Wall -Jfortran")
	set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	set (CMAKE_Fortran_FLAGS_RELEASE "-g -O2")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    # Intel has special flags on Windows
    # ihrasko: here I tried to translate all flags from Linux to Windows flags
    if (CMAKE_SYSTEM_NAME MATCHES "Windows")
        set (CMAKE_Fortran_FLAGS "")
        set (CMAKE_Fortran_FLAGS_DEBUG "/Z7 /Od")
        set (CMAKE_Fortran_FLAGS_RELEASE "/Z7 /O2")
        set (CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}")
        if(ENABLE_64BIT_INTEGERS)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} /4I8"
                )
        endif()
    else ()
        set (CMAKE_Fortran_FLAGS "")
        set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
        set (CMAKE_Fortran_FLAGS_RELEASE "-g -O2")
        set (CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -shared-intel")
        if(ENABLE_64BIT_INTEGERS)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -i8"
                )
        endif()
    endif ()
endif ()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    # remove -rdynamic flag offensive for PGI Fortran in static linking
    if (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
        list(REMOVE_ITEM CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-rdynamic")
        message(STATUS "xcfun: removed -rdynamic flag due to pgf90")
    endif()
endif()

include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(src)
if (XCFUN_ENABLE_TESTS)
	enable_testing()
	add_subdirectory(fortran)
	add_subdirectory(test)
endif()

#Miro: Doxygen searching is in main CMake-in-DIRAC, prevent thus configure error
# because ${DOXYGEN_EXECUTABLE} tends to be filled with wrong executable!
set(DOXYGEN_EXECUTABLE "")
find_package(Doxygen QUIET)
message("xcfun DOXYGEN_FOUND ?  ${DOXYGEN_FOUND}")

if (DOXYGEN_FOUND)
	add_custom_target(docs
		COMMAND ${DOXYGEN_EXECUTABLE}
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
	)

find_program(VALGRIND_EXECUTABLE valgrind)
if (VALGRIND_EXECUTABLE)
	set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
	set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full")
	#set(MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_BINARY_DIR}/valgrind.supp)
endif()

include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "XCFun")
set(CPACK_PACKAGE_VENDOR "Ulf Ekstrom")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
	"XCFun ${MRCPP_VERSION_MAJOR}.${MRCPP_VERSION_MINOR}")
include(CPack)
