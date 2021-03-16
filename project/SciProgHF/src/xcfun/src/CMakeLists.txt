
include_directories(taylor functionals)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/include)

set (BASE_LIB_SRCS 
	xcfun.cpp
	evaluators.cpp
	functional.cpp
	parameters.cpp
	potential.cpp
	fortran.c
	)

set(FUNCTIONAL_SRCS
    functionals/beckex.cpp
    functionals/brx.cpp
    functionals/ktx.cpp
#    functionals/lb94.cpp
    functionals/ldaerfc.cpp
    functionals/ldaerfc_jt.cpp
    functionals/ldaerfx.cpp
    functionals/list_of_functionals.cpp
    functionals/lypc.cpp
    functionals/m05c.cpp
    functionals/m05x.cpp
    functionals/m05x2c.cpp
    functionals/m05x2x.cpp
    functionals/m06c.cpp
    functionals/m06hfx.cpp
    functionals/m06lc.cpp
    functionals/m06lx.cpp
    functionals/m06x.cpp
    functionals/m06x2c.cpp
    functionals/m06x2x.cpp
    functionals/optx.cpp
    functionals/p86c.cpp
    functionals/pbec.cpp
    functionals/pbex.cpp
    functionals/pw86x.cpp
    functionals/pw91k.cpp
    functionals/pw91x.cpp
    functionals/pw92c.cpp
    functionals/revpbex.cpp
    functionals/revtpssc.cpp
    functionals/revtpssx.cpp
    functionals/rpbex.cpp
    functionals/slaterx.cpp
    functionals/spbec.cpp
    functionals/tfk.cpp
    functionals/tpssc.cpp
    functionals/tpssx.cpp
    functionals/vwn5c.cpp
    )

#miro: these variables are passed from DIRAC CMake to xcfun submodule
if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL AND ENABLE_64BIT_INTEGERS)
    SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libxcfun.a submodule library, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

add_library(xcfun
    ${BASE_LIB_SRCS}
    ${FUNCTIONAL_SRCS}
    )

add_custom_command (
	OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/include
		#${CMAKE_CURRENT_BINARY_DIR}/fortran # cmake bug
	COMMAND ${CMAKE_COMMAND} -E make_directory 
		${CMAKE_CURRENT_BINARY_DIR}/include
	COMMAND ${CMAKE_COMMAND} -E make_directory 
		${CMAKE_CURRENT_BINARY_DIR}/fortran
	)

install(TARGETS xcfun ARCHIVE DESTINATION lib)
