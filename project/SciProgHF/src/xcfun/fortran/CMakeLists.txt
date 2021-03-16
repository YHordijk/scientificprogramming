
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(example_fortran 
	example.F90
	xcfun_module.F90
	xcfun_autogen.F90
	)

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#   miro: necessary for PGI
    target_link_libraries(example_fortran stdc++)
    message(STATUS "added libstdc++ to linker flags of xcfun due to PGI ")
endif()

set_target_properties(example_fortran 
	PROPERTIES LINKER_LANGUAGE Fortran
	)

target_link_libraries(example_fortran xcfun m)

# workaround for strange problem with cilkrts lib
# intel does not ship the static version of cilkrts
# and for some reason it gets detected and included in these variables
# breaking the static linking buildup with intel compilers
foreach(_lang C CXX)
    if(CMAKE_${_lang}_IMPLICIT_LINK_LIBRARIES)
        list(REMOVE_ITEM CMAKE_${_lang}_IMPLICIT_LINK_LIBRARIES "cilkrts")
    endif()
endforeach()

