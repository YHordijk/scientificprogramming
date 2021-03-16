option(ENABLE_OLD_LINKER "Enable static linking with older linkers" OFF)

if(ENABLE_STATIC_LINKING)
    if(ENABLE_OLD_LINKER)
       set(_no_dyn_export "")
    else()
       # Miro: the "-Wl,--no-export-dynamic" flag is only for newer
       # ld-linkers, fix for older linkers
       set(_no_dyn_export "-Wl,--no-export-dynamic")
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Bstatic ${_no_dyn_export}")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static ${_no_dyn_export}")
    endif()
    unset(_no_dyn_export)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
endif()
