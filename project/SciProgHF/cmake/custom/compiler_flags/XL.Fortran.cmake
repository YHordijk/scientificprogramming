if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qzerosize -qextname -qsuppress=cmpmsg -qtbtable -qport=mod -qflag=w:w -DVAR_XLF")
    # -qport=mod to allow MOD(I,J) with I,J different types, is needed for relccsd (i*4 and i*8 mixed).
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g")
# miro: needed for linking C/C++ objects with XLF compiled sources
    add_definitions(-DVAR_XLF -DVAR_MFDS)
endif()
