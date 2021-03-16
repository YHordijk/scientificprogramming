option(ENABLE_MATLAB_LOG "Log interesting matrices in matlab format (very verbose!)" OFF)

if(ENABLE_MATLAB_LOG)
    add_definitions(-DMOD_MATLAB_LOG)
endif()
