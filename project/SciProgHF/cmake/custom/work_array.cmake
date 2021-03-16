set(WORK_MEM_WORDS "64000000" CACHE STRING "Work memory in words")

add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})
