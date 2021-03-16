set(EXTERNAL_LIBS
    ${EXTERNAL_LIBS}
    ${MATH_LIBS}
    )

if(NOT EXPLICIT_LIBS STREQUAL "off")
    set(EXTERNAL_LIBS
        ${EXTERNAL_LIBS}
        ${EXPLICIT_LIBS}
        )
endif()
