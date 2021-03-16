# - This module looks for Sphinx
# Find the Sphinx documentation generator
#
# This modules defines
#  SPHINX_EXECUTABLE
#  SPHINX_FOUND

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
# Copyright 2009-2011 Peter Colberg
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file COPYING-CMAKE-SCRIPTS for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_program(SPHINX_EXECUTABLE NAMES sphinx-build
  HINTS
  $ENV{SPHINX_DIR}
  PATH_SUFFIXES bin
  DOC "Sphinx documentation generator"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx DEFAULT_MSG
  SPHINX_EXECUTABLE
)

mark_as_advanced(
  SPHINX_EXECUTABLE
)

find_package(Sphinx QUIET)

if(SPHINX_FOUND)
    add_custom_target(
        html
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/doc/preprocess_sphinx ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR} html
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_custom_target(
        latex
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/doc/preprocess_sphinx ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR} latex
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_custom_target(
        slides
        COMMAND sphinx-build -b slides ${PROJECT_SOURCE_DIR}/doc/workshop-slides ${PROJECT_BINARY_DIR}/html_slides
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
else()
    add_custom_target(
        html
        COMMAND echo error: please install Sphinx, see http://diracprogram.org/doc/master/programmers/documentation_howto.html
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
endif()
