# Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
# Applications (ASGA). All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of ASGA nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#     http://www.ring-team.org
#
#     RING Project
#     Ecole Nationale Superieure de Geologie - GeoRessources
#     2 Rue du Doyen Marcel Roubault - TSA 70605
#     54518 VANDOEUVRE-LES-NANCY
#     FRANCE

set(TEST_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set(TEST_DATA_INPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/data")
set(TEST_DATA_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/output")
set(CTEST_OUTPUT_ON_FAILURE ON)

# Create the directory for the outputs of the tests
# It is output in the configuration file as well as the test input data directory
file(MAKE_DIRECTORY ${TEST_DATA_OUTPUT_DIRECTORY})

# Configure a header file to pass some of the CMake settings
# to the source code
set(test_config_file_in ${TEST_DIRECTORY}/ringmesh_tests_config.h.in)
set(test_config_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_tests_config.h)
configure_file(${test_config_file_in} ${test_config_file})

add_subdirectory(basic)
add_subdirectory(geogram_extension)
add_subdirectory(geomodel/builder)
add_subdirectory(geomodel/core)
add_subdirectory(geomodel/tools)
add_subdirectory(io)
add_subdirectory(mesh)
add_subdirectory(tetrahedralize)
