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

# CMakeList for RINGMesh tutorials

set(TUTORIALS_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/tutorials)
set(TUTORIALS_SOURCE_DIRECTORY ${TUTORIALS_DIRECTORY}/src)
set(TUTORIALS_DATA_INPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/data)
set(TUTORIALS_DATA_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tutorials/output)

# Create the directory for the outputs of the tests
# It is output in the configuration file as well as the input data for tutorials
file(MAKE_DIRECTORY ${TUTORIALS_DATA_OUTPUT_DIRECTORY})

# Configure a header file to pass some of the CMake settings
# to the source code
set(tutorial_config_file_in ${TUTORIALS_SOURCE_DIRECTORY}/ringmesh_tutorials_config.h.in)
set(tutorial_config_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_tutorials_config.h)
configure_file(${tutorial_config_file_in} ${tutorial_config_file})

# Get the name of the tutorial files
file(GLOB TUTORIALS "${TUTORIALS_SOURCE_DIRECTORY}/*.cpp")
# ringmesh_files is defined in the root RINGMesh CMakeLists.txt.
# This line is for clang utilities.
set(ringmesh_files ${ringmesh_files} ${TUTORIALS} PARENT_SCOPE)

# For all the test files found
foreach(tutorialsourcefile ${TUTORIALS})
    add_ringmesh_tutorial(${tutorialsourcefile} geomodel_tools io)
endforeach(tutorialsourcefile ${TUTORIALS})


