# Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#  Contacts:
#     Arnaud.Botella@univ-lorraine.fr
#     Antoine.Mazuyer@univ-lorraine.fr
#     Jeanne.Pellerin@wias-berlin.de
# 
#     http://www.ring-team.org
# 
#     RING Project
#     Ecole Nationale Superieure de Geologie - Georessources
#     2 Rue du Doyen Marcel Roubault - TSA 70605
#     54518 VANDOEUVRE-LES-NANCY
#     FRANCE

# CMakeLists for RINGMesh tests
# This file is included by the root CMakeList of the RINGMesh project 

set(TEST_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/src)

# Create the directory for the outputs of the tests
# It is output in the configuration file as well as the test input data directory 
set(test_output_directory ${PROJECT_BINARY_DIR}/tests/output)
file(MAKE_DIRECTORY ${test_output_directory}) 

# Configure a header file to pass some of the CMake settings
# to the source code
set(test_config_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_tests_config.h)
configure_file(
  "${TEST_DIRECTORY}/ringmesh_tests_config.h.in"
   ${test_config_file}
)

# The libraries with which the test are supposed to link
set(RINGMESH_LIBS geogram ringmesh)

# Get the name of the test files 
file(GLOB TESTS RELATIVE ${TEST_DIRECTORY} "${TEST_DIRECTORY}/*.cpp")

# For all the test files found
foreach(testsourcefile ${TESTS})
    # Get the name of the file without extension
    get_filename_component(testname ${testsourcefile} NAME_WE)
    
    # Add the target for that file
    add_executable(${testname} ${TEST_DIRECTORY}/${testsourcefile})
    
    # Link RINGMesh Make sure YourLib is linked to each app
    target_link_libraries(${testname} geogram RINGMesh)
    
    add_dependencies(${testname} RINGMesh)
    
    # Add the project to a folder of projects for the tests 
    set_property(TARGET ${testname} PROPERTY FOLDER "Tests")
endforeach(testsourcefile ${TESTS})