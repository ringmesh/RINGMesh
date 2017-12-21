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

cmake_minimum_required(VERSION 3.1)

find_program(GCOV_EXECUTABLE gcov)
if(NOT GCOV_EXECUTABLE)
	message(FATAL_ERROR "gcov not found! Aborting...")
endif()

set(GCOV_OUTPUT ${PROJECT_ROOT}/gcov-outputs)
file(MAKE_DIRECTORY ${GCOV_OUTPUT})

# Get the coverage data.
file(GLOB_RECURSE GCDA_FILES "${COV_PATH}/*.gcda")
message("GCDA files:")

# Get a list of all the object directories needed by gcov
# (The directories the .gcda files and .o files are found in)
# and run gcov on those.
foreach(GCDA ${GCDA_FILES})
	message("Process: ${GCDA}")
	message("------------------------------------------------------------------------------")
	get_filename_component(GCDA_DIR ${GCDA} PATH)
	execute_process(
		COMMAND ${GCOV_EXECUTABLE} -p -o ${GCDA_DIR} ${GCDA}
		WORKING_DIRECTORY ${GCOV_OUTPUT}
	)
endforeach()
