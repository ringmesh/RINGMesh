# Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

# additional target to perform doxygen run, requires doxygen
set(DOXYGEN_MINIMUM_VERSION 1.6.0)
find_package(Doxygen ${DOXYGEN_MINIMUM_VERSION} QUIET)
if(NOT DOXYGEN_FOUND)
    message(STATUS "Doxygen >= ${DOXYGEN_MINIMUM_VERSION} not found, cannot generate documentation")
else()
    message(STATUS "Configuring RINGMesh with Doxygen")
    # Generate a target for generating a specific documentation type
    function(add_doc_target doc_type)
        set(doc_output_dir ${CMAKE_CURRENT_BINARY_DIR})
        
        configure_file(${PROJECT_SOURCE_DIR}/doc/${doc_type}.dox.in ${doc_type}.dox)
        
        set(doc_target doc-${doc_type})
        
        add_custom_target(
            ${doc_target}
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${doc_output_dir}/${doc_type}
            COMMAND ${DOXYGEN_EXECUTABLE} ${doc_type}.dox # > ${doc_type}.log 2>&1
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generate ${doc_type} documentation for ${PROJECT_NAME}"
        )
        
        set_property(
            DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${doc_output_dir}/${doc_type}
        )
    endfunction()

    # Add a custom documentation targets
    add_doc_target(devkit)
endif()

# additional target to perform clang-format run, requires clang-format
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
    message(STATUS "Configuring RINGMesh with clang-format")
    file(GLOB_RECURSE ringmesh_files
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*/*.h
        ${PROJECT_SOURCE_DIR}/tests/*/*.cpp
        ${PROJECT_SOURCE_DIR}/doc/tutorials/src/*.cpp
    )
    add_custom_target(
        format
        COMMAND ${CLANG_FORMAT}
        -style=file
        -i
        ${ringmesh_files}
    )
endif()
