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

target_sources(${target_name}
    PRIVATE
        "${lib_source_dir}/common.cpp"
        "${lib_source_dir}/geogram_extension.cpp"
        "${lib_source_dir}/geogram_mesh.cpp"
    PRIVATE # Could be PUBLIC from CMake 3.3
        "${lib_include_dir}/common.h"
        "${lib_include_dir}/geogram_extension.h"
        "${lib_include_dir}/geogram_mesh.h"
        "${lib_include_dir}/geogram_mesh_builder.h"
)

set( GEOGRAM_EXTENSION_RINGMESH_LIB_INTERNAL_DEPENDENCIES
    basic
    mesh
)
#set(MANDATORY_RINGMESH_LIBRARIES ${MANDATORY_RINGMESH_LIBRARIES} ${GEOGRAM_EXTENSION_RINGMESH_LIB_INTERNAL_DEPENDENCIES} PARENT_SCOPE)
foreach(ringmesh_lib ${GEOGRAM_EXTENSION_RINGMESH_LIB_INTERNAL_DEPENDENCIES})
    if(NOT ${ringmesh_lib} IN_LIST ${MANDATORY_RINGMESH_LIBRARIES})
        set(MANDATORY_RINGMESH_LIBRARIES ${MANDATORY_RINGMESH_LIBRARIES} ${ringmesh_lib} PARENT_SCOPE)
    endif()
endforeach()
set(WANTED_RINGMESH_LIBRARIES ${WANTED_RINGMESH_LIBRARIES} ${target_name} PARENT_SCOPE)

target_link_libraries(${target_name} 
    PUBLIC 
        ${GEOGRAM_EXTENSION_RINGMESH_LIB_INTERNAL_DEPENDENCIES}
)
