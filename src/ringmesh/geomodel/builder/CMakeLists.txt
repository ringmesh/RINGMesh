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

target_sources(${target_name}
    PRIVATE
        "${lib_source_dir}/common.cpp"
        "${lib_source_dir}/geomodel_builder_2d_from_3d.cpp"
        "${lib_source_dir}/geomodel_builder_access.cpp"
        "${lib_source_dir}/geomodel_builder_from_mesh.cpp"
        "${lib_source_dir}/geomodel_builder_geology.cpp"
        "${lib_source_dir}/geomodel_builder_geometry.cpp"
        "${lib_source_dir}/geomodel_builder_remove.cpp"
        "${lib_source_dir}/geomodel_builder_topology.cpp"
        "${lib_source_dir}/geomodel_builder.cpp"
        "${lib_source_dir}/stratigraphic_column_builder.cpp"
    PRIVATE # Could be PUBLIC from CMake 3.3
        "${lib_include_dir}/common.h"
        "${lib_include_dir}/geomodel_builder_2d_from_3d.h"
        "${lib_include_dir}/geomodel_builder_access.h"
        "${lib_include_dir}/geomodel_builder_from_mesh.h"
        "${lib_include_dir}/geomodel_builder_geology.h"
        "${lib_include_dir}/geomodel_builder_geometry.h"
        "${lib_include_dir}/geomodel_builder_remove.h"
        "${lib_include_dir}/geomodel_builder_topology.h"
        "${lib_include_dir}/geomodel_builder.h"
        "${lib_include_dir}/stratigraphic_column_builder.h"
)

target_link_libraries(${target_name} 
    PUBLIC 
        basic 
        mesh 
        geomodel_core 
    PRIVATE 
        tinyxml2
)
