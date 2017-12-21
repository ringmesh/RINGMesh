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
        "${lib_source_dir}/tetgen_mesher.cpp"
        "${lib_source_dir}/tetra_gen.cpp"
    PRIVATE # Could be PUBLIC from CMake 3.3
        "${lib_include_dir}/common.h"
        "${lib_include_dir}/tetgen_mesher.h"
        "${lib_include_dir}/tetra_gen.h"
)

target_link_libraries(${target_name} PUBLIC geomodel_builder)

if(MG_TETRA)
    target_include_directories(${target_name} PRIVATE ${MG_TETRA}/include)

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(PLATFORM _64)
    else(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(PLATFORM)
    endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    
    if(WIN32)
        link_directories(${MG_TETRA}/lib/WIN7${PLATFORM}_VC10)
    else(WIN32)
        link_directories(${MG_TETRA}/lib/Linux${PLATFORM})
    endif(WIN32)
    
    target_link_libraries(${target_name} LINK_PRIVATE meshgems mg-tetra meshgems_stubs)
endif()

