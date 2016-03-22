/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/common.h>

#include <geogram/mesh/mesh_io.h>

/*!
 * @file ringmesh_surface_convert/main.cpp
 * @brief Executable to High level functions on GeoModel
 * @author Gautier Laurent
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;


    try {
        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        std::string ext_mesh = argv[1] ;
        std::string ext_geogram = "test_geogram.geogram" ;

        GEO::Mesh M ;
        GEO::mesh_load( ext_mesh, M ) ;
//        M.vertices.create_vertices(3) ;
//        M.vertices.point( 0 ) = vec3( 0, 0, 0 ) ;
//        M.vertices.point( 1 ) = vec3( 0, 1, 0 ) ;
//        M.vertices.point( 2 ) = vec3( 0, 0, 1 ) ;
//        M.facets.create_triangle( 0, 1, 2 ) ;
        M.show_stats() ;

        GEO::mesh_save( M, ext_geogram ) ;
        GEO::Mesh M2 ;
        GEO::mesh_load( ext_geogram, M2 ) ;
        M2.show_stats() ;
        
    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
