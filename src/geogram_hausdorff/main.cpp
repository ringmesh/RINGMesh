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

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_distance.h>
#include <geogram/mesh/mesh_distance.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/command_line_args.h>

/*!
 * @author Benjamin Chauvin
 */

int main( int argc, char** argv )
{
    try {

        GEO::initialize() ;
        // From RINGMesh::configure_geogram
        GEO::CmdLine::import_arg_group( "sys" ) ;
        GEO::CmdLine::set_arg( "sys:assert", "abort" ) ;
        GEO::CmdLine::set_arg( "sys:FPE", true ) ;
        GEO::CmdLine::import_arg_group( "algo" ) ;
        GEO::CmdLine::set_arg( "algo:predicates", "exact" ) ;
        GEO::CmdLine::import_arg_group( "log" ) ;
        GEO::CmdLine::set_arg( "sys:use_doubles", true ) ;

        if( argc != 4 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        std::string first_mesh_path = argv[1] ;
        GEO::Mesh first_mesh ;
        GEO::mesh_load( first_mesh_path, first_mesh ) ;

        std::string second_mesh_path = argv[2] ;
        GEO::Mesh second_mesh ;
        GEO::mesh_load( second_mesh_path, second_mesh ) ;

        double sampling_distance = GEO::String::to_double( argv[3] ) ;

        GEO::Logger::div( "Distance between 2 meshes" ) ;
        double one_way = GEO::mesh_one_sided_Hausdorff_distance( first_mesh,
            second_mesh, sampling_distance ) ;
        GEO::Logger::out( "Hausdorff one way" ) << one_way << std::endl ;
        double other_way = GEO::mesh_one_sided_Hausdorff_distance( second_mesh,
            first_mesh, sampling_distance ) ;
        GEO::Logger::out( "Hausdorff other way" ) << other_way << std::endl ;
        GEO::Logger::out( "Hausdorff sym = max" )
            << GEO::geo_max( one_way, other_way ) << std::endl ;

    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
