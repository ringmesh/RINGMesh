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

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/io.h>

/*!
 * @author Benjamin Chauvin
 */

namespace {
    using namespace RINGMesh ;

}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshDistance" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshDistance !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "distance" ) ;

        if( argc == 1 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        std::vector< std::string > filenames ;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1 ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        std::string input_geomodel_name = GEO::CmdLine::get_arg( "in:geomodel" ) ;
        if( input_geomodel_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in in:geomodel" ) ;
        }
        GeoModel geomodel ;
        geomodel_load( geomodel, input_geomodel_name ) ;
        geomodel.mesh.facets.test_and_initialize() ;
        Mesh copy_mesh(geomodel,3,false) ;
        geomodel.mesh.copy_mesh( copy_mesh ) ;
        GEO::Mesh copy_geomodel ;
        copy_geomodel.copy(copy_mesh.gfx_mesh()); /// @todo names are not teally good.

        std::string geomodel_to_compare_name = GEO::CmdLine::get_arg(
            "distance:geomodel" ) ;
        if( geomodel_to_compare_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in distance:geomodel" ) ;
        }
        GeoModel geomodel_to_compare ;
        geomodel_load( geomodel_to_compare, geomodel_to_compare_name ) ;
        geomodel_to_compare.mesh.facets.test_and_initialize() ;

        Mesh copy_mesh_to_compare(geomodel_to_compare,3,false) ;
        geomodel_to_compare.mesh.copy_mesh( copy_mesh_to_compare ) ;
        GEO::Mesh copy_geomodel_to_compare ;
        copy_geomodel_to_compare.copy(copy_mesh_to_compare.gfx_mesh());

        double sampling_distance = GEO::CmdLine::get_arg_double(
            "distance:sampling_distance" ) ;
        if( sampling_distance <= epsilon ) {
            throw RINGMeshException( "I/O",
                "Sampling distance cannot be negative." ) ;
        }

        /// @todo for now it is just the Hausdorff distance but we can imagine
        /// that several algo of distances may be used to compute the distance
        /// between 2 geomodels.
        GEO::Logger::div( "Distance between 2 geomodels" ) ;
        double one_way = GEO::mesh_one_sided_Hausdorff_distance( copy_geomodel,
            copy_geomodel_to_compare, sampling_distance ) ;
        GEO::Logger::out( "Hausdorff one way" ) << one_way << std::endl ;
        double other_way = GEO::mesh_one_sided_Hausdorff_distance(
            copy_geomodel_to_compare, copy_geomodel, sampling_distance ) ;
        GEO::Logger::out( "Hausdorff other way" ) << other_way << std::endl ;
        GEO::Logger::out( "Hausdorff sym = max" )
            << GEO::geo_max( one_way, other_way ) << std::endl ;
        ringmesh_assert( std::abs(GEO::geo_max( one_way, other_way ) -
                GEO::mesh_symmetric_Hausdorff_distance( copy_geomodel,
                    copy_geomodel_to_compare, sampling_distance ) ) < epsilon ) ;

    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
