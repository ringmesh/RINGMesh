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

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/io.h>
#include <ringmesh/geo_model_api.h>

/*!
 * @author Benjamin Chauvin
 */

namespace {
    using namespace RINGMesh ;

    vec3 extract_coords_from_string( const std::string& coords_in_string )
    {
        std::vector< std::string > split_coords ;
        split_coords.reserve( 3 ) ;
        GEO::String::split_string( coords_in_string, ' ', split_coords, true ) ;
        if( split_coords.size() != 3 ) {
            throw RINGMeshException( "I/O",
                "Vector" + coords_in_string + "has not exactly 3 components" ) ;
        }
        vec3 coords_vec ;
        for( index_t split_coords_itr = 0; split_coords_itr < 3;
            ++split_coords_itr ) {
            coords_vec[split_coords_itr] = GEO::String::to_double(
                split_coords[split_coords_itr] ) ;
        }
        return coords_vec ;
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshRotateGeoModel" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshRotateGeoModel !"
            << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "rotation" ) ;
        CmdLine::import_arg_group( "out" ) ;

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

        std::string rotation_origin_string = GEO::CmdLine::get_arg(
            "rotation:origin" ) ;
        vec3 rotation_origin_vec = extract_coords_from_string(
            rotation_origin_string ) ;

        std::string rotation_axis_string = GEO::CmdLine::get_arg( "rotation:axis" ) ;
        vec3 rotation_axis_vec = extract_coords_from_string( rotation_axis_string ) ;

        double rotation_angle = GEO::CmdLine::get_arg_double( "rotation:angle" ) ;
        std::string rotation_unit = GEO::CmdLine::get_arg( "rotation:unit" ) ;
        bool is_deg ;
        if( rotation_unit == "deg" ) {
            is_deg = true ;
        } else if( rotation_unit == "rad" ) {
            is_deg = false ;
        } else {
            throw RINGMeshException( "I/O", "Unknown angle unit " + rotation_unit ) ;
        }

        rotate( geomodel, rotation_origin_vec, rotation_axis_vec, rotation_angle,
            is_deg ) ;

        std::string output_geomodel_name = GEO::CmdLine::get_arg( "out:geomodel" ) ;
        if( output_geomodel_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in out:geomodel" ) ;
        }
        geomodel_save( geomodel, output_geomodel_name ) ;

    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
