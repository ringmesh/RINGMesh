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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_repair.h>
#include <ringmesh/io/io.h>

/*!
 * @author Benjamin Chauvin
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        Logger::div( "RINGMeshRepair" ) ;
        Logger::out( "" ) << "Welcome to RINGMeshRepair !" << std::endl ;
        Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        Logger::out( "" ) << "Benjamin Chauvin <Benjamin.chauvin@univ-lorraine.fr> "
            << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "repair" ) ;

        if( argc == 1 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        std::vector< std::string > filenames ;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1 ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        std::string model_name = GEO::CmdLine::get_arg( "in:geomodel" ) ;
        if( model_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in in:geomodel" ) ;
        }
        GeoModel geomodel ;
        geomodel_load( geomodel, model_name ) ;

        GeoModelRepair geo_model_repair( geomodel ) ;
        if( GEO::CmdLine::get_arg_bool( "repair:basic" ) ) {
            geo_model_repair.geo_model_mesh_repair() ;
        }

        if( GEO::CmdLine::get_arg_bool( "repair:line_boundary" ) ) {
            geo_model_repair.repair_line_boundary_vertex_order() ;
        }

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
