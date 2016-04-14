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
#include <ringmesh/geo_model_editor.h>

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

        GEO::Logger::div( "RINGMeshRemoveGME" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshDuplicateFaults !"
            << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "remove_gme" ) ;
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

        std::string gme_type_name = GEO::CmdLine::get_arg( "remove_gme:type" ) ;
        GME::TYPE gme_type = GME::NO_TYPE ;
        for( index_t type_itr = 0; type_itr < GME::NO_TYPE; ++type_itr ) {
            if( GME::type_name( static_cast< GME::TYPE >( type_itr ) )
                == gme_type_name ) {
                gme_type = static_cast< GME::TYPE >( type_itr ) ;
                break ;
            }
        }

        if( gme_type >= GME::NO_TYPE ) {
            throw RINGMeshException( "I/O",
                "Unknown GeoModelElement type " + gme_type_name ) ;
        }

        index_t gme_id = GEO::CmdLine::get_arg_uint( "remove_gme:id" ) ;
        if( gme_id >= geomodel.nb_elements( gme_type ) ) {
            throw RINGMeshException( "I/O",
                "Index superior to number of elements of type " + gme_type_name ) ;
        }

        std::set< GME::gme_t > to_delete ;
        to_delete.insert( GME::gme_t( gme_type, gme_id ) ) ;

        GeoModelEditor editor( geomodel ) ;
        editor.remove_elements_and_dependencies( to_delete ) ;

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
