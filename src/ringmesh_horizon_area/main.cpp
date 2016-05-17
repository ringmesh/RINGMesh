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
#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/io.h>

/*!
 * @author Benjamin Chauvin
 */

/// @todo for now it is for horizons but it could be easily generalized for all
/// the kinds of interfaces.
int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshHorizonArea" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshHorizonArea !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;

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

        const char separator = ',' ;
        const char new_line = '\n' ;

        std::ofstream out( "horizon_areas.csv" ) ;
        out.precision( 16 ) ;
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: horizon_areas.csv" ) ;
        }

        out << "Horizon name" ;
        out << separator ;
        out << "Horizon id or surface id" ;
        out << separator ;
        out << "area (m2)" ;
        out << new_line ;
        out << new_line ;

        for( index_t interface_itr = 0; interface_itr < geomodel.nb_interfaces();
            ++interface_itr ) {
            if( !GME::is_stratigraphic_limit(
                geomodel.one_interface( interface_itr ).geological_feature() ) ) {
                continue ;
            }

            const GME& cur_horizon = geomodel.one_interface( interface_itr ) ;
            out << cur_horizon.name() ;
            out << separator ;
            out << cur_horizon.index() ;
            out << new_line ;

            double total_interface_area = 0. ;
            for( index_t child_itr = 0; child_itr < cur_horizon.nb_children();
                ++child_itr ) {
                const Surface& cur_surface = geomodel.surface(
                    cur_horizon.child( child_itr ).index() ) ;
                double cur_area = GEO::Geom::mesh_area( cur_surface.mesh() ) ;
                total_interface_area += cur_area ;

                out << separator ; // first column is left empty
                out << cur_surface.index() ;
                out << separator ;
                out << cur_area ;
                out << new_line ;
            }
            out << separator ;
            out << "all" ;
            out << separator ;
            out << total_interface_area ;
            out << new_line ;
            out << new_line ;
        }

    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
