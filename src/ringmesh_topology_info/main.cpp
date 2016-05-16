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

/// @todo NOT FINIESHED AT ALL, FOR NOW IT IS FOR GETTING CONTACT AND CORNER from
/// the intersection of two interfaces.
namespace {
    using namespace RINGMesh ;

    index_t get_contact_id_from_two_interfaces(
        const GeoModel& geomodel,
        index_t interface_one_id,
        index_t interface_two_id )
    {
        index_t contact_id = NO_ID ;
        for( index_t contact_itr = 0; contact_itr < geomodel.nb_contacts();
            ++contact_itr ) {
            const GME& cur_contact = geomodel.contact( contact_itr ) ;
            const GME& first_line = cur_contact.child( 0 ) ;
            bool interface_one_found = false ;
            bool interface_two_found = false ;
            for( index_t surf_in_boun_itr = 0;
                surf_in_boun_itr < first_line.nb_in_boundary();
                ++surf_in_boun_itr ) {
                const GME& cur_surf = first_line.in_boundary( surf_in_boun_itr ) ;
                if( cur_surf.parent().index() == interface_one_id ) {
                    interface_one_found = true ;
                } else if( cur_surf.parent().index() == interface_two_id ) {
                    interface_two_found = true ;
                }

                if( interface_one_found && interface_two_found ) {
                    contact_id = contact_itr ;
                    break ;
                }
            }
            if( contact_id != NO_ID ) {
                break ;
            }
        }
        return contact_id ;
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshTopologyInfo" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshTopologyInfo !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "topology_info" ) ;

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

        index_t interface_one_id = GEO::CmdLine::get_arg_uint(
            "topology_info:interface_one_id" ) ;
        if( interface_one_id >= geomodel.nb_interfaces() ) {
            throw RINGMeshException( "I/O",
                "Index superior to number of interfaces " ) ;
        }

        GEO::Logger::div( "Get contacts with interfaces" ) ;

        for( index_t interface_itr = 0; interface_itr < geomodel.nb_interfaces();
            ++interface_itr ) {
            // For now I just want contact between a fault and the horizons
            // So skip contact with model boundaries.
            if( !GME::is_fault(
                geomodel.one_interface( interface_itr ).geological_feature() ) ) {
                continue ;
            }
            // In my case it should not happen since I just want contact between a faut and the horizons
            if( interface_itr == interface_one_id ) {
                continue ;
            }
            index_t contact_id = get_contact_id_from_two_interfaces( geomodel,
                interface_one_id, interface_itr ) ;
            if( contact_id == NO_ID ) {
                continue ;
            }
            GEO::Logger::out( "Interface" ) << interface_itr << std::endl ;
            GEO::Logger::out( "Interface" )
                << geomodel.one_interface( interface_itr ).name() << std::endl ;
            GEO::Logger::out( "Contact" ) << contact_id << std::endl ;

            const GME& cur_contact = geomodel.contact( contact_id ) ;
            for( index_t line_itr = 0; line_itr < cur_contact.nb_children();
                ++line_itr ) {
                const GME& cur_line = cur_contact.child( line_itr ) ;
                for( index_t corner_itr = 0; corner_itr < cur_line.nb_boundaries();
                    ++corner_itr ) {
                    GEO::Logger::out( "Corner" )
                        << cur_line.boundary( corner_itr ).index() << std::endl ;
                }
            }
        }

        GEO::Logger::div( "Get fault names and ids" ) ;
        for( index_t interface_itr = 0; interface_itr < geomodel.nb_interfaces();
            ++interface_itr ) {
            if( !GME::is_fault(
                geomodel.one_interface( interface_itr ).geological_feature() ) ) {
                continue ;
            }
            GEO::Logger::out( "Interface" ) << interface_itr << std::endl ;
            GEO::Logger::out( "Interface" )
                << geomodel.one_interface( interface_itr ).name() << std::endl ;
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
