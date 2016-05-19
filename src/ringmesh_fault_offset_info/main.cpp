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

/*!
 * @author Benjamin Chauvin
 */

/// @todo NOT FINIESHED AT ALL, FOR NOW IT IS FOR GETTING CONTACT AND CORNER from
/// the intersection of two interfaces.
namespace {
    using namespace RINGMesh ;

    // Same code than in the main of ringmesh_topology_info
    // It may happen that even in the case of duplicated fault model there are
    // several contacts between a fault and a horizon (case where another fault
    // has a line (and so a contact) in common with the 2 other interfaces).
    // That should not happen too often.
    // return copy of a very small vector...
    std::vector< index_t > get_contact_id_from_two_interfaces(
        const GeoModel& geomodel,
        index_t interface_one_id,
        index_t interface_two_id )
    {
        std::vector< index_t > contacts ;
        contacts.reserve( 3 ) ; // In most cases it should be just 1
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
                    contacts.push_back( contact_itr ) ;
                    break ;
                }
            }
        }
        return contacts ;
    }

    /// @todo that only works with models with duplicated fault network.
    index_t find_other_fault_mirror( const GME& one_fault_side )
    {
        ringmesh_assert( one_fault_side.type() == GME::INTERFACE ) ;
        ringmesh_assert( GME::is_fault(one_fault_side.geological_feature()) ) ;

        std::string cur_name = one_fault_side.name() ;
        std::string other_name = cur_name ;
        if( GEO::String::string_ends_with( cur_name, "_side_minus" ) ) {
            other_name.erase( other_name.end() - 11, other_name.end() ) ;
            other_name += "_side_plus" ;
        } else {
            ringmesh_assert( GEO::String::string_ends_with(cur_name, "_side_plus") ) ;
            other_name.erase( other_name.end() - 10, other_name.end() ) ;
            other_name += "_side_minus" ;
        }

        for( index_t interface_itr = 0;
            interface_itr < one_fault_side.model().nb_interfaces();
            ++interface_itr ) {
            // several asserts may be set
            if( one_fault_side.model().one_interface( interface_itr ).name()
                == other_name ) {
                return interface_itr ;
            }
        }

        return NO_ID ;
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshFaultOffsetInfo" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshFaultOffsetInfo !"
            << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "fault_offset_info" ) ;

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

        std::string axis = GEO::CmdLine::get_arg( "fault_offset_info:axis" ) ;
        if( axis != "x" && axis != "y" ) {
            throw RINGMeshException( "FaultOffset",
                "Axis must be x or y, not " + axis ) ;
        }
        const double value = GEO::CmdLine::get_arg_double(
            "fault_offset_info:value" ) ;
        const double tolerance = GEO::CmdLine::get_arg_double(
            "fault_offset_info:tolerance" ) ;

        GEO::Logger::div( "Get fault offset informations" ) ;
        std::vector< bool > already_treated( geomodel.nb_interfaces(), false ) ;

        const char separator = ',' ;
        const char new_line = '\n' ;

        std::ofstream out( "fault_offset_informations.csv" ) ;
        out.precision( 16 ) ;
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: fault_offset_informations.csv" ) ;
        }

        out << "Contact (interface names and indices)" ;
        out << separator ;
        out << "Upper corner (id)" ;
        out << separator ;
        out << "Lower corner (id)" ;
        out << separator ;
        out << "Throw" ;
        out << separator ;
        out << "Heave" ;
        out << separator ;
        out << "Dip" ;
        out << new_line ;

        for( index_t interface_itr = 0; interface_itr < geomodel.nb_interfaces();
            ++interface_itr ) {

            if( already_treated[interface_itr] ) {
                continue ;
            }
            already_treated[interface_itr] = true ;

            const GME& cur_interface = geomodel.one_interface( interface_itr ) ;
            if( !GME::is_fault( cur_interface.geological_feature() ) ) {
                continue ;
            }

            index_t other_side_id = find_other_fault_mirror( cur_interface ) ;
            if( other_side_id == NO_ID ) {
                continue ;
            }
            ringmesh_assert( !already_treated[other_side_id] ) ;
            already_treated[other_side_id] = true ;

            /// @todo dirty code with interface_itr2... split into small functions
            for( index_t interface_itr2 = 0;
                interface_itr2 < geomodel.nb_interfaces(); ++interface_itr2 ) {
                const GME& cur_interface2 = geomodel.one_interface(
                    interface_itr2 ) ;
                if( !GME::is_stratigraphic_limit(
                    cur_interface2.geological_feature() ) ) {
                    continue ;
                }

                std::vector< index_t > first_contact_ids =
                    get_contact_id_from_two_interfaces( geomodel, interface_itr2,
                        interface_itr ) ;
                if( first_contact_ids.empty() ) {
                    continue ;
                }
                std::vector< index_t > second_contact_ids =
                    get_contact_id_from_two_interfaces( geomodel, interface_itr2,
                        other_side_id ) ;
                if( second_contact_ids.empty() ) {
                    continue ;
                }

                bool found = false ;
                vec3 first_wanted_corner_vec ;
                index_t first_wanted_corner_index = NO_ID ;
                for( index_t first_contact_ids_itr = 0;
                    first_contact_ids_itr < first_contact_ids.size();
                    ++first_contact_ids_itr ) {
                    const GME& cur_first_contact = geomodel.contact(
                        first_contact_ids[first_contact_ids_itr] ) ;
                    for( index_t line_itr = 0;
                        line_itr < cur_first_contact.nb_children(); ++line_itr ) {
                        const GME& cur_line = cur_first_contact.child( line_itr ) ;
                        for( index_t corner_itr = 0;
                            corner_itr < cur_line.nb_boundaries(); ++corner_itr ) {
                            const Corner& cur_corner = geomodel.corner(
                                cur_line.boundary( corner_itr ).index() ) ;
                            double to_check ;
                            if( axis == "x" ) {
                                to_check = cur_corner.vertex( 0 ).x ;
                            } else {
                                ringmesh_assert( axis == "y" ) ;
                                to_check = cur_corner.vertex( 0 ).y ;
                            }

                            if( std::abs( to_check - value ) < tolerance ) {
                                first_wanted_corner_vec = cur_corner.vertex( 0 ) ;
                                first_wanted_corner_index = cur_corner.index() ;
                                found = true ;
                                break ;
                            }
                        }
                        if( found ) {
                            break ;
                        }
                    }
                    if( found ) {
                        break ;
                    }
                }

                if( !found ) {
                    continue ;
                }

                /// @todo dirty copy paste: do a function
                found = false ;
                vec3 second_wanted_corner_vec ;
                index_t second_wanted_corner_index = NO_ID ;
                for( index_t second_contact_ids_itr = 0;
                    second_contact_ids_itr < second_contact_ids.size();
                    ++second_contact_ids_itr ) {
                    const GME& cur_second_contact = geomodel.contact(
                        second_contact_ids[second_contact_ids_itr] ) ;
                    for( index_t line_itr = 0;
                        line_itr < cur_second_contact.nb_children(); ++line_itr ) {
                        const GME& cur_line = cur_second_contact.child( line_itr ) ;
                        for( index_t corner_itr = 0;
                            corner_itr < cur_line.nb_boundaries(); ++corner_itr ) {
                            const Corner& cur_corner = geomodel.corner(
                                cur_line.boundary( corner_itr ).index() ) ;
                            double to_check ;
                            if( axis == "x" ) {
                                to_check = cur_corner.vertex( 0 ).x ;
                            } else {
                                ringmesh_assert( axis == "y" ) ;
                                to_check = cur_corner.vertex( 0 ).y ;
                            }

                            if( std::abs( to_check - value ) < tolerance ) {
                                second_wanted_corner_vec = cur_corner.vertex( 0 ) ;
                                second_wanted_corner_index = cur_corner.index() ;
                                found = true ;
                                break ;
                            }
                        }
                        if( found ) {
                            break ;
                        }
                    }
                    if( found ) {
                        break ;
                    }
                }

                if( !found ) {
                    continue ;
                }
                ringmesh_assert( first_wanted_corner_index != NO_ID ) ;
                ringmesh_assert( second_wanted_corner_index != NO_ID ) ;
                ringmesh_assert( first_wanted_corner_index != second_wanted_corner_index ) ;

                const vec3 dip_vec = first_wanted_corner_vec
                    - second_wanted_corner_vec ;

                // Display the info
                // Interfaces (with indices) responsible of the fault contact
                out
                    << "Between " + cur_interface.name() + " ("
                        + GEO::String::to_string( cur_interface.index() ) + ")" + "/"
                        + geomodel.one_interface( other_side_id ).name() + " ("
                        + GEO::String::to_string( other_side_id ) + ")" + " and "
                        + cur_interface2.name() + " ("
                        + GEO::String::to_string( cur_interface2.index() ) + ")" ;
                out << separator ;
                // Upper corner (index)
                out
                    << (
                        first_wanted_corner_vec.z > second_wanted_corner_vec.z ?
                            first_wanted_corner_index : second_wanted_corner_index ) ;
                out << separator ;
                // Lower corner (index)
                out
                    << (
                        first_wanted_corner_vec.z < second_wanted_corner_vec.z ?
                            first_wanted_corner_index : second_wanted_corner_index ) ;
                out << separator ;
                // Throw
                out << std::abs( dip_vec.z ) ;
                out << separator ;
                // Heave
                out
                    << (
                        ( axis == "x" ) ?
                            std::abs( dip_vec.y ) : std::abs( dip_vec.x ) ) ;
                out << separator ;
                // Dip
                out << dip_vec.length() ;
                out << new_line ;
            }
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
