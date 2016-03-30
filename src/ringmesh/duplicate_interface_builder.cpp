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

#include <ringmesh/duplicate_interface_builder.h>

#include <geogram/mesh/mesh_io.h>
#include <ringmesh/geometry.h>

/*!
 * @file ringmesh/duplicate_interface_builder.cpp
 * @brief Class to duplicate GeoModel Interface to
 * enable sliding along them (faults) and unconformal
 * mesh generation.
 * @author Benjamin Chauvin
 */

namespace RINGMesh {

    DuplicateInterfaceBuilder::DuplicateInterfaceBuilder( GeoModel& model )
        : GeoModelBuilder( model )
    {
    }

    DuplicateInterfaceBuilder::~DuplicateInterfaceBuilder()
    {

    }

    void DuplicateInterfaceBuilder::duplicate_interface(
        index_t interface_id_to_duplicate )
    {
        ringmesh_assert(interface_id_to_duplicate < model_.nb_interfaces()) ;

        const GeoModelElement& interface_to_duplicate = model_.one_interface(
            interface_id_to_duplicate ) ;

        index_t interface_to_duplicate_nb_children =
            interface_to_duplicate.nb_children() ;

        create_elements( GME::SURFACE, interface_to_duplicate_nb_children ) ;

        for( index_t interface_child_itr = 0;
            interface_child_itr < interface_to_duplicate.nb_children();
            ++interface_child_itr ) {
            const GeoModelElement& cur_child = interface_to_duplicate.child(
                interface_child_itr ) ;
            ringmesh_assert( cur_child.type() == GME::SURFACE ) ;
        }
    }

    index_t find_local_boundary_id(
        const GeoModelElement& reg,
        const GeoModelElement& surf )
    {
        ringmesh_assert(reg.type()==GME::REGION) ;
        ringmesh_assert(surf.type()==GME::SURFACE) ;
        for( index_t boundary_i = 0; boundary_i < reg.nb_boundaries();
            ++boundary_i ) {
            if( reg.boundary( boundary_i ).index() == surf.index() ) {
                return boundary_i ;
            }
        }

        return NO_ID ;
    }

    void DuplicateInterfaceBuilder::get_new_surfaces(
        index_t interface_id_to_duplicate )
    {
        ringmesh_assert(interface_id_to_duplicate < model_.nb_interfaces()) ;

        const GeoModelElement& interface_to_duplicate = model_.one_interface(
            interface_id_to_duplicate ) ;

        const index_t interface_to_duplicate_nb_children =
            interface_to_duplicate.nb_children() ;
        ringmesh_assert(interface_to_duplicate_nb_children >= 1) ;

        //////////////////////////////////////////
        // copy paste from void GeoModelEditor::remove_elements( const std::set< gme_t >& elements )
        std::vector< std::vector< index_t > > to_erase_by_type ;
        to_erase_by_type.reserve( GME::NO_TYPE ) ;
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
            to_erase_by_type.push_back(
                std::vector< index_t >(
                    model_.nb_elements( static_cast< GME::TYPE >( i ) ), 0 ) ) ;
        }
        // Delete of the interface (will be replaced by a custom interface with
        // side informations)
        to_erase_by_type[GME::INTERFACE][interface_id_to_duplicate] = NO_ID ;
        //////////////////////////////////////////

        // minus = false, plus = true
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_minus ;
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_plus ;
        GME::gme_t interface_minus_gme_t = create_element( GME::INTERFACE ) ;
        GME::gme_t interface_plus_gme_t = create_element( GME::INTERFACE ) ;
        to_erase_by_type[GME::INTERFACE].push_back( 0 ) ;
        to_erase_by_type[GME::INTERFACE].push_back( 0 ) ;

        // Find for each region, what surfaces are in boundary.
        for( index_t interface_child_itr = 0;
            interface_child_itr < interface_to_duplicate.nb_children();
            ++interface_child_itr ) {
            const GeoModelElement& cur_child = interface_to_duplicate.child(
                interface_child_itr ) ;
            ringmesh_assert( cur_child.type() == GME::SURFACE ) ;
            to_erase_by_type[GME::SURFACE][cur_child.index()] = NO_ID ;

            const index_t nb_in_boundary_cur_child = cur_child.nb_in_boundary() ;
            ringmesh_assert( nb_in_boundary_cur_child == 1 || nb_in_boundary_cur_child == 2 ) ;
            if( nb_in_boundary_cur_child == 2 ) {
                ringmesh_assert( !cur_child.is_on_voi() ) ;
                const GeoModelElement& cur_in_boundary = cur_child.in_boundary( 0 ) ;
                ringmesh_assert( cur_in_boundary.type() == GME::REGION ) ;
                const Region& cur_reg =
                    dynamic_cast< const Region& >( cur_in_boundary ) ;

                const GeoModelElement& cur_in_boundary2 = cur_child.in_boundary(
                    1 ) ;
                ringmesh_assert( cur_in_boundary2.type() == GME::REGION ) ;
                const Region& cur_reg2 =
                    dynamic_cast< const Region& >( cur_in_boundary2 ) ;

                if( cur_in_boundary.index() == cur_in_boundary2.index() ) {
                    // if same region there is 2 in_boundarie even if they are the same
                    // The surface is internal and on the both side there is the
                    // same region. This surface is duplicated.
                    surfaces_boundary_regions_side_plus[cur_in_boundary.index()].push_back(
                        cur_child.index() ) ;
                    surfaces_boundary_regions_side_minus[cur_in_boundary.index()].push_back(
                        cur_child.index() ) ;
                } else {
                    index_t local_boundary_id = find_local_boundary_id(
                        cur_in_boundary, cur_child ) ;

#ifdef RINGMESH_DEBUG
                    index_t local_boundary_id2 = find_local_boundary_id(
                        cur_in_boundary2, cur_child ) ;
                    ringmesh_assert(
                        cur_reg.side( local_boundary_id )
                        != cur_reg2.side( local_boundary_id2 ) ) ;
#endif

                    if( cur_reg.side( local_boundary_id ) ) {
                        surfaces_boundary_regions_side_plus[cur_in_boundary.index()].push_back(
                            cur_child.index() ) ;
                        surfaces_boundary_regions_side_minus[cur_in_boundary2.index()].push_back(
                            cur_child.index() ) ;
                    } else {
                        surfaces_boundary_regions_side_minus[cur_in_boundary.index()].push_back(
                            cur_child.index() ) ;
                        surfaces_boundary_regions_side_plus[cur_in_boundary2.index()].push_back(
                            cur_child.index() ) ;
                    }
                }
            } else {
                ringmesh_assert( nb_in_boundary_cur_child == 1 ) ;
                if( cur_child.is_on_voi() ) {
                    const GeoModelElement& cur_in_boundary = cur_child.in_boundary(
                        0 ) ;
                    ringmesh_assert( cur_in_boundary.type() == GME::REGION ) ;
                    const Region& cur_reg =
                        dynamic_cast< const Region& >( cur_in_boundary ) ;

                    index_t local_boundary_id = find_local_boundary_id(
                        cur_in_boundary, cur_child ) ;
                    if( cur_reg.side( local_boundary_id ) ) {
                        surfaces_boundary_regions_side_plus[cur_in_boundary.index()].push_back(
                            cur_child.index() ) ;
                    } else {
                        surfaces_boundary_regions_side_minus[cur_in_boundary.index()].push_back(
                            cur_child.index() ) ;
                    }
                }
            }
        }

        build_merged_and_bad_lines( surfaces_boundary_regions_side_plus, "_plus",
            to_erase_by_type, interface_plus_gme_t ) ;
        build_merged_and_bad_lines( surfaces_boundary_regions_side_minus, "_minus",
            to_erase_by_type, interface_minus_gme_t ) ;

        delete_elements( to_erase_by_type ) ;
    }

    void DuplicateInterfaceBuilder::build_merged_and_bad_lines(
        const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
        const std::string& side_name,
        std::vector< std::vector< index_t > >& to_erase_by_type,
        const GME::gme_t& sided_interface_gme_t )
    {
        for( std::map< index_t, std::vector< index_t > >::const_iterator map_itr =
            surfaces_boundary_regions.begin();
            map_itr != surfaces_boundary_regions.end(); ++map_itr ) {

            GEO::Mesh new_surface_mesh ;
            // first = line index in geomodel, second = count.
            std::map< index_t, index_t > all_surface_lines ; // TODO better to handle that with boolean?

            index_t region_index = map_itr->first ;
            std::vector< vec3 > all_points ;
            for( std::vector< index_t >::const_iterator surf_itr =
                map_itr->second.begin(); surf_itr != map_itr->second.end();
                ++surf_itr ) {
                index_t surf_id = *surf_itr ;

                const Surface& cur_surf = model_.surface( surf_id ) ;

                for( index_t vertex_surf_i = 0;
                    vertex_surf_i < cur_surf.nb_vertices(); ++vertex_surf_i ) {
                    all_points.push_back( cur_surf.vertex( vertex_surf_i ) ) ;
                }
            }

            MakeUnique make_unique_surf( all_points ) ;
            make_unique_surf.unique() ;
            std::vector< vec3 > facet_points ;
            make_unique_surf.unique_points( facet_points ) ;
            const std::vector< index_t >& unique_id = make_unique_surf.indices() ;
            index_t offset_vertices = 0 ;
            std::vector< index_t > facet_indices ;
            std::vector< index_t > facet_ptr ;
            index_t count_facet_vertices = 0 ;
            facet_ptr.push_back( count_facet_vertices ) ;
            for( std::vector< index_t >::const_iterator surf_itr =
                map_itr->second.begin(); surf_itr != map_itr->second.end();
                ++surf_itr ) {
                index_t surf_id = *surf_itr ;

                const Surface& cur_surf = model_.surface( surf_id ) ;
                const GEO::Mesh& cur_surf_mesh = cur_surf.mesh() ;

                // Add current surface to merged surface
                for( index_t facet_itr = 0; facet_itr < cur_surf_mesh.facets.nb();
                    ++facet_itr ) {
                    for( index_t point_i = 0;
                        point_i < cur_surf.nb_vertices_in_facet( facet_itr );
                        ++point_i ) {

                        index_t index = cur_surf.surf_vertex_id( facet_itr,
                            point_i ) ;
                        facet_indices.push_back(
                            unique_id[index + offset_vertices] ) ;

                    }
                    count_facet_vertices += cur_surf.nb_vertices_in_facet(
                        facet_itr ) ;
                    facet_ptr.push_back( count_facet_vertices ) ;
                }

                // Update the lines in common
                for( index_t line_itr = 0; line_itr < cur_surf.nb_boundaries();
                    ++line_itr ) {
                    const GeoModelElement& cur_line_gme = cur_surf.boundary(
                        line_itr ) ;
                    ringmesh_assert( cur_line_gme.type() == GME::LINE ) ;

                    const Line& ll = dynamic_cast< const Line& >( cur_line_gme ) ;
                    GEO::mesh_save( ll.mesh(),
                        "line_" + GEO::String::to_string( line_itr ) + ".meshb" ) ;

                    if( all_surface_lines.find( cur_line_gme.index() )
                        == all_surface_lines.end() ) {
                        all_surface_lines[cur_line_gme.index()] = 0 ; // initialization
                    }
                    ++all_surface_lines[cur_line_gme.index()] ;
                }
                offset_vertices += cur_surf.nb_vertices() ;
            }

            // Create RINGMesh::Surface and fill it.
            GME::gme_t new_surface_gme_t = create_element( GME::SURFACE ) ;
            set_surface_geometry( new_surface_gme_t.index, facet_points,
                facet_indices, facet_ptr ) ;
            set_element_parent( new_surface_gme_t, sided_interface_gme_t ) ;
            add_element_child( sided_interface_gme_t, new_surface_gme_t ) ;
            add_element_in_boundary( new_surface_gme_t,
                GME::gme_t( GME::REGION, region_index ) ) ;
            bool side = ( side_name == "_plus" ) ;
            add_element_boundary( GME::gme_t( GME::REGION, region_index ),
                new_surface_gme_t, side ) ;
            // no child
            // boundaries = all the lines of the previous surfaces excepted the
            // one in common
            // parent = the new interface
            // in_boundaries = as the interface is duplicated there is only one single region as in_boundary
            to_erase_by_type[GME::SURFACE].push_back( 0 ) ;
            GEO::mesh_save( new_surface_mesh,
                "merged_surf_reg_" + GEO::String::to_string( region_index )
                    + side_name + ".meshb" ) ;

            // Lines not boundary of the final merged surface
            for( std::map< index_t, index_t >::iterator all_surface_lines_itr =
                all_surface_lines.begin();
                all_surface_lines_itr != all_surface_lines.end();
                ++all_surface_lines_itr ) {

                if( all_surface_lines_itr->second > 1 ) {
                    std::string name = "line_not_kept_merged_surf_" ;
                    name += GEO::String::to_string( region_index ) ;
                    name += "_index_" ;
                    name += GEO::String::to_string( all_surface_lines_itr->first ) ;
                    name += side_name ;
                    name += ".meshb" ;
                    GEO::mesh_save(
                        model_.line( all_surface_lines_itr->first ).mesh(), name ) ;
                    /*to_erase_by_type[GME::LINE][all_surface_lines_itr->first] =
                     NO_ID ;*/
                } else {
                    GME::gme_t line_gme_t( GME::LINE,
                        all_surface_lines_itr->first ) ;
                    add_element_boundary( new_surface_gme_t, line_gme_t ) ;
                    add_element_in_boundary( line_gme_t, new_surface_gme_t ) ;
                }
            }
        }
    }
}
