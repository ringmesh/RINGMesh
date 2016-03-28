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

    void DuplicateInterfaceBuilder::get_new_surfaces(
        index_t interface_id_to_duplicate ) const
    {
        ringmesh_assert(interface_id_to_duplicate < model_.nb_interfaces()) ;

        const GeoModelElement& interface_to_duplicate = model_.one_interface(
            interface_id_to_duplicate ) ;

        const index_t interface_to_duplicate_nb_children =
            interface_to_duplicate.nb_children() ;
        ringmesh_assert(interface_to_duplicate_nb_children >= 1) ;

        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions ;

        // Find for each region, what surfaces are in boundary.
        for( index_t interface_child_itr = 0;
            interface_child_itr < interface_to_duplicate.nb_children();
            ++interface_child_itr ) {
            const GeoModelElement& cur_child = interface_to_duplicate.child(
                interface_child_itr ) ;
            ringmesh_assert( cur_child.type() == GME::SURFACE ) ;

            const index_t nb_in_boundary_cur_child = cur_child.nb_in_boundary() ;
            ringmesh_assert( nb_in_boundary_cur_child == 1 || nb_in_boundary_cur_child == 2 ) ;
            for( index_t in_boundary_itr = 0;
                in_boundary_itr < nb_in_boundary_cur_child; ++in_boundary_itr ) {

                const GeoModelElement& cur_in_boundary = cur_child.in_boundary(
                    in_boundary_itr ) ;
                ringmesh_assert( cur_in_boundary.type() == GME::REGION ) ;

                surfaces_boundary_regions[cur_in_boundary.index()].push_back(
                    cur_child.index() ) ;
            }
        }

        for( std::map< index_t, std::vector< index_t > >::iterator map_itr =
            surfaces_boundary_regions.begin();
            map_itr != surfaces_boundary_regions.end(); ++map_itr ) {

            GEO::Mesh new_surface_mesh ;
            // first = line index in geomodel, second = count.
            std::map< index_t, index_t > all_surface_lines ; // TODO better to handle that with boolean?
            // first = contact index in geomodel, second = vector of line indexes.
            std::map< index_t, std::vector< index_t > > contact_to_lines ;

            index_t region_index = map_itr->first ;
            for( std::vector< index_t >::iterator surf_itr = map_itr->second.begin();
                surf_itr != map_itr->second.end(); ++surf_itr ) {
                index_t surf_id = *surf_itr ;

                const Surface& cur_surf = model_.surface( surf_id ) ;
                const GEO::Mesh& cur_surf_mesh = cur_surf.mesh() ;

                // Add current surface to merged surface
                for( index_t facet_itr = 0; facet_itr < cur_surf_mesh.facets.nb();
                    ++facet_itr ) {
                    index_t one = find_or_create_vertex_facet( cur_surf_mesh,
                        facet_itr, 0, new_surface_mesh ) ;
                    index_t two = find_or_create_vertex_facet( cur_surf_mesh,
                        facet_itr, 1, new_surface_mesh ) ;
                    index_t three = find_or_create_vertex_facet( cur_surf_mesh,
                        facet_itr, 2, new_surface_mesh ) ;
                    new_surface_mesh.facets.create_triangle( one, two, three ) ;
                }

                // Update the lines in common and to merge
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
                        DEBUG("initialization") ;
                        all_surface_lines[cur_line_gme.index()] = 0 ; // initialization
                    }
                    ++all_surface_lines[cur_line_gme.index()] ;
                    const GeoModelElement& cur_parent_gme = cur_line_gme.parent() ;
                    ringmesh_assert( cur_parent_gme.type() == GME::CONTACT ) ;
                    if( contact_to_lines.find( cur_parent_gme.index() )
                        == contact_to_lines.end() ) {
                        contact_to_lines[cur_parent_gme.index()] = std::vector<
                            index_t >() ; // TODO necessary???
                    }
                    contact_to_lines[cur_parent_gme.index()].push_back(
                        cur_line_gme.index() ) ;
                }
            }

            GEO::mesh_save( new_surface_mesh,
                "merged_surf_reg_" + GEO::String::to_string( region_index )
                    + ".meshb" ) ;

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
                    name += ".meshb" ;
                    GEO::mesh_save(
                        model_.line( all_surface_lines_itr->first ).mesh(), name ) ;
                    std::map< index_t, std::vector< index_t > >::iterator bad_contact =
                        contact_to_lines.find(
                            model_.line( all_surface_lines_itr->first ).parent().index() ) ;
                    if( bad_contact != contact_to_lines.end() ) {
                        contact_to_lines.erase( bad_contact ) ;
                    }
                }
            }

            // All the lines good and to merged (bad have been deleted in the previous
            // for loop.
            index_t oo = 0 ;
            for( std::map< index_t, std::vector< index_t > >::iterator contact_to_lines_itr =
                contact_to_lines.begin();
                contact_to_lines_itr != contact_to_lines.end();
                ++contact_to_lines_itr ) {

                index_t cur_contact_id = contact_to_lines_itr->first ;
                GEO::Mesh new_line_mesh ;
                for( std::vector< index_t >::iterator contact_line_itr =
                    contact_to_lines_itr->second.begin();
                    contact_line_itr != contact_to_lines_itr->second.end();
                    ++contact_line_itr ) {
                    index_t cur_line = *contact_line_itr ;
                    const GEO::Mesh& cur_line_mesh = model_.line( cur_line ).mesh() ;

                    for( index_t edge_itr = 0; edge_itr < cur_line_mesh.edges.nb();
                        ++edge_itr ) {
                        index_t one = find_or_create_vertex_edge( cur_line_mesh,
                            edge_itr, 0, new_line_mesh ) ;
                        index_t two = find_or_create_vertex_edge( cur_line_mesh,
                            edge_itr, 1, new_line_mesh ) ;
                        new_surface_mesh.edges.create_edge( one, two ) ;
                    }
                }

                std::string name = "merged_line_contact_" ;
                name += GEO::String::to_string( cur_contact_id ) ;
                name += "_nb_" ;
                name += GEO::String::to_string( oo ) ;
                name += ".meshb" ;
                GEO::mesh_save( new_line_mesh, name ) ;
                ++oo ;
            }

#if 0
            // To complete the region, add the other surfaces of the region to
            // the newly built merged surface. The surfaces before the merge must
            // not be added.
            for( index_t surf_reg_itr = 0 ;
                surf_reg_itr < model_.region( region_index ).nb_boundaries() ;
                ++surf_reg_itr ) {

                const GeoModelElement& cur_elt =
                model_.region( region_index ).boundary( surf_reg_itr ) ;
                ringmesh_assert( cur_elt.type()==GME::SURFACE ) ;
                std::vector< index_t >::iterator found = std::find(
                    map_itr->second.begin(), map_itr->second.end(),
                    cur_elt.index() ) ;
                // One surface used for previous merging. Avoided.
                // THESE SURFACES MUST BE REMOVED.
                if( found != map_itr->second.end() ) {
                    continue ;
                }

                // Add the facets for debugging for now. But I found all the
                // surfaces of the regions.
                const Surface& cur_surf2 = dynamic_cast< const Surface& >( cur_elt ) ;
                const GEO::Mesh& cur_surf_mesh2 = cur_surf2.mesh() ;
                for( index_t facet_itr = 0 ; facet_itr < cur_surf_mesh2.facets.nb() ;
                    ++facet_itr ) {
                    index_t one = find_or_create_vertex( cur_surf_mesh2, facet_itr,
                        0, new_mesh ) ;
                    index_t two = find_or_create_vertex( cur_surf_mesh2, facet_itr,
                        1, new_mesh ) ;
                    index_t three = find_or_create_vertex( cur_surf_mesh2, facet_itr,
                        2, new_mesh ) ;
                    new_mesh.facets.create_triangle( one, two, three ) ;
                }
            }
            GEO::mesh_save( new_mesh,
                "surf_reg_" + GEO::String::to_string( region_index ) + ".meshb" ) ;
#endif
        }

    }

    index_t DuplicateInterfaceBuilder::find_or_create_vertex_facet(
        const GEO::Mesh& cur_surf_mesh,
        index_t facet_itr,
        index_t v,
        GEO::Mesh& new_mesh ) const
    {
        const GEO::MeshFacets& cur_surf_mesh_facets = cur_surf_mesh.facets ;
        const GEO::MeshVertices& cur_surf_mesh_verticess = cur_surf_mesh.vertices ;
        ColocaterANN ann( new_mesh, ColocaterANN::VERTICES ) ;
        std::vector< index_t > colocated ;
        const vec3& cur_point = cur_surf_mesh.vertices.point(
            cur_surf_mesh_facets.vertex( facet_itr, v ) ) ;
        if( ann.get_colocated( cur_point, colocated ) ) {
            ringmesh_assert( colocated.size() == 1 ) ;
            return colocated[0] ;
        }
        return new_mesh.vertices.create_vertex( cur_point.data() ) ;
    }

    index_t DuplicateInterfaceBuilder::find_or_create_vertex_edge(
        const GEO::Mesh& cur_line_mesh,
        index_t edge_itr,
        index_t v,
        GEO::Mesh& new_mesh ) const
    {
        const GEO::MeshEdges& cur_line_mesh_edges = cur_line_mesh.edges ;
        const GEO::MeshVertices& cur_line_mesh_vertices = cur_line_mesh.vertices ;
        ColocaterANN ann( new_mesh, ColocaterANN::VERTICES ) ;
        std::vector< index_t > colocated ;
        const vec3& cur_point = cur_line_mesh.vertices.point(
            cur_line_mesh_edges.vertex( edge_itr, v ) ) ;
        if( ann.get_colocated( cur_point, colocated ) ) {
            ringmesh_assert( colocated.size() == 1 ) ;
            return colocated[0] ;
        }
        return new_mesh.vertices.create_vertex( cur_point.data() ) ;
    }
}
