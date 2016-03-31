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
#include <ringmesh/algorithm.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_geometry.h>

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

        translate_new_interface_by_epsilon_to_avoid_colocation( interface_plus_gme_t,
            to_erase_by_type ) ;
        translate_new_interface_by_epsilon_to_avoid_colocation(
            interface_minus_gme_t, to_erase_by_type ) ;

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
                } else {
                    GME::gme_t line_gme_t( GME::LINE,
                        all_surface_lines_itr->first ) ;
                    add_element_boundary( new_surface_gme_t, line_gme_t ) ;
                    add_element_in_boundary( line_gme_t, new_surface_gme_t ) ;
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::translate_new_interface_by_epsilon_to_avoid_colocation(
        const GME::gme_t& interface_gme_t,
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        // Initialization of the mapping to know which GME of the interface to move.
        std::vector< std::vector< bool > > gme_to_move ;
        std::vector< std::vector< index_t > > gme_in_interface ;
        fill_info_gme_interfation_motion( interface_gme_t, to_erase_by_type,
            gme_to_move, gme_in_interface ) ;
        apply_translation_on_gme_to_move( interface_gme_t, to_erase_by_type,
            gme_to_move, gme_in_interface ) ;
    }

    void DuplicateInterfaceBuilder::fill_info_gme_interfation_motion(
        const GME::gme_t& interface_gme_t,
        const std::vector< std::vector< index_t > >& to_erase_by_type,
        std::vector< std::vector< bool > >& gme_to_move,
        std::vector< std::vector< index_t > >& gme_in_interface )
    {
        gme_to_move.resize( GME::REGION + 1 ) ; // 4 = Corner, Line, Surface, Region (same as enum GME::TYPE)
        gme_in_interface.resize( GME::REGION + 1 ) ; // 4 = Corner, Line, Surface, Region (same as enum GME::TYPE)

        index_t new_nb_surfaces = 0 ;
        for( index_t child_itr = 0;
            child_itr < model_.one_interface( interface_gme_t.index ).nb_children();
            ++child_itr ) {
            const GeoModelElement& child_gme = model_.one_interface(
                interface_gme_t.index ).child( child_itr ) ;
            index_t child_id = child_gme.index() ;
            ringmesh_assert( child_gme.type()==GME::SURFACE ) ;
            if( to_erase_by_type[GME::SURFACE][child_id] != NO_ID ) {
                ringmesh_assert( to_erase_by_type[GME::SURFACE][child_id] == 0 ) ;
                ++new_nb_surfaces ;
                gme_in_interface[GME::SURFACE].push_back( child_id ) ;
                ringmesh_assert(child_gme.nb_in_boundary()==1) ;
                index_t in_boundary_id = child_gme.in_boundary( 0 ).index() ;
                // In theory the surfaces have no in_boundary (region) in common.
                ringmesh_assert(std::find(gme_in_interface[GME::REGION].begin(),
                        gme_in_interface[GME::REGION].end(),in_boundary_id)==gme_in_interface[GME::REGION].end()) ;
                gme_in_interface[GME::REGION].push_back( in_boundary_id ) ;

                for( index_t line_itr = 0; line_itr < child_gme.nb_boundaries();
                    ++line_itr ) {
                    const GeoModelElement& line_gme = child_gme.boundary(
                        line_itr ) ;
                    ringmesh_assert( line_gme.type() == GME::LINE ) ;
                    // For now no line is removed.
                    ringmesh_assert(to_erase_by_type[GME::LINE][line_gme.index()]==0) ;
                    ringmesh_assert(to_erase_by_type[GME::LINE][line_gme.index()]!=NO_ID) ;

                    bool is_line_somewhere_else = false ;
                    for( index_t child_itr2 = 0;
                        child_itr2
                            < model_.one_interface( interface_gme_t.index ).nb_children();
                        ++child_itr2 ) {
                        if( child_itr == child_itr2 ) {
                            continue ;
                        }
                        const GeoModelElement& child_gme2 = model_.one_interface(
                            interface_gme_t.index ).child( child_itr2 ) ;
                        for( index_t line_itr2 = 0;
                            line_itr2 < child_gme2.nb_boundaries(); ++line_itr2 ) {

                            if( child_gme2.boundary( line_itr2 ).index()
                                == line_gme.index() ) {
                                is_line_somewhere_else = true ;
                                break ;
                            }
                        }
                        if( is_line_somewhere_else ) {
                            break ;
                        }
                    }

                    if( !is_line_somewhere_else ) {
                        // As the line is not shared by another surface of the
                        // interface, it cannot be already in the vector.
                        ringmesh_assert(std::find(gme_in_interface[GME::LINE].begin(),
                                gme_in_interface[GME::LINE].end(),line_gme.index())==gme_in_interface[GME::LINE].end()) ;
                        gme_in_interface[GME::LINE].push_back( line_gme.index() ) ;
                        gme_to_move[GME::LINE].push_back( false ) ;
                    } else {
                        // Not added twice.
                        if( std::find( gme_in_interface[GME::LINE].begin(),
                            gme_in_interface[GME::LINE].end(), line_gme.index() )
                            == gme_in_interface[GME::LINE].end() ) {
                            gme_in_interface[GME::LINE].push_back(
                                line_gme.index() ) ;
                            gme_to_move[GME::LINE].push_back( true ) ;
                        }
                    }

                    /// @TODO I have a lot of std::find which is very dirty.
                    /// @TODO When it works, find a better way.
                    /// @todo see RINGMesh::contains in algorithm.h.

                    for( index_t corner_itr = 0;
                        corner_itr < line_gme.nb_boundaries(); ++corner_itr ) {
                        const GeoModelElement& corner_gme = line_gme.boundary(
                            corner_itr ) ;
                        index_t corner_id = corner_gme.index() ;
                        ringmesh_assert( corner_gme.type() == GME::CORNER ) ;
                        if( std::find( gme_in_interface[GME::CORNER].begin(),
                            gme_in_interface[GME::CORNER].end(), corner_id )
                            == gme_in_interface[GME::CORNER].end() ) {
                            gme_in_interface[GME::CORNER].push_back( corner_id ) ;
                            gme_to_move[GME::CORNER].push_back( false ) ; // for now no motion of the corners
                        }
                    }
                }
            }
        }
        ringmesh_assert( new_nb_surfaces != 0 ) ;

        // Surfaces and regions move
        gme_to_move[GME::SURFACE].resize( new_nb_surfaces, true ) ;
        // Each surface on the duplicated interface has only one single region has in_boundary.
        gme_to_move[GME::REGION].resize( new_nb_surfaces, true ) ;
    }

    void DuplicateInterfaceBuilder::apply_translation_on_gme_to_move(
        const GME::gme_t& interface_gme_t,
        const std::vector< std::vector< index_t > >& to_erase_by_type,
        const std::vector< std::vector< bool > >& gme_to_move,
        const std::vector< std::vector< index_t > >& gme_in_interface )
    {
        const GeoModelElement& interface_gme = model_.one_interface(
            interface_gme_t.index ) ;
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
            ringmesh_assert(cur_child.type() == GME::SURFACE) ;
            const Surface& cur_surface = model_.surface( cur_child.index() ) ;
            GEO::Mesh& cur_surf_mesh = cur_surface.mesh() ;
            GEO::compute_normals( cur_surf_mesh ) ;
        }

        // Clear to take into account the new gme in the geomodel.
        model_.mesh.vertices.clear() ;
        const GeoModelMeshVertices& gmmv = model_.mesh.vertices ;
        // In this function we iterate on all the nodes of the interface
        // by iterating on all the nodes of all the children (surfaces).
        // Some nodes are shared between among these surfaces (lines in common).
        // To avoid to iterate twice on these nodes, this vector stores if
        // the corresponding node in the GeoModelMesh is moved (for a node in
        // common between 2 surfaces, there are 2 nodes colocated (one by surface)
        // and so just one in the GeoModelMesh.
        std::vector< bool > has_moved( gmmv.nb(), false ) ; /// @todo treated is a better than moved since motion is not mandatory

        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
            ringmesh_assert(cur_child.type() == GME::SURFACE) ;
            const Surface& cur_surface = model_.surface( cur_child.index() ) ; // avoid dynamic_cast of cur_child
            for( index_t surf_vertex_itr = 0;
                surf_vertex_itr < cur_surface.nb_vertices(); ++surf_vertex_itr ) {
                index_t vertex_id_in_gmm = cur_surface.model_vertex_id(
                    surf_vertex_itr ) ;
                if( has_moved[vertex_id_in_gmm] ) {
                    continue ;
                }

                has_moved[vertex_id_in_gmm] = true ;
                // Gets all the GME with a vertex colocated to the one of vertex_id_in_gmm
                const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
                    vertex_id_in_gmm ) ;

                // The idea of these commented lines was the store separately
                // the corners, the lines, the surfaces and the regions and then
                // check the motion by order of priority. If the corner cannot
                // move so no motion. If it may move so motion. So the motion is
                // imposed by the corner if any. If no corner so the motion is imposed
                // by the lines if any...
                // For now I do not do that. For now if at least one gme forbids
                // the motion so no motion. We will see if it works.
                /*std::vector< std::vector< index_t > > only_kept_gme ;
                 only_kept_gme.resize( GME::REGION ) ;

                 for( index_t type_itr = 0; type_itr < GME::REGION; ++type_itr ) {
                 only_kept_gme[type_itr].reserve( gme_vertices.size() ) ; // a little large but no implicit resize
                 // needed to ensure if there are both corner and lines.... to evaluate motion priority
                 ringmesh_assert(only_kept_gme[type_itr].size()==0) ;
                 }*/

                bool to_move = true ;
                std::vector< GMEVertex > only_kept_gme_vertices ;
                only_kept_gme_vertices.reserve( gme_vertices.size() ) ; // a little large but no implicit resize
                for( index_t gme_vertex_itr = 0;
                    gme_vertex_itr < gme_vertices.size(); ++gme_vertex_itr ) {
                    const GME::gme_t& cur_gme_t = gme_vertices[gme_vertex_itr].gme_id ;
                    if( to_erase_by_type[cur_gme_t.type][cur_gme_t.index]
                        == NO_ID ) {
                        continue ; // It is an old element to remove.
                    }
                    ringmesh_assert(
                        to_erase_by_type[cur_gme_t.type][cur_gme_t.index] == 0 ) ;
                    //only_kept_gme[cur_gme_t.type].push_back( cur_gme_t.index ) ;

                    index_t pos = RINGMesh::find( gme_in_interface[cur_gme_t.type],
                        cur_gme_t.index ) ;
                    if( pos == NO_ID ) {
                        // It is a gme not inside the interface but which has a
                        // node with the same location (example: an horizon in
                        // contact with a fault interface).
                        continue ;
                    }
                    only_kept_gme_vertices.push_back(
                        gme_vertices[gme_vertex_itr] ) ;
                    if( !gme_to_move[cur_gme_t.type][pos] ) {
                        to_move = false ;
                        break ;
                    }
                }

                if( !to_move ) {
                    continue ;
                }

                // Finds normal and side of displacement
                vec3 displacement ;
                for( index_t gme_vertex_itr = 0;
                    gme_vertex_itr < only_kept_gme_vertices.size();
                    ++gme_vertex_itr ) {
                    const GME::gme_t& cur_gme_t =
                        only_kept_gme_vertices[gme_vertex_itr].gme_id ;
                    if( cur_gme_t.type == GME::SURFACE ) {
                        const Surface& cur_surf = model_.surface( cur_gme_t.index ) ;
                        // only one side for the sided interface
                        ringmesh_assert( cur_surf.nb_in_boundary() == 1 ) ;
                        const GeoModelElement& in_boun = cur_surf.in_boundary( 0 ) ;
                        ringmesh_assert( in_boun.type() == GME::REGION ) ;
                        const Region& cur_reg = model_.region( in_boun.index() ) ;
                        index_t local_surf_id = find_local_boundary_id( cur_reg,
                            cur_surf ) ;
                        bool side = cur_reg.side( local_surf_id ) ;
                        const vec3& normal = GEO::Geom::mesh_vertex_normal(
                            cur_surf.mesh(),
                            only_kept_gme_vertices[gme_vertex_itr].v_id ) ;
                        ringmesh_assert( std::abs(normal.length() -1.)<epsilon ) ;
                        if( side ) {
                            displacement = normal ;
                        } else {
                            displacement = -1 * normal ;
                        }
                        displacement *= 1.5 * 200 ;
                        break ;
                    }
                }

                // move all the gme.
                for( index_t gme_vertex_itr = 0;
                    gme_vertex_itr < only_kept_gme_vertices.size();
                    ++gme_vertex_itr ) {

                    // Move only_kept_gme_vertices[gme_vertex_itr].gme_t
                    // at vertex only_kept_gme_vertices[gme_vertex_itr].v_id
                    GEO::Mesh& cur_mesh = model_.mesh_element(
                        only_kept_gme_vertices[gme_vertex_itr].gme_id ).mesh() ;
                    vec3& cur_pos = cur_mesh.vertices.point(
                        only_kept_gme_vertices[gme_vertex_itr].v_id ) ;
                    cur_pos += displacement ;
                }
            }
        }
    }
}
