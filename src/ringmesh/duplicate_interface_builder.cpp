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

    void fill_vect_with_NO_ID( std::vector< index_t >& to_fill )
    {
        std::fill( to_fill.begin(), to_fill.end(), NO_ID ) ;
    }

    void DuplicateInterfaceBuilder::duplicate_fault_network()
    {
        //////////////////////////////////////////
        // copy paste from void GeoModelEditor::remove_elements( const std::set< gme_t >& elements )
        std::vector< std::vector< index_t > > to_erase_by_type ;
        to_erase_by_type.reserve( GME::NO_TYPE ) ;
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
            to_erase_by_type.push_back(
                std::vector< index_t >(
                    model_.nb_elements( static_cast< GME::TYPE >( i ) ), 0 ) ) ;
        }
        //////////////////////////////////////////

        const index_t nb_initial_interfaces = model_.nb_interfaces() ;
        // Loop to nb_initial_interfaces. model_.nb_interfaces() cannot be inside
        // the for statement since the number of interfaces will increase during
        // the duplication.
        for( index_t interface_itr = 0; interface_itr < nb_initial_interfaces;
            ++interface_itr ) {
            const GeoModelElement& cur_interface = model_.one_interface(
                interface_itr ) ;
            if( !GeoModelElement::is_fault( cur_interface.geological_feature() ) ) {
                continue ;
            }

            // Delete of the interface (will be replaced by a custom interface with
            // side informations)
            to_erase_by_type[GME::INTERFACE][cur_interface.index()] = NO_ID ;
            get_new_surfaces( cur_interface, to_erase_by_type ) ;
        }

        /// @todo a clear of the geomodelmesh vertices may be necessary somewhere

        // compute translation vectors
        compute_translation_vectors_duplicated_fault_network( nb_initial_interfaces,
            to_erase_by_type ) ;
        // set no translation on fault real extension (only on fault ending inside
        // the model).
        // apply translation
        translate_duplicated_fault_network( to_erase_by_type ) ;

        // Put here for new.
        /// @todo if all the lines are removed, is it still necessary to fill the new
        /// interface children with them? They should be recomputed with
        /// build_lines_and_corners_from_surfaces.
        fill_vect_with_NO_ID( to_erase_by_type[GME::CORNER] ) ;
        fill_vect_with_NO_ID( to_erase_by_type[GME::LINE] ) ;
        fill_vect_with_NO_ID( to_erase_by_type[GME::CONTACT] ) ;

        model_.mesh.vertices.clear() ;
        model_.mesh.vertices.test_and_initialize() ;

        DEBUG( model_.nb_regions() ) ;

        delete_elements( to_erase_by_type ) ;

        /*model_.mesh.vertices.clear() ;
         model_.mesh.vertices.test_and_initialize() ;
         build_model_from_surfaces() ;
         build_contacts() ;
         return;*/

        DEBUG( model_.nb_regions() ) ;

        model_.mesh.vertices.clear() ;
        model_.mesh.vertices.test_and_initialize() ;

        build_lines_and_corners_from_surfaces() ;
        DEBUG( model_.nb_regions() ) ;

        model_.mesh.vertices.clear() ;
        model_.mesh.vertices.test_and_initialize() ;

        /*build_brep_regions_from_surfaces() ;
         DEBUG( model_.nb_regions() ) ;
         model_.mesh.vertices.clear() ;
         model_.mesh.vertices.test_and_initialize() ;*/
        fill_elements_boundaries( GME::SURFACE ) ;

        build_contacts() ;
        DEBUG( model_.nb_regions() ) ;
        model_.mesh.vertices.clear() ;
        model_.mesh.vertices.test_and_initialize() ;

        // To debug a tmp result
        for( index_t corner_itr = 0; corner_itr < model_.nb_corners();
            ++corner_itr ) {
            const Corner& cur_corner = model_.corner( corner_itr ) ;
            const GEO::Mesh& cur_mesh = cur_corner.mesh() ;
            GEO::mesh_save( cur_mesh,
                "mesh_corner_" + GEO::String::to_string( corner_itr ) + ".meshb" ) ;
        }

        for( index_t line_itr = 0; line_itr < model_.nb_lines(); ++line_itr ) {
            const Line& cur_line = model_.line( line_itr ) ;
            const GEO::Mesh& cur_mesh = cur_line.mesh() ;
            GEO::mesh_save( cur_mesh,
                "mesh_line_" + GEO::String::to_string( line_itr ) + ".meshb" ) ;
        }

        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            const Surface& cur_surf = model_.surface( surf_itr ) ;
            const GEO::Mesh& cur_mesh = cur_surf.mesh() ;
            GEO::mesh_save( cur_mesh,
                "mesh_surf_" + GEO::String::to_string( surf_itr ) + ".meshb" ) ;
        }
        // To debug a tmp result
    }

    void DuplicateInterfaceBuilder::get_new_surfaces(
        const GeoModelElement& interface_to_duplicate,
        std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        // minus = false, plus = true
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_minus ;
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_plus ;
        GME::gme_t interface_minus_gme_t = create_element( GME::INTERFACE ) ;
        const_cast< GeoModelElement& >( model_.one_interface(
            interface_minus_gme_t.index ) ).set_name(
            interface_to_duplicate.name() + "_side_minus" ) ;
        const_cast< GeoModelElement& >( model_.one_interface(
            interface_minus_gme_t.index ) ).set_geological_feature(
            interface_to_duplicate.geological_feature() ) ;
        GME::gme_t interface_plus_gme_t = create_element( GME::INTERFACE ) ;
        const_cast< GeoModelElement& >( model_.one_interface(
            interface_plus_gme_t.index ) ).set_name(
            interface_to_duplicate.name() + "_side_plus" ) ;
        const_cast< GeoModelElement& >( model_.one_interface(
            interface_plus_gme_t.index ) ).set_geological_feature(
            interface_to_duplicate.geological_feature() ) ;
        to_erase_by_type[GME::INTERFACE].push_back( 0 ) ;
        to_erase_by_type[GME::INTERFACE].push_back( 0 ) ;

        const index_t interface_to_duplicate_nb_children =
            interface_to_duplicate.nb_children() ;
        ringmesh_assert(interface_to_duplicate_nb_children >= 1) ;
        // Find for each region, what surfaces are in boundary.
        for( index_t interface_child_itr = 0;
            interface_child_itr < interface_to_duplicate_nb_children;
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
                ringmesh_assert( cur_child.is_on_voi() ) ;
                const GeoModelElement& cur_in_boundary = cur_child.in_boundary( 0 ) ;
                ringmesh_assert( cur_in_boundary.type() == GME::REGION ) ;
                const Region& cur_reg =
                    dynamic_cast< const Region& >( cur_in_boundary ) ;

                index_t local_boundary_id = find_local_boundary_id( cur_in_boundary,
                    cur_child ) ;
                if( cur_reg.side( local_boundary_id ) ) {
                    surfaces_boundary_regions_side_plus[cur_in_boundary.index()].push_back(
                        cur_child.index() ) ;
                } else {
                    surfaces_boundary_regions_side_minus[cur_in_boundary.index()].push_back(
                        cur_child.index() ) ;
                }
            }
        }

        build_merged_and_bad_lines( surfaces_boundary_regions_side_plus, "_plus",
            to_erase_by_type, interface_plus_gme_t ) ;
        build_merged_and_bad_lines( surfaces_boundary_regions_side_minus, "_minus",
            to_erase_by_type, interface_minus_gme_t ) ;
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

    void DuplicateInterfaceBuilder::initialize_translation_attributes(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        for( index_t reg_itr = 0; reg_itr < model_.nb_regions(); ++reg_itr ) {
            const Region& reg = model_.region( reg_itr ) ;
            if( !reg.is_meshed() ) {
                continue ;
            }
            ringmesh_assert(to_erase_by_type[GME::REGION][reg_itr]!=NO_ID) ;
            ringmesh_assert(to_erase_by_type[GME::REGION][reg_itr]==0) ;
            GEO::Mesh& reg_mesh = reg.mesh() ;
            GEO::AttributesManager& att_mgr = reg_mesh.vertices.attributes() ;
            GEO::Attribute< double > translation_att_x( att_mgr,
                "translation_attr_x" ) ;
            translation_att_x.fill( 0. ) ;
            GEO::Attribute< double > translation_att_y( att_mgr,
                "translation_attr_y" ) ;
            translation_att_y.fill( 0. ) ;
            GEO::Attribute< double > translation_att_z( att_mgr,
                "translation_attr_z" ) ;
            translation_att_z.fill( 0. ) ;
        }

        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            const Surface& surf = model_.surface( surf_itr ) ;

            if( to_erase_by_type[GME::SURFACE][surf_itr] == NO_ID ) {
                continue ;
            }

            GEO::Mesh& surf_mesh = surf.mesh() ;
            GEO::AttributesManager& att_mgr = surf_mesh.vertices.attributes() ;
            GEO::Attribute< double > translation_att_x( att_mgr,
                "translation_attr_x" ) ;
            translation_att_x.fill( 0. ) ;
            GEO::Attribute< double > translation_att_y( att_mgr,
                "translation_attr_y" ) ;
            translation_att_y.fill( 0. ) ;
            GEO::Attribute< double > translation_att_z( att_mgr,
                "translation_attr_z" ) ;
            translation_att_z.fill( 0. ) ;
        }
    }

    void DuplicateInterfaceBuilder::save_normals_on_one_new_interface(
        const std::vector< std::vector< index_t > >& to_erase_by_type,
        const GeoModelElement& interface_gme ) const
    {
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
            ringmesh_assert(cur_child.type() == GME::SURFACE) ;
            // As the loop begins at the first new interface, no surface
            // met in this loop should be to delete.
            ringmesh_assert( to_erase_by_type[GME::SURFACE][cur_child.index()] == 0 ) ;
            ringmesh_assert( to_erase_by_type[GME::SURFACE][cur_child.index()] != NO_ID ) ;
            const Surface& cur_surface = model_.surface( cur_child.index() ) ;
            GEO::Mesh& cur_surf_mesh = cur_surface.mesh() ;
            // GEO::compute_normals cannot be used because the dimension
            // of the vertices from 3 to 6 and that provokes a problem
            // of copying in GeoModelMeshVertices::initialize
            // with GEO::Memory::copy( mesh_.vertices.point_ptr( count ),
            // E.vertex( 0 ).data(), 3 * E.nb_vertices() * sizeof(double) ) ;
            // 3 means vertices of dimension 3 and not another dimension.
            GEO::AttributesManager& att_mgr = cur_surf_mesh.vertices.attributes() ;
            GEO::Attribute< double > normal_att_x( att_mgr, "normal_attr_x" ) ;
            GEO::Attribute< double > normal_att_y( att_mgr, "normal_attr_y" ) ;
            GEO::Attribute< double > normal_att_z( att_mgr, "normal_attr_z" ) ;
            normal_att_x.fill( 0. ) ;
            normal_att_y.fill( 0. ) ;
            normal_att_z.fill( 0. ) ;
            // begin copy paste from GEO::compute_normals
            for( index_t f = 0; f < cur_surf_mesh.facets.nb(); f++ ) {
                vec3 N = GEO::Geom::mesh_facet_normal( cur_surf_mesh, f ) ;
                for( index_t corner = cur_surf_mesh.facets.corners_begin( f );
                    corner < cur_surf_mesh.facets.corners_end( f ); corner++ ) {
                    index_t v = cur_surf_mesh.facet_corners.vertex( corner ) ;
                    normal_att_x[v] += N.x ;
                    normal_att_y[v] += N.y ;
                    normal_att_z[v] += N.z ;
                }
            }
            for( index_t i = 0; i < cur_surf_mesh.vertices.nb(); i++ ) {
                vec3 cur_normal( normal_att_x[i], normal_att_y[i],
                    normal_att_z[i] ) ;
                cur_normal = normalize( cur_normal ) ;
                normal_att_x[i] = cur_normal.x ;
                normal_att_y[i] = cur_normal.y ;
                normal_att_z[i] = cur_normal.z ;
            }
            // end copy paste from GEO::compute_normals
        }
    }

    vec3 DuplicateInterfaceBuilder::get_local_translation_vector(
        const Surface& surface,
        index_t vertex_id_in_surface ) const
    {
        vec3 displacement ;
        // only one side for the sided interface
        ringmesh_assert( surface.nb_in_boundary() == 1 ) ;
        const GeoModelElement& in_boun = surface.in_boundary( 0 ) ;
        ringmesh_assert( in_boun.type() == GME::REGION ) ;
        const Region& cur_reg = model_.region( in_boun.index() ) ;
        index_t local_surf_id = find_local_boundary_id( cur_reg, surface ) ;
        bool side = cur_reg.side( local_surf_id ) ;

        GEO::AttributesManager& att_mgr = surface.mesh().vertices.attributes() ;
        GEO::Attribute< double > normal_att_x( att_mgr, "normal_attr_x" ) ;
        GEO::Attribute< double > normal_att_y( att_mgr, "normal_attr_y" ) ;
        GEO::Attribute< double > normal_att_z( att_mgr, "normal_attr_z" ) ;
        vec3 normal( normal_att_x[vertex_id_in_surface],
            normal_att_y[vertex_id_in_surface],
            normal_att_z[vertex_id_in_surface] ) ;

        ringmesh_assert( std::abs(normal.length() -1.)<epsilon ) ;
        if( side ) {
            displacement = normal ;
        } else {
            displacement = -1 * normal ;
        }
        displacement *= 1.5 * epsilon ;
        return displacement ;
    }

    void DuplicateInterfaceBuilder::compute_translation_vectors_duplicated_fault_network(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        initialize_translation_attributes( to_erase_by_type ) ;

        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr < model_.nb_interfaces(); ++new_interface_itr ) {
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] != NO_ID ) ;
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] == 0 ) ;

            const GeoModelElement& interface_gme = model_.one_interface(
                new_interface_itr ) ;

            save_normals_on_one_new_interface( to_erase_by_type, interface_gme ) ;

            // Clear to take into account the new gme in the geomodel.
            model_.mesh.vertices.clear() ;
            model_.mesh.vertices.test_and_initialize() ; // not done in model_vertex_id

            const GeoModelMeshVertices& gmmv = model_.mesh.vertices ;

            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
                ringmesh_assert(cur_child.type() == GME::SURFACE) ;
                ringmesh_assert(to_erase_by_type[GME::SURFACE][cur_child.index()] != NO_ID) ;
                const Surface& cur_surface = model_.surface( cur_child.index() ) ; // avoid dynamic_cast of cur_child

                ringmesh_assert(cur_surface.nb_in_boundary()==1) ;
                ringmesh_assert(cur_surface.in_boundary(0).type()==GME::REGION) ;
                const index_t surf_in_boun_reg =
                    cur_surface.in_boundary( 0 ).index() ;
                for( index_t surf_vertex_itr = 0;
                    surf_vertex_itr < cur_surface.nb_vertices();
                    ++surf_vertex_itr ) {
                    const index_t vertex_id_in_gmm = cur_surface.model_vertex_id(
                        surf_vertex_itr ) ;

                    /// @todo potentially for the vertices in common between 2 surfaces
                    /// the translation is applied twice. Does not seem to be
                    // a problem. To check.

                    // Gets all the GME with a vertex colocated to the one of vertex_id_in_gmm
                    const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
                        vertex_id_in_gmm ) ;
                    ringmesh_assert(std::find(gme_vertices.begin(),gme_vertices.end(),
                            GMEVertex(GME::gme_t(GME::SURFACE,cur_child.index()), surf_vertex_itr))
                        != gme_vertices.end()) ;

                    const vec3 local_translation_vector =
                        get_local_translation_vector( cur_surface,
                            surf_vertex_itr ) ;

                    for( index_t gme_vertex_itr = 0;
                        gme_vertex_itr < gme_vertices.size(); ++gme_vertex_itr ) {
                        const GME::gme_t& cur_gme_t =
                            gme_vertices[gme_vertex_itr].gme_id ;
                        if( to_erase_by_type[cur_gme_t.type][cur_gme_t.index]
                            == NO_ID ) {
                            continue ; // It is an old element to remove.
                        }
                        ringmesh_assert(
                            to_erase_by_type[cur_gme_t.type][cur_gme_t.index] == 0 ) ;

                        if( cur_gme_t.type != GME::SURFACE
                            && cur_gme_t.type != GME::REGION ) {
                            continue ;
                        }

                        if( cur_gme_t.type == GME::REGION ) {
                            if( !model_.region( cur_gme_t.index ).is_meshed() ) {
                                continue ;
                            }
                            if( cur_gme_t.index != surf_in_boun_reg ) {
                                // Region on the other side of the fault
                                continue ;
                            }
                            store_displacement_in_gme(
                                model_.region( cur_gme_t.index ),
                                gme_vertices[gme_vertex_itr].v_id,
                                local_translation_vector ) ;
                        } else {
                            ringmesh_assert(cur_gme_t.type == GME::SURFACE) ;
                            const Surface& cur_gme_surf = model_.surface(
                                cur_gme_t.index ) ;
                            bool ok = false ;
                            for( index_t cur_gme_surf_in_boun_itr = 0;
                                cur_gme_surf_in_boun_itr
                                    < cur_gme_surf.nb_in_boundary();
                                ++cur_gme_surf_in_boun_itr ) {
                                const GeoModelElement& cur_gme_surf_in_boun =
                                    cur_gme_surf.in_boundary(
                                        cur_gme_surf_in_boun_itr ) ;
                                ringmesh_assert(
                                    cur_gme_surf_in_boun.type() == GME::REGION ) ;
                                if( cur_gme_surf_in_boun.index()
                                    == surf_in_boun_reg ) {
                                    ok = true ;
                                    break ;
                                }
                            }

                            if( !ok ) {
                                continue ;
                            }
                            store_displacement_in_gme( cur_gme_surf,
                                gme_vertices[gme_vertex_itr].v_id,
                                local_translation_vector ) ;
                        }

                    }
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::store_displacement_in_gme(
        const GeoModelMeshElement& gmme,
        index_t vertex_id_in_gmme,
        const vec3& translation ) const
    {
        GEO::AttributesManager& att_mgr = gmme.mesh().vertices.attributes() ;
        GEO::Attribute< double > translation_att_x( att_mgr, "translation_attr_x" ) ;
        translation_att_x[vertex_id_in_gmme] += translation.x ;
        GEO::Attribute< double > translation_att_y( att_mgr, "translation_attr_y" ) ;
        translation_att_y[vertex_id_in_gmme] += translation.y ;
        GEO::Attribute< double > translation_att_z( att_mgr, "translation_attr_z" ) ;
        translation_att_z[vertex_id_in_gmme] += translation.z ;
    }

    void DuplicateInterfaceBuilder::translate_duplicated_fault_network(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        for( index_t reg_itr = 0; reg_itr < model_.nb_regions(); ++reg_itr ) {
            const Region& reg = model_.region( reg_itr ) ;
            if( !reg.is_meshed() ) {
                continue ;
            }
            ringmesh_assert(to_erase_by_type[GME::REGION][reg_itr]!=NO_ID) ;
            ringmesh_assert(to_erase_by_type[GME::REGION][reg_itr]==0) ;
            GEO::Mesh& reg_mesh = reg.mesh() ;
            GEO::AttributesManager& att_mgr = reg_mesh.vertices.attributes() ;
            GEO::Attribute< double > translation_att_x( att_mgr,
                "translation_attr_x" ) ;
            GEO::Attribute< double > translation_att_y( att_mgr,
                "translation_attr_y" ) ;
            GEO::Attribute< double > translation_att_z( att_mgr,
                "translation_attr_z" ) ;
            for( index_t vertex_itr = 0; vertex_itr < reg_mesh.vertices.nb();
                ++vertex_itr ) {
                reg_mesh.vertices.point( vertex_itr ).x +=
                    translation_att_x[vertex_itr] ;
                reg_mesh.vertices.point( vertex_itr ).y +=
                    translation_att_y[vertex_itr] ;
                reg_mesh.vertices.point( vertex_itr ).z +=
                    translation_att_z[vertex_itr] ;
            }
        }

        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            const Surface& surf = model_.surface( surf_itr ) ;

            if( to_erase_by_type[GME::SURFACE][surf_itr] == NO_ID ) {
                continue ;
            }

            GEO::Mesh& surf_mesh = surf.mesh() ;
            GEO::AttributesManager& att_mgr = surf_mesh.vertices.attributes() ;
            GEO::Attribute< double > translation_att_x( att_mgr,
                "translation_attr_x" ) ;
            GEO::Attribute< double > translation_att_y( att_mgr,
                "translation_attr_y" ) ;
            GEO::Attribute< double > translation_att_z( att_mgr,
                "translation_attr_z" ) ;
            for( index_t vertex_itr = 0; vertex_itr < surf_mesh.vertices.nb();
                ++vertex_itr ) {
                surf_mesh.vertices.point( vertex_itr ).x +=
                    translation_att_x[vertex_itr] ;
                surf_mesh.vertices.point( vertex_itr ).y +=
                    translation_att_y[vertex_itr] ;
                surf_mesh.vertices.point( vertex_itr ).z +=
                    translation_att_z[vertex_itr] ;
            }
        }
    }
}
