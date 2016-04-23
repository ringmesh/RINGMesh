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
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_topology.h>

/*!
 * @file ringmesh/duplicate_interface_builder.cpp
 * @brief Class to duplicate GeoModel Interface to
 * enable sliding along them (faults) and unconformal
 * mesh generation.
 * @author Benjamin Chauvin
 */

namespace RINGMesh {

    GEO::Mesh to_debug66 ;

    DuplicateInterfaceBuilder::DuplicateInterfaceBuilder( GeoModel& model )
        : GeoModelBuilder( model ), all_meshed_( true ), gme_vertices_links_()
    {
    }

    DuplicateInterfaceBuilder::~DuplicateInterfaceBuilder()
    {
        for( index_t gme_vertices_links_itr = 0;
            gme_vertices_links_itr < gme_vertices_links_.size();
            ++gme_vertices_links_itr ) {
            ringmesh_assert( gme_vertices_links_[gme_vertices_links_itr] ) ;
            delete gme_vertices_links_[gme_vertices_links_itr] ;
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

    index_t find_second_local_boundary_id(
        const GeoModelElement& reg,
        const GeoModelElement& surf )
    {
        ringmesh_assert(reg.type()==GME::REGION) ;
        ringmesh_assert(surf.type()==GME::SURFACE) ;
        bool is_second = false ;
        for( index_t boundary_i = 0; boundary_i < reg.nb_boundaries();
            ++boundary_i ) {
            if( reg.boundary( boundary_i ).index() == surf.index() ) {
                if( !is_second ) {
                    is_second = true ;
                } else {
                    return boundary_i ;
                }
            }
        }

        return NO_ID ;
    }

    void fill_vect_with_NO_ID( std::vector< index_t >& to_fill )
    {
        std::fill( to_fill.begin(), to_fill.end(), NO_ID ) ;
    }

    void DuplicateInterfaceBuilder::homogenize_normal_orientation_surface_all_interfaces()
    {
        std::vector< index_t > surfaces_to_inverse_normals ;
        surfaces_to_inverse_normals.reserve( model_.nb_surfaces() ) ;
        for( index_t interface_itr = 0; interface_itr < model_.nb_interfaces();
            ++interface_itr ) {
            const GeoModelElement& cur_interface = model_.one_interface(
                interface_itr ) ;
            if( !GeoModelElement::is_fault( cur_interface.geological_feature() ) ) {
                continue ;
            }
            save_normals_on_one_old_interface( cur_interface ) ;
            homogenize_normal_orientation_surface_one_interface( cur_interface,
                surfaces_to_inverse_normals ) ;
        }

        for( index_t to_inverse_itr = 0;
            to_inverse_itr < surfaces_to_inverse_normals.size(); ++to_inverse_itr ) {
            const Surface& surface = model_.surface(
                surfaces_to_inverse_normals[to_inverse_itr] ) ;
            GEO::invert_normals( surface.mesh() ) ;
        }
        recompute_geomodel_mesh() ;
    }

    void DuplicateInterfaceBuilder::homogenize_normal_orientation_surface_one_interface(
        const GeoModelElement& fault_interface,
        std::vector< index_t >& surfaces_to_inverse_normals )
    {
        ringmesh_assert( fault_interface.type() == GME::INTERFACE ) ;
        ringmesh_assert( GeoModelElement::is_fault( fault_interface.geological_feature() ) ) ;
        ringmesh_assert( fault_interface.nb_children() >= 1 ) ;
        if( fault_interface.nb_children() == 1 ) {
            return ;
        }

        std::vector< bool > already_seen( model_.nb_surfaces(), false ) ;

        const Surface& first_child = model_.surface(
            fault_interface.child( 0 ).index() ) ;
        already_seen[first_child.index()] = true ;
        /// @todo for now I assume that all the surfaces of the fault interface
        /// form one connected component. What happens in case of a fault cut into
        /// 2 distinct parts by another fault??? A trick would be to define each part
        /// as a different fault interface...
        homogenize_surfaces_around_surface( fault_interface, first_child,
            already_seen, surfaces_to_inverse_normals ) ;
    }

    void DuplicateInterfaceBuilder::homogenize_surfaces_around_surface(
        const GeoModelElement& fault_interface,
        const Surface& first_child,
        std::vector< bool >& already_seen,
        std::vector< index_t >& surfaces_to_inverse_normals )
    {
        for( index_t line_boundary_itr = 0;
            line_boundary_itr < first_child.nb_boundaries(); ++line_boundary_itr ) {
            const GeoModelElement& cur_line_boun = first_child.boundary(
                line_boundary_itr ) ;
            ringmesh_assert(cur_line_boun.type() == GME::LINE) ;
            for( index_t surf_in_boun_itr = 0;
                surf_in_boun_itr < cur_line_boun.nb_in_boundary();
                ++surf_in_boun_itr ) {
                const GeoModelElement& cur_in_boun = cur_line_boun.in_boundary(
                    surf_in_boun_itr ) ;
                ringmesh_assert( cur_in_boun.type() == GME::SURFACE ) ;
                const Surface& cur_surf_in_boun = model_.surface(
                    cur_in_boun.index() ) ;
                if( cur_surf_in_boun.index() == first_child.index() ) {
                    continue ;
                }
                if( !does_surface_belong_to_interface( cur_surf_in_boun,
                    fault_interface ) ) {
                    continue ;
                }

                if( already_seen[cur_in_boun.index()] ) {
                    continue ;
                }

                const Line& cur_line = model_.line( cur_line_boun.index() ) ;
                ringmesh_assert(cur_line.nb_vertices()>0) ;
                index_t first_line_vertex_id_in_gmm = cur_line.model_vertex_id( 0 ) ;
                const std::vector< GMEVertex >& gme_vertices =
                    model_.mesh.vertices.gme_vertices(
                        first_line_vertex_id_in_gmm ) ;
                index_t v_id_in_first_surf = NO_ID ;
                bool v_id_in_second_surf = NO_ID ;
                for( index_t gme_vertex_itr = 0;
                    gme_vertex_itr < gme_vertices.size(); ++gme_vertex_itr ) {
                    const GMEVertex& cur_gme_vertex = gme_vertices[gme_vertex_itr] ;
                    if( cur_gme_vertex.gme_id.type != GME::SURFACE ) {
                        continue ;
                    }
                    if( cur_gme_vertex.gme_id.index == first_child.index() ) {
                        v_id_in_first_surf = cur_gme_vertex.v_id ;
                    } else if( cur_gme_vertex.gme_id.index
                        == cur_surf_in_boun.index() ) {
                        v_id_in_second_surf = cur_gme_vertex.v_id ;
                    } else {
                        continue ;
                    }

                    if( v_id_in_first_surf != NO_ID
                        && v_id_in_second_surf != NO_ID ) {
                        break ;
                    }
                }

                ringmesh_assert(
                    v_id_in_first_surf != NO_ID && v_id_in_second_surf != NO_ID ) ;

                vec3 first_normal = get_normal_on_surface_vertex( first_child,
                    v_id_in_first_surf ) ;
                vec3 second_normal = get_normal_on_surface_vertex( cur_surf_in_boun,
                    v_id_in_second_surf ) ;
                // As the normals are computed on different surfaces, they cannot
                // be exactly the same. However they should be very close (0.1 here is arbitrary).
                ringmesh_assert( std::abs( first_normal.length() - 1 ) < epsilon ) ;
                ringmesh_assert( std::abs( second_normal.length() - 1 ) < epsilon ) ;
                double dot_product = GEO::dot( first_normal, second_normal ) ;
                ringmesh_assert( std::abs( std::abs( dot_product ) -1 ) < 1e-1 ) ;
                if( dot_product < 0 ) {
                    surfaces_to_inverse_normals.push_back(
                        cur_surf_in_boun.index() ) ;
                    inverse_normal_attribute_one_surface( cur_surf_in_boun ) ;
                    update_region_polarity( cur_surf_in_boun ) ;
                }
                already_seen[cur_surf_in_boun.index()] = true ;
                homogenize_surfaces_around_surface( fault_interface,
                    cur_surf_in_boun, already_seen, surfaces_to_inverse_normals ) ;
                break ;
            }
        }
    }

    void DuplicateInterfaceBuilder::update_region_polarity( const Surface& surface )
    {
        index_t nb_in_boundaries = surface.nb_in_boundary() ;
        ringmesh_assert( nb_in_boundaries == 1 || nb_in_boundaries == 2 ) ;
        if( nb_in_boundaries == 1 ) {
            const GeoModelElement& reg_gme = surface.in_boundary( 0 ) ;
            index_t boundary_id = find_local_boundary_id( reg_gme, surface ) ;
            Region& reg = const_cast< Region& >( model_.region( reg_gme.index() ) ) ;
            bool cur_side = reg.side( boundary_id ) ;
            set_region_side( reg, boundary_id, !cur_side ) ;

            // Add to universe (other side of the surface)
            boundary_id = find_local_boundary_id( model_.universe(), surface ) ;
            ringmesh_assert(model_.universe().side(boundary_id) != cur_side) ;
            set_region_side( const_cast< Region& >( model_.universe() ), boundary_id,
                cur_side ) ;
        } else {
            const GeoModelElement& reg_gme1 = surface.in_boundary( 0 ) ;
            index_t boundary_id1 = find_local_boundary_id( reg_gme1, surface ) ;
            Region& reg1 = const_cast< Region& >( model_.region( reg_gme1.index() ) ) ;
            bool cur_side1 = reg1.side( boundary_id1 ) ;

            const GeoModelElement& reg_gme2 = surface.in_boundary( 1 ) ;
            index_t boundary_id2 = find_local_boundary_id( reg_gme2, surface ) ;
            Region& reg2 = const_cast< Region& >( model_.region( reg_gme2.index() ) ) ;
            bool cur_side2 = reg2.side( boundary_id2 ) ;
            // to check. if reg1 == reg2 no side in particular.

            if( reg1.index() != reg2.index() ) {
                ringmesh_assert( cur_side1 != cur_side2 ) ;
                set_region_side( reg1, boundary_id1, !cur_side1 ) ;
                set_region_side( reg2, boundary_id2, !cur_side2 ) ;
            } else {
                // Same region on the both side of the interface.
                boundary_id2 = find_second_local_boundary_id( reg1, surface ) ;
                ringmesh_assert( boundary_id2 != NO_ID ) ;
                cur_side2 = reg1.side( boundary_id2 ) ;
                ringmesh_assert( cur_side1 != cur_side2 ) ;
                set_region_side( reg1, boundary_id1, !cur_side1 ) ;
                set_region_side( reg1, boundary_id2, !cur_side2 ) ;
            }
        }
    }

    vec3 DuplicateInterfaceBuilder::get_normal_on_surface_vertex(
        const Surface& surface,
        index_t vertex_id_on_surface ) const
    {
        ringmesh_assert( vertex_id_on_surface < surface.nb_vertices() ) ;
        GEO::AttributesManager& att_mgr = surface.mesh().vertices.attributes() ;
        GEO::Attribute< double > normal_att_x( att_mgr, "normal_attr_x" ) ;
        GEO::Attribute< double > normal_att_y( att_mgr, "normal_attr_y" ) ;
        GEO::Attribute< double > normal_att_z( att_mgr, "normal_attr_z" ) ;
        return vec3( normal_att_x[vertex_id_on_surface],
            normal_att_y[vertex_id_on_surface], normal_att_z[vertex_id_on_surface] ) ;
    }

    void DuplicateInterfaceBuilder::inverse_normal_attribute_one_surface(
        const Surface& surface ) const
    {
        GEO::AttributesManager& att_mgr = surface.mesh().vertices.attributes() ;
        GEO::Attribute< double > normal_att_x( att_mgr, "normal_attr_x" ) ;
        GEO::Attribute< double > normal_att_y( att_mgr, "normal_attr_y" ) ;
        GEO::Attribute< double > normal_att_z( att_mgr, "normal_attr_z" ) ;
        for( index_t vertex_id_on_surface = 0;
            vertex_id_on_surface < surface.nb_vertices(); ++vertex_id_on_surface ) {
            normal_att_x[vertex_id_on_surface] *= -1 ;
            normal_att_y[vertex_id_on_surface] *= -1 ;
            normal_att_z[vertex_id_on_surface] *= -1 ;
        }
    }

    /// @todo COPY PASTE OF save_normals_on_one_new_interface
    /// @todo TO REFACTORE
    void DuplicateInterfaceBuilder::save_normals_on_one_old_interface(
        const GeoModelElement& interface_gme ) const
    {
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
            ringmesh_assert(cur_child.type() == GME::SURFACE) ;
            // As the loop begins at the first new interface, no surface
            // met in this loop should be to delete.
            const Surface& cur_surface = model_.surface( cur_child.index() ) ;
            save_normal_on_one_surface( cur_surface ) ;
        }
    }

    /// @todo could be inside an unamed namespace.
    /// @todo could a more generic function (is_child_of ?)
    /// @todo such function exists?
    bool DuplicateInterfaceBuilder::does_surface_belong_to_interface(
        const Surface& surface,
        const GeoModelElement& interface ) const
    {
        ringmesh_assert( interface.type() == GME::INTERFACE ) ;
        for( index_t child_itr = 0; child_itr < interface.nb_children();
            ++child_itr ) {
            const GeoModelElement& cur_child = interface.child( child_itr ) ;
            ringmesh_assert( cur_child.type() == GME::SURFACE ) ;
            if( cur_child.index() == surface.index() ) {
                return true ;
            }
        }
        return false ;
    }

    void DuplicateInterfaceBuilder::duplicate_fault_network()
    {
        if( model_.nb_regions() == 0 ) {
            throw RINGMeshException( "Duplication",
                "The model must contain at least one region." ) ;
        }
        all_meshed_ = model_.region( 0 ).is_meshed() ;
        for( index_t reg_itr = 1; reg_itr < model_.nb_regions(); ++reg_itr ) {
            if( model_.region( reg_itr ).is_meshed() != all_meshed_ ) {
                throw RINGMeshException( "Duplication",
                    "The regions must be all meshed or all unmeshed." ) ;
            }
        }

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

        homogenize_normal_orientation_surface_all_interfaces() ;

        index_t nb_faults = 0 ;
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
            ++nb_faults ;

            // Delete of the interface (will be replaced by a custom interface with
            // side informations)
            to_erase_by_type[GME::INTERFACE][cur_interface.index()] = NO_ID ;
            get_new_surfaces( cur_interface, to_erase_by_type ) ;
        }

        if( nb_faults == 0 ) {
            std::string message = "There is no fault inside the model.\n" ;
            message += "Assign your fault interfaces to fault geological features." ;
            throw std::runtime_error( message ) ;
        }

//        recompute_geomodel_mesh() ; /// @todo necessary???

        // The new surfaces has no internal borders made by the intersections
        // with other new surfaces. Some of these may cut some of these
        // new surfaces by two (two connected components). In this last case,
        // two new surfaces are built to replace the new uncut surface.
        /*mutual_cut_between_new_merged_surfaces( nb_initial_interfaces,
         to_erase_by_type ) */;

        // compute translation vectors
        if( all_meshed_ ) {
            ringmesh_assert( model_.nb_interfaces() - nb_initial_interfaces >= 2 ) ;
            for( index_t new_interface_itr = nb_initial_interfaces;
                new_interface_itr < model_.nb_interfaces(); ++new_interface_itr ) {
                ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] != NO_ID ) ;
                ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] == 0 ) ;

                const GeoModelElement& interface_gme = model_.one_interface(
                    new_interface_itr ) ;
                save_normals_on_one_new_interface( to_erase_by_type,
                    interface_gme ) ;
            }
            define_global_motion_relation( to_erase_by_type ) ;
            compute_translation_vectors_duplicated_fault_network_surfaces_and_regions(
                nb_initial_interfaces, to_erase_by_type ) ;
        } else {
            // Only for models with surfaces (region parts to remove?)
            compute_translation_vectors_duplicated_fault_network(
                nb_initial_interfaces, to_erase_by_type ) ;
        }
        GEO::mesh_save( to_debug66, "to_debug66.mesh" ) ;
        // set no translation on fault real extension (only on fault ending inside
        // the model).
        set_no_displacement_on_fault_real_extension( to_erase_by_type ) ;
        // apply translation
        translate_duplicated_fault_network( to_erase_by_type ) ;

        // Put here for new.
        /// @todo if all the lines are removed, is it still necessary to fill the new
        /// interface children with them? They should be recomputed with
        /// build_lines_and_corners_from_surfaces.
        fill_vect_with_NO_ID( to_erase_by_type[GME::CORNER] ) ;
        fill_vect_with_NO_ID( to_erase_by_type[GME::LINE] ) ;
        fill_vect_with_NO_ID( to_erase_by_type[GME::CONTACT] ) ;

        recompute_geomodel_mesh() ;
        delete_elements( to_erase_by_type ) ;
        recompute_geomodel_mesh() ;
        build_lines_and_corners_from_surfaces() ;
        complete_element_connectivity() ;
        build_contacts() ;
    }

    void DuplicateInterfaceBuilder::define_global_motion_relation(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        initialize_gme_vertices_links( to_erase_by_type ) ;
        do_define_motion_relation( to_erase_by_type ) ;
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation_on_not_voi_surface(
        const Surface& cur_surface,
        index_t surf_facet_itr,
        const Region& reg1,
        const std::vector< index_t >& colocated_facets_reg1,
        const vec3& facet_bary,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1,
        std::vector< bool >& surf_vertex_visited )
    {
        ringmesh_assert( cur_surface.nb_in_boundary() == 2 ) ;
        const GME::gme_t& reg2_gme_t = cur_surface.in_boundary_gme( 1 ) ;
        ringmesh_assert( reg2_gme_t.type == GME::REGION ) ;
        const Region& reg2 = model_.region( reg2_gme_t.index ) ;
        ringmesh_assert( reg2.is_meshed() ) ;
        ColocaterANN reg2_ann( reg2.mesh(), ColocaterANN::FACETS ) ;
        GEO::Attribute< index_t > id_in_link_vector_reg2(
            reg2.mesh().vertices.attributes(), "id_in_link_vector" ) ;

        // Working on horizon not voi
        ringmesh_assert( reg1.index() != reg2.index() ) ;
        ringmesh_assert( colocated_facets_reg1.size() == 1 ) ;
        std::vector< index_t > colocated_facets_reg2 ;
        colocated_facets_reg2.reserve( 1 ) ;
        reg2_ann.get_colocated( facet_bary, colocated_facets_reg2 ) ;
        ringmesh_assert( colocated_facets_reg2.size() == 1 ) ;

        for( index_t v_in_surf_facet = 0;
            v_in_surf_facet < cur_surface.nb_vertices_in_facet( surf_facet_itr );
            ++v_in_surf_facet ) {
            index_t v_id_in_surf = cur_surface.surf_vertex_id( surf_facet_itr,
                v_in_surf_facet ) ;
            if( surf_vertex_visited[v_id_in_surf] ) {
                continue ;
            }

            const index_t surf_v_id_in_gmm = cur_surface.model_vertex_id(
                v_id_in_surf ) ;

            surf_vertex_visited[v_id_in_surf] = true ;
            index_t v_id_in_reg1 = NO_ID ;
            index_t v_id_in_reg2 = NO_ID ;

            for( index_t v_in_reg1_itr = 0;
                v_in_reg1_itr
                    < reg1.mesh().facets.nb_vertices( colocated_facets_reg1[0] );
                ++v_in_reg1_itr ) {
                index_t reg1_v_id_in_gmme = reg1.mesh().facets.vertex(
                    colocated_facets_reg1[0], v_in_reg1_itr ) ;
                if( reg1.model_vertex_id( reg1_v_id_in_gmme ) == surf_v_id_in_gmm ) {
                    v_id_in_reg1 = reg1_v_id_in_gmme ;
                    break ;
                }
            }
            ringmesh_assert( v_id_in_reg1 != NO_ID ) ;

            for( index_t v_in_reg2_itr = 0;
                v_in_reg2_itr
                    < reg2.mesh().facets.nb_vertices( colocated_facets_reg2[0] );
                ++v_in_reg2_itr ) {
                index_t reg2_v_id_in_gmme = reg2.mesh().facets.vertex(
                    colocated_facets_reg2[0], v_in_reg2_itr ) ;
                if( reg2.model_vertex_id( reg2_v_id_in_gmme ) == surf_v_id_in_gmm ) {
                    v_id_in_reg2 = reg2_v_id_in_gmme ;
                    break ;
                }
            }
            ringmesh_assert( v_id_in_reg2 != NO_ID ) ;

            // LINKING
            index_t link_id_surf = id_in_link_vector_surf[v_id_in_surf] ;
            index_t link_id_reg1 = id_in_link_vector_reg1[v_id_in_reg1] ;
            index_t link_id_reg2 = id_in_link_vector_reg2[v_id_in_reg2] ;
            link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg1 ) ;
            link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg2 ) ;
        }
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation_on_voi_surface(
        const Surface& cur_surface,
        index_t surf_facet_itr,
        const Region& reg1,
        const std::vector< index_t >& colocated_facets_reg1,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1,
        std::vector< bool >& surf_vertex_visited )
    {
        if( colocated_facets_reg1.size() == 1 ) {
            do_define_motion_relation_on_voi_fault_crossing_or_model_boundary_or_voi_horizon(
                cur_surface, surf_facet_itr, reg1, colocated_facets_reg1,
                id_in_link_vector_surf, id_in_link_vector_reg1,
                surf_vertex_visited ) ;
        } else {
            do_define_motion_relation_on_voi_fault_not_crossing( cur_surface,
                surf_facet_itr, reg1, colocated_facets_reg1, id_in_link_vector_surf,
                id_in_link_vector_reg1, surf_vertex_visited ) ;
        }
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation_on_voi_fault_crossing_or_model_boundary_or_voi_horizon(
        const Surface& cur_surface,
        index_t surf_facet_itr,
        const Region& reg1,
        const std::vector< index_t >& colocated_facets_reg1,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1,
        std::vector< bool >& surf_vertex_visited )
    {
        ringmesh_assert( colocated_facets_reg1.size() == 1 ) ;
        // It is not the same region on the other side of the fault.
        // Or it is a model boundary so defining universe, no cell
        // on the other side.

        for( index_t v_in_surf_facet = 0;
            v_in_surf_facet < cur_surface.nb_vertices_in_facet( surf_facet_itr );
            ++v_in_surf_facet ) {
            index_t v_id_in_surf = cur_surface.surf_vertex_id( surf_facet_itr,
                v_in_surf_facet ) ;
            if( surf_vertex_visited[v_id_in_surf] ) {
                continue ;
            }

            const index_t surf_v_id_in_gmm = cur_surface.model_vertex_id(
                v_id_in_surf ) ;

            surf_vertex_visited[v_id_in_surf] = true ;
            index_t v_id_in_reg1 = NO_ID ;

            for( index_t v_in_reg1_itr = 0;
                v_in_reg1_itr
                    < reg1.mesh().facets.nb_vertices( colocated_facets_reg1[0] );
                ++v_in_reg1_itr ) {
                index_t reg1_v_id_in_gmme = reg1.mesh().facets.vertex(
                    colocated_facets_reg1[0], v_in_reg1_itr ) ;
                if( reg1.model_vertex_id( reg1_v_id_in_gmme ) == surf_v_id_in_gmm ) {
                    v_id_in_reg1 = reg1_v_id_in_gmme ;
                    break ;
                }
            }

            ringmesh_assert( v_id_in_reg1 != NO_ID ) ;

            // LINKING
            index_t link_id_surf = id_in_link_vector_surf[v_id_in_surf] ;
            index_t link_id_reg1 = id_in_link_vector_reg1[v_id_in_reg1] ;
            link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg1 ) ;
        }
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation_on_voi_fault_not_crossing(
        const Surface& cur_surface,
        index_t surf_facet_itr,
        const Region& reg1,
        const std::vector< index_t >& colocated_facets_reg1,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1,
        std::vector< bool >& surf_vertex_visited )
    {
        // It is a fault with the region on the other side.
        ringmesh_assert( GME::is_fault(
                cur_surface.parent().geological_feature()) ) ;
        ringmesh_assert( colocated_facets_reg1.size() == 2 ) ;

        // Try to find which region facet is realy on the surface.
        index_t v0_id_in_surf = cur_surface.surf_vertex_id( surf_facet_itr, 0 ) ;
        const index_t surf_v0_id_in_gmm = cur_surface.model_vertex_id(
            v0_id_in_surf ) ;

        index_t v0_id_in_reg1 = NO_ID ;
        for( index_t v_in_reg1_itr = 0;
            v_in_reg1_itr
                < reg1.mesh().facets.nb_vertices( colocated_facets_reg1[0] );
            ++v_in_reg1_itr ) {

            index_t reg1_v0_id_in_gmme = reg1.mesh().facets.vertex(
                colocated_facets_reg1[0], v_in_reg1_itr ) ;
            if( reg1.model_vertex_id( reg1_v0_id_in_gmme ) == surf_v0_id_in_gmm ) {
                v0_id_in_reg1 = reg1_v0_id_in_gmme ;
                break ;
            }
        }
        ringmesh_assert( v0_id_in_reg1 != NO_ID ) ;

        const vec3 local_translation_normal = get_local_translation_normal(
            cur_surface, v0_id_in_surf ) ;
        bool is_first = is_region_on_right_side_of_sided_interface( reg1.index(),
            local_translation_normal, v0_id_in_reg1,
            cur_surface.vertex( v0_id_in_surf ) ) ;

        index_t right_colocated_facet =
            is_first ? colocated_facets_reg1[0] : colocated_facets_reg1[1] ;

        for( index_t v_in_surf_facet = 0;
            v_in_surf_facet < cur_surface.nb_vertices_in_facet( surf_facet_itr );
            ++v_in_surf_facet ) {
            index_t v_id_in_surf = cur_surface.surf_vertex_id( surf_facet_itr,
                v_in_surf_facet ) ;
            if( surf_vertex_visited[v_id_in_surf] ) {
                continue ;
            }

            const index_t surf_v_id_in_gmm = cur_surface.model_vertex_id(
                v_id_in_surf ) ;

            surf_vertex_visited[v_id_in_surf] = true ;
            index_t v_id_in_reg1 = NO_ID ;

            for( index_t v_in_reg1_itr = 0;
                v_in_reg1_itr
                    < reg1.mesh().facets.nb_vertices( right_colocated_facet );
                ++v_in_reg1_itr ) {
                index_t reg1_v_id_in_gmme = reg1.mesh().facets.vertex(
                    right_colocated_facet, v_in_reg1_itr ) ;
                if( reg1.model_vertex_id( reg1_v_id_in_gmme ) == surf_v_id_in_gmm ) {
                    v_id_in_reg1 = reg1_v_id_in_gmme ;
                    break ;
                }
            }

            ringmesh_assert( v_id_in_reg1 != NO_ID ) ;

            // LINKING
            index_t link_id_surf = id_in_link_vector_surf[v_id_in_surf] ;
            index_t link_id_reg1 = id_in_link_vector_reg1[v_id_in_reg1] ;
            link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg1 ) ;
        }
    }

    void DuplicateInterfaceBuilder::link_surf_vertex_id_to_reg_vertex_id(
        index_t link_id_surf,
        index_t link_id_reg )
    {
        gme_vertices_links_[link_id_surf]->add_linked_gme_vertex( link_id_reg ) ;
        gme_vertices_links_[link_id_reg]->add_linked_gme_vertex( link_id_surf ) ;
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        recompute_geomodel_mesh() ; /// @todo check if it is really necessary

        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            if( to_erase_by_type[GME::SURFACE][surf_itr] == NO_ID ) {
                continue ;
            }
            ringmesh_assert( to_erase_by_type[GME::SURFACE][surf_itr] == 0 ) ;

            // For each vertex of the surf, I must find all the region vertex
            // ids which are linked to it, on one side for a fault
            // (and model boundary) and on the both sides for horizons.
            // In theory, on the both sides of a horizon there are two different
            // regions, else it is not a horizon.
            const Surface& cur_surface = model_.surface( surf_itr ) ;
            GEO::Attribute< index_t > id_in_link_vector_surf(
                cur_surface.mesh().vertices.attributes(), "id_in_link_vector" ) ;

            const index_t nb_in_boundaries = cur_surface.nb_in_boundary() ;

            bool is_horizon_not_voi = GME::is_stratigraphic_limit(
                cur_surface.parent().geological_feature() )
                && !cur_surface.is_on_voi() ;

            ringmesh_assert( nb_in_boundaries == 1 || nb_in_boundaries == 2 ) ;

            const GME::gme_t& reg1_gme_t = cur_surface.in_boundary_gme( 0 ) ;
            ringmesh_assert( reg1_gme_t.type == GME::REGION ) ;
            const Region& reg1 = model_.region( reg1_gme_t.index ) ;
            ringmesh_assert( reg1.is_meshed() ) ;
            ColocaterANN reg1_ann( reg1.mesh(), ColocaterANN::FACETS ) ;
            GEO::Attribute< index_t > id_in_link_vector_reg1(
                reg1.mesh().vertices.attributes(), "id_in_link_vector" ) ;

            std::vector< bool > surf_vertex_visited( cur_surface.nb_cells(),
                false ) ;

            for( index_t surf_facet_itr = 0; surf_facet_itr < cur_surface.nb_cells();
                ++surf_facet_itr ) {

                const vec3 facet_bary = GEO::Geom::mesh_facet_center(
                    cur_surface.mesh(), surf_facet_itr ) ;
                std::vector< index_t > colocated_facets_reg1 ;
                colocated_facets_reg1.reserve( 2 ) ;
                reg1_ann.get_colocated( facet_bary, colocated_facets_reg1 ) ;

                ringmesh_assert( colocated_facets_reg1.size() == 1 ||
                    colocated_facets_reg1.size() == 2 ) ;
                if( is_horizon_not_voi ) {
                    ringmesh_assert( nb_in_boundaries == 2 ) ;
                    ringmesh_assert( colocated_facets_reg1.size() == 1 ) ;
                    do_define_motion_relation_on_not_voi_surface( cur_surface,
                        surf_facet_itr, reg1, colocated_facets_reg1, facet_bary,
                        id_in_link_vector_surf, id_in_link_vector_reg1,
                        surf_vertex_visited ) ;
                } else {
                    ringmesh_assert( nb_in_boundaries == 1 ) ;
                    do_define_motion_relation_on_voi_surface( cur_surface,
                        surf_facet_itr, reg1, colocated_facets_reg1,
                        id_in_link_vector_surf, id_in_link_vector_reg1,
                        surf_vertex_visited ) ;
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::initialize_gme_vertices_links(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        index_t count = 0 ;
        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            if( to_erase_by_type[GME::SURFACE][surf_itr] == NO_ID ) {
                continue ;
            }
            ringmesh_assert( to_erase_by_type[GME::SURFACE][surf_itr] == 0 ) ;
            count += model_.surface( surf_itr ).nb_vertices() ;
        }

        if( all_meshed_ ) {

            for( index_t reg_itr = 0; reg_itr < model_.nb_regions(); ++reg_itr ) {
                if( to_erase_by_type[GME::REGION][reg_itr] == NO_ID ) {
                    continue ;
                }
                ringmesh_assert( to_erase_by_type[GME::REGION][reg_itr] == 0 ) ;
                ringmesh_assert(model_.region( reg_itr ).is_meshed()) ;
                count += model_.region( reg_itr ).nb_vertices() ;
            }
        }

        gme_vertices_links_.reserve( count ) ;

        for( index_t surf_itr = 0; surf_itr < model_.nb_surfaces(); ++surf_itr ) {
            if( to_erase_by_type[GME::SURFACE][surf_itr] == NO_ID ) {
                continue ;
            }
            ringmesh_assert( to_erase_by_type[GME::SURFACE][surf_itr] == 0 ) ;

            GEO::Attribute< index_t > id_in_link_vector(
                model_.surface( surf_itr ).mesh().vertices.attributes(),
                "id_in_link_vector" ) ;

            for( index_t v_itr = 0; v_itr < model_.surface( surf_itr ).nb_vertices();
                ++v_itr ) {
                id_in_link_vector[v_itr] = gme_vertices_links_.size() ;
                gme_vertices_links_.push_back(
                    new GMEVertexLink(
                        GMEVertex( GME::gme_t( GME::SURFACE, surf_itr ), v_itr ),
                        model_, gme_vertices_links_ ) ) ;
            }
        }

        if( all_meshed_ ) {

            for( index_t reg_itr = 0; reg_itr < model_.nb_regions(); ++reg_itr ) {
                if( to_erase_by_type[GME::REGION][reg_itr] == NO_ID ) {
                    continue ;
                }
                ringmesh_assert( to_erase_by_type[GME::REGION][reg_itr] == 0 ) ;
                ringmesh_assert(model_.region( reg_itr ).is_meshed()) ;
                GEO::Attribute< index_t > id_in_link_vector(
                    model_.region( reg_itr ).mesh().vertices.attributes(),
                    "id_in_link_vector" ) ;
                for( index_t v_itr = 0;
                    v_itr < model_.region( reg_itr ).nb_vertices(); ++v_itr ) {
                    id_in_link_vector[v_itr] = gme_vertices_links_.size() ;
                    gme_vertices_links_.push_back(
                        new GMEVertexLink(
                            GMEVertex( GME::gme_t( GME::REGION, reg_itr ), v_itr ),
                            model_, gme_vertices_links_ ) ) ;
                }
            }
        }

        ringmesh_assert( gme_vertices_links_.size() == gme_vertices_links_.capacity() ) ;
        ringmesh_assert( gme_vertices_links_.size() == count ) ;
    }

    DuplicateInterfaceBuilder::GMEVertexLink::GMEVertexLink(
        const GMEVertex& gme_vertex,
        const GeoModel& model,
        const std::vector< GMEVertexLink* >& gme_vertices_links )
        :
            has_moved_( false ),
            model_( model ),
            gme_vertex_( gme_vertex ),
            gme_vertices_links_( gme_vertices_links )
    {

    }

    void DuplicateInterfaceBuilder::GMEVertexLink::displace(
        const vec3& displacement_vector )
    {
//#ifdef RINGMESH_DEBUG
        if( gme_vertex_.gme_id.type == GME::SURFACE ) {
            if( linked_gme_vertices_.size() >= 3 ) {
                DEBUG( "more than 2" ) ;
            }
            ringmesh_assert( linked_gme_vertices_.size() <= 2 ) ;
        } else {
            ringmesh_assert(gme_vertex_.gme_id.type == GME::REGION) ;
            if( linked_gme_vertices_.size() >= 4 ) {
                DEBUG( "more than 3" ) ;
                to_debug66.vertices.create_vertex(
                    model_.region( gme_vertex_.gme_id.index ).vertex(
                        gme_vertex_.v_id ).data() ) ;
            }

            ringmesh_assert( linked_gme_vertices_.size() <= 3 ) ;
        }

        for( index_t link_itr = 0; link_itr < linked_gme_vertices_.size();
            ++link_itr ) {
            for( index_t link_itr2 = 0; link_itr2 < linked_gme_vertices_.size();
                ++link_itr2 ) {
                if( link_itr == link_itr2 ) {
                    continue ;
                }
                if( linked_gme_vertices_[link_itr]
                    == linked_gme_vertices_[link_itr2] ) {
                    DEBUG( "several time" ) ;
                }
                ringmesh_assert( linked_gme_vertices_[link_itr] != linked_gme_vertices_[link_itr2] ) ;
            }
        }
//#endif

        if( has_moved_ ) {
            return ;
        }
        has_moved_ = true ;

        GEO::AttributesManager& att_mgr =
            model_.mesh_element( gme_vertex_.gme_id ).mesh().vertices.attributes() ;
        GEO::Attribute< double > translation_att_x( att_mgr, "translation_attr_x" ) ;
        translation_att_x[gme_vertex_.v_id] += displacement_vector.x ;
        GEO::Attribute< double > translation_att_y( att_mgr, "translation_attr_y" ) ;
        translation_att_y[gme_vertex_.v_id] += displacement_vector.y ;
        GEO::Attribute< double > translation_att_z( att_mgr, "translation_attr_z" ) ;
        translation_att_z[gme_vertex_.v_id] += displacement_vector.z ;

        for( index_t link_itr = 0; link_itr < linked_gme_vertices_.size();
            ++link_itr ) {
            gme_vertices_links_[linked_gme_vertices_[link_itr]]->displace(
                displacement_vector ) ;
        }
        has_moved_ = false ;
    }

    void DuplicateInterfaceBuilder::mutual_cut_between_new_merged_surfaces(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        ringmesh_assert( model_.nb_interfaces() - first_new_interface_index >= 2 ) ;
        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr < model_.nb_interfaces(); ++new_interface_itr ) {
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] != NO_ID ) ;
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] == 0 ) ;

            const GeoModelElement& interface_gme = model_.one_interface(
                new_interface_itr ) ;
            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GME::gme_t& cur_gme_t =
                    interface_gme.child( child_itr ).gme_id() ;
                ringmesh_assert( cur_gme_t.type == GME::SURFACE ) ;
                const Surface& cur_surface = model_.surface( cur_gme_t.index ) ;
                // ========= bad copy paste from geo model repair
                std::set< index_t > cutting_lines ;
                for( index_t l = 0; l < cur_surface.nb_boundaries(); ++l ) {
                    const Line& L = model_.line(
                        cur_surface.boundary_gme( l ).index ) ;
                    if( /*to_remove.count( L.gme_id() ) == 0 &&*/L.is_inside_border(
                        cur_surface ) ) {
                        cutting_lines.insert( L.index() ) ;
                    }
                }

                for( std::set< index_t >::iterator it = cutting_lines.begin();
                    it != cutting_lines.end(); ++it ) {
                    // Force the recomputing of the model vertices
                    // before performing the cut.
//                    model_.mesh.vertices.clear() ;
                    disconnect_surface_facets_along_line_edges(
                        const_cast< Surface& >( cur_surface ), model_.line( *it ) ) ;
                }

                for( std::set< index_t >::iterator it = cutting_lines.begin();
                    it != cutting_lines.end(); ++it ) {
                    exit( 0 ) ;
                    // Force the recomputing of the model vertices
                    // before performing the cut.
//                    model_.mesh.vertices.clear() ;
                    duplicate_surface_vertices_along_line_benjamin(
                        const_cast< Surface& >( cur_surface ), model_.line( *it ) ) ;
                    cur_surface.mesh().vertices.remove_isolated() ;
                }
                // ========= bad copy paste from geo model repair
            }
        }
    }

    void DuplicateInterfaceBuilder::get_new_surfaces(
        const GeoModelElement& interface_to_duplicate,
        std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        // minus = false, plus = true
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_minus ;
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_plus ;
        GME::gme_t interface_minus_gme_t = create_element( GME::INTERFACE ) ;
        set_element_name( interface_minus_gme_t,
            interface_to_duplicate.name() + "_side_minus" ) ;
        set_element_geol_feature( interface_minus_gme_t,
            interface_to_duplicate.geological_feature() ) ;
        GME::gme_t interface_plus_gme_t = create_element( GME::INTERFACE ) ;
        set_element_name( interface_plus_gme_t,
            interface_to_duplicate.name() + "_side_plus" ) ;
        set_element_geol_feature( interface_plus_gme_t,
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

                /// @todo it seems that this if statement is contained in the else
                /// case. To check and simplify if necessary.
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

        build_merged_surfaces( surfaces_boundary_regions_side_plus, "_plus",
            to_erase_by_type, interface_plus_gme_t, interface_to_duplicate ) ;
        build_merged_surfaces( surfaces_boundary_regions_side_minus, "_minus",
            to_erase_by_type, interface_minus_gme_t, interface_to_duplicate ) ;
    }

    void DuplicateInterfaceBuilder::build_merged_surfaces(
        const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
        const std::string& side_name,
        std::vector< std::vector< index_t > >& to_erase_by_type,
        const GME::gme_t& sided_interface_gme_t,
        const GeoModelElement& interface_to_duplicate )
    {
        for( std::map< index_t, std::vector< index_t > >::const_iterator map_itr =
            surfaces_boundary_regions.begin();
            map_itr != surfaces_boundary_regions.end(); ++map_itr ) {

            // first = line index in geomodel, second = count.
            std::map< index_t, index_t > all_surface_lines ;

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
                /// @todo this line part is necessary if used for the mutural cut
                /// else delete it.
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
            /*set_element_parent( new_surface_gme_t, sided_interface_gme_t ) ;
             add_element_child( sided_interface_gme_t, new_surface_gme_t ) ;
             */
            // Boundary information is necessary for get_local_translation_normal
            add_element_in_boundary( new_surface_gme_t,
                GME::gme_t( GME::REGION, region_index ) ) ;
            bool side = ( side_name == "_plus" ) ;
            add_element_boundary( GME::gme_t( GME::REGION, region_index ),
                new_surface_gme_t, side ) ;
            /*
             // Add to universe (other side of the surface)
             add_element_boundary( model_.universe().gme_id(), new_surface_gme_t,
             !side ) ;*/

//            to_erase_by_type[GME::SURFACE].push_back( 0 ) ;
            to_erase_by_type[GME::SURFACE].push_back( NO_ID ) ;

            add_fake_internal_boudnary_lines_to_merged_surface( all_surface_lines,
                side_name, sided_interface_gme_t, interface_to_duplicate,
                new_surface_gme_t, to_erase_by_type, region_index ) ;
        }
    }

    void DuplicateInterfaceBuilder::add_fake_internal_boudnary_lines_to_merged_surface(
        const std::map< index_t, index_t >& all_surface_lines,
        const std::string& side_name,
        const GME::gme_t& sided_interface_gme_t,
        const GeoModelElement& interface_to_duplicate,
        const GME::gme_t& new_surface_gme_t,
        std::vector< std::vector< index_t > >& to_erase_by_type,
        index_t region_index )
    {
        recompute_geomodel_mesh() ; // to take into account the new surface in gme_vertices.
        save_normal_on_one_surface( model_.surface( new_surface_gme_t.index ) ) ;
        // I DO THAT HERE BUT MAYBE IT IS SIMPLIER TO DO THAT AFTER HAVING
        // CREATED ALL THE NEW INTERFACES.
        // Add all the lines of the old surfaces, which served to build the
        // new merged surface, to the new merged surface. These lines will be
        // removed, they are here only to ease the mutual cut when the new
        // interfaces are all built (see mutual_cut_between_new_merged_surfaces)
        // IN FACT I ADD ONLY THE LINES RESPONSIBLE OF INTERNAL BORDER OR CUTTING
        // THE SURFACE INTO TWO SURFACES.
        for( std::map< index_t, index_t >::const_iterator all_surface_lines_itr =
            all_surface_lines.begin();
            all_surface_lines_itr != all_surface_lines.end();
            ++all_surface_lines_itr ) {
            if( all_surface_lines_itr->second == 1 ) {
                // Line on the border of the new surface (not internal border)
                continue ;
            }
            ringmesh_assert( all_surface_lines_itr->second != 0 ) ;
            const Line& cur_line = model_.line( all_surface_lines_itr->first ) ;
            // As the line is not on the border, its number of in boundaries
            // superior to 1 strictly.
            ringmesh_assert( cur_line.nb_in_boundary() > 1 ) ;
            bool good_line = false ;
            for( index_t in_boun_itr = 0; in_boun_itr < cur_line.nb_in_boundary();
                ++in_boun_itr ) {
                const GeoModelElement& cur_in_boun_gme = cur_line.in_boundary(
                    in_boun_itr ) ;
                ringmesh_assert( cur_in_boun_gme.type() == GME::SURFACE ) ;

                if( !cur_in_boun_gme.has_parent() ) {
                    continue ;
                }

                if( does_surface_belong_to_interface(
                    model_.surface( cur_in_boun_gme.index() ),
                    model_.one_interface( sided_interface_gme_t.index ) ) ) {
                    continue ;
                }

                // minus side is done before the plus side
                if( side_name == "_minus" ) {

                    ringmesh_assert(sided_interface_gme_t.index + 1 < model_.nb_interfaces()) ;
                    const GeoModelElement& plus_side_gme = model_.one_interface(
                        sided_interface_gme_t.index + 1 ) ;

                    if( does_surface_belong_to_interface(
                        model_.surface( cur_in_boun_gme.index() ),
                        plus_side_gme ) ) {
                        continue ;
                    }
                }

                // Check if the surface is not from the old interface to
                // duplicate.
                if( does_surface_belong_to_interface(
                    model_.surface( cur_in_boun_gme.index() ),
                    interface_to_duplicate ) ) {
                    continue ;
                }

                GME::GEOL_FEATURE parent_geol_feature =
                    cur_in_boun_gme.parent().geological_feature() ;
                if( !GME::is_fault( parent_geol_feature ) ) {
                    continue ;
                }

                // Check if the found fault is on the right side
                ringmesh_assert( cur_line.nb_vertices() > 0 ) ;
                index_t first_vertex_id_in_gmm = cur_line.model_vertex_id( 0 ) ;
                const std::vector< GMEVertex >& gme_vertices =
                    model_.mesh.vertices.gme_vertices( first_vertex_id_in_gmm ) ;
                ringmesh_assert( !gme_vertices.empty() ) ; /// @todo I think that this assert may be more restrictive
                index_t vertex_id_in_new_surface = NO_ID ;
                index_t vertex_id_in_found_surface = NO_ID ;
                for( index_t gme_vertices_itr = 0;
                    gme_vertices_itr < gme_vertices.size(); ++gme_vertices_itr ) {
                    const GMEVertex& cur_gme_vertex = gme_vertices[gme_vertices_itr] ;
                    if( cur_gme_vertex.gme_id.type != GME::SURFACE ) {
                        continue ;
                    }
                    if( cur_gme_vertex.gme_id.index == new_surface_gme_t.index ) {
                        vertex_id_in_new_surface = cur_gme_vertex.v_id ;
                    } else if( cur_gme_vertex.gme_id.index
                        == cur_in_boun_gme.index() ) {
                        vertex_id_in_found_surface = cur_gme_vertex.v_id ;
                    }
                    if( vertex_id_in_new_surface != NO_ID
                        && vertex_id_in_found_surface != NO_ID ) {
                        break ;
                    }
                }
                ringmesh_assert(
                    vertex_id_in_new_surface != NO_ID
                    && vertex_id_in_found_surface != NO_ID ) ;

                // On the found surface, the normal has been computed before
                // for the homogenization of normals.

                const vec3 local_translation_normal = get_local_translation_normal(
                    model_.surface( new_surface_gme_t.index ),
                    vertex_id_in_new_surface ) ;

                vec3 vertex_pos = model_.surface( new_surface_gme_t.index ).vertex(
                    vertex_id_in_new_surface ) ;
                if( !is_surface_on_right_side_of_sided_interface(
                    cur_in_boun_gme.index(), local_translation_normal,
                    vertex_id_in_found_surface, vertex_pos ) ) {
                    continue ;
                }

                add_element_boundary( new_surface_gme_t, cur_line.gme_id() ) ;
                // Adds twice in boundary for internal border.
                add_element_in_boundary( cur_line.gme_id(), new_surface_gme_t ) ;
                add_element_in_boundary( cur_line.gme_id(), new_surface_gme_t ) ;

                good_line = true ;
                break ;
            }
            if( !good_line ) {
                continue ;
            }
        }

        split_merged_surface( new_surface_gme_t, side_name, sided_interface_gme_t,
            to_erase_by_type, region_index ) ;
    }

    void DuplicateInterfaceBuilder::split_merged_surface(
        const GME::gme_t& new_surface_gme_t,
        const std::string& side_name,
        const GME::gme_t& sided_interface_gme_t,
        std::vector< std::vector< index_t > >& to_erase_by_type,
        index_t region_index )
    {
        const Surface& cur_surface = model_.surface( new_surface_gme_t.index ) ;
        // ========= bad copy paste from geo model repair
        std::set< index_t > cutting_lines ;
        for( index_t l = 0; l < cur_surface.nb_boundaries(); ++l ) {
            const Line& L = model_.line( cur_surface.boundary_gme( l ).index ) ;
            if( /*to_remove.count( L.gme_id() ) == 0 &&*/L.is_inside_border(
                cur_surface ) ) {
                cutting_lines.insert( L.index() ) ;
            }
        }

        for( std::set< index_t >::iterator it = cutting_lines.begin();
            it != cutting_lines.end(); ++it ) {
            // Force the recomputing of the model vertices
            // before performing the cut.
            //                    model_.mesh.vertices.clear() ;
            disconnect_surface_facets_along_line_edges(
                const_cast< Surface& >( cur_surface ), model_.line( *it ) ) ;
        }

        ringmesh_assert( new_surface_gme_t.type == GME::SURFACE ) ;
        const GEO::Mesh& surface_mesh =
            model_.surface( new_surface_gme_t.index ).mesh() ;
        GEO::vector< index_t > components ;
        index_t nb_connected_components = GEO::get_connected_components(
            surface_mesh, components ) ;
        if( nb_connected_components == 1 ) {
            /// @todo HANDLE THE INTERNAL BORDER
            set_element_parent( new_surface_gme_t, sided_interface_gme_t ) ;
            add_element_child( sided_interface_gme_t, new_surface_gme_t ) ;
            // boundary informations are defined in build_merged_surfaces
            // expected for the universe.
            /*add_element_in_boundary( new_surface_gme_t,
             GME::gme_t( GME::REGION, region_index ) ) ;*/
            bool side = ( side_name == "_plus" ) ;
            /*add_element_boundary( GME::gme_t( GME::REGION, region_index ),
             new_surface_gme_t, side ) ;*/

            // Add to universe (other side of the surface)
            add_element_boundary( model_.universe().gme_id(), new_surface_gme_t,
                !side ) ;

            to_erase_by_type[GME::SURFACE][new_surface_gme_t.index] = 0 ;
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                cut_surface_by_line(
                    const_cast< Surface& >( model_.surface( new_surface_gme_t.index ) ),
                    model_.line( *it ) ) ;
            }
            return ;
        }
        ringmesh_assert( nb_connected_components != 0 ) ;

        // We create a new surface for each connected component
        // (since a RINGMesh::Surface has only one connected component).
        std::vector< std::vector< vec3 > > all_points( nb_connected_components,
            std::vector< vec3 >() ) ;
        for( index_t all_points_itr = 0; all_points_itr < all_points.size();
            ++all_points_itr ) {
            all_points[all_points_itr].reserve( surface_mesh.vertices.nb() ) ; // a little big but to avoid copy
        }
        for( index_t components_itr = 0; components_itr < components.size();
            ++components_itr ) {
            ringmesh_assert( components[components_itr] < all_points.size() ) ;
            for( index_t v_in_f_itr = 0;
                v_in_f_itr < surface_mesh.facets.nb_vertices( components_itr );
                ++v_in_f_itr ) {
                index_t v_id = surface_mesh.facets.vertex( components_itr,
                    v_in_f_itr ) ;
                all_points[components[components_itr]].push_back(
                    surface_mesh.vertices.point( v_id ) ) ;
            }
        }

        for( index_t all_points_itr = 0; all_points_itr < all_points.size();
            ++all_points_itr ) {
            MakeUnique make_unique_surf( all_points[all_points_itr] ) ;
            make_unique_surf.unique() ;
            std::vector< vec3 > facet_points ;
            make_unique_surf.unique_points( facet_points ) ;
            const std::vector< index_t >& unique_id = make_unique_surf.indices() ;
            std::vector< index_t > facet_indices ;
            std::vector< index_t > facet_ptr ;
            index_t count_facet_vertices = 0 ;
            facet_ptr.push_back( count_facet_vertices ) ;
            index_t offset = 0 ;

            for( index_t components_itr = 0; components_itr < components.size();
                ++components_itr ) {
                ringmesh_assert( components[components_itr] < all_points.size() ) ;
                if( components[components_itr] == all_points_itr ) {
                    for( index_t v_in_f_itr = 0;
                        v_in_f_itr
                            < surface_mesh.facets.nb_vertices( components_itr );
                        ++v_in_f_itr ) {
                        facet_indices.push_back( unique_id[offset] ) ;
                        ++offset ;
                    }
                    count_facet_vertices += surface_mesh.facets.nb_vertices(
                        components_itr ) ;
                    facet_ptr.push_back( count_facet_vertices ) ;
                }
            }
            GME::gme_t new_new_surface_gme_t = create_element( GME::SURFACE ) ;
            set_surface_geometry( new_new_surface_gme_t.index, facet_points,
                facet_indices, facet_ptr ) ;
            set_element_parent( new_new_surface_gme_t, sided_interface_gme_t ) ;
            add_element_child( sided_interface_gme_t, new_new_surface_gme_t ) ;
            add_element_in_boundary( new_new_surface_gme_t,
                GME::gme_t( GME::REGION, region_index ) ) ;
            bool side = ( side_name == "_plus" ) ;
            add_element_boundary( GME::gme_t( GME::REGION, region_index ),
                new_new_surface_gme_t, side ) ;

            // Add to universe (other side of the surface)
            add_element_boundary( model_.universe().gme_id(), new_new_surface_gme_t,
                !side ) ;
            to_erase_by_type[GME::SURFACE].push_back( 0 ) ;

#ifdef RINGMESH_DEBUG
            // In theory there is no isolated vertex
            index_t previous =
                model_.surface( new_new_surface_gme_t.index ).mesh().vertices.nb() ;
            model_.surface( new_new_surface_gme_t.index ).mesh().vertices.remove_isolated() ;
            ringmesh_assert(previous==model_.surface(new_new_surface_gme_t.index).mesh().vertices.nb()) ;
#endif
        }
        /// @todo HANDLE THE INTERNAL BORDER
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
            save_normal_on_one_surface( cur_surface ) ;
        }
    }

    void DuplicateInterfaceBuilder::save_normal_on_one_surface(
        const Surface& surface ) const
    {
        GEO::Mesh& cur_surf_mesh = surface.mesh() ;
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
            vec3 cur_normal( normal_att_x[i], normal_att_y[i], normal_att_z[i] ) ;
            cur_normal = normalize( cur_normal ) ;
            normal_att_x[i] = cur_normal.x ;
            normal_att_y[i] = cur_normal.y ;
            normal_att_z[i] = cur_normal.z ;
        }
        // end copy paste from GEO::compute_normals
    }

    vec3 DuplicateInterfaceBuilder::get_local_translation_normal(
        const Surface& surface,
        index_t vertex_id_in_surface ) const
    {
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
        if( !side ) {
            normal *= -1 ;
        }
        return normal ;
    }

    vec3 DuplicateInterfaceBuilder::get_local_translation_vector(
        const vec3& normal ) const
    {
        vec3 displacement = normal * 1.5 * 1 ;
        return displacement ;
    }

    void DuplicateInterfaceBuilder::compute_translation_vectors_duplicated_fault_network_surfaces_and_regions(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        initialize_translation_attributes( to_erase_by_type ) ;

        ringmesh_assert( model_.nb_interfaces() - first_new_interface_index >= 2 ) ;
        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr < model_.nb_interfaces(); ++new_interface_itr ) {
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] != NO_ID ) ;
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] == 0 ) ;

            const GeoModelElement& interface_gme = model_.one_interface(
                new_interface_itr ) ;
//            save_normals_on_one_new_interface( to_erase_by_type, interface_gme ) ;

            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
                ringmesh_assert(cur_child.type() == GME::SURFACE) ;
                ringmesh_assert(to_erase_by_type[GME::SURFACE][cur_child.index()] != NO_ID) ;
                const Surface& cur_surface = model_.surface( cur_child.index() ) ; // avoid dynamic_cast of cur_child

                ringmesh_assert(cur_surface.nb_in_boundary()==1) ;
                ringmesh_assert(cur_surface.in_boundary(0).type()==GME::REGION) ;
                GEO::Attribute< index_t > id_in_link_vector(
                    model_.surface( cur_surface.index() ).mesh().vertices.attributes(),
                    "id_in_link_vector" ) ;

                for( index_t surf_vertex_itr = 0;
                    surf_vertex_itr < cur_surface.nb_vertices();
                    ++surf_vertex_itr ) {

                    const vec3 local_translation_normal =
                        get_local_translation_normal( cur_surface,
                            surf_vertex_itr ) ;
                    const vec3 local_translation_vector =
                        get_local_translation_vector( local_translation_normal ) ;

                    index_t id_in_link = id_in_link_vector[surf_vertex_itr] ;

                    ringmesh_assert( gme_vertices_links_[id_in_link] != nil ) ;
                    gme_vertices_links_[id_in_link]->displace(
                        local_translation_vector ) ;
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::compute_translation_vectors_duplicated_fault_network(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {
        initialize_translation_attributes( to_erase_by_type ) ;

        int step_to_other_side = 1 ;
        ringmesh_assert( model_.nb_interfaces() - first_new_interface_index >= 2 ) ;
        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr < model_.nb_interfaces(); ++new_interface_itr ) {
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] != NO_ID ) ;
            ringmesh_assert( to_erase_by_type[GME::INTERFACE][new_interface_itr] == 0 ) ;

            const GeoModelElement& interface_gme = model_.one_interface(
                new_interface_itr ) ;
            ringmesh_assert( new_interface_itr + step_to_other_side >=first_new_interface_index ) ;
            ringmesh_assert( new_interface_itr + step_to_other_side <model_.nb_interfaces()) ;
            const GeoModelElement& other_side_interface_gme = model_.one_interface(
                new_interface_itr + step_to_other_side ) ;
            step_to_other_side *= -1 ;

            save_normals_on_one_new_interface( to_erase_by_type, interface_gme ) ;

            // Clear to take into account the new gme in the geomodel.
            recompute_geomodel_mesh() ; // not done in model_vertex_id

            const GeoModelMeshVertices& gmmv = model_.mesh.vertices ;

            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GeoModelElement& cur_child = interface_gme.child( child_itr ) ;
                ringmesh_assert(cur_child.type() == GME::SURFACE) ;
                ringmesh_assert(to_erase_by_type[GME::SURFACE][cur_child.index()] != NO_ID) ;
                const Surface& cur_surface = model_.surface( cur_child.index() ) ; // avoid dynamic_cast of cur_child

                ringmesh_assert(cur_surface.nb_in_boundary()==1) ;
                ringmesh_assert(cur_surface.in_boundary(0).type()==GME::REGION) ;
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

                    const vec3 local_translation_normal =
                        get_local_translation_normal( cur_surface,
                            surf_vertex_itr ) ;
                    const vec3 local_translation_vector =
                        get_local_translation_vector( local_translation_normal ) ;

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

                        if( is_surface_or_region_on_the_right_side_of_the_fault(
                            cur_gme_t, local_translation_normal,
                            gme_vertices[gme_vertex_itr].v_id,
                            model_.mesh.vertices.vertex( vertex_id_in_gmm ),
                            interface_gme, other_side_interface_gme ) ) {
                            store_displacement_in_gme(
                                model_.mesh_element( cur_gme_t ),
                                gme_vertices[gme_vertex_itr].v_id,
                                local_translation_vector ) ;
                        }
                    }
                }
            }
        }
    }

    bool DuplicateInterfaceBuilder::is_surface_or_region_on_the_right_side_of_the_fault(
        const GME::gme_t& cur_gme_t,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_gmme,
        const vec3& vertex_pos,
        const GeoModelElement& interface_gme,
        const GeoModelElement& other_side_interface_gme ) const
    {
        if( cur_gme_t.type == GME::REGION ) {
            if( !model_.region( cur_gme_t.index ).is_meshed() ) {
                return false ;
            }
            if( !is_region_on_right_side_of_sided_interface( cur_gme_t.index,
                normal_on_vertex_interface, vertex_id_in_gmme, vertex_pos ) ) {
                // Region on the other side of the fault
                return false ;
            }
        } else {
            ringmesh_assert(cur_gme_t.type == GME::SURFACE) ;
            ringmesh_assert( interface_gme.type() == GME::INTERFACE ) ;
            for( index_t interface__child_itr = 0;
                interface__child_itr < interface_gme.nb_children();
                ++interface__child_itr ) {
                const GeoModelElement& cur_child = interface_gme.child(
                    interface__child_itr ) ;
                ringmesh_assert( cur_child.type() == GME::SURFACE ) ;
                if( cur_child.index() == cur_gme_t.index ) {
                    return true ;
                }
            }

            ringmesh_assert( other_side_interface_gme.type() == GME::INTERFACE ) ;
            for( index_t other_side_interface__child_itr = 0;
                other_side_interface__child_itr
                    < other_side_interface_gme.nb_children();
                ++other_side_interface__child_itr ) {
                const GeoModelElement& cur_other_side_child =
                    other_side_interface_gme.child(
                        other_side_interface__child_itr ) ;
                ringmesh_assert( cur_other_side_child.type() == GME::SURFACE ) ;
                if( cur_other_side_child.index() == cur_gme_t.index ) {
                    return false ;
                }
            }

            if( !is_surface_on_right_side_of_sided_interface( cur_gme_t.index,
                normal_on_vertex_interface, vertex_id_in_gmme, vertex_pos ) ) {
                return false ;
            }
        }
        return true ;
    }

    bool DuplicateInterfaceBuilder::is_region_on_right_side_of_sided_interface(
        index_t region_to_check_id,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_region,
        const vec3& vertex_pos ) const
    {
        ringmesh_assert(region_to_check_id<model_.nb_regions()) ;

        const Region& region_to_check = model_.region( region_to_check_id ) ;
        std::vector< index_t > cells_around ;
        cells_around.reserve( 10 ) ;
        region_to_check.cells_around_vertex( vertex_id_in_region, cells_around,
            false ) ;
        ringmesh_assert( !cells_around.empty() ) ;

        vec3 region_to_check_mean_normal_on_vertex( 0., 0., 0. ) ;
        for( index_t cells_around_itr = 0; cells_around_itr < cells_around.size();
            ++cells_around_itr ) {

            index_t cur_cell_id_in_region = cells_around[cells_around_itr] ;
            vec3 cur_cell_barycenter = region_to_check.cell_barycenter(
                cur_cell_id_in_region ) ;
            vec3 p = cur_cell_barycenter - vertex_pos ;
            region_to_check_mean_normal_on_vertex += p ;
        }
        ringmesh_assert(
            std::abs(region_to_check_mean_normal_on_vertex.x ) > epsilon ||
            std::abs(region_to_check_mean_normal_on_vertex.y ) > epsilon ||
            std::abs(region_to_check_mean_normal_on_vertex.z ) > epsilon ) ;

        if( GEO::dot( normal_on_vertex_interface,
            region_to_check_mean_normal_on_vertex ) > epsilon ) {
            return true ;
        }
        return false ;
    }

    bool DuplicateInterfaceBuilder::is_surface_on_right_side_of_sided_interface(
        index_t surface_to_check_id,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_surface,
        const vec3& vertex_pos ) const
    {
        ringmesh_assert(surface_to_check_id<model_.nb_surfaces()) ;

        const Surface& surface_to_check = model_.surface( surface_to_check_id ) ;
        ringmesh_assert(surface_to_check.nb_in_boundary()==1 || surface_to_check.nb_in_boundary()==2) ;

        std::vector< index_t > facets_around ;
        facets_around.reserve( 10 ) ;
        surface_to_check.facets_around_vertex( vertex_id_in_surface, facets_around,
            false ) ;
        ringmesh_assert( !facets_around.empty() ) ;

        vec3 surf_to_check_mean_normal_on_vertex( 0., 0., 0. ) ;
        for( index_t facets_around_itr = 0; facets_around_itr < facets_around.size();
            ++facets_around_itr ) {
            index_t cur_facet_id_in_surf = facets_around[facets_around_itr] ;
            vec3 cur_facet_barycenter = surface_to_check.facet_barycenter(
                cur_facet_id_in_surf ) ;
            vec3 p = cur_facet_barycenter - vertex_pos ;
            surf_to_check_mean_normal_on_vertex += p ;
        }
        ringmesh_assert(
            std::abs(surf_to_check_mean_normal_on_vertex.x ) > epsilon ||
            std::abs(surf_to_check_mean_normal_on_vertex.y ) > epsilon ||
            std::abs(surf_to_check_mean_normal_on_vertex.z ) > epsilon ) ;

        if( GEO::dot( 50 * normal_on_vertex_interface, /// @todo remove this 50
        surf_to_check_mean_normal_on_vertex ) > epsilon ) {
            return true ;
        }
        return false ;
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

    void DuplicateInterfaceBuilder::set_no_displacement_on_fault_real_extension(
        const std::vector< std::vector< index_t > >& to_erase_by_type )
    {

        for( index_t line_itr = 0; line_itr < model_.nb_lines(); ++line_itr ) {
            const Line& cur_line = model_.line( line_itr ) ;
            if( cur_line.nb_in_boundary() != 1 ) {
                continue ;
            }
            // cur_line.nb_in_boundary() == 1 means fault extension
            /// @todo put assert to check that it is indeed a fault
            for( index_t line_vertex_itr = 0;
                line_vertex_itr < cur_line.nb_vertices(); ++line_vertex_itr ) {

                index_t vertex_id_in_gmm = cur_line.model_vertex_id(
                    line_vertex_itr ) ;
                const std::vector< GMEVertex >& gme_vertices =
                    model_.mesh.vertices.gme_vertices( vertex_id_in_gmm ) ;
                for( index_t gme_vertex_itr = 0;
                    gme_vertex_itr < gme_vertices.size(); ++gme_vertex_itr ) {

                    const GMEVertex& cur_gme_vertex = gme_vertices[gme_vertex_itr] ;
                    if( cur_gme_vertex.gme_id.type != GME::SURFACE
                        && cur_gme_vertex.gme_id.type != GME::REGION ) {
                        continue ;
                    }
                    if( to_erase_by_type[cur_gme_vertex.gme_id.type][cur_gme_vertex.gme_id.index]
                        == NO_ID ) {
                        continue ;
                    }
                    ringmesh_assert( to_erase_by_type[cur_gme_vertex.gme_id.type][cur_gme_vertex.gme_id.index] == 0 ) ;
                    const GeoModelMeshElement& gmme = model_.mesh_element(
                        cur_gme_vertex.gme_id ) ;
                    GEO::AttributesManager& att_mgr =
                        gmme.mesh().vertices.attributes() ;
                    GEO::Attribute< double > translation_att_x( att_mgr,
                        "translation_attr_x" ) ;
                    translation_att_x[cur_gme_vertex.v_id] = 0 ;
                    GEO::Attribute< double > translation_att_y( att_mgr,
                        "translation_attr_y" ) ;
                    translation_att_y[cur_gme_vertex.v_id] = 0 ;
                    GEO::Attribute< double > translation_att_z( att_mgr,
                        "translation_attr_z" ) ;
                    translation_att_z[cur_gme_vertex.v_id] = 0 ;
                }
            }
        }
    }
}
