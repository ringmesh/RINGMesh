/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geomodel/duplicate_fntk_builder.h>

#include <geogram/mesh/mesh_io.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/algorithm.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_topology.h>

/*!
 * @file ringmesh/geomodel/duplicate_fntk_builder.cpp
 * @brief Class to duplicate GeoModel Interface to
 * enable sliding along them (faults) and unconformal
 * mesh generation.
 * @author Benjamin Chauvin
 */

namespace {
    using namespace RINGMesh;

    /// @todo VERY BAD COPY PASTER FROM THE FUNCTION IN GEO_MODEL_BUILDER.CPP
    /*!
     * Finds a facet and its edge index that are colocalised with an edge
     * defined by its two geomodel vertex indices
     * @param[in] ann a ColocatorANN of the Surface \p surface using the keyword FACETS
     * @param[in] surface the surface where to find the facet
     * @param[in] geomodel_v0 the first geomodel vertex index of the edge
     * @param[in] geomodel_v1 the second geomodel vertex index of the edge
     * @param[out] f the found facet index
     * @param[out] e the found edge index
     * @return True if the facet and the edge indices are found
     * @todo RENAME these parameters and split in smaller functions !! [JP]
     */
    bool find_facet_and_edge(
        const NNSearch3D& nn_search,
        const Surface3D& surface,
        index_t geomodel_v0,
        index_t geomodel_v1,
        index_t& f,
        index_t& e )
    {
        const GeoModelMeshVertices3D& gmmv = surface.geomodel().mesh.vertices;
        // This is bad ! One level of abstraction is far far away
        const vec3& v0 = gmmv.vertex( geomodel_v0 );
        const vec3& v1 = gmmv.vertex( geomodel_v1 );
        vec3 v_bary = 0.5 * ( v0 + v1 );

        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() );
        std::vector< index_t > neighbors;
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() );
            neighbors.resize( cur_neighbor );
            neighbors = nn_search.get_neighbors( v_bary, cur_neighbor );
            // nb_neighbors can be less than cur_neighbor.
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                f = neighbors[i];
                for( index_t j = 0; j < surface.nb_mesh_element_vertices( f );
                    j++ ) {
                    if( gmmv.geomodel_vertex_id( surface.gmme(),
                        ElementLocalVertex( f, j ) ) == geomodel_v0 ) {
                        index_t j_next =
                            surface.low_level_mesh_storage().next_polygon_vertex(
                                ElementLocalVertex( f, j ) ).local_vertex_id_;
                        if( gmmv.geomodel_vertex_id( surface.gmme(),
                            ElementLocalVertex( f, j_next ) ) == geomodel_v1 ) {
                            e = j;
                            return true;
                        }
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor );

        f = NO_ID;
        e = NO_ID;
        return false;
    }

    /// TODO can be more generic and take a gme_t + an index_t (type deduced from first gme_t.type)
    index_t find_local_boundary_id( const Region3D& reg, const Surface3D& surf )
    {
        for( index_t boundary_i = 0; boundary_i < reg.nb_boundaries();
            ++boundary_i ) {
            if( reg.boundary( boundary_i ).index() == surf.index() ) {
                return boundary_i;
            }
        }

        return NO_ID;
    }

    index_t find_local_boundary_id( const Universe3D& reg, const Surface3D& surf )
    {
        for( index_t boundary_i = 0; boundary_i < reg.nb_boundaries();
            ++boundary_i ) {
            if( reg.boundary_gmme( boundary_i ).index() == surf.index() ) {
                return boundary_i;
            }
        }

        return NO_ID;
    }

    index_t find_second_local_boundary_id(
        const Region3D& reg,
        const Surface3D& surf )
    {
        bool is_second = false;
        for( index_t boundary_i = 0; boundary_i < reg.nb_boundaries();
            ++boundary_i ) {
            if( reg.boundary( boundary_i ).index() == surf.index() ) {
                if( !is_second ) {
                    is_second = true;
                } else {
                    return boundary_i;
                }
            }
        }

        return NO_ID;
    }

    void fill_vect_with_NO_ID( std::vector< index_t >& to_fill )
    {
        std::fill( to_fill.begin(), to_fill.end(), NO_ID );
    }
}

namespace RINGMesh {

    /// TODO should we check if these attributes already exist?
    const std::string DuplicateInterfaceBuilder::translation_attribute_name_ =
        "translation_vector_for_fault_duplication";
    const std::string DuplicateInterfaceBuilder::normal_attribute_name_ =
        "normal_at_vertex";

    DuplicateInterfaceBuilder::DuplicateInterfaceBuilder( GeoModel3D& geomodel )
        : GeoModelBuilder3D( geomodel ), all_meshed_( true ), gme_vertices_links_()
    {
        fill_entity_type_to_index_map();
    }

    DuplicateInterfaceBuilder::~DuplicateInterfaceBuilder()
    {
        for( GMEVertexLink*& gme_vertices_links_itr : gme_vertices_links_ ) {
            ringmesh_assert( gme_vertices_links_itr );
            delete gme_vertices_links_itr;
        }
        /// TODO delete the created attributes?
    }

    void DuplicateInterfaceBuilder::duplicate_fault_network( bool gap )
    {
        DEBUG( geomodel_.epsilon() );
        check_geomodel_validity_for_duplication();

        std::vector< std::vector< index_t > > geological_entities_to_erase;
        std::vector< std::vector< index_t > > mesh_entities_to_erase;
        const MeshEntityTypeManager3D& mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeologicalTypeManager& geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;

        geological_entities_to_erase.resize(
            geological_entity_type_manager.nb_geological_entity_types() );
        mesh_entities_to_erase.resize(
            mesh_entity_type_manager.nb_mesh_entity_types() );

        for( auto geomodel_mesh_entity_index : range(
            mesh_entity_type_manager.nb_mesh_entity_types() ) ) {
            mesh_entities_to_erase.emplace_back(
                std::vector< index_t >(
                    geomodel_.nb_mesh_entities(
                        MeshEntityType(
                            mesh_entity_type_manager.mesh_entity_types()[geomodel_mesh_entity_index] ) ),
                    0 ) );
        }

        for( auto geomodel_geological_entity_index : range(
            geological_entity_type_manager.nb_geological_entity_types() ) ) {
            geological_entities_to_erase.emplace_back(
                std::vector< index_t >(
                    geomodel_.nb_geological_entities(
                        geological_entity_type_manager.geological_entity_type(
                            geomodel_geological_entity_index ) ), 0 ) );
        }

        homogenize_normal_orientation_surface_all_interfaces();
        const index_t nb_initial_interfaces = geomodel_.nb_geological_entities(
            Interface3D::type_name_static() );
        build_new_fault_surfaces( mesh_entities_to_erase,
            geological_entities_to_erase );

        add_hole_between_faults( mesh_entities_to_erase,
            geological_entities_to_erase, nb_initial_interfaces );

        // Put here for new.
        /// @todo if all the lines are removed, is it still necessary to fill the new
        /// interface children with them? They should be recomputed with
        /// build_lines_and_corners_from_surfaces.
        flag_corners_lines_contacts_to_be_deleted( mesh_entities_to_erase,
            geological_entities_to_erase );
        delete_old_entities( mesh_entities_to_erase, geological_entities_to_erase );
        rebuild_valid_geomodel();
        if( !gap ) {
            remove_gap();
        }
    }

    void DuplicateInterfaceBuilder::check_geomodel_validity_for_duplication()
    {
        if( geomodel_.nb_regions() == 0 ) {
            throw RINGMeshException( "Fault Network Duplication",
                "The geomodel must contain at least one region." );
        }
        all_meshed_ = geomodel_.region( 0 ).is_meshed();
        for( index_t reg_itr = 1; reg_itr < geomodel_.nb_regions(); ++reg_itr ) {
            if( geomodel_.region( reg_itr ).is_meshed() != all_meshed_ ) {
                throw RINGMeshException( "Fault Network Duplication",
                    "The regions must be all meshed or all unmeshed." );
            }
        }

        for( index_t interface_itr = 0;
            interface_itr
                < geomodel_.nb_geological_entities( Interface3D::type_name_static() );
            ++interface_itr ) {
            const GeoModelGeologicalEntity3D& cur_interface =
                geomodel_.geological_entity( Interface3D::type_name_static(),
                    interface_itr );
            if( GeoModelGeologicalEntity3D::is_fault(
                cur_interface.geological_feature() ) ) {
                return;
            }
        }
        std::string message = "There is no fault inside the geomodel.\n";
        message += "Assign your fault interfaces to fault geological features.";
        throw RINGMeshException( "Fault Network Duplication", message );
    }

    void DuplicateInterfaceBuilder::homogenize_normal_orientation_surface_all_interfaces()
    {
        std::vector< index_t > surfaces_to_inverse_normals;
        surfaces_to_inverse_normals.reserve( geomodel_.nb_surfaces() );
        for( index_t interface_itr = 0;
            interface_itr
                < geomodel_.nb_geological_entities( Interface3D::type_name_static() );
            ++interface_itr ) {
            const GeoModelGeologicalEntity3D& cur_interface = interface(
                interface_itr );
            if( !GeoModelGeologicalEntity3D::is_fault(
                cur_interface.geological_feature() ) ) {
                continue;
            }
            save_normals_on_one_old_interface( interface_itr );
            homogenize_normal_orientation_surface_one_interface( interface_itr,
                surfaces_to_inverse_normals );
        }

        invert_normals_of_surface_list( surfaces_to_inverse_normals );
        geometry.clear_geomodel_mesh(); // TODO check if rebuild is needed (as before: recompute geomodel mesh)
    }

    void DuplicateInterfaceBuilder::invert_normals_of_surface_list(
        std::vector< index_t >& surfaces_to_inverse_normals )
    {
        for( index_t to_inverse_itr = 0;
            to_inverse_itr < surfaces_to_inverse_normals.size(); ++to_inverse_itr ) {
            geometry.invert_surface_normals(
                surfaces_to_inverse_normals[to_inverse_itr] );
        }
    }

    void DuplicateInterfaceBuilder::homogenize_normal_orientation_surface_one_interface(
        index_t fault_interface_id,
        std::vector< index_t >& surfaces_to_inverse_normals )
    {
        const GeoModelGeologicalEntity3D& fault_interface = interface(
            fault_interface_id );
        ringmesh_assert(
            GeoModelGeologicalEntity3D::is_fault(
                fault_interface.geological_feature() ) );
        ringmesh_assert( fault_interface.nb_children() >= 1 );
        if( fault_interface.nb_children() == 1 ) {
            return;
        }

        std::vector< bool > already_seen( geomodel_.nb_surfaces(), false );

        const Surface3D& first_child = geomodel_.surface(
            fault_interface.child( 0 ).index() );
        already_seen[first_child.index()] = true;
        /// @todo for now I assume that all the surfaces of the fault interface
        /// form one connected component. What happens in case of a fault cut into
        /// 2 distinct parts by another fault??? A trick would be to define each part
        /// as a different fault interface...
        homogenize_surfaces_around_surface( fault_interface_id, first_child,
            already_seen, surfaces_to_inverse_normals );
    }

    const GeoModelGeologicalEntity3D& DuplicateInterfaceBuilder::interface(
        index_t interface_id ) const
    {
        return geomodel_.geological_entity( Interface3D::type_name_static(),
            interface_id );
    }

    void DuplicateInterfaceBuilder::homogenize_surfaces_around_surface(
        index_t fault_interface_id,
        const Surface3D& first_child,
        std::vector< bool >& already_seen,
        std::vector< index_t >& surfaces_to_inverse_normals )
    {
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
        for( index_t line_boundary_itr = 0;
            line_boundary_itr < first_child.nb_boundaries(); ++line_boundary_itr ) {
            const GeoModelMeshEntity3D& cur_line_boun = first_child.boundary(
                line_boundary_itr );
            ringmesh_assert(
                cur_line_boun.type_name() == Line3D::type_name_static() );
            for( index_t surf_incident_entity_itr = 0;
                surf_incident_entity_itr < cur_line_boun.nb_incident_entities();
                ++surf_incident_entity_itr ) {
                const GeoModelMeshEntity3D& cur_incident_entity =
                    cur_line_boun.incident_entity( surf_incident_entity_itr );
                ringmesh_assert(
                    cur_incident_entity.type_name()
                        == Surface3D::type_name_static() );
                const Surface3D& cur_surf_incident_entity = geomodel_.surface(
                    cur_incident_entity.index() );
                if( cur_surf_incident_entity.index() == first_child.index() ) {
                    continue;
                }
                if( !does_surface_belong_to_interface(
                    cur_surf_incident_entity.index(), fault_interface_id ) ) {
                    continue;
                }

                if( already_seen[cur_incident_entity.index()] ) {
                    continue;
                }

                const Line3D& cur_line = geomodel_.line( cur_line_boun.index() );
                ringmesh_assert( cur_line.nb_vertices() > 0 );
                index_t first_line_vertex_id_in_gmm = gmmv.geomodel_vertex_id(
                    cur_line.gmme(), 0 );
                const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
                    first_line_vertex_id_in_gmm );
                index_t v_id_in_first_surf = NO_ID;
                bool v_id_in_second_surf = NO_ID;
                for( const GMEVertex& gme_vertex_itr : gme_vertices ) {
                    if( gme_vertex_itr.gmme.type()
                        != Surface3D::type_name_static() ) {
                        continue;
                    }
                    if( gme_vertex_itr.gmme.index() == first_child.index() ) {
                        v_id_in_first_surf = gme_vertex_itr.v_index;
                    } else if( gme_vertex_itr.gmme.index()
                        == cur_surf_incident_entity.index() ) {
                        v_id_in_second_surf = gme_vertex_itr.v_index;
                    } else {
                        continue;
                    }

                    if( v_id_in_first_surf != NO_ID
                        && v_id_in_second_surf != NO_ID ) {
                        break;
                    }
                }

                ringmesh_assert(
                    v_id_in_first_surf != NO_ID && v_id_in_second_surf != NO_ID );

                vec3 first_normal = get_normal_on_surface_vertex( first_child,
                    v_id_in_first_surf );
                vec3 second_normal = get_normal_on_surface_vertex(
                    cur_surf_incident_entity, v_id_in_second_surf );
                // As the normals are computed on different surfaces, they cannot
                // be exactly the same. However they should be very close (0.3 here is arbitrary).
                ringmesh_assert(
                    std::abs( first_normal.length() - 1 ) < global_epsilon );
                ringmesh_assert(
                    std::abs( second_normal.length() - 1 ) < global_epsilon );
                double dot_product = GEO::dot( first_normal, second_normal );
                ringmesh_assert( std::abs( std::abs( dot_product ) - 1 ) < 3e-1 ); //TODO to reput
                if( dot_product < 0 ) {
                    surfaces_to_inverse_normals.push_back(
                        cur_surf_incident_entity.index() );
                    inverse_normal_attribute_one_surface( cur_surf_incident_entity );
                    update_region_polarity( cur_surf_incident_entity );
                }
                already_seen[cur_surf_incident_entity.index()] = true;
                homogenize_surfaces_around_surface( fault_interface_id,
                    cur_surf_incident_entity, already_seen,
                    surfaces_to_inverse_normals );
                break;
            }
        }
    }

    void DuplicateInterfaceBuilder::update_region_polarity(
        const Surface3D& surface )
    {
        index_t nb_incident_entities = surface.nb_incident_entities();
        ringmesh_assert( nb_incident_entities == 1 || nb_incident_entities == 2 );
        if( nb_incident_entities == 1 ) {
            const Region3D& reg_gme =
                dynamic_cast< const Region3D& >( surface.incident_entity( 0 ) );
            index_t boundary_id = find_local_boundary_id( reg_gme, surface );
            Region3D& reg = const_cast< Region3D& >( geomodel_.region(
                reg_gme.index() ) );
            bool cur_side = reg.side( boundary_id );
            // Just the sign is changed, not the boundary. Assert to be sure that
            // it is well the case.
            ringmesh_assert(
                reg_gme.boundary_gmme( boundary_id ) == surface.gmme() );
            topology.set_mesh_entity_boundary( reg_gme.gmme(), boundary_id,
                surface.gmme().index(), !cur_side );
//            set_boundary_sign( reg, boundary_id, !cur_side ) ;

            // Add to universe (other side of the surface)
            boundary_id = find_local_boundary_id( geomodel_.universe(), surface );
            ringmesh_assert( geomodel_.universe().side( boundary_id ) != cur_side );
            ringmesh_assert(
                geomodel_.universe().boundary_gmme( boundary_id )
                    == surface.gmme() );
            topology.set_universe_boundary( boundary_id, surface.gmme().index(),
                cur_side );
//            set_universe_boundary_sign( boundary_id, cur_side ) ;
        } else {
            const GeoModelMeshEntity3D& reg_gme1 = surface.incident_entity( 0 );
            index_t boundary_id1 = find_local_boundary_id(
                dynamic_cast< const Region3D& >( reg_gme1 ), surface );
            Region3D& reg1 = const_cast< Region3D& >( geomodel_.region(
                reg_gme1.index() ) );
            bool cur_side1 = reg1.side( boundary_id1 );

            const GeoModelMeshEntity3D& reg_gme2 = surface.incident_entity( 1 );
            index_t boundary_id2 = find_local_boundary_id(
                dynamic_cast< const Region3D& >( reg_gme2 ), surface );
            Region3D& reg2 = const_cast< Region3D& >( geomodel_.region(
                reg_gme2.index() ) );
            bool cur_side2 = reg2.side( boundary_id2 );
            // to check. if reg1 == reg2 no side in particular.

            if( reg1.index() != reg2.index() ) {
                ringmesh_assert( cur_side1 != cur_side2 );
                // Just the sign is changed, not the boundary. Assert to be sure that
                // it is well the case.
                ringmesh_assert(
                    reg1.boundary_gmme( boundary_id1 ) == surface.gmme() );
                topology.set_mesh_entity_boundary( reg1.gmme(), boundary_id1,
                    surface.gmme().index(), !cur_side1 );
                ringmesh_assert(
                    reg2.boundary_gmme( boundary_id2 ) == surface.gmme() );
                topology.set_mesh_entity_boundary( reg2.gmme(), boundary_id2,
                    surface.gmme().index(), !cur_side2 );
            } else {
                // Same region on the both side of the interface.
                boundary_id2 = find_second_local_boundary_id( reg1, surface );
                ringmesh_assert( boundary_id2 != NO_ID );
                cur_side2 = reg1.side( boundary_id2 );
                ringmesh_assert( cur_side1 != cur_side2 );
                // Just the sign is changed, not the boundary. Assert to be sure that
                // it is well the case.
                ringmesh_assert(
                    reg1.boundary_gmme( boundary_id1 ) == surface.gmme() );
                ringmesh_assert(
                    reg1.boundary_gmme( boundary_id2 ) == surface.gmme() );
                topology.set_mesh_entity_boundary( reg1.gmme(), boundary_id1,
                    surface.gmme().index(), !cur_side1 );
                topology.set_mesh_entity_boundary( reg1.gmme(), boundary_id2,
                    surface.gmme().index(), !cur_side2 );
            }
        }
    }

    vec3 DuplicateInterfaceBuilder::get_normal_on_surface_vertex(
        const Surface3D& surface,
        index_t vertex_id_on_surface ) const
    {
        ringmesh_assert( vertex_id_on_surface < surface.nb_vertices() );
        GEO::AttributesManager& att_mgr = surface.vertex_attribute_manager();
        GEO::Attribute< vec3 > normal_attr( att_mgr, normal_attribute_name_ );
        return normal_attr[vertex_id_on_surface];
    }

    void DuplicateInterfaceBuilder::inverse_normal_attribute_one_surface(
        const Surface3D& surface ) const
    {
        GEO::AttributesManager& att_mgr = surface.vertex_attribute_manager();
        GEO::Attribute< double > normal_attr( att_mgr, normal_attribute_name_ );
        for( index_t vertex_id_on_surface = 0;
            vertex_id_on_surface < surface.nb_vertices(); ++vertex_id_on_surface ) {
            normal_attr[vertex_id_on_surface] *= -1;
        }
    }

    /// @todo COPY PASTE OF save_normals_on_one_new_interface
    /// @todo TO REFACTORE
    void DuplicateInterfaceBuilder::save_normals_on_one_old_interface(
        index_t interface_id ) const
    {
        const GeoModelGeologicalEntity3D& interface_gme = interface( interface_id );
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelMeshEntity3D& cur_child = interface_gme.child( child_itr );
            ringmesh_assert(
                cur_child.mesh_entity_type() == Surface3D::type_name_static() );
            // As the loop begins at the first new interface, no surface
            // met in this loop should be to delete.
            const Surface3D& cur_surface = geomodel_.surface( cur_child.index() );
            save_normal_on_one_surface( cur_surface );
        }
    }

    /// @todo could be inside an unamed namespace.
    /// @todo could a more generic function (is_child_of ?)
    /// @todo such function exists?
    bool DuplicateInterfaceBuilder::does_surface_belong_to_interface(
        index_t surface_id,
        index_t interface_id ) const
    {
        const GeoModelGeologicalEntity3D& interface_gme = interface( interface_id );
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelMeshEntity3D& cur_child = interface_gme.child( child_itr );
            ringmesh_assert(
                cur_child.mesh_entity_type() == Surface3D::type_name_static() );
            if( cur_child.index() == surface_id ) {
                return true;
            }
        }
        return false;
    }

    void DuplicateInterfaceBuilder::build_new_fault_surfaces(
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        std::vector< std::vector< index_t > >& geological_entities_to_erase )
    {
        const index_t nb_initial_interfaces = geomodel_.nb_geological_entities(
            Interface3D::type_name_static() );
        // Loop to nb_initial_interfaces. geomodel_.nb_interfaces() cannot be inside
        // the for statement since the number of interfaces will increase during
        // the duplication.
        for( index_t interface_itr = 0; interface_itr < nb_initial_interfaces;
            ++interface_itr ) {
            const GeoModelGeologicalEntity3D& cur_interface = interface(
                interface_itr );
            if( !GeoModelGeologicalEntity3D::is_fault(
                cur_interface.geological_feature() ) ) {
                continue;
            }

            // Delete of the interface (will be replaced by a custom interface with
            // side informations)

            geological_entities_to_erase[geomodel_.entity_type_manager().geological_entity_manager.geological_entity_type_index(
                Interface3D::type_name_static() )][cur_interface.index()] = NO_ID;
            get_new_surfaces( interface_itr, mesh_entities_to_erase,
                geological_entities_to_erase );
        }
    }

    void DuplicateInterfaceBuilder::flag_corners_lines_contacts_to_be_deleted(
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        std::vector< std::vector< index_t > >& geological_entities_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;

        fill_vect_with_NO_ID(
            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Corner3D::type_name_static() )] );
        fill_vect_with_NO_ID(
            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Line3D::type_name_static() )] );
        fill_vect_with_NO_ID(
            geological_entities_erase[geological_entity_type_manager.geological_entity_type_index(
                Contact3D::type_name_static() )] );
    }

    void DuplicateInterfaceBuilder::add_hole_between_faults(
        const std::vector< std::vector< index_t > >& geological_entities_to_erase,
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        index_t nb_initial_interfaces )
    {
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;
        // compute translation vectors
        if( all_meshed_ ) {
            ringmesh_assert(
                geomodel_.nb_geological_entities( Interface3D::type_name_static() )
                    - nb_initial_interfaces >= 2 );
            for( index_t new_interface_itr = nb_initial_interfaces;
                new_interface_itr
                    < geomodel_.nb_geological_entities(
                        Interface3D::type_name_static() ); ++new_interface_itr ) {

#ifdef RINGMESH_DEBUG
                // Only the new interfaces with no child are removed.
                // Such interfaces have no child because there are entirely in
                // the boundary of the geomodel (no surface not voi).
                if( geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                    Interface3D::type_name_static() )][new_interface_itr]
                    == NO_ID ) {
                    const GeoModelGeologicalEntity3D& interface_gme =
                        geomodel_.geological_entity( Interface3D::type_name_static(),
                            new_interface_itr );
                    ringmesh_assert( interface_gme.nb_children() == 0 );
                } else {
                    ringmesh_assert(
                        geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                            Interface3D::type_name_static() )][new_interface_itr]
                            == 0 );
                }
#endif

                save_normals_on_one_new_interface( mesh_entities_to_erase,
                    new_interface_itr );
            }
            DEBUG( "define_global_motion_relation" );
            define_global_motion_relation( mesh_entities_to_erase );
            DEBUG(
                "compute_translation_vectors_duplicated_fault_network_surfaces_and_regions" );
            compute_translation_vectors_duplicated_fault_network_surfaces_and_regions(
                nb_initial_interfaces, mesh_entities_to_erase,
                geological_entities_to_erase );
        } else {
            // Only for models with surfaces (region parts to remove?)
            compute_translation_vectors_duplicated_fault_network(
                nb_initial_interfaces, mesh_entities_to_erase,
                geological_entities_to_erase );
        }
        // set no translation on fault real extension (only on fault ending inside
        // the geomodel).
        DEBUG( "set_no_displacement_on_fault_real_extension" );
        set_no_displacement_on_fault_real_extension( mesh_entities_to_erase );
        // apply translation
        DEBUG( "translate_duplicated_fault_network" );
        translate_duplicated_fault_network( mesh_entities_to_erase );
    }

    void DuplicateInterfaceBuilder::delete_old_entities(
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        std::vector< std::vector< index_t > >& geological_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;

        geometry.clear_geomodel_mesh(); // TODO check if rebuild is needed (as before: recompute geomodel mesh)

        std::set< gmme_id > to_delete_mesh_entities;
        std::set< gmge_id > to_delete_geological_entities;
        for( index_t i = 0; i < mesh_entities_to_erase.size(); ++i ) {
            for( index_t j = 0; j < mesh_entities_to_erase[i].size(); ++j ) {
                if( mesh_entities_to_erase[i][j] == NO_ID ) {
                    to_delete_mesh_entities.insert(
                        gmme_id( mesh_entity_type_manager.mesh_entity_types()[i],
                            j ) );
                }
            }
        }

        for( index_t i = 0; i < geological_entities_to_erase.size(); ++i ) {
            for( index_t j = 0; j < geological_entities_to_erase[i].size(); ++j ) {
                if( geological_entities_to_erase[i][j] == NO_ID ) {
                    to_delete_geological_entities.insert(
                        gmge_id(
                            geological_entity_type_manager.geological_entity_type(
                                i ), j ) );
                }
            }
        }

        removal.remove_mesh_entities( to_delete_mesh_entities );
//        remove_geological_entities( to_delete_geological_entities ) ; // TODO check if needed, I think yes.
    }

    void DuplicateInterfaceBuilder::rebuild_valid_geomodel()
    {
        //todo WTF
//        geometry.clear_geomodel_mesh(); // TODO check if rebuild is needed (as before: recompute geomodel mesh)
//        from_surfaces.build_lines_and_corners_from_surfaces();
////        topology.complete_entity_connectivity();
//        geology.build_contacts();
    }

    void DuplicateInterfaceBuilder::define_global_motion_relation(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        DEBUG( "initialize_gme_vertices_links" );
        initialize_gme_vertices_links( mesh_entities_to_erase );
        DEBUG( "do_define_motion_relation" );
        do_define_motion_relation( mesh_entities_to_erase );
    }

    /// @todo split this function into several smaller functions
    void DuplicateInterfaceBuilder::do_define_motion_relation_on_not_voi_surface(
        const Surface3D& cur_surface,
        index_t surf_facet_itr,
        index_t v_id_in_surf,
        const Region3D& reg1,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1 )
    {
        ringmesh_assert( cur_surface.nb_incident_entities() == 2 );
        ringmesh_assert(
            reg1.index() == cur_surface.incident_entity_gmme( 0 ).index() );
        const gmme_id& reg2_gmme_t = cur_surface.incident_entity_gmme( 1 );
        ringmesh_assert( reg2_gmme_t.type() == Region3D::type_name_static() );
        const Region3D& reg2 = geomodel_.region( reg2_gmme_t.index() );
        ringmesh_assert( reg2.is_meshed() );
        GEO::Attribute< index_t > id_in_link_vector_reg2(
            reg2.vertex_attribute_manager(), "id_in_link_vector" );

        // Working on horizon not voi
        ringmesh_assert( reg1.index() != reg2.index() );
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;

        const index_t surf_v_id_in_gmm = gmmv.geomodel_vertex_id( cur_surface.gmme(),
            v_id_in_surf );
        const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
            surf_v_id_in_gmm );
        std::vector< index_t > found_gmev_reg1;
        found_gmev_reg1.reserve( gme_vertices.size() );
        std::vector< index_t > found_gmev_reg2;
        found_gmev_reg2.reserve( gme_vertices.size() );
        for( index_t gme_vertices_itr = 0; gme_vertices_itr < gme_vertices.size();
            ++gme_vertices_itr ) {
            if( gme_vertices[gme_vertices_itr].gmme.type()
                != Region3D::type_name_static() ) {
                continue;
            }
            if( gme_vertices[gme_vertices_itr].gmme.index() == reg1.index() ) {
                found_gmev_reg1.push_back( gme_vertices_itr );
            } else if( gme_vertices[gme_vertices_itr].gmme.index()
                == reg2.index() ) {
                found_gmev_reg2.push_back( gme_vertices_itr );
            }
        }
        ringmesh_assert( !found_gmev_reg1.empty() );
        ringmesh_assert( !found_gmev_reg2.empty() );
        index_t v_id_in_reg1 =
            find_region_vertex_id_from_surface_facet_among_colocated_points_in_one_region(
                gme_vertices, found_gmev_reg1, cur_surface, surf_facet_itr, reg1,
                surf_v_id_in_gmm );
        index_t v_id_in_reg2 =
            find_region_vertex_id_from_surface_facet_among_colocated_points_in_one_region(
                gme_vertices, found_gmev_reg2, cur_surface, surf_facet_itr, reg2,
                surf_v_id_in_gmm );
        ringmesh_assert( v_id_in_reg1 != NO_ID );
        ringmesh_assert( v_id_in_reg2 != NO_ID );

        // LINKING
        index_t link_id_surf = id_in_link_vector_surf[v_id_in_surf];
        index_t link_id_reg1 = id_in_link_vector_reg1[v_id_in_reg1];
        index_t link_id_reg2 = id_in_link_vector_reg2[v_id_in_reg2];
        link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg1 );
        link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg2 );
    }

    std::set< std::vector< index_t > > DuplicateInterfaceBuilder::get_region_cell_facet_colocated_to_surface_facet_from_hints(
        const Region3D& reg,
        const Surface3D& cur_surface,
        index_t surf_facet_itr,
        const std::vector< index_t >& found_gmev_reg ) const
    {
        // set is necessary to avoid doublons on the surface border of the dead end fault.
        std::set< std::vector< index_t > > found_cells_and_cell_facets;
        /// TODO surface_facet_barycenter parameter of this method instead of cur_surface
        /// and surf_facet_itr???
        vec3 surface_facet_barycenter = cur_surface.mesh_element_barycenter(
            surf_facet_itr );
        for( index_t cur_found_gmev_reg : found_gmev_reg ) {
            std::vector< index_t > cells_around = reg.cells_around_vertex(
                cur_found_gmev_reg, NO_ID );
/// TODO if border_only option in reg.cells_around_vertex, put an assert,
/// else write a if and continue for loop if not in border...
            for( index_t cur_cell : cells_around ) {
                for( index_t cell_facet_i = 0;
                    cell_facet_i < reg.nb_cell_facets( cur_cell ); ++cell_facet_i ) {
                    if( !reg.is_cell_facet_on_border(
                        CellLocalFacet( cur_cell, cell_facet_i ) ) ) {
                        continue;
                    }
                    vec3 cur_cell_facet_barycenter =
                        reg.low_level_mesh_storage().cell_facet_barycenter(
                            CellLocalFacet( cur_cell, cell_facet_i ) );
                    vec3 diff = cur_cell_facet_barycenter - surface_facet_barycenter;
                    if( std::abs( diff.x ) < geomodel_.epsilon()
                        && std::abs( diff.y ) < geomodel_.epsilon()
                        && std::abs( diff.z ) < geomodel_.epsilon() ) {
                        std::vector< index_t > new_couple( 2 );
                        new_couple[0] = cur_cell;
                        new_couple[1] = cell_facet_i;
                        found_cells_and_cell_facets.insert( new_couple );
                    }
                }
            }
        }
        return found_cells_and_cell_facets;
    }

    index_t DuplicateInterfaceBuilder::find_region_vertex_id_from_surface_facet_among_colocated_points_in_one_region(
        const std::vector< GMEVertex >& gme_vertices,
        const std::vector< index_t >& found_gmev_reg,
        const Surface3D& cur_surface,
        index_t surf_facet_itr,
        const Region3D& reg,
        index_t surf_v_id_in_gmm ) const
    {
        if( found_gmev_reg.size() == 1 ) {
// It is the most common case. NNSearch is not used here to reduce
// the computational time.
            return gme_vertices[found_gmev_reg[0]].v_index;
        } else {
// Happens when the current horizon terminates on a fault which ends
// within Region reg...
// No choice. Test of the facet colocated. More time consuming.
            /*const vec3 facet_bary = cur_surface.mesh_element_barycenter(
             surf_facet_itr );*/
            std::set< std::vector< index_t > > colocated_facets_reg =
                get_region_cell_facet_colocated_to_surface_facet_from_hints( reg,
                    cur_surface, surf_facet_itr, found_gmev_reg );
            /*reg_nn_searches_[reg.index()]->get_neighbors( facet_bary,
             geomodel_.epsilon() );*/
            ringmesh_assert( colocated_facets_reg.size() == 1 );
            return find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
                reg, *colocated_facets_reg.begin(), surf_v_id_in_gmm );
        }
    }

    /// @todo split this function into several smaller functions
    void DuplicateInterfaceBuilder::do_define_motion_relation_on_voi_surface(
        const Surface3D& cur_surface,
        index_t surf_facet_itr,
        index_t v_id_in_surf,
        const Region3D& reg1,
        GEO::Attribute< index_t >& id_in_link_vector_surf,
        GEO::Attribute< index_t >& id_in_link_vector_reg1 )
    {
        ringmesh_assert( cur_surface.nb_incident_entities() == 1 );
        ringmesh_assert( cur_surface.incident_entity( 0 ).index() == reg1.index() );
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
        const index_t surf_v_id_in_gmm = gmmv.geomodel_vertex_id( cur_surface.gmme(),
            v_id_in_surf );
        const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
            surf_v_id_in_gmm );
        std::vector< index_t > found_regions;
        found_regions.reserve( gme_vertices.size() );
        for( index_t gme_vertices_itr = 0; gme_vertices_itr < gme_vertices.size();
            ++gme_vertices_itr ) {
            if( gme_vertices[gme_vertices_itr].gmme.type()
                != Region3D::type_name_static() ) {
                continue;
            }
            if( gme_vertices[gme_vertices_itr].gmme.index() != reg1.index() ) {
                continue;
            }
            found_regions.push_back( gme_vertices_itr );
        }
        index_t v_id_in_reg1 = NO_ID;
        /// Actually found_regions corresponds to the number of vertex of a region
        /// colocated with v_id_in_surf.
        if( found_regions.size() == 1 ) {
/// To check but for me (03 july 2017) this case corresponds to the case
/// of voi horizon or voi boundary or voi fault which has no other side
/// (only a + or a - side).
            v_id_in_reg1 = gme_vertices[found_regions[0]].v_index;
        } else {
// Configuration in which there is a fault which ends inside the geomodel.
// The current surface may be a voi boundary, a voi horizon or a fault. ((03 july 2017) not so sure about that).
// No choice => use of ColocaterAnn on the current facet (more time
// consuming).
            /*const vec3 facet_bary = cur_surface.mesh_element_barycenter(
             surf_facet_itr );*/
            /*std::vector< index_t > colocated_facets_reg1 =
             reg_nn_searches_[reg1.index()]->get_neighbors( facet_bary,
             geomodel_.epsilon() );*/
            std::set< std::vector< index_t > > colocated_facets_reg1 =
                get_region_cell_facet_colocated_to_surface_facet_from_hints( reg1,
                    cur_surface, surf_facet_itr, found_regions );
            if( colocated_facets_reg1.size() == 1 ) { /// (03 july 2017) Does this case really happen??? We are in a voi with only a vertex found...
// Case of a voi horizon or voi boundary
                v_id_in_reg1 =
                    find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
                        reg1, *colocated_facets_reg1.begin(), surf_v_id_in_gmm );
            } else {
// Case of a fault ending inside the geomodel
                ringmesh_assert( colocated_facets_reg1.size() == 2 );
// It is a fault with the region on the other side.
                ringmesh_assert(
                    GeoModelGeologicalEntity3D::is_fault(
                        cur_surface.parent( Interface3D::type_name_static() ).geological_feature() ) );

// Try to find which region facet is really on the surface.
                v_id_in_reg1 =
                    find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
                        reg1, *colocated_facets_reg1.begin(), surf_v_id_in_gmm );
                ringmesh_assert( v_id_in_reg1 != NO_ID );

                const vec3 local_translation_normal = get_local_translation_normal(
                    cur_surface, v_id_in_surf );
                bool is_first = is_region_on_right_side_of_sided_interface(
                    reg1.index(), local_translation_normal, v_id_in_reg1,
                    cur_surface.vertex( v_id_in_surf ) );

                if( !is_first ) {
                    v_id_in_reg1 = NO_ID;
                    v_id_in_reg1 =
                        find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
                            reg1, *( ++( colocated_facets_reg1.begin() ) ),
                            surf_v_id_in_gmm );
                    ringmesh_assert( v_id_in_reg1 != NO_ID );
                }
            }
        }
        ringmesh_assert( v_id_in_reg1 != NO_ID );
        // LINKING
        index_t link_id_surf = id_in_link_vector_surf[v_id_in_surf];
        index_t link_id_reg1 = id_in_link_vector_reg1[v_id_in_reg1];
        link_surf_vertex_id_to_reg_vertex_id( link_id_surf, link_id_reg1 );
    }

    index_t DuplicateInterfaceBuilder::find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
        const Region3D& reg,
        const std::vector< index_t >& reg_facet_cell,
        index_t surf_v_id_in_gmm ) const
    {
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
        for( index_t v_in_reg_itr = 0;
            v_in_reg_itr
                < reg.nb_cell_facet_vertices( reg_facet_cell[0], reg_facet_cell[1] );
            ++v_in_reg_itr ) {
            index_t reg_v_id_in_gmme = reg.cell_facet_vertex_index(
                reg_facet_cell[0], reg_facet_cell[1], v_in_reg_itr );
            if( gmmv.geomodel_vertex_id( reg.gmme(), reg_v_id_in_gmme )
                == surf_v_id_in_gmm ) {
                return reg_v_id_in_gmme;
            }
        }
        return NO_ID;
    }

    void DuplicateInterfaceBuilder::link_surf_vertex_id_to_reg_vertex_id(
        index_t link_id_surf,
        index_t link_id_reg )
    {
        gme_vertices_links_[link_id_surf]->add_linked_gme_vertex( link_id_reg );
        gme_vertices_links_[link_id_reg]->add_linked_gme_vertex( link_id_surf );
    }

    void DuplicateInterfaceBuilder::do_define_motion_relation(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        geometry.clear_geomodel_mesh(); /// @todo check if it is really necessary // TODO check if rebuild is needed (as before: recompute geomodel mesh)

        for( index_t surf_itr = 0; surf_itr < geomodel_.nb_surfaces(); ++surf_itr ) {
            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][surf_itr] == NO_ID ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Surface3D::type_name_static() )][surf_itr] == 0 );

// For each vertex of the surf, I must find all the region vertex
// ids which are linked to it, on one side for a fault
// (and geomodel boundary) and on the both sides for horizons.
// In theory, on the both sides of a horizon there are two different
// regions, else it is not a horizon.
            const Surface3D& cur_surface = geomodel_.surface( surf_itr );
            GEO::Attribute< index_t > id_in_link_vector_surf(
                cur_surface.vertex_attribute_manager(), "id_in_link_vector" );

            const index_t nb_incident_entities = cur_surface.nb_incident_entities();

            bool is_horizon_not_voi =
                GeoModelGeologicalEntity3D::is_stratigraphic_limit(
                    cur_surface.parent( Interface3D::type_name_static() ).geological_feature() )
                    && !cur_surface.is_on_voi();

            ringmesh_assert(
                nb_incident_entities == 1 || nb_incident_entities == 2 );

            const gmme_id& reg1_gme_t = cur_surface.incident_entity_gmme( 0 );
            ringmesh_assert( reg1_gme_t.type() == Region3D::type_name_static() );
            const Region3D& reg1 = geomodel_.region( reg1_gme_t.index() );
            ringmesh_assert( reg1.is_meshed() );
            GEO::Attribute< index_t > id_in_link_vector_reg1(
                reg1.vertex_attribute_manager(), "id_in_link_vector" );

            std::vector< bool > surf_vertex_visited( cur_surface.nb_vertices(),
                false );

            for( index_t surf_facet_itr = 0;
                surf_facet_itr < cur_surface.nb_mesh_elements(); ++surf_facet_itr ) {

                for( index_t v_in_surf_facet = 0;
                    v_in_surf_facet
                        < cur_surface.nb_mesh_element_vertices( surf_facet_itr );
                    ++v_in_surf_facet ) {
                    index_t v_id_in_surf = cur_surface.mesh_element_vertex_index(
                        ElementLocalVertex( surf_facet_itr, v_in_surf_facet ) );
                    if( surf_vertex_visited[v_id_in_surf] ) {
                        continue;
                    }
                    surf_vertex_visited[v_id_in_surf] = true;

                    if( is_horizon_not_voi ) {
                        ringmesh_assert( nb_incident_entities == 2 );
                        do_define_motion_relation_on_not_voi_surface( cur_surface,
                            surf_facet_itr, v_id_in_surf, reg1,
                            id_in_link_vector_surf, id_in_link_vector_reg1 );
                    } else {
                        ringmesh_assert( nb_incident_entities == 1 );
                        do_define_motion_relation_on_voi_surface( cur_surface,
                            surf_facet_itr, v_id_in_surf, reg1,
                            id_in_link_vector_surf, id_in_link_vector_reg1 );
                    }
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::initialize_gme_vertices_links(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        index_t count = 0;
        for( index_t surf_itr = 0; surf_itr < geomodel_.nb_surfaces(); ++surf_itr ) {
            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][surf_itr] == NO_ID ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Surface3D::type_name_static() )][surf_itr] == 0 );
            count += geomodel_.surface( surf_itr ).nb_vertices();
        }

        if( all_meshed_ ) {

            for( index_t reg_itr = 0; reg_itr < geomodel_.nb_regions(); ++reg_itr ) {
                if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] == NO_ID ) {
                    continue;
                }
                ringmesh_assert(
                    mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                        Region3D::type_name_static() )][reg_itr] == 0 );
                ringmesh_assert( geomodel_.region( reg_itr ).is_meshed() );
                count += geomodel_.region( reg_itr ).nb_vertices();
            }
        }

        gme_vertices_links_.reserve( count );

        for( index_t surf_itr = 0; surf_itr < geomodel_.nb_surfaces(); ++surf_itr ) {
            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][surf_itr] == NO_ID ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Surface3D::type_name_static() )][surf_itr] == 0 );

            GEO::Attribute< index_t > id_in_link_vector(
                geomodel_.surface( surf_itr ).vertex_attribute_manager(),
                "id_in_link_vector" );

            for( index_t v_itr = 0;
                v_itr < geomodel_.surface( surf_itr ).nb_vertices(); ++v_itr ) {
                id_in_link_vector[v_itr] = gme_vertices_links_.size();
                gme_vertices_links_.push_back(
                    new GMEVertexLink(
                        GMEVertex(
                            gmme_id( Surface3D::type_name_static(), surf_itr ),
                            v_itr ), geomodel_, gme_vertices_links_ ) );
            }
        }

        if( all_meshed_ ) {

            for( index_t reg_itr = 0; reg_itr < geomodel_.nb_regions(); ++reg_itr ) {
                if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] == NO_ID ) {
                    continue;
                }
                ringmesh_assert(
                    mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                        Region3D::type_name_static() )][reg_itr] == 0 );
                ringmesh_assert( geomodel_.region( reg_itr ).is_meshed() );
                GEO::Attribute< index_t > id_in_link_vector(
                    geomodel_.region( reg_itr ).vertex_attribute_manager(),
                    "id_in_link_vector" );
                for( index_t v_itr = 0;
                    v_itr < geomodel_.region( reg_itr ).nb_vertices(); ++v_itr ) {
                    id_in_link_vector[v_itr] = gme_vertices_links_.size();
                    gme_vertices_links_.push_back(
                        new GMEVertexLink(
                            GMEVertex(
                                gmme_id( Region3D::type_name_static(), reg_itr ),
                                v_itr ), geomodel_, gme_vertices_links_ ) );
                }
            }
        }

        ringmesh_assert(
            gme_vertices_links_.size() == gme_vertices_links_.capacity() );
        ringmesh_assert( gme_vertices_links_.size() == count );
    }

    DuplicateInterfaceBuilder::GMEVertexLink::GMEVertexLink(
        const GMEVertex& gme_vertex,
        const GeoModel3D& geomodel,
        const std::vector< GMEVertexLink* >& gme_vertices_links )
        :
            has_moved_( false ),
            geomodel_( geomodel ),
            gme_vertex_( gme_vertex ),
            gme_vertices_links_( gme_vertices_links )
    {

    }

    void DuplicateInterfaceBuilder::GMEVertexLink::displace(
        const vec3& displacement_vector )
    {
#ifdef RINGMESH_DEBUG
        if( gme_vertex_.gmme.type() == Surface3D::type_name_static() ) {
            ringmesh_assert( linked_gme_vertices_.size() <= 2 );
        } else {
            ringmesh_assert(
                gme_vertex_.gmme.type() == Region3D::type_name_static() );
            ringmesh_assert( linked_gme_vertices_.size() <= 5 ); // before it was 3 but 5 may exist, see internal_border_in_fault geomodel
        }

        for( index_t link_itr = 0; link_itr < linked_gme_vertices_.size();
            ++link_itr ) {
            for( index_t link_itr2 = 0; link_itr2 < linked_gme_vertices_.size();
                ++link_itr2 ) {
                if( link_itr == link_itr2 ) {
                    continue;
                }
                if( linked_gme_vertices_[link_itr]
                    == linked_gme_vertices_[link_itr2] ) {
                    DEBUG( "several time" );
                }
                ringmesh_assert(
                    linked_gme_vertices_[link_itr]
                        != linked_gme_vertices_[link_itr2] );
            }
        }
#endif

        if( has_moved_ ) {
            return;
        }
        has_moved_ = true;

        GEO::AttributesManager& att_mgr =
            geomodel_.mesh_entity( gme_vertex_.gmme ).vertex_attribute_manager();
        GEO::Attribute< vec3 > translation_attr( att_mgr,
            translation_attribute_name_ );
        translation_attr[gme_vertex_.v_index] += displacement_vector;

        for( index_t link_itr = 0; link_itr < linked_gme_vertices_.size();
            ++link_itr ) {
            gme_vertices_links_[linked_gme_vertices_[link_itr]]->displace(
                displacement_vector );
        }
        has_moved_ = false;
    }

    void DuplicateInterfaceBuilder::get_new_surfaces(
        index_t interface_to_duplicate_id,
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        std::vector< std::vector< index_t > >& geological_entities_to_erase )
    {
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeoModelGeologicalEntity3D& interface_to_duplicate = interface(
            interface_to_duplicate_id );
        // minus = false, plus = true
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_minus;
        std::map< index_t, std::vector< index_t > > surfaces_boundary_regions_side_plus;
        gmge_id interface_minus_gme_t = geology.create_geological_entity(
            Interface3D::type_name_static() );
        info.set_geological_entity_name( interface_minus_gme_t,
            interface_to_duplicate.name() + "_side_minus" );
        geology.set_geological_entity_geol_feature( interface_minus_gme_t,
            interface_to_duplicate.geological_feature() );
        gmge_id interface_plus_gme_t = geology.create_geological_entity(
            Interface3D::type_name_static() );
        info.set_geological_entity_name( interface_plus_gme_t,
            interface_to_duplicate.name() + "_side_plus" );
        geology.set_geological_entity_geol_feature( interface_plus_gme_t,
            interface_to_duplicate.geological_feature() );
        geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
            Interface3D::type_name_static() )].push_back( 0 );
        geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
            Interface3D::type_name_static() )].push_back( 0 );

        const index_t interface_to_duplicate_nb_children =
            interface_to_duplicate.nb_children();
        ringmesh_assert( interface_to_duplicate_nb_children >= 1 );
        // Find for each region, what surfaces are in boundary.
        for( index_t interface_child_itr = 0;
            interface_child_itr < interface_to_duplicate_nb_children;
            ++interface_child_itr ) {
            const GeoModelMeshEntity3D& cur_child = interface_to_duplicate.child(
                interface_child_itr );
            ringmesh_assert(
                cur_child.type_name() == Surface3D::type_name_static() );
            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][cur_child.index()] = NO_ID;

            const index_t nb_incident_entities_cur_child =
                cur_child.nb_incident_entities();
            ringmesh_assert(
                nb_incident_entities_cur_child == 1
                    || nb_incident_entities_cur_child == 2 );
            if( nb_incident_entities_cur_child == 2 ) {
                ringmesh_assert( !cur_child.is_on_voi() );
                const GeoModelMeshEntity3D& cur_incident_entity =
                    cur_child.incident_entity( 0 );
                ringmesh_assert(
                    cur_incident_entity.type_name()
                        == Region3D::type_name_static() );
                const Region3D& cur_reg =
                    dynamic_cast< const Region3D& >( cur_incident_entity );

                const GeoModelMeshEntity3D& cur_incident_entity2 =
                    cur_child.incident_entity( 1 );
                ringmesh_assert(
                    cur_incident_entity2.type_name()
                        == Region3D::type_name_static() );
                const Region3D& cur_reg2 =
                    dynamic_cast< const Region3D& >( cur_incident_entity2 );

/// @todo it seems that this if statement is contained in the else
/// case. To check and simplify if necessary.
                if( cur_incident_entity.index() == cur_incident_entity2.index() ) {
// if same region there is 2 incident_entities even if they are the same
// The surface is internal and on the both side there is the
// same region. This surface is duplicated.
                    surfaces_boundary_regions_side_plus[cur_incident_entity.index()].push_back(
                        cur_child.index() );
                    surfaces_boundary_regions_side_minus[cur_incident_entity.index()].push_back(
                        cur_child.index() );
                } else {
                    index_t local_boundary_id = find_local_boundary_id( cur_reg,
                        geomodel_.surface( cur_child.index() ) );

#ifdef RINGMESH_DEBUG
                    index_t local_boundary_id2 = find_local_boundary_id( cur_reg2,
                        geomodel_.surface( cur_child.index() ) );
                    ringmesh_assert(
                        cur_reg.side( local_boundary_id )
                            != cur_reg2.side( local_boundary_id2 ) );
#endif

                    if( cur_reg.side( local_boundary_id ) ) {
                        surfaces_boundary_regions_side_plus[cur_incident_entity.index()].push_back(
                            cur_child.index() );
                        surfaces_boundary_regions_side_minus[cur_incident_entity2.index()].push_back(
                            cur_child.index() );
                    } else {
                        surfaces_boundary_regions_side_minus[cur_incident_entity.index()].push_back(
                            cur_child.index() );
                        surfaces_boundary_regions_side_plus[cur_incident_entity2.index()].push_back(
                            cur_child.index() );
                    }
                }
            } else {
                ringmesh_assert( nb_incident_entities_cur_child == 1 );
                ringmesh_assert( cur_child.is_on_voi() );
                const GeoModelMeshEntity3D& cur_incident_entity =
                    cur_child.incident_entity( 0 );
                ringmesh_assert(
                    cur_incident_entity.type_name()
                        == Region3D::type_name_static() );
                const Region3D& cur_reg =
                    dynamic_cast< const Region3D& >( cur_incident_entity );

                index_t local_boundary_id = find_local_boundary_id( cur_reg,
                    geomodel_.surface( cur_child.index() ) );
                if( cur_reg.side( local_boundary_id ) ) {
                    surfaces_boundary_regions_side_plus[cur_incident_entity.index()].push_back(
                        cur_child.index() );
                } else {
                    surfaces_boundary_regions_side_minus[cur_incident_entity.index()].push_back(
                        cur_child.index() );
                }
            }
        }

        build_merged_surfaces( surfaces_boundary_regions_side_plus, "_plus",
            mesh_entities_to_erase, geological_entities_to_erase,
            interface_plus_gme_t.index(), interface_to_duplicate_id );
        build_merged_surfaces( surfaces_boundary_regions_side_minus, "_minus",
            mesh_entities_to_erase, geological_entities_to_erase,
            interface_minus_gme_t.index(), interface_to_duplicate_id );
    }

    void DuplicateInterfaceBuilder::build_merged_surfaces(
        const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
        const std::string& side_name,
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        std::vector< std::vector< index_t > >& geological_entities_to_erase,
        index_t sided_interface_id,
        index_t interface_to_duplicate_id )
    {
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;

        if( surfaces_boundary_regions.empty() ) {
// May happen when a fault is entirely a geomodel boundary
            geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                Interface3D::type_name_static() )][sided_interface_id] = NO_ID;
            return;
        }
        for( std::map< index_t, std::vector< index_t > >::const_iterator map_itr =
            surfaces_boundary_regions.begin();
            map_itr != surfaces_boundary_regions.end(); ++map_itr ) {

// first = line index in geomodel, second = count.
            std::map< index_t, index_t > all_surface_lines;

            index_t region_index = map_itr->first;
            std::vector< vec3 > all_points;
            for( std::vector< index_t >::const_iterator surf_itr =
                map_itr->second.begin(); surf_itr != map_itr->second.end();
                ++surf_itr ) {
                index_t surf_id = *surf_itr;

                const Surface3D& cur_surf = geomodel_.surface( surf_id );

                for( index_t vertex_surf_i = 0;
                    vertex_surf_i < cur_surf.nb_vertices(); ++vertex_surf_i ) {
                    all_points.push_back( cur_surf.vertex( vertex_surf_i ) );
                }
            }

            NNSearch3D make_unique_surf( all_points );
            std::vector< index_t > unique_id;
            std::vector< vec3 > facet_points;
            index_t not_used = 0;
            std::tie( not_used, unique_id, facet_points ) =
                make_unique_surf.get_colocated_index_mapping_and_unique_points(
                    geomodel_.epsilon() );

            index_t offset_vertices = 0;
            std::vector< index_t > facet_indices;
            std::vector< index_t > facet_ptr;
            index_t count_facet_vertices = 0;
            facet_ptr.push_back( count_facet_vertices );
            for( std::vector< index_t >::const_iterator surf_itr =
                map_itr->second.begin(); surf_itr != map_itr->second.end();
                ++surf_itr ) {
                index_t surf_id = *surf_itr;

                const Surface3D& cur_surf = geomodel_.surface( surf_id );

// Add current surface to merged surface
                for( index_t facet_itr = 0; facet_itr < cur_surf.nb_mesh_elements();
                    ++facet_itr ) {
                    for( index_t point_i = 0;
                        point_i < cur_surf.nb_mesh_element_vertices( facet_itr );
                        ++point_i ) {

                        index_t index = cur_surf.mesh_element_vertex_index(
                            ElementLocalVertex( facet_itr, point_i ) );
                        facet_indices.push_back(
                            unique_id[index + offset_vertices] );

                    }
                    count_facet_vertices += cur_surf.nb_mesh_element_vertices(
                        facet_itr );
                    facet_ptr.push_back( count_facet_vertices );
                }

// Update the lines in common
/// @todo this line part is necessary if used for the mutural cut
/// else delete it.
                for( index_t line_itr = 0; line_itr < cur_surf.nb_boundaries();
                    ++line_itr ) {
                    const GeoModelMeshEntity3D& cur_line_gme = cur_surf.boundary(
                        line_itr );
                    ringmesh_assert(
                        cur_line_gme.mesh_entity_type()
                            == Line3D::type_name_static() );

                    if( all_surface_lines.find( cur_line_gme.index() )
                        == all_surface_lines.end() ) {
                        all_surface_lines[cur_line_gme.index()] = 0; // initialization
                    }
                    ++all_surface_lines[cur_line_gme.index()];
                }

                offset_vertices += cur_surf.nb_vertices();
            }

// Create RINGMesh::Surface and fill it.
            gmme_id new_surface_gme_t = topology.create_mesh_entity< Surface >();
/// TODO weird that the GEO::vector cannot be put as parameter
/// of set_surface_geometry since it is a std::vector... to see.
            geometry.set_surface_geometry( new_surface_gme_t.index(), facet_points,
                facet_indices, facet_ptr );
            /*set_entity_parent( new_surface_gme_t, sided_interface_gme_t ) ;
             add_entity_child( sided_interface_gme_t, new_surface_gme_t ) ;
             */
// Boundary information is necessary for get_local_translation_normal
            bool side = ( side_name == "_plus" );
            topology.add_mesh_entity_boundary_relation(
                gmme_id( Region3D::type_name_static(), region_index ),
                new_surface_gme_t, side );
            /*
             // Add to universe (other side of the surface)
             add_entity_boundary( geomodel_.universe().gme_id(), new_surface_gme_t,
             !side ) ;*/

//            to_erase_by_type[GME::SURFACE].push_back( 0 ) ;
            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )].push_back( NO_ID );

            add_fake_internal_boundary_lines_to_merged_surface( all_surface_lines,
                side_name, sided_interface_id, interface_to_duplicate_id,
                new_surface_gme_t.index(), mesh_entities_to_erase, region_index );
        }
    }

    void DuplicateInterfaceBuilder::add_fake_internal_boundary_lines_to_merged_surface(
        const std::map< index_t, index_t >& all_surface_lines,
        const std::string& side_name,
        index_t sided_interface_id,
        index_t interface_to_duplicate_id,
        index_t new_surface_id,
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        index_t region_index )
    {
        geometry.clear_geomodel_mesh(); // to take into account the new surface in gme_vertices. // TODO check if rebuild is needed (as before: recompute geomodel mesh)
        save_normal_on_one_surface( geomodel_.surface( new_surface_id ) );
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
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
                continue;
            }
            ringmesh_assert( all_surface_lines_itr->second != 0 );
            const Line3D& cur_line = geomodel_.line( all_surface_lines_itr->first );
// As the line is not on the border, its number of in boundaries
// superior to 1 strictly.
            ringmesh_assert( cur_line.nb_incident_entities() > 1 );
            bool good_line = false;
            for( index_t incident_ent_itr = 0;
                incident_ent_itr < cur_line.nb_incident_entities();
                ++incident_ent_itr ) {
                const GeoModelMeshEntity3D& cur_incident_ent_gme =
                    cur_line.incident_entity( incident_ent_itr );
                ringmesh_assert(
                    cur_incident_ent_gme.type_name()
                        == Surface3D::type_name_static() );

                if( !cur_incident_ent_gme.has_parent() ) {
                    continue;
                }

                if( does_surface_belong_to_interface( cur_incident_ent_gme.index(),
                    sided_interface_id ) ) {
                    continue;
                }

// minus side is done before the plus side
                if( side_name == "_minus" ) {

                    ringmesh_assert(
                        sided_interface_id + 1
                            < geomodel_.nb_geological_entities(
                                Interface3D::type_name_static() ) );
                    const GeoModelGeologicalEntity3D& plus_side_gme =
                        geomodel_.geological_entity( Interface3D::type_name_static(),
                            sided_interface_id + 1 );

                    if( does_surface_belong_to_interface(
                        cur_incident_ent_gme.index(), plus_side_gme.index() ) ) {
                        continue;
                    }
                }

// Check if the surface is not from the old interface to
// duplicate.
                if( does_surface_belong_to_interface( cur_incident_ent_gme.index(),
                    interface_to_duplicate_id ) ) {
                    continue;
                }

                GeoModelGeologicalEntity3D::GEOL_FEATURE parent_geol_feature =
                    cur_incident_ent_gme.parent( Interface3D::type_name_static() ).geological_feature();
                if( !GeoModelGeologicalEntity3D::is_fault( parent_geol_feature ) ) {
                    continue;
                }

// Check if the found fault is on the right side
                ringmesh_assert( cur_line.nb_vertices() > 0 );
                index_t first_vertex_id_in_gmm = gmmv.geomodel_vertex_id(
                    cur_line.gmme(), 0 );
                const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
                    first_vertex_id_in_gmm );
                ringmesh_assert( !gme_vertices.empty() ); /// @todo I think that this assert may be more restrictive
                index_t vertex_id_in_new_surface = NO_ID;
                index_t vertex_id_in_found_surface = NO_ID;
                for( index_t gme_vertices_itr = 0;
                    gme_vertices_itr < gme_vertices.size(); ++gme_vertices_itr ) {
                    const GMEVertex& cur_gme_vertex = gme_vertices[gme_vertices_itr];
                    if( cur_gme_vertex.gmme.type()
                        != Surface3D::type_name_static() ) {
                        continue;
                    }
                    if( cur_gme_vertex.gmme.index() == new_surface_id ) {
                        vertex_id_in_new_surface = cur_gme_vertex.v_index;
                    } else if( cur_gme_vertex.gmme.index()
                        == cur_incident_ent_gme.index() ) {
                        vertex_id_in_found_surface = cur_gme_vertex.v_index;
                    }
                    if( vertex_id_in_new_surface != NO_ID
                        && vertex_id_in_found_surface != NO_ID ) {
                        break;
                    }
                }
                ringmesh_assert(
                    vertex_id_in_new_surface != NO_ID
                        && vertex_id_in_found_surface != NO_ID );

// On the found surface, the normal has been computed before
// for the homogenization of normals.

                const vec3 local_translation_normal = get_local_translation_normal(
                    geomodel_.surface( new_surface_id ), vertex_id_in_new_surface );

                vec3 vertex_pos = geomodel_.surface( new_surface_id ).vertex(
                    vertex_id_in_new_surface );
                if( !is_surface_on_right_side_of_sided_interface(
                    cur_incident_ent_gme.index(), local_translation_normal,
                    vertex_id_in_found_surface, vertex_pos ) ) {
                    continue;
                }

// Adds twice in boundary for internal border.
                topology.add_mesh_entity_boundary_relation(
                    gmme_id( Surface3D::type_name_static(), new_surface_id ),
                    cur_line.gmme() );
                topology.add_mesh_entity_boundary_relation(
                    gmme_id( Surface3D::type_name_static(), new_surface_id ),
                    cur_line.gmme() );

                good_line = true;
                break;
            }
            if( !good_line ) {
                continue;
            }
        }

        split_merged_surface( new_surface_id, side_name, sided_interface_id,
            mesh_entities_to_erase, region_index );
    }

    void DuplicateInterfaceBuilder::split_merged_surface(
        index_t new_surface_id,
        const std::string& side_name,
        index_t sided_interface_id,
        std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        index_t region_index )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const Surface3D& cur_surface = geomodel_.surface( new_surface_id );
        const gmge_id& sided_interface_gme_t =
            interface( sided_interface_id ).gmge();
        // ========= bad copy paste from geo geomodel repair
        std::set< index_t > cutting_lines;
        for( index_t l = 0; l < cur_surface.nb_boundaries(); ++l ) {
            const Line3D& L = geomodel_.line(
                cur_surface.boundary_gmme( l ).index() );
            if( /*to_remove.count( L.gme_id() ) == 0 &&*/L.is_inside_border(
                cur_surface ) ) {
                cutting_lines.insert( L.index() );
            }
        }

        for( std::set< index_t >::iterator it = cutting_lines.begin();
            it != cutting_lines.end(); ++it ) {
// Force the recomputing of the geomodel vertices
// before performing the cut.
//                    geomodel_.mesh.vertices.clear() ;
            geometry.disconnect_surface_polygons_along_line_edges(
                cur_surface.index(), *it );
        }

        // cutting_lines std::set contains only the lines to get the different
        // connected components and not the other lines (such as borders of the entire
        // surface).

        GEO::vector< index_t > components;
        /// TODO to avoid to use the low_level_mesh_storage we may create
        /// a Mesh2D, work on it and then assign the meshes to the new surfaces?
        /// To discuss.
        /// But we need to find the meshtype to create (from existing surface type?).
        index_t nb_connected_components =
            cur_surface.low_level_mesh_storage().get_connected_components(
                components );
        if( nb_connected_components == 1 ) {
            DEBUG( "ONE CONNECTED COMPONENT" );
            ringmesh_assert( cur_surface.nb_parents() == 0 );
            geology.add_parent_children_relation( sided_interface_gme_t,
                cur_surface.gmme() );
// boundary informations are defined in build_merged_surfaces
// expected for the universe.
            /*add_entity_in_boundary( new_surface_gme_t,
             gme_t( GME::REGION, region_index ) ) ;*/
            bool side = ( side_name == "_plus" );
            /*add_entity_boundary( gme_t( GME::REGION, region_index ),
             new_surface_gme_t, side ) ;*/

// Add to universe (other side of the surface)
            topology.add_universe_boundary( new_surface_id, !side );

            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][new_surface_id] = 0;

// Cut for internal border. Some of the cutting lines are on the
// surface border or outside the surface. THAT MAY PROVOKE SOME
// PROBLEMS? TO CHECK!
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                DEBUG( "RECUT FOR INTERNAL BORDERS" );
                DEBUG( new_surface_id );
                geometry.cut_surface_by_line( new_surface_id, *it );
            }
            if( !cutting_lines.empty() ) {
/// TODO TO READD. COMMENTED TO COMPILE BUT WE MUST HANDLE THAT!
//const_cast< GEO::MeshVertices& >( geomodel_.surface( new_surface_id ).gfx_mesh().vertices ).remove_isolated();
            }
            return;
        }
        DEBUG( "MORE THAN ONE CONNECTED COMPONENT" );
        ringmesh_assert( nb_connected_components != 0 );

        // We create a new surface for each connected component
        // (since a RINGMesh::Surface has only one connected component).
        std::vector< std::vector< vec3 > > all_points( nb_connected_components,
            std::vector< vec3 >() );
        for( index_t all_points_itr = 0; all_points_itr < all_points.size();
            ++all_points_itr ) {
            all_points[all_points_itr].reserve( cur_surface.nb_vertices() ); // a little big but to avoid copy
        }
        for( index_t components_itr = 0; components_itr < components.size();
            ++components_itr ) {
            ringmesh_assert( components[components_itr] < all_points.size() );
            for( index_t v_in_f_itr = 0;
                v_in_f_itr < cur_surface.nb_mesh_element_vertices( components_itr );
                ++v_in_f_itr ) {
                index_t v_id = cur_surface.mesh_element_vertex_index(
                    ElementLocalVertex( components_itr, v_in_f_itr ) );
                all_points[components[components_itr]].push_back(
                    cur_surface.vertex( v_id ) );
            }
        }

        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
        for( index_t all_points_itr = 0; all_points_itr < all_points.size();
            ++all_points_itr ) {
            NNSearch3D make_unique_surf( all_points[all_points_itr] );
            std::vector< index_t > unique_id;
            std::vector< vec3 > facet_points;
            index_t not_used = 0;
            std::tie( not_used, unique_id, facet_points ) =
                make_unique_surf.get_colocated_index_mapping_and_unique_points(
                    geomodel_.epsilon() );

            std::vector< index_t > facet_indices;
            std::vector< index_t > facet_ptr;
            index_t count_facet_vertices = 0;
            facet_ptr.push_back( count_facet_vertices );
            index_t offset = 0;

            for( index_t components_itr = 0; components_itr < components.size();
                ++components_itr ) {
                ringmesh_assert( components[components_itr] < all_points.size() );
                if( components[components_itr] == all_points_itr ) {
                    for( index_t v_in_f_itr = 0;
                        v_in_f_itr
                            < cur_surface.nb_mesh_element_vertices( components_itr );
                        ++v_in_f_itr ) {
                        facet_indices.push_back( unique_id[offset] );
                        ++offset;
                    }
                    count_facet_vertices += cur_surface.nb_mesh_element_vertices(
                        components_itr );
                    facet_ptr.push_back( count_facet_vertices );
                }
            }
            gmme_id new_new_surface_gme_t = topology.create_mesh_entity< Surface >();
/// TODO weird that the GEO::vector cannot be put as parameter
/// of set_surface_geometry since it is a std::vector... to see.
            geometry.set_surface_geometry( new_new_surface_gme_t.index(),
                facet_points, facet_indices, facet_ptr );

// THE LINE BELOW IS REPLACED BY THE NEXT 2 LINES.
//set_entity_parent( new_new_surface_gme_t, sided_interface_gme_t ) ;
            ringmesh_assert(
                geomodel_.surface( new_new_surface_gme_t.index() ).nb_parents()
                    == 0 );
            geology.add_parent_children_relation( sided_interface_gme_t,
                new_new_surface_gme_t );
            bool side = ( side_name == "_plus" );
            topology.add_mesh_entity_boundary_relation(
                gmme_id( Region3D::type_name_static(), region_index ),
                new_new_surface_gme_t, side );

// Add to universe (other side of the surface)
            topology.add_universe_boundary( new_new_surface_gme_t.index(), !side );
            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )].push_back( 0 );

#ifdef RINGMESH_DEBUG
// In theory there is no isolated vertex
            index_t previous =
                geomodel_.surface( new_new_surface_gme_t.index() ).nb_vertices();
            geometry.delete_mesh_entity_isolated_vertices( new_new_surface_gme_t );
            ringmesh_assert(
                previous
                    == geomodel_.surface( new_new_surface_gme_t.index() ).nb_vertices() );
#endif
            geometry.clear_geomodel_mesh(); // TODO check if rebuild is needed (as before: recompute geomodel mesh)
// HANDLE THE INTERNAL BORDER
            Surface3D& new_new_surf = const_cast< Surface3D& >( geomodel_.surface(
                new_new_surface_gme_t.index() ) );
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                const Line3D& cur_cutting_line = geomodel_.line( *it );

                index_t facet_index = NO_ID;
                index_t edge_index = NO_ID;
                if( !find_facet_and_edge( new_new_surf.polygon_nn_search(),
                    new_new_surf,
                    gmmv.geomodel_vertex_id( cur_cutting_line.gmme(), 0 ),
                    gmmv.geomodel_vertex_id( cur_cutting_line.gmme(), 1 ),
                    facet_index, edge_index ) ) {
                    continue;
                }

                if( new_new_surf.polygon_adjacent_index(
                    PolygonLocalEdge( facet_index, edge_index ) ) == NO_ID ) {
                    continue;
                }

                DEBUG( "RECUT FOR INTERNAL BORDERS 2" );
                DEBUG( new_new_surf.index() );
                geometry.cut_surface_by_line( new_new_surf.index(),
                    cur_cutting_line.index() );
                /*if( !cutting_lines.empty() ) {
                 const_cast< GEO::MeshVertices& >( geomodel_.surface(
                 new_new_surf.index() ).gfx_mesh().vertices ).remove_isolated() ;
                 }*/
            }
            if( !cutting_lines.empty() ) {
/// TODO TO READD. COMMENTED TO COMPILE BUT WE MUST HANDLE THAT!
                /*const_cast< GEO::MeshVertices& >( geomodel_.surface(
                 new_new_surf.index() ).gfx_mesh().vertices ).remove_isolated();*/
            }
        }
    }

    void DuplicateInterfaceBuilder::initialize_translation_attributes(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        for( index_t reg_itr = 0; reg_itr < geomodel_.nb_regions(); ++reg_itr ) {
            const Region3D& reg = geomodel_.region( reg_itr );
            if( !reg.is_meshed() ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] != NO_ID );
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] == 0 );
            GEO::AttributesManager& att_mgr = reg.vertex_attribute_manager();
            GEO::Attribute< vec3 > translation_attr( att_mgr,
                translation_attribute_name_ );
            translation_attr.fill( vec3( 0., 0., 0. ) );
        }

        for( index_t surf_itr = 0; surf_itr < geomodel_.nb_surfaces(); ++surf_itr ) {
            const Surface3D& surf = geomodel_.surface( surf_itr );

            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][surf_itr] == NO_ID ) {
                continue;
            }

            GEO::AttributesManager& att_mgr = surf.vertex_attribute_manager();
            GEO::Attribute< vec3 > translation_attr( att_mgr,
                translation_attribute_name_ );
            translation_attr.fill( vec3( 0., 0., 0. ) );
        }
    }

    /// TODO that should not be necessary anymore since there is a new functionality
    /// to compute normals.... to see.
    void DuplicateInterfaceBuilder::save_normals_on_one_new_interface(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        index_t interface_id ) const
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeoModelGeologicalEntity3D& interface_gme = interface( interface_id );
        for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
            ++child_itr ) {
            const GeoModelMeshEntity3D& cur_child = interface_gme.child( child_itr );
            ringmesh_assert(
                cur_child.mesh_entity_type() == Surface3D::type_name_static() );
// As the loop begins at the first new interface, no surface
// met in this loop should be to delete.
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Surface3D::type_name_static() )][cur_child.index()] == 0 );
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Surface3D::type_name_static() )][cur_child.index()] != NO_ID );
            const Surface3D& cur_surface = geomodel_.surface( cur_child.index() );
            save_normal_on_one_surface( cur_surface );
        }
    }

    void DuplicateInterfaceBuilder::save_normal_on_one_surface(
        const Surface3D& surface ) const
    {
        GEO::AttributesManager& att_mgr = surface.vertex_attribute_manager();
        GEO::Attribute< vec3 > normal_attr( att_mgr, normal_attribute_name_ );
        for( index_t v_i = 0; v_i < surface.nb_vertices(); ++v_i ) {
            normal_attr[v_i] = surface.low_level_mesh_storage().normal_at_vertex(
                v_i );
        }
    }

    vec3 DuplicateInterfaceBuilder::get_local_translation_normal(
        const Surface3D& surface,
        index_t vertex_id_in_surface ) const
    {
        // only one side for the sided interface
        ringmesh_assert( surface.nb_incident_entities() == 1 );
        const GeoModelMeshEntity3D& incident_ent = surface.incident_entity( 0 );
        ringmesh_assert( incident_ent.type_name() == Region3D::type_name_static() );
        const Region3D& cur_reg = geomodel_.region( incident_ent.index() );
        index_t local_surf_id = find_local_boundary_id( cur_reg, surface );
        bool side = cur_reg.side( local_surf_id );

        GEO::AttributesManager& att_mgr = surface.vertex_attribute_manager();
        GEO::Attribute< vec3 > normal_attr( att_mgr, normal_attribute_name_ );
        vec3 normal = normal_attr[vertex_id_in_surface];

        ringmesh_assert( std::abs( normal.length() - 1. ) < global_epsilon );
        if( !side ) {
            normal *= -1;
        }
        return normal;
    }

    vec3 DuplicateInterfaceBuilder::get_local_translation_vector(
        const vec3& normal ) const
    {
        vec3 displacement = normal * 1.5 * geomodel_.epsilon(); //global_epsilon ;
        return displacement;
    }

    void DuplicateInterfaceBuilder::compute_translation_vectors_duplicated_fault_network_surfaces_and_regions(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        const std::vector< std::vector< index_t > >& geological_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;
        initialize_translation_attributes( mesh_entities_to_erase );

        ringmesh_assert(
            geomodel_.nb_geological_entities( Interface3D::type_name_static() )
                - first_new_interface_index >= 2 );
        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr
                < geomodel_.nb_geological_entities( Interface3D::type_name_static() );
            ++new_interface_itr ) {
#ifdef RINGMESH_DEBUG
// Only the new interfaces with no child are removed.
// Such interfaces have no child because there are entirely in
// the boundary of the geomodel (no surface not voi).
            if( geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                Interface3D::type_name_static() )][new_interface_itr] == NO_ID ) {
                const GeoModelGeologicalEntity3D& interface_gme =
                    geomodel_.geological_entity( Interface3D::type_name_static(),
                        new_interface_itr );
                ringmesh_assert( interface_gme.nb_children() == 0 );
            } else {
                ringmesh_assert(
                    geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                        Interface3D::type_name_static() )][new_interface_itr] == 0 );
            }
#endif

            const GeoModelGeologicalEntity3D& interface_gme =
                geomodel_.geological_entity( Interface3D::type_name_static(),
                    new_interface_itr );
//            save_normals_on_one_new_interface( to_erase_by_type, interface_gme ) ;

            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GeoModelMeshEntity3D& cur_child = interface_gme.child(
                    child_itr );
                ringmesh_assert(
                    cur_child.mesh_entity_type() == Surface3D::type_name_static() );
                ringmesh_assert(
                    mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                        Surface3D::type_name_static() )][cur_child.index()]
                        != NO_ID );
                const Surface3D& cur_surface = geomodel_.surface(
                    cur_child.index() ); // avoid dynamic_cast of cur_child

                ringmesh_assert( cur_surface.nb_incident_entities() == 1 );
                ringmesh_assert(
                    cur_surface.incident_entity( 0 ).type_name()
                        == Region3D::type_name_static() );
                GEO::Attribute< index_t > id_in_link_vector(
                    geomodel_.surface( cur_surface.index() ).vertex_attribute_manager(),
                    "id_in_link_vector" );

                for( index_t surf_vertex_itr = 0;
                    surf_vertex_itr < cur_surface.nb_vertices();
                    ++surf_vertex_itr ) {

                    const vec3 local_translation_normal =
                        get_local_translation_normal( cur_surface, surf_vertex_itr );
                    const vec3 local_translation_vector =
                        get_local_translation_vector( local_translation_normal );

                    index_t id_in_link = id_in_link_vector[surf_vertex_itr];

                    ringmesh_assert( gme_vertices_links_[id_in_link] != nullptr );
                    gme_vertices_links_[id_in_link]->displace(
                        local_translation_vector );
                }
            }
        }
    }

    void DuplicateInterfaceBuilder::compute_translation_vectors_duplicated_fault_network(
        index_t first_new_interface_index,
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase,
        const std::vector< std::vector< index_t > >& geological_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const GeologicalTypeManager geological_entity_type_manager =
            geomodel_.entity_type_manager().geological_entity_manager;

        initialize_translation_attributes( mesh_entities_to_erase );

        int step_to_other_side = 1;
        ringmesh_assert(
            geomodel_.nb_geological_entities( Interface3D::type_name_static() )
                - first_new_interface_index >= 2 );
        for( index_t new_interface_itr = first_new_interface_index;
            new_interface_itr
                < geomodel_.nb_geological_entities( Interface3D::type_name_static() );
            ++new_interface_itr ) {
            ringmesh_assert(
                geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                    Interface3D::type_name_static() )][new_interface_itr] != NO_ID );
            ringmesh_assert(
                geological_entities_to_erase[geological_entity_type_manager.geological_entity_type_index(
                    Interface3D::type_name_static() )][new_interface_itr] == 0 );

            const GeoModelGeologicalEntity3D& interface_gme =
                geomodel_.geological_entity( Interface3D::type_name_static(),
                    new_interface_itr );
            ringmesh_assert(
                new_interface_itr + step_to_other_side
                    >= first_new_interface_index );
            ringmesh_assert(
                new_interface_itr + step_to_other_side
                    < geomodel_.nb_geological_entities(
                        Interface3D::type_name_static() ) );
/// @bug in the case that there is no other side (fault = geomodel boundary)
/// this method (trick) does not work.
            const GeoModelGeologicalEntity3D& other_side_interface_gme =
                geomodel_.geological_entity( Interface3D::type_name_static(),
                    new_interface_itr + step_to_other_side );
            step_to_other_side *= -1;

            save_normals_on_one_new_interface( mesh_entities_to_erase,
                interface_gme.index() );

// Clear to take into account the new gme in the geomodel.
            geometry.clear_geomodel_mesh(); // not done in geomodel_vertex_id // TODO check if rebuild is needed (as before: recompute geomodel mesh)

            const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;

            for( index_t child_itr = 0; child_itr < interface_gme.nb_children();
                ++child_itr ) {
                const GeoModelMeshEntity3D& cur_child = interface_gme.child(
                    child_itr );
                ringmesh_assert(
                    cur_child.mesh_entity_type() == Surface3D::type_name_static() );
                ringmesh_assert(
                    mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                        Surface3D::type_name_static() )][cur_child.index()]
                        != NO_ID );
                const Surface3D& cur_surface = geomodel_.surface(
                    cur_child.index() ); // avoid dynamic_cast of cur_child

                ringmesh_assert( cur_surface.nb_incident_entities() == 1 );
                ringmesh_assert(
                    cur_surface.incident_entity( 0 ).type_name()
                        == Region3D::type_name_static() );
                for( index_t surf_vertex_itr = 0;
                    surf_vertex_itr < cur_surface.nb_vertices();
                    ++surf_vertex_itr ) {
                    const index_t vertex_id_in_gmm = gmmv.geomodel_vertex_id(
                        cur_surface.gmme(), surf_vertex_itr );

/// @todo potentially for the vertices in common between 2 surfaces
/// the translation is applied twice. Does not seem to be
// a problem. To check.

// Gets all the GME with a vertex colocated to the one of vertex_id_in_gmm
                    const std::vector< GMEVertex >& gme_vertices = gmmv.gme_vertices(
                        vertex_id_in_gmm );
                    ringmesh_assert(
                        std::find( gme_vertices.begin(), gme_vertices.end(),
                            GMEVertex(
                                gmme_id( Surface3D::type_name_static(),
                                    cur_child.index() ), surf_vertex_itr ) )
                            != gme_vertices.end() );

                    const vec3 local_translation_normal =
                        get_local_translation_normal( cur_surface, surf_vertex_itr );
                    const vec3 local_translation_vector =
                        get_local_translation_vector( local_translation_normal );

                    for( index_t gme_vertex_itr = 0;
                        gme_vertex_itr < gme_vertices.size(); ++gme_vertex_itr ) {
                        const gmme_id& cur_gmme_id =
                            gme_vertices[gme_vertex_itr].gmme;
                        if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                            cur_gmme_id.type() )][cur_gmme_id.index()] == NO_ID ) {
                            continue; // It is an old entity to remove.
                        }
                        ringmesh_assert(
                            mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                                cur_gmme_id.type() )][cur_gmme_id.index()] == 0 );

                        if( cur_gmme_id.type() != Surface3D::type_name_static()
                            && cur_gmme_id.type() != Region3D::type_name_static() ) {
                            continue;
                        }

                        if( is_surface_or_region_on_the_right_side_of_the_fault(
                            cur_gmme_id, local_translation_normal,
                            gme_vertices[gme_vertex_itr].v_index,
                            geomodel_.mesh.vertices.vertex( vertex_id_in_gmm ),
                            interface_gme.index(),
                            other_side_interface_gme.index() ) ) {
                            store_displacement_in_gme(
                                geomodel_.mesh_entity( cur_gmme_id ),
                                gme_vertices[gme_vertex_itr].v_index,
                                local_translation_vector );
                        }
                    }
                }
            }
        }
    }

    bool DuplicateInterfaceBuilder::is_surface_or_region_on_the_right_side_of_the_fault(
        const gmme_id& cur_gmme_t,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_gmme,
        const vec3& vertex_pos,
        index_t interface_id,
        index_t other_side_interface_id ) const
    {
        const GeoModelGeologicalEntity3D& interface_gme = interface( interface_id );
        const GeoModelGeologicalEntity3D& other_side_interface_gme = interface(
            other_side_interface_id );
        if( cur_gmme_t.type() == Region3D::type_name_static() ) {
            if( !geomodel_.region( cur_gmme_t.index() ).is_meshed() ) {
                return false;
            }
            if( !is_region_on_right_side_of_sided_interface( cur_gmme_t.index(),
                normal_on_vertex_interface, vertex_id_in_gmme, vertex_pos ) ) {
// Region on the other side of the fault
                return false;
            }
        } else {
            ringmesh_assert( cur_gmme_t.type() == Surface3D::type_name_static() );
            for( index_t interface__child_itr = 0;
                interface__child_itr < interface_gme.nb_children();
                ++interface__child_itr ) {
                const GeoModelMeshEntity3D& cur_child = interface_gme.child(
                    interface__child_itr );
                ringmesh_assert(
                    cur_child.mesh_entity_type() == Surface3D::type_name_static() );
                if( cur_child.index() == cur_gmme_t.index() ) {
                    return true;
                }
            }

            for( index_t other_side_interface__child_itr = 0;
                other_side_interface__child_itr
                    < other_side_interface_gme.nb_children();
                ++other_side_interface__child_itr ) {
                const GeoModelMeshEntity3D& cur_other_side_child =
                    other_side_interface_gme.child(
                        other_side_interface__child_itr );
                ringmesh_assert(
                    cur_other_side_child.mesh_entity_type()
                        == Surface3D::type_name_static() );
                if( cur_other_side_child.index() == cur_gmme_t.index() ) {
                    return false;
                }
            }

            if( !is_surface_on_right_side_of_sided_interface( cur_gmme_t.index(),
                normal_on_vertex_interface, vertex_id_in_gmme, vertex_pos ) ) {
                return false;
            }
        }
        return true;
    }

    bool DuplicateInterfaceBuilder::is_region_on_right_side_of_sided_interface(
        index_t region_to_check_id,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_region,
        const vec3& vertex_pos ) const
    {
        ringmesh_assert( region_to_check_id < geomodel_.nb_regions() );

        const Region3D& region_to_check = geomodel_.region( region_to_check_id );
        std::vector< index_t > cells_around = region_to_check.cells_around_vertex(
            vertex_id_in_region, NO_ID );
//            false ) ;
        ringmesh_assert( !cells_around.empty() );

        vec3 region_to_check_mean_normal_on_vertex( 0., 0., 0. );
        for( index_t cells_around_itr = 0; cells_around_itr < cells_around.size();
            ++cells_around_itr ) {

            index_t cur_cell_id_in_region = cells_around[cells_around_itr];
            vec3 cur_cell_barycenter = region_to_check.mesh_element_barycenter(
                cur_cell_id_in_region );
            vec3 p = cur_cell_barycenter - vertex_pos;
            region_to_check_mean_normal_on_vertex += p;
        }
        ringmesh_assert(
            std::abs( region_to_check_mean_normal_on_vertex.x ) > global_epsilon
                || std::abs( region_to_check_mean_normal_on_vertex.y )
                    > global_epsilon
                || std::abs( region_to_check_mean_normal_on_vertex.z )
                    > global_epsilon );

        if( GEO::dot( normal_on_vertex_interface,
            region_to_check_mean_normal_on_vertex ) > global_epsilon ) {
            return true;
        }
        return false;
    }

    bool DuplicateInterfaceBuilder::is_surface_on_right_side_of_sided_interface(
        index_t surface_to_check_id,
        const vec3& normal_on_vertex_interface,
        index_t vertex_id_in_surface,
        const vec3& vertex_pos ) const
    {
        ringmesh_assert( surface_to_check_id < geomodel_.nb_surfaces() );

        const Surface3D& surface_to_check = geomodel_.surface( surface_to_check_id );
        ringmesh_assert(
            surface_to_check.nb_incident_entities() == 1
                || surface_to_check.nb_incident_entities() == 2 );

        std::vector< index_t > polygons_around =
            surface_to_check.low_level_mesh_storage().polygons_around_vertex( vertex_id_in_surface, false, NO_ID );
        ringmesh_assert( !polygons_around.empty() );

        vec3 surf_to_check_mean_normal_on_vertex( 0., 0., 0. );
        for( index_t facets_around_itr = 0;
            facets_around_itr < polygons_around.size(); ++facets_around_itr ) {
            index_t cur_facet_id_in_surf = polygons_around[facets_around_itr];
            vec3 cur_facet_barycenter = surface_to_check.mesh_element_barycenter(
                cur_facet_id_in_surf );
            vec3 p = cur_facet_barycenter - vertex_pos;
            surf_to_check_mean_normal_on_vertex += p;
        }
        ringmesh_assert(
            std::abs( surf_to_check_mean_normal_on_vertex.x ) > global_epsilon
                || std::abs( surf_to_check_mean_normal_on_vertex.y ) > global_epsilon
                || std::abs( surf_to_check_mean_normal_on_vertex.z )
                    > global_epsilon );

        if( GEO::dot( 50 * normal_on_vertex_interface, /// @todo remove this 50
        surf_to_check_mean_normal_on_vertex ) > global_epsilon ) {
            return true;
        }
        return false;
    }

    void DuplicateInterfaceBuilder::store_displacement_in_gme(
        const GeoModelMeshEntity3D& gmme,
        index_t vertex_id_in_gmme,
        const vec3& translation ) const
    {
        GEO::AttributesManager& att_mgr = gmme.vertex_attribute_manager();
        GEO::Attribute< vec3 > translation_attr( att_mgr,
            translation_attribute_name_ );
        translation_attr[vertex_id_in_gmme] += translation;
    }

    void DuplicateInterfaceBuilder::translate_duplicated_fault_network(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        for( index_t reg_itr = 0; reg_itr < geomodel_.nb_regions(); ++reg_itr ) {
            const Region3D& reg = geomodel_.region( reg_itr );
            if( !reg.is_meshed() ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] != NO_ID );
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    Region3D::type_name_static() )][reg_itr] == 0 );
            GEO::AttributesManager& att_mgr = reg.vertex_attribute_manager();
            GEO::Attribute< vec3 > translation_attr( att_mgr,
                translation_attribute_name_ );
            for( index_t vertex_itr = 0; vertex_itr < reg.nb_vertices();
                ++vertex_itr ) {

                vec3 new_vertex_coords = reg.vertex( vertex_itr );
                new_vertex_coords += translation_attr[vertex_itr];
                geometry.set_mesh_entity_vertex( reg.gmme(), vertex_itr,
                    new_vertex_coords, false );
            }
        }

        for( index_t surf_itr = 0; surf_itr < geomodel_.nb_surfaces(); ++surf_itr ) {
            const Surface3D& surf = geomodel_.surface( surf_itr );

            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                Surface3D::type_name_static() )][surf_itr] == NO_ID ) {
                continue;
            }

            GEO::AttributesManager& att_mgr = surf.vertex_attribute_manager();
            GEO::Attribute< vec3 > translation_attr( att_mgr,
                translation_attribute_name_ );
            for( index_t vertex_itr = 0; vertex_itr < surf.nb_vertices();
                ++vertex_itr ) {

                vec3 new_vertex_coords = surf.vertex( vertex_itr );
                new_vertex_coords += translation_attr[vertex_itr];
                geometry.set_mesh_entity_vertex( surf.gmme(), vertex_itr,
                    new_vertex_coords, false );
            }
        }
    }

    void DuplicateInterfaceBuilder::set_no_displacement_on_fault_real_extension(
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const GeoModelMeshVertices3D& gmmv = geomodel_.mesh.vertices;
        for( index_t line_itr = 0; line_itr < geomodel_.nb_lines(); ++line_itr ) {
            const Line3D& cur_line = geomodel_.line( line_itr );
            if( cur_line.nb_incident_entities() != 1 ) {
                continue;
            }
// cur_line.nb_incident_entity() == 1 means fault extension
            ringmesh_assert(
                cur_line.incident_entity( 0 ).type_name()
                    == Surface3D::type_name_static() );
            ringmesh_assert(
                GeoModelGeologicalEntity3D::is_fault(
                    cur_line.incident_entity( 0 ).parent(
                        Interface3D::type_name_static() ).geological_feature() ) );

            ringmesh_assert( cur_line.nb_vertices() >= 2 );
// Vertices not corner
            for( index_t line_vertex_itr = 1;
                line_vertex_itr < cur_line.nb_vertices() - 1; ++line_vertex_itr ) {

                index_t vertex_id_in_gmm = gmmv.geomodel_vertex_id( cur_line.gmme(),
                    line_vertex_itr );
                set_no_displacement_on_gme_sharing_vertex( vertex_id_in_gmm,
                    mesh_entities_to_erase );
            }
// Line corners
            const Corner3D& first_corner = geomodel_.corner(
                cur_line.boundary_gmme( 0 ).index() );
            if( !displace_corner( first_corner, cur_line ) ) {
                index_t vertex_id_in_gmm = gmmv.geomodel_vertex_id(
                    first_corner.gmme(), 0 );
                set_no_displacement_on_gme_sharing_vertex( vertex_id_in_gmm,
                    mesh_entities_to_erase );
            }
            const Corner3D& second_corner = geomodel_.corner(
                cur_line.boundary_gmme( 1 ).index() );
            if( !displace_corner( second_corner, cur_line ) ) {
                index_t vertex_id_in_gmm = gmmv.geomodel_vertex_id(
                    second_corner.gmme(), 0 );
                set_no_displacement_on_gme_sharing_vertex( vertex_id_in_gmm,
                    mesh_entities_to_erase );
            }
        }
    }

    bool DuplicateInterfaceBuilder::displace_corner(
        const Corner3D& corner,
        const Line3D& line_one_incident_entity ) const
    {
        ringmesh_assert( line_one_incident_entity.nb_incident_entities() == 1 );
        ringmesh_assert(
            GeoModelGeologicalEntity3D::is_fault(
                line_one_incident_entity.incident_entity( 0 ).parent(
                    Interface3D::type_name_static() ).geological_feature() ) );
        ringmesh_assert(
            line_one_incident_entity.boundary( 0 ).index() == corner.index()
                || line_one_incident_entity.boundary( 1 ).index()
                    == corner.index() );
        for( index_t corner_incident_entities_itr = 0;
            corner_incident_entities_itr < corner.nb_incident_entities();
            ++corner_incident_entities_itr ) {
            const GeoModelMeshEntity3D& cur_line_gme = corner.incident_entity(
                corner_incident_entities_itr );
            ringmesh_assert(
                cur_line_gme.type_name() == Line3D::type_name_static() );
            if( cur_line_gme.index() == line_one_incident_entity.index() ) {
                continue;
            }
            for( index_t surface_incident_ent_itr = 0;
                surface_incident_ent_itr < cur_line_gme.nb_incident_entities();
                ++surface_incident_ent_itr ) {

                if( !GeoModelGeologicalEntity3D::is_fault(
                    cur_line_gme.incident_entity( surface_incident_ent_itr ).parent(
                        Interface3D::type_name_static() ).geological_feature() ) ) {
                    continue;
                }

                if( cur_line_gme.incident_entity( surface_incident_ent_itr ).parent(
                    Interface3D::type_name_static() ).index()
                    == line_one_incident_entity.incident_entity( 0 ).parent(
                        Interface3D::type_name_static() ).index() ) {
                    continue;
                }
                return true;
            }
        }
        return false;
    }

    void DuplicateInterfaceBuilder::set_no_displacement_on_gme_sharing_vertex(
        index_t vertex_id_in_gmm,
        const std::vector< std::vector< index_t > >& mesh_entities_to_erase )
    {
        const MeshEntityTypeManager3D mesh_entity_type_manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const std::vector< GMEVertex >& gme_vertices =
            geomodel_.mesh.vertices.gme_vertices( vertex_id_in_gmm );

        for( index_t gme_vertex_itr = 0; gme_vertex_itr < gme_vertices.size();
            ++gme_vertex_itr ) {

            const GMEVertex& cur_gme_vertex = gme_vertices[gme_vertex_itr];
            if( cur_gme_vertex.gmme.type() != Surface3D::type_name_static()
                && cur_gme_vertex.gmme.type() != Region3D::type_name_static() ) {
                continue;
            }
            if( mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                cur_gme_vertex.gmme.type() )][cur_gme_vertex.gmme.index()]
                == NO_ID ) {
                continue;
            }
            ringmesh_assert(
                mesh_entities_to_erase[mesh_entity_type_manager.mesh_entity_type_index(
                    cur_gme_vertex.gmme.type() )][cur_gme_vertex.gmme.index()]
                    == 0 );
            const GeoModelMeshEntity3D& gmme = geomodel_.mesh_entity(
                cur_gme_vertex.gmme );
            GEO::AttributesManager& att_mgr = gmme.vertex_attribute_manager();
            GEO::Attribute< vec3 > translation_attr( att_mgr,
                translation_attribute_name_ );
            translation_attr[cur_gme_vertex.v_index] = vec3( 0, 0, 0 );
        }
    }

    void DuplicateInterfaceBuilder::remove_gap()
    {
        // TODO
    }
}
