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

#include <ringmesh/geomodel/geomodel_builder_repair.h>

#include <ringmesh/geomodel/geomodel_builder.h>

/*!
 * @file Implementation of repair function of the surfaces of a GeoModel
 * @author Jeanne Pellerin
 * @author Pierre Anquez
 */

namespace RINGMesh {

    GeoModelBuilderRepair::GeoModelBuilderRepair(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderRepair::repair( RepairMode repair_mode )
    {
        switch( repair_mode ) {
            case ALL:
                geomodel_mesh_repair();
                break;
            case BASIC:
                builder_.end_geomodel();
                break;
            case COLOCATED_VERTICES:
                remove_colocated_entity_vertices_and_update_gm();
                break;
            case DEGENERATE_FACETS_EDGES:
                remove_degenerate_polygons_and_edges_and_update_gm();
                break;
            case LINE_BOUNDARY_ORDER:
                repair_line_boundary_vertex_order();
                break;
            case CONTACTS:
                build_contacts();
                break;
            default:
                ringmesh_assert_not_reached;
        }
    }

    void GeoModelBuilderRepair::geomodel_mesh_repair()
    {
        // Remove colocated vertices in each entity
        remove_colocated_entity_vertices_and_update_gm();

        // Basic mesh repair for surfaces and lines
        remove_degenerate_polygons_and_edges_and_update_gm();

        // Proper reordering of line boundaries
        repair_line_boundary_vertex_order();

        // This is basic requirement ! no_colocated geomodel vertices !
        // So remove them if there are any
        geomodel_.mesh.remove_colocated_vertices();

        // Builds the contacts
        build_contacts();

        builder_.end_geomodel();
    }

    void GeoModelBuilderRepair::remove_colocated_entity_vertices_and_update_gm()
    {
        std::set< gmme_id > empty_mesh_entities;
        std::set< gmge_id > empty_geological_entities;

        remove_colocated_entity_vertices( empty_mesh_entities );
        if( !empty_mesh_entities.empty() ) {
            builder_.topology.get_dependent_entities( empty_mesh_entities,
                empty_geological_entities );
            builder_.removal.remove_mesh_entities( empty_mesh_entities );
        }
    }

    void GeoModelBuilderRepair::remove_degenerate_polygons_and_edges_and_update_gm()
    {
        std::set< gmme_id > empty_mesh_entities;
        remove_degenerate_polygons_and_edges( empty_mesh_entities );
        /// TODO when it will works,
        /// use GeoModelBuilderRemoval::remove_entities_and_dependencies
        if( !empty_mesh_entities.empty() ) {
            builder_.removal.remove_mesh_entities( empty_mesh_entities );
        }

        // This is basic requirement ! no_colocated geomodel vertices !
        // So remove them if there are any
        geomodel_.mesh.remove_colocated_vertices();

        builder_.end_geomodel();
    }

    void GeoModelBuilderRepair::repair_line_boundary_vertex_order()
    {
        for( index_t line_itr = 0; line_itr < geomodel_.nb_lines(); ++line_itr ) {
            const Line& cur_line = geomodel_.line( line_itr );
            if( !cur_line.is_first_corner_first_vertex() ) {
                const index_t first_boundary_index = cur_line.boundary( 0 ).index();
                builder_.topology.set_mesh_entity_boundary( cur_line.gmme(), 0,
                    cur_line.boundary_gmme( 1 ).index() );
                builder_.topology.set_mesh_entity_boundary( cur_line.gmme(), 1,
                    first_boundary_index );
            }
        }
    }

    bool GeoModelBuilderRepair::polygon_is_degenerate(
        const Surface& surface,
        index_t polygon,
        std::vector< index_t >& colocated_vertices )
    {
        index_t nb_vertices = surface.nb_mesh_element_vertices( polygon );
        if( nb_vertices != 3 ) {
            index_t* vertices = (index_t*) alloca( nb_vertices * sizeof(index_t) );
            for( index_t lv = 0; lv < nb_vertices; ++lv ) {
                vertices[lv] = colocated_vertices[surface.mesh_element_vertex_index( polygon,
                    lv )];
            }
            std::sort( vertices, vertices + nb_vertices );
            return std::unique( vertices, vertices + nb_vertices )
                != vertices + nb_vertices;
        }
        index_t v1 = colocated_vertices[surface.mesh_element_vertex_index( polygon, 0 )];
        index_t v2 = colocated_vertices[surface.mesh_element_vertex_index( polygon, 1 )];
        index_t v3 = colocated_vertices[surface.mesh_element_vertex_index( polygon, 2 )];
        return v1 == v2 || v2 == v3 || v3 == v1;
    }

    std::vector< index_t > GeoModelBuilderRepair::surface_detect_degenerate_polygons(
        const Surface& surface,
        std::vector< index_t >& colocated_vertices )
    {
        std::vector< index_t > f_is_degenerate( surface.nb_mesh_elements() );
        for( index_t p = 0; p < surface.nb_mesh_elements(); ++p ) {
            f_is_degenerate[p] = polygon_is_degenerate( surface, p, colocated_vertices );
        }
        return f_is_degenerate;
    }

    index_t GeoModelBuilderRepair::detect_degenerate_polygons( const Surface& surface )
    {
        std::vector< index_t > colocated;
        const NNSearch& nn_search = surface.vertex_nn_search();
        std::tie( std::ignore, colocated ) = nn_search.get_colocated_index_mapping(
            geomodel_.epsilon() );

        std::vector< index_t > degenerate = surface_detect_degenerate_polygons( surface,
            colocated );
        return static_cast< index_t >( std::count( degenerate.begin(),
            degenerate.end(), 1 ) );
    }

    std::vector< bool > GeoModelBuilderRepair::line_detect_degenerate_edges(
        const Line& line,
        std::vector< index_t >& colocated_vertices )
    {
        std::vector< bool > e_is_degenerate( line.nb_mesh_elements() );
        for( index_t e = 0; e < line.nb_mesh_elements(); ++e ) {
            e_is_degenerate[e] = edge_is_degenerate( line, e, colocated_vertices );
        }
        return e_is_degenerate;
    }

    index_t GeoModelBuilderRepair::repair_line_mesh( const Line& line )
    {
        std::vector< index_t > colocated;
        const NNSearch& nn_search = line.vertex_nn_search();
        std::tie( std::ignore, colocated ) = nn_search.get_colocated_index_mapping(
            geomodel_.epsilon() );

        std::vector< bool > degenerate = line_detect_degenerate_edges( line,
            colocated );
        index_t nb = static_cast< index_t >( std::count( degenerate.begin(),
            degenerate.end(), 1 ) );
        /// We have a problem if some vertices are left isolated
        /// If we remove them here we can kill all indices correspondances
        builder_.geometry.delete_line_edges( line.index(), degenerate, false );
        return nb;
    }

    void GeoModelBuilderRepair::remove_degenerate_polygons_and_edges(
        std::set< gmme_id >& to_remove )
    {
        to_remove.clear();
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            const Line& line = geomodel_.line( i );
            index_t nb = repair_line_mesh( line );
            if( nb > 0 ) {
                Logger::out( "GeoModel", nb, " degenerated edges removed in LINE ",
                    i );
                // If the Line is set it to remove
                if( geomodel_.line( i ).nb_mesh_elements() == 0 ) {
                    to_remove.insert( geomodel_.line( i ).gmme() );
                }
            }
        }
        // The builder might be needed

        double epsilon_sq = geomodel_.epsilon() * geomodel_.epsilon();
        for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
            const Surface& surface = geomodel_.surface( i );
            index_t nb = detect_degenerate_polygons( surface );
            /// @todo Check if that cannot be simplified
            if( nb > 0 ) {
                // If there are some degenerated polygons
                // Using repair function of geogram
                // Warning - This triangulates the mesh
                if( surface.nb_vertices() > 0 ) {
                    // Colocated vertices must be processed before
                    // MESH_REPAIR_DUP_F 2 ;
                    GEO::MeshRepairMode mode =
                        static_cast< GEO::MeshRepairMode >( 2 );
                    std::unique_ptr< SurfaceMeshBuilder > builder =
                        builder_.geometry.create_surface_builder( i );
                    builder->mesh_repair( mode, 0.0 );

                    // This might create some small components - remove them
                    builder->remove_small_connected_components( epsilon_sq, 3 );

                    // Alright, this is a bit of an overkill [JP]
                    if( surface.nb_vertices() > 0 ) {
                        builder->mesh_repair( mode, 0.0 );
                    }
                }
                if( surface.nb_vertices() == 0 || surface.nb_mesh_elements() == 0 ) {
                    to_remove.insert( geomodel_.surface( i ).gmme() );
                }
            }
        }
    }

    std::set< index_t > GeoModelBuilderRepair::vertices_on_inside_boundary(
        const gmme_id& E_id )
    {
        std::set< index_t > vertices;
        if( E_id.type() == Corner::type_name_static() ) {
            return vertices;
        }
        const GeoModelMeshEntity& E = geomodel_.mesh_entity( E_id );
        if( E_id.type() == Line::type_name_static() ) {
            if( E.boundary( 0 ).is_inside_border( E ) ) {
                vertices.insert( E.nb_vertices() - 1 );
            }
            return vertices;
        }
        std::vector< const GeoModelMeshEntity* > inside_border;
        for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
            if( E.boundary( i ).is_inside_border( E ) ) {
                inside_border.push_back(
                    dynamic_cast< const GeoModelMeshEntity* >( &E.boundary( i ) ) );
            }
        }
        if( !inside_border.empty() ) {
            // We want to get the indices of the vertices in E
            // that are colocated with those of the inside boundary
            // We assume that the geomodel vertices are not computed
            const NNSearch& nn_search = E.vertex_nn_search();

            for( const GeoModelMeshEntity*& entity : inside_border ) {
                for( index_t v = 0; v < entity->nb_vertices(); ++v ) {
                    std::vector< index_t > colocated_indices =
                        nn_search.get_neighbors( entity->vertex( v ),
                            geomodel_.epsilon() );
                    if( colocated_indices.size() > 1 ) {
                        std::sort( colocated_indices.begin(),
                            colocated_indices.end() );
                        // Add colocated vertices except one to the duplicated vertices set
                        vertices.insert( colocated_indices.begin() + 1,
                            colocated_indices.end() );
                    }
                }
            }
        }
        return vertices;
    }

    void GeoModelBuilderRepair::remove_colocated_entity_vertices(
        std::set< gmme_id >& to_remove )
    {
        to_remove.clear();
        // For all Lines and Surfaces
        const MeshEntityType types[2] = { Line::type_name_static(),
                                          Surface::type_name_static() };
        for( index_t t = 0; t < 2; ++t ) {
            const MeshEntityType& T = types[t];

            for( index_t e = 0; e < geomodel_.nb_mesh_entities( T ); ++e ) {
                gmme_id entity_id( T, e );
                const GeoModelMeshEntity& E = geomodel_.mesh_entity( entity_id );

                const NNSearch& kdtree = E.vertex_nn_search();
                std::vector< index_t > colocated;
                std::tie( std::ignore, colocated ) =
                    kdtree.get_colocated_index_mapping( geomodel_.epsilon() );

                // Get the vertices to delete
                std::set< index_t > inside_border = vertices_on_inside_boundary(
                    entity_id );

                std::vector< bool > to_delete( colocated.size(), false );
                index_t nb_todelete = 0;
                for( index_t v = 0; v < colocated.size(); ++v ) {
                    if( colocated[v] == v
                        || inside_border.find( v ) != inside_border.end() ) {
                        // This point is kept
                        // No colocated or on an inside boundary
                    } else {
                        // The point is to remove
                        to_delete[v] = true;
                        nb_todelete++;
                    }
                }

                if( nb_todelete == 0 ) {
                    // Nothing to do there
                    continue;
                } else if( nb_todelete == E.nb_vertices() ) {
                    // The complete entity should be removed
                    to_remove.insert( E.gmme() );
                    continue;
                } else {
                    if( t == 1 ) {
                        std::unique_ptr< SurfaceMeshBuilder > builder =
                            builder_.geometry.create_surface_builder( e );
                        for( index_t p_itr = 0; p_itr < E.nb_mesh_elements();
                            p_itr++ ) {
                            for( index_t fpv_itr = 0;
                                fpv_itr < E.nb_mesh_element_vertices( p_itr );
                                fpv_itr++ ) {
                                builder->set_polygon_vertex( p_itr, fpv_itr,
                                    colocated[E.mesh_element_vertex_index( p_itr,
                                        fpv_itr )] );
                            }
                        }
                        builder->delete_vertices( to_delete );
                        Logger::out( "Repair", nb_todelete,
                            " colocated vertices deleted in ", entity_id );

                    } else if( t == 0 ) {
                        std::unique_ptr< LineMeshBuilder > builder =
                            builder_.geometry.create_line_builder( e );
                        for( index_t e_itr = 0; e_itr < E.nb_mesh_elements();
                            e_itr++ ) {
                            builder->set_edge_vertex( e_itr, 0,
                                colocated[E.mesh_element_vertex_index( e_itr, 0 )] );
                            builder->set_edge_vertex( e_itr, 1,
                                colocated[E.mesh_element_vertex_index( e_itr, 1 )] );
                        }
                        builder->delete_vertices( to_delete );
                        Logger::out( "Repair", nb_todelete,
                            " colocated vertices deleted in ", entity_id );
                    } else {
                        ringmesh_assert_not_reached;
                    }
                }
            }
        }
    }

    bool GeoModelBuilderRepair::edge_is_degenerate(
        const Line& line,
        index_t edge,
        const std::vector< index_t >& colocated_vertices )
    {
        index_t v1 = colocated_vertices[line.mesh_element_vertex_index( edge, 0 )];
        index_t v2 = colocated_vertices[line.mesh_element_vertex_index( edge, 1 )];
        return v1 == v2;
    }

    void GeoModelBuilderRepair::build_contacts()
    {
        builder_.geology.build_contacts();
    }
}
