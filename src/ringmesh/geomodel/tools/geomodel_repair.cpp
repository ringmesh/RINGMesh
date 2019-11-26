/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <array>

#include <geogram/basic/algorithm.h>

#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_repair.h>

#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>

/*!
 * @file Implementation of repair function of the surfaces of a GeoModel
 * @author Jeanne Pellerin
 * @author Pierre Anquez
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    class GeoModelRepair
    {
        ringmesh_disable_copy_and_move( GeoModelRepair );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        GeoModelRepair( GeoModel< DIMENSION >& geomodel )
            : builder_( geomodel ), geomodel_( geomodel )
        {
        }

        /*!
         * @param[in] repair_mode repair mode to apply.
         */
        void repair( RepairMode repair_mode )
        {
            switch( repair_mode )
            {
            case RepairMode::ALL:
                geomodel_mesh_repair();
                break;
            case RepairMode::BASIC:
                builder_.end_geomodel();
                break;
            case RepairMode::COLOCATED_VERTICES:
                remove_colocated_entity_vertices_and_update_gm();
                break;
            case RepairMode::DEGENERATE_POLYGONS_EDGES:
                remove_degenerate_polygons_and_edges_and_update_gm();
                break;
            case RepairMode::LINE_BOUNDARY_ORDER:
                repair_line_boundary_vertex_order();
                break;
            case RepairMode::CONTACTS:
                build_contacts();
                break;
            case RepairMode::ISOLATED_VERTICES:
                remove_isolated_vertices();
                break;
            default:
                ringmesh_assert_not_reached;
            }
            geomodel_.mesh.vertices.clear();
        }

        ~GeoModelRepair() = default;

    private:
        /*!
         * All implemented repair for a GeoModel.
         */
        void geomodel_mesh_repair()
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

            // Remove isolated vertices on mesh entities
            remove_isolated_vertices();

            builder_.end_geomodel();
        }
        /*!
         * Remove the colocated vertices in all the GeoModelMeshEntities within
         * the GeoModel. GeoModelMeshEntities without any vertex anymore
         * (after the removal of the vertices) are removed off the GeoModel.
         */
        void remove_colocated_entity_vertices_and_update_gm()
        {
            std::set< gmme_id > empty_mesh_entities;
            std::set< gmge_id > empty_geological_entities;

            remove_colocated_entity_vertices( empty_mesh_entities );
            if( !empty_mesh_entities.empty() )
            {
                builder_.topology.get_dependent_entities(
                    empty_mesh_entities, empty_geological_entities );
                builder_.remove.remove_mesh_entities( empty_mesh_entities );
            }
        }

        /*!
         * Remove the degenerated polygons in all the Surfaces and all the
         * degenerate edges in all the Lines within
         * the GeoModel. Degeneration is due to colocated vertices.
         * Surfaces and Lines without any vertex anymore
         * (after the removal of the vertices) are removed off the GeoModel.
         */
        void remove_degenerate_polygons_and_edges_and_update_gm()
        {
            std::set< gmme_id > empty_mesh_entities;
            remove_degenerate_polygons_and_edges( empty_mesh_entities );
            if( !empty_mesh_entities.empty() )
            {
                builder_.remove.remove_mesh_entities( empty_mesh_entities );
            }

            // This is basic requirement ! no_colocated geomodel vertices !
            // So remove them if there are any
            geomodel_.mesh.remove_colocated_vertices();

            builder_.end_geomodel();
        }

        /*!
         * @brief For all the lines in the geomodel, switch line boundaries
         * if the way of their indices does not follow the way of the vertex
         * indices.
         */
        void repair_line_boundary_vertex_order()
        {
            for( const auto& line : geomodel_.lines() )
            {
                if( !line.is_first_corner_first_vertex() )
                {
                    const auto first_boundary_index =
                        line.boundary( 0 ).index();
                    builder_.topology.set_line_corner_boundary(
                        line.gmme().index(), 0,
                        line.boundary_gmme( 1 ).index() );
                    builder_.topology.set_line_corner_boundary(
                        line.gmme().index(), 1, first_boundary_index );
                }
            }
        }

        /*!
         * @brief remove isolated vertices on GeoModelMeshEntities
         */
        void remove_isolated_vertices();
        void remove_isolated_vertices_base()
        {
            for( const auto& line : geomodel_.lines() )
            {
                remove_isolated_vertices_on_mesh_entity( line );
            }
        }

        /*!
         * @brief remove isolated vertices on a GeoModelMeshEntity
         * @param[in] geomodel_mesh_entity The GeoModelMeshEntity to repair
         */
        void remove_isolated_vertices_on_mesh_entity(
            const GeoModelMeshEntity< DIMENSION >& geomodel_mesh_entity )
        {
            std::vector< bool > vertices_to_delete(
                geomodel_mesh_entity.nb_vertices(), true );
            for( auto mesh_element_index :
                range( geomodel_mesh_entity.nb_mesh_elements() ) )
            {
                for( auto vertex :
                    range( geomodel_mesh_entity.nb_mesh_element_vertices(
                        mesh_element_index ) ) )
                {
                    vertices_to_delete
                        [geomodel_mesh_entity.mesh_element_vertex_index(
                            { mesh_element_index, vertex } )] = false;
                }
            }
            builder_.geometry.delete_mesh_entity_vertices(
                geomodel_mesh_entity.gmme(), vertices_to_delete );
        }

        /*!
         * @brief Detect and remove degenerate edges in a \param line.
         * @return the number of degenerate edges that have been removed from
         * the line.
         */
        index_t repair_line_mesh( const Line< DIMENSION >& line )
        {
            std::vector< index_t > colocated;
            const auto& nn_search = line.vertex_nn_search();
            std::tie( std::ignore, colocated ) =
                nn_search.get_colocated_index_mapping( geomodel_.epsilon() );

            auto degenerate = line_detect_degenerate_edges( line, colocated );
            auto nb = static_cast< index_t >(
                std::count( degenerate.begin(), degenerate.end(), 1 ) );
            /// We have a problem if some vertices are left isolated
            /// If we remove them here we can kill all index correspondences
            builder_.geometry.delete_line_edges(
                line.index(), degenerate, false );
            return nb;
        }

        /*!
         * @return a vector of boolean. Element i of this vector corresponds
         * to the edge i of the line. If the element is true, the edge is
         * degenerated.
         */
        std::vector< bool > line_detect_degenerate_edges(
            const Line< DIMENSION >& line,
            std::vector< index_t >& colocated_vertices )
        {
            std::vector< bool > e_is_degenerate( line.nb_mesh_elements() );
            for( auto e : range( line.nb_mesh_elements() ) )
            {
                e_is_degenerate[e] =
                    edge_is_degenerate( line, e, colocated_vertices );
            }
            return e_is_degenerate;
        }

        /*!
         * \note Copied and modified from geogram\mesh\mesh_repair.cpp.
         * @return a vector of boolean. Element i of this vector corresponds
         * to the facet i of the surface. If the element is true, the facet is
         * degenerated.
         */
        std::vector< index_t > surface_detect_degenerate_polygons(
            const Surface< DIMENSION >& surface,
            std::vector< index_t >& colocated_vertices )
        {
            std::vector< index_t > f_is_degenerate(
                surface.nb_mesh_elements() );
            for( auto p : range( surface.nb_mesh_elements() ) )
            {
                f_is_degenerate[p] =
                    polygon_is_degenerate( surface, p, colocated_vertices );
            }
            return f_is_degenerate;
        }

        /*!
         * \note Copied and modified from geogram\mesh\mesh_repair.cpp
         *
         * @brief Tests whether a polygon is degenerate.
         * @param[in] surface the Surface that the polygon belongs to
         * @param[in] polygon_id the index of the polygon in \p S
         * @param[out] colocated_vertices contains the found colocated vertices
         * in \p f if any.
         * \return true if polygon \p f has duplicated vertices,
         *  false otherwise
         */
        bool polygon_is_degenerate( const Surface< DIMENSION >& surface,
            index_t polygon_id,
            std::vector< index_t >& colocated_vertices )
        {
            auto nb_vertices = surface.nb_mesh_element_vertices( polygon_id );
            if( nb_vertices != 3 )
            {
                std::vector< index_t > vertices( nb_vertices );
                for( auto v : range( nb_vertices ) )
                {
                    vertices[v] =
                        colocated_vertices[surface.mesh_element_vertex_index(
                            ElementLocalVertex( polygon_id, v ) )];
                }
                GEO::sort_unique( vertices );
                return vertices.size() != nb_vertices;
            }
            auto v1 = colocated_vertices[surface.mesh_element_vertex_index(
                { polygon_id, 0 } )];
            auto v2 = colocated_vertices[surface.mesh_element_vertex_index(
                { polygon_id, 1 } )];
            auto v3 = colocated_vertices[surface.mesh_element_vertex_index(
                { polygon_id, 2 } )];
            return v1 == v2 || v2 == v3 || v3 == v1;
        }

        /*!
         * @brief Detect and remove degenerated polygons in a Surface
         * @param[in,out] surface Surface to check for potential degenerate
         * polygons.
         * @return the number of degenerate polygons in \p surface.
         */
        index_t detect_degenerate_polygons(
            const Surface< DIMENSION >& surface )
        {
            std::vector< index_t > colocated;
            const auto& nn_search = surface.vertex_nn_search();
            std::tie( std::ignore, colocated ) =
                nn_search.get_colocated_index_mapping( geomodel_.epsilon() );

            auto degenerate =
                surface_detect_degenerate_polygons( surface, colocated );
            return static_cast< index_t >(
                std::count( degenerate.begin(), degenerate.end(), 1 ) );
        }

        /*!
         * @brief Remove degenerate polygons and edges from the Surface
         *        and Line of the geomodel.
         * @param[out] to_remove gmme_t of the entities (Surface and Line)
         * of the geomodel that are empty once degenerate entities are removed
         * @pre Colocated vertices have already been removed
         */
        void remove_degenerate_polygons_and_edges(
            std::set< gmme_id >& to_remove )
        {
            to_remove.clear();
            for( const auto& line : geomodel_.lines() )
            {
                auto nb = repair_line_mesh( line );
                if( nb > 0 )
                {
                    Logger::out( "Repair", nb, " degenerated edges removed in ",
                        line.gmme() );
                    // If the Line is set it to remove
                    if( line.nb_mesh_elements() == 0 )
                    {
                        to_remove.insert( line.gmme() );
                    }
                }
            }
            double epsilon_sq = geomodel_.epsilon() * geomodel_.epsilon();
            for( const auto& surface : geomodel_.surfaces() )
            {
                auto nb = detect_degenerate_polygons( surface );
                /// @todo Check if that cannot be simplified
                if( nb > 0 )
                {
                    if( surface.nb_vertices() > 0 )
                    {
                        auto builder = builder_.geometry.create_surface_builder(
                            surface.index() );
                        remove_duplicated_or_degenerated_polygons(
                            surface.mesh(), *builder );
                        remove_small_connected_components(
                            surface.mesh(), *builder, epsilon_sq, 3 );
                    }
                    if( surface.nb_vertices() == 0
                        || surface.nb_mesh_elements() == 0 )
                    {
                        to_remove.insert( surface.gmme() );
                    }
                }
            }
        }

        bool polygon_is_degenerate(
            const SurfaceMesh< DIMENSION >& surface, index_t polygon_id )
        {
            if( surface.polygon_area( polygon_id ) < geomodel_.epsilon2() )
            {
                return true;
            }

            auto min_length = geomodel_.epsilon();
            for( auto c : range( surface.nb_polygon_vertices( polygon_id ) ) )
            {
                if( surface.polygon_edge_length( { polygon_id, c } )
                    < min_length )
                {
                    return false;
                }
            }
            return false;
        }

        void detect_bad_facets( const SurfaceMesh< DIMENSION >& surface,
            std::vector< bool >& remove_polygon )
        {
            const auto& polygon_search = surface.polygon_nn_search();
            index_t nb_duplicates;
            std::vector< index_t > mapping;
            std::tie( nb_duplicates, mapping ) =
                polygon_search.get_colocated_index_mapping(
                    geomodel_.epsilon() );
            for( auto p : range( surface.nb_polygons() ) )
            {
                if( mapping[p] != p )
                {
                    remove_polygon[p] = true;
                    // Check if duplicated polygons are adjacent
                    for( auto v : range( surface.nb_polygon_vertices( p ) ) )
                    {
                        if( surface.polygon_adjacent( { p, v } ) == mapping[p] )
                        {
                            // If the duplicated polygons are adjacent, the
                            // shared
                            // edges will become a non-manifold edge.
                            // The two polygons should be removed.
                            remove_polygon[mapping[p]] = true;
                            break;
                        }
                    }
                }
            }

            index_t nb_degenerate = 0;
            for( auto p : range( surface.nb_polygons() ) )
            {
                if( !remove_polygon[p] && polygon_is_degenerate( surface, p ) )
                {
                    nb_degenerate++;
                    remove_polygon[p] = true;
                }
            }
            if( nb_duplicates != 0 || nb_degenerate != 0 )
            {
                Logger::out( "Repair", "Detected ", nb_duplicates,
                    " duplicate and ", nb_degenerate, " degenerate facets." );
            }
        }

        void remove_duplicated_or_degenerated_polygons(
            const SurfaceMesh< DIMENSION >& surface,
            SurfaceMeshBuilder< DIMENSION >& builder )
        {
            std::vector< bool > remove_polygon( surface.nb_polygons(), false );
            detect_bad_facets( surface, remove_polygon );
            builder.delete_polygons( remove_polygon, false );
            for( auto p : range( surface.nb_polygons() ) )
            {
                for( auto v : range( surface.nb_polygon_vertices( p ) ) )
                {
                    builder.set_polygon_adjacent( { p, v }, NO_ID );
                }
            }
            builder.connect_polygons();
        }

        /*!
         * \brief Removes the connected components that have an area
         *  smaller than a given threshold.
         * \param[in] min_area the connected components with an
         *  area smaller than this threshold are removed
         * \param[in] min_polygons the connected components with
         *  less than \param min_polygons polygons are removed
         */
        void remove_small_connected_components(
            const SurfaceMesh< DIMENSION >& surface,
            SurfaceMeshBuilder< DIMENSION >& builder,
            double min_area,
            index_t min_polygons )
        {
            std::vector< index_t > components;
            index_t nb_components;
            std::tie( nb_components, components ) =
                surface.connected_components();
            if( nb_components == 0 )
            {
                return;
            }
            std::vector< double > comp_area( nb_components, 0.0 );
            std::vector< index_t > comp_polygons( nb_components, 0 );
            for( auto p : range( surface.nb_polygons() ) )
            {
                comp_area[components[p]] += surface.polygon_area( p );
                ++comp_polygons[components[p]];
            }

            std::vector< bool > polygon_to_delete(
                surface.nb_polygons(), false );
            for( auto p : range( surface.nb_polygons() ) )
            {
                auto component = components[p];
                if( comp_area[component] < min_area
                    || comp_polygons[component] < min_polygons )
                {
                    polygon_to_delete[p] = true;
                }
            }
            builder.delete_polygons( polygon_to_delete, true );
        }

        /*!
         * @brief Remove colocated vertices of the geomodel.
         * @param[out] to_remove gmme_t of the entities of the geomodel that
         *  are empty once degenerate entities are removed
         */
        void remove_colocated_entity_vertices( std::set< gmme_id >& to_remove )
        {
            to_remove.clear();
            // For all Lines and Surfaces
            std::array< const MeshEntityType, 2 > types{
                { Line< DIMENSION >::type_name_static(),
                    Surface< DIMENSION >::type_name_static() }
            };
            for( const auto& type : types )
            {
                for( auto e : range( geomodel_.nb_mesh_entities( type ) ) )
                {
                    gmme_id entity_id{ type, e };
                    const auto& E = geomodel_.mesh_entity( entity_id );

                    const auto& kdtree = E.vertex_nn_search();
                    std::vector< index_t > colocated;
                    std::tie( std::ignore, colocated ) =
                        kdtree.get_colocated_index_mapping(
                            geomodel_.epsilon() );

                    // Get the vertices to delete
                    auto inside_border =
                        vertices_on_inside_boundary( entity_id );

                    std::vector< bool > to_delete( colocated.size(), false );
                    index_t nb_todelete{ 0 };
                    for( auto v : range( colocated.size() ) )
                    {
                        if( colocated[v] == v
                            || inside_border.find( v ) != inside_border.end() )
                        {
                            // This point is kept
                            // No colocated or on an inside boundary
                        }
                        else
                        {
                            // The point is to remove
                            to_delete[v] = true;
                            nb_todelete++;
                        }
                    }

                    if( nb_todelete == 0 )
                    {
                        // Nothing to do there
                        continue;
                    }
                    if( nb_todelete == E.nb_vertices() )
                    {
                        // The complete entity should be removed
                        to_remove.insert( E.gmme() );
                        continue;
                    }
                    if( type == Surface< DIMENSION >::type_name_static() )
                    {
                        auto builder =
                            builder_.geometry.create_surface_builder( e );
                        for( auto p_itr : range( E.nb_mesh_elements() ) )
                        {
                            for( auto fpv_itr :
                                range( E.nb_mesh_element_vertices( p_itr ) ) )
                            {
                                builder->set_polygon_vertex( { p_itr, fpv_itr },
                                    colocated[E.mesh_element_vertex_index(
                                        { p_itr, fpv_itr } )] );
                            }
                        }
                        builder->delete_vertices( to_delete );
                        Logger::out( "Repair", nb_todelete,
                            " colocated vertices deleted in ", entity_id );
                    }
                    else if( type == Line< DIMENSION >::type_name_static() )
                    {
                        auto builder =
                            builder_.geometry.create_line_builder( e );
                        for( auto e_itr : range( E.nb_mesh_elements() ) )
                        {
                            builder->set_edge_vertex( { e_itr, 0 },
                                colocated[E.mesh_element_vertex_index(
                                    { e_itr, 0 } )] );
                            builder->set_edge_vertex( { e_itr, 1 },
                                colocated[E.mesh_element_vertex_index(
                                    { e_itr, 1 } )] );
                        }
                        builder->delete_vertices( to_delete );
                        Logger::out( "Repair", nb_todelete,
                            " colocated vertices deleted in ", entity_id );
                    }
                    else
                    {
                        ringmesh_assert_not_reached;
                    }
                }
            }
        }

        /*!
         * Get the indices of the duplicated vertices that are on an inside
         * border.
         * Only the vertex with the biggest index are added.
         * @param[in] E_id GeoModelMeshEntity to check.
         * @return vector of the vertex indexes on an inside boundary.
         */
        std::set< index_t > vertices_on_inside_boundary( const gmme_id& E_id )
        {
            std::set< index_t > vertices;
            if( E_id.type() == Corner< DIMENSION >::type_name_static() )
            {
                return vertices;
            }
            const auto& mesh_entity = geomodel_.mesh_entity( E_id );
            if( E_id.type() == Line< DIMENSION >::type_name_static() )
            {
                if( mesh_entity.boundary( 0 ).is_inside_border( mesh_entity ) )
                {
                    vertices.insert( mesh_entity.nb_vertices() - 1 );
                }
                return vertices;
            }
            std::vector< const GeoModelMeshEntity< DIMENSION >* > inside_border;
            for( auto i : range( mesh_entity.nb_boundaries() ) )
            {
                if( mesh_entity.boundary( i ).is_inside_border( mesh_entity ) )
                {
                    inside_border.push_back(
                        dynamic_cast< const GeoModelMeshEntity< DIMENSION >* >(
                            &mesh_entity.boundary( i ) ) );
                }
            }
            if( !inside_border.empty() )
            {
                // We want to get the indices of the vertices in E
                // that are colocated with those of the inside boundary
                // We assume that the geomodel vertices are not computed
                const auto& nn_search = mesh_entity.vertex_nn_search();

                for( const auto& entity : inside_border )
                {
                    for( auto v : range( entity->nb_vertices() ) )
                    {
                        auto colocated_indices = nn_search.get_neighbors(
                            entity->vertex( v ), geomodel_.epsilon() );
                        if( colocated_indices.size() > 1 )
                        {
                            std::sort( colocated_indices.begin(),
                                colocated_indices.end() );
                            // Add colocated vertices except one to the
                            // duplicated
                            // vertices set
                            vertices.insert( colocated_indices.begin() + 1,
                                colocated_indices.end() );
                        }
                    }
                }
            }
            return vertices;
        }

        /*!
         * @brief Checks if an edge is degenerate.
         *
         * An edge is degenerate if both vertices are colocated.
         *
         * @param[in] line Line to check the edge \p edge.
         * @param[in] edge edge index in Line \p line.
         * @param[in] colocated_vertices contains the colocated mapping of the
         * Line.
         * @return true if the edge is degenerate. Else false.
         */
        bool edge_is_degenerate( const Line< DIMENSION >& line,
            index_t edge,
            const std::vector< index_t >& colocated_vertices )
        {
            auto v1 = colocated_vertices[line.mesh_element_vertex_index(
                { edge, 0 } )];
            auto v2 = colocated_vertices[line.mesh_element_vertex_index(
                { edge, 1 } )];
            return v1 == v2;
        }

        void build_contacts()
        {
            builder_.geology.build_contacts();
        }

    private:
        GeoModelBuilder< DIMENSION > builder_;
        GeoModel< DIMENSION >& geomodel_;
    };

    template <>
    void GeoModelRepair< 3 >::remove_isolated_vertices()
    {
        remove_isolated_vertices_base();
        for( const auto& surface : geomodel_.surfaces() )
        {
            remove_isolated_vertices_on_mesh_entity( surface );
        }
        for( const auto& region : geomodel_.regions() )
        {
            if( region.is_meshed() )
            {
                remove_isolated_vertices_on_mesh_entity( region );
            }
        }
    }

    template <>
    void GeoModelRepair< 2 >::remove_isolated_vertices()
    {
        remove_isolated_vertices_base();
        for( const auto& surface : geomodel_.surfaces() )
        {
            if( surface.is_meshed() )
            {
                remove_isolated_vertices_on_mesh_entity( surface );
            }
        }
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    void repair_geomodel(
        GeoModel< DIMENSION >& geomodel, RepairMode repair_mode )
    {
        GeoModelRepair< DIMENSION > repairer( geomodel );
        repairer.repair( repair_mode );
    }

    template void geomodel_tools_api repair_geomodel( GeoModel2D&, RepairMode );

    template void geomodel_tools_api repair_geomodel( GeoModel3D&, RepairMode );
} // namespace RINGMesh
