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

#include <ringmesh/geomodel/geomodel_builder_geometry.h>

#include <ringmesh/geomodel/geomodel_builder.h>

namespace {
    using namespace RINGMesh;

    void get_internal_borders(
        const GeoModelMeshEntity& entity,
        std::set< index_t >& internal_borders )
    {
        for( index_t i = 0; i < entity.nb_boundaries(); ++i ) {
            const GeoModelMeshEntity& border = entity.boundary( i );
            if( border.is_inside_border( entity ) ) {
                internal_borders.insert( border.index() );
            }
        }
    }

    bool inexact_equal( const vec3& v1, const vec3& v2, double epsilon )
    {
        return length( v2 - v1 ) < epsilon;
    }

    /*!
     * Finds a polygon and its edge index that are colocalised with an edge
     * defined by its two geomodel vertex indices
     * @param[in] surface the surface where to find the polygon
     * @param[in] v0 the first vertex of the edge
     * @param[in] v1 the second vertex of the edge
     * @param[out] f the found polygon index
     * @param[out] e the found edge index
     * @return True if the polygon and the edge indices are found
     */
    bool find_polygon_from_edge_vertices(
        const Surface& surface,
        const vec3& v0,
        const vec3& v1,
        index_t& p,
        index_t& e )
    {
        vec3 v_bary = 0.5 * ( v0 + v1 );
        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() );
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() );
            std::vector< index_t > neighbors =
                surface.polygon_nn_search().get_neighbors( v_bary, cur_neighbor );
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                p = neighbors[i];
                for( index_t j = 0; j < surface.nb_mesh_element_vertices( p );
                    j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( p, j ), v0,
                        surface.geomodel().epsilon() ) ) {
                        index_t j_next = surface.next_polygon_vertex_index( p, j );
                        if( inexact_equal( surface.mesh_element_vertex( p, j_next ),
                            v1, surface.geomodel().epsilon() ) ) {
                            e = j;
                            return true;
                        }
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor );

        p = NO_ID;
        e = NO_ID;
        return false;
    }

    bool are_cell_facet_and_polygon_equal(
        const Region& region,
        index_t cell,
        index_t cell_facet,
        const Surface& surface,
        index_t polygon )
    {
        index_t nb_cell_facet_vertices = region.nb_cell_facet_vertices( cell,
            cell_facet );
        index_t nb_polygon_vertices = surface.nb_mesh_element_vertices( polygon );
        if( nb_cell_facet_vertices != nb_polygon_vertices ) {
            return false;
        }
        vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell,
            cell_facet );
        vec3 polygon_barycenter = surface.mesh_element_barycenter( polygon );
        return inexact_equal( cell_facet_barycenter, polygon_barycenter,
            region.geomodel().epsilon() );
    }

    bool find_cell_facet_from_polygon(
        const Region& region,
        const Surface& surface,
        index_t polygon,
        index_t& cell,
        index_t& cell_facet )
    {
        vec3 v_bary = surface.mesh_element_barycenter( polygon );
        index_t nb_neighbors = std::min( index_t( 5 ), region.nb_mesh_elements() );
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, region.nb_mesh_elements() );
            std::vector< index_t > neighbors = region.cell_nn_search().get_neighbors(
                v_bary, cur_neighbor );
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                cell = neighbors[i];
                for( cell_facet = 0; cell_facet < region.nb_cell_facets( cell );
                    cell_facet++ ) {
                    if( are_cell_facet_and_polygon_equal( region, cell, cell_facet,
                        surface, polygon ) ) {
                        return true;
                    }
                }
            }
        } while( region.nb_mesh_elements() != cur_neighbor );

        cell = NO_ID;
        cell_facet = NO_ID;
        return false;
    }

    bool find_polygon_from_vertex(
        const Surface& surface,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() );
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() );
            std::vector< index_t > neighbors =
                surface.polygon_nn_search().get_neighbors( v, cur_neighbor );
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i];
                for( index_t j = 0;
                    j < surface.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( element_id, j ),
                        v, surface.geomodel().epsilon() ) ) {
                        vertex_id = surface.mesh_element_vertex_index( element_id,
                            j );
                        return true;
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor );

        element_id = NO_ID;
        vertex_id = NO_ID;
        return false;
    }

    bool find_cell_from_vertex(
        const Region& entity,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        index_t nb_neighbors = std::min( index_t( 5 ), entity.nb_mesh_elements() );
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, entity.nb_mesh_elements() );
            std::vector< index_t > neighbors = entity.cell_nn_search().get_neighbors(
                v, cur_neighbor );
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i];
                for( index_t j = 0;
                    j < entity.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( entity.mesh_element_vertex( element_id, j ),
                        v, entity.geomodel().epsilon() ) ) {
                        vertex_id = entity.mesh_element_vertex_index( element_id,
                            j );
                        return true;
                    }
                }
            }
        } while( entity.nb_mesh_elements() != cur_neighbor );

        element_id = NO_ID;
        vertex_id = NO_ID;
        return false;
    }

    index_t edge_index_from_polygon_and_edge_vertex_indices(
        const Surface& surface,
        index_t p,
        const vec3& v0,
        const vec3& v1 )
    {
        for( index_t v = 0; v < surface.nb_mesh_element_vertices( p ); v++ ) {
            if( !inexact_equal( surface.mesh_element_vertex( p, v ), v0,
                surface.geomodel().epsilon() ) ) continue;
            index_t prev_v = surface.prev_polygon_vertex_index( p, v );
            index_t next_v = surface.next_polygon_vertex_index( p, v );
            if( inexact_equal( surface.mesh_element_vertex( p, prev_v ), v1,
                surface.geomodel().epsilon() ) ) {
                return prev_v;
            } else if( inexact_equal( surface.mesh_element_vertex( p, next_v ), v1,
                surface.geomodel().epsilon() ) ) {
                return v;
            }
        }
        return NO_ID;
    }

    index_t cell_facet_index_from_cell_and_polygon(
        const Region& region,
        index_t cell,
        const Surface& surface,
        index_t polygon )
    {
        vec3 polygon_barycenter = surface.mesh_element_barycenter( polygon );
        for( index_t f = 0; f < region.nb_cell_facets( cell ); f++ ) {
            vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell, f );
            if( inexact_equal( cell_facet_barycenter, polygon_barycenter,
                surface.geomodel().epsilon() ) ) {
                return f;
            }
        }
        return NO_ID;
    }

    void check_and_initialize_corner_vertex( GeoModel& geomodel, index_t corner_id )
    {
        if( geomodel.corner( corner_id ).nb_vertices() == 0 ) {
            GeoModelBuilder builder( geomodel );
            builder.geometry.create_mesh_entity_vertices(
                gmme_id( Corner::type_name_static(), corner_id ), 1 );
        }
    }
}

namespace RINGMesh {
    GeoModelBuilderGeometry::GeoModelBuilderGeometry(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderGeometry::recompute_geomodel_mesh()
    {
        geomodel_.mesh.vertices.clear();
        geomodel_.mesh.vertices.test_and_initialize();
    }

    void GeoModelBuilderGeometry::copy_meshes( const GeoModel& geomodel )
    {
        copy_meshes( geomodel, Corner::type_name_static() );
        copy_meshes( geomodel, Line::type_name_static() );
        copy_meshes( geomodel, Surface::type_name_static() );
        copy_meshes( geomodel, Region::type_name_static() );
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertex(
        const gmme_id& t,
        index_t v,
        const vec3& point,
        bool update )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( t );
        GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices;
        ringmesh_assert( v < E.nb_vertices() );
        if( update ) {
            geomodel_vertices.update_point(
                geomodel_vertices.geomodel_vertex_id( E.gmme(), v ), point );
        } else {
            GeoModelMeshEntityAccess gmme_access( E );
            std::unique_ptr< MeshBaseBuilder > builder =
                MeshBaseBuilder::create_builder( *gmme_access.modifiable_mesh() );
            builder->set_vertex( v, point );
        }
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertex(
        const gmme_id& entity_id,
        index_t v,
        index_t geomodel_vertex )
    {
        GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices;
        set_mesh_entity_vertex( entity_id, v,
            geomodel_vertices.vertex( geomodel_vertex ), false );

        ringmesh_assert( v < geomodel_.mesh_entity( entity_id ).nb_vertices() );
        geomodel_vertices.update_vertex_mapping( entity_id, v, geomodel_vertex );
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertices(
        const gmme_id& id,
        const std::vector< vec3 >& points,
        bool clear )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( id );
        GeoModelMeshEntityAccess gmme_access( E );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true );
        }
        if( !points.empty() ) {
            index_t nb_points = static_cast< index_t >( points.size() );
            index_t start = builder->create_vertices( nb_points );
            for( index_t v = 0; v < nb_points; v++ ) {
                builder->set_vertex( start + v, points[v] );
            }
        }
    }

    index_t GeoModelBuilderGeometry::create_mesh_entity_vertices(
        const gmme_id& entity_id,
        index_t nb_vertices )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        GeoModelMeshEntityAccess gmme_access( E );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        return builder->create_vertices( nb_vertices );
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertices(
        const gmme_id& entity_id,
        const std::vector< index_t >& geomodel_vertices,
        bool clear )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        GeoModelMeshEntityAccess gmme_access( E );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true );
        }
        index_t nb_model_vertices =
            static_cast< index_t >( geomodel_vertices.size() );
        index_t start = builder->create_vertices( nb_model_vertices );
        for( index_t v = 0; v < nb_model_vertices; v++ ) {
            set_mesh_entity_vertex( entity_id, start + v, geomodel_vertices[v] );
        }
    }

    void GeoModelBuilderGeometry::set_corner( index_t corner_id, const vec3& point )
    {
        check_and_initialize_corner_vertex( geomodel_, corner_id );
        set_mesh_entity_vertex( gmme_id( Corner::type_name_static(), corner_id ), 0,
            point, false );
    }

    void GeoModelBuilderGeometry::set_line(
        index_t line_id,
        const std::vector< vec3 >& vertices )
    {
        set_mesh_entity_vertices( gmme_id( Line::type_name_static(), line_id ),
            vertices, true );

        Line& line = dynamic_cast< Line& >( geomodel_access_.modifiable_mesh_entity(
            gmme_id( Line::type_name_static(), line_id ) ) );
        std::unique_ptr< LineMeshBuilder > builder = create_line_builder( line_id );
        for( index_t e = 1; e < line.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e );
        }
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< vec3 >& surface_vertices,
        const std::vector< index_t >& surface_polygons,
        const std::vector< index_t >& surface_polygon_ptr )
    {
        set_mesh_entity_vertices( gmme_id( Surface::type_name_static(), surface_id ),
            surface_vertices, true );
        set_surface_geometry( surface_id, surface_polygons, surface_polygon_ptr );
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& polygons,
        const std::vector< index_t >& polygon_ptr )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        builder->create_polygons( polygons, polygon_ptr );
        compute_surface_adjacencies( surface_id );
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        builder->assign_triangle_mesh( triangle_vertices );
        compute_surface_adjacencies( surface_id );
    }

    void GeoModelBuilderGeometry::set_region_geometry(
        index_t region_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& tetras )
    {
        set_mesh_entity_vertices( gmme_id( Region::type_name_static(), region_id ),
            points, true );
        assign_region_tet_mesh( region_id, tetras );
    }

    void GeoModelBuilderGeometry::set_corner(
        index_t corner_id,
        index_t geomodel_vertex_id )
    {
        check_and_initialize_corner_vertex( geomodel_, corner_id );
        set_mesh_entity_vertex( gmme_id( Corner::type_name_static(), corner_id ), 0,
            geomodel_vertex_id );
    }

    void GeoModelBuilderGeometry::set_line(
        index_t line_id,
        const std::vector< index_t >& unique_vertices )
    {
        bool clear_vertices = false;
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gmme_id( Line::type_name_static(), line_id ) );

        ringmesh_assert( E.nb_vertices() == 0 );
        // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_mesh_entity_vertices( E.gmme(), unique_vertices, clear_vertices );

        std::unique_ptr< LineMeshBuilder > builder = create_line_builder( line_id );
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e );
        }
    }

    void GeoModelBuilderGeometry::set_surface_element_geometry(
        index_t surface_id,
        index_t polygon_id,
        const std::vector< index_t >& corners )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        for( index_t polygon_vertex = 0; polygon_vertex < corners.size();
            polygon_vertex++ ) {
            builder->set_polygon_vertex( polygon_id, polygon_vertex,
                corners[polygon_vertex] );
        }
    }

    void GeoModelBuilderGeometry::set_region_element_geometry(
        index_t region_id,
        index_t cell_id,
        const std::vector< index_t >& corners )
    {
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        for( index_t cell_vertex = 0; cell_vertex < corners.size(); cell_vertex++ ) {
            builder->set_cell_vertex( cell_id, cell_vertex, corners[cell_vertex] );
        }
    }

    index_t GeoModelBuilderGeometry::create_surface_polygon(
        index_t surface_id,
        const std::vector< index_t >& vertex_indices )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        return builder->create_polygon( vertex_indices );
    }

    index_t GeoModelBuilderGeometry::create_region_cells(
        index_t region_id,
        GEO::MeshCellType type,
        index_t nb_cells )
    {
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        return builder->create_cells( nb_cells, type );
    }

    index_t GeoModelBuilderGeometry::create_region_cell(
        index_t region_id,
        GEO::MeshCellType type,
        const std::vector< index_t >& vertex_indices )
    {
        index_t cell_id = create_region_cells( region_id, type, 1 );
        set_region_element_geometry( region_id, cell_id, vertex_indices );
        return cell_id;
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_mesh( const gmme_id& E_id )
    {
        GeoModelMeshEntityAccess gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->clear( true, false );
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_isolated_vertices(
        const gmme_id& E_id )
    {
        if( geomodel_.entity_type_manager().mesh_entity_manager.is_line(
            E_id.type() ) ) {
            std::unique_ptr< LineMeshBuilder > builder = create_line_builder(
                E_id.index() );
            builder->remove_isolated_vertices();
        } else if( geomodel_.entity_type_manager().mesh_entity_manager.is_surface(
            E_id.type() ) ) {
            std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
                E_id.index() );
            builder->remove_isolated_vertices();
        } else if( geomodel_.entity_type_manager().mesh_entity_manager.is_region(
            E_id.type() ) ) {
            std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
                E_id.index() );
            builder->remove_isolated_vertices();
        } else if( geomodel_.entity_type_manager().mesh_entity_manager.is_corner(
            E_id.type() ) ) {
            std::unique_ptr< PointSetMeshBuilder > builder = create_corner_builder(
                E_id.index() );
            builder->remove_isolated_vertices();
        } else {
            ringmesh_assert_not_reached;
        }
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_vertices(
        const gmme_id& E_id,
        const std::vector< bool >& to_delete )
    {
        GeoModelMeshEntityAccess gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->delete_vertices( to_delete );
    }

    void GeoModelBuilderGeometry::delete_corner_vertex( index_t corner_id )
    {
        gmme_id corner( Corner::type_name_static(), corner_id );
        std::vector< bool > to_delete;
        to_delete.push_back( true );
        delete_mesh_entity_vertices( corner, to_delete );
    }
    void GeoModelBuilderGeometry::delete_line_edges(
        index_t line_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        std::unique_ptr< LineMeshBuilder > builder = create_line_builder( line_id );
        builder->delete_edges( to_delete, remove_isolated_vertices );
    }
    void GeoModelBuilderGeometry::delete_surface_polygons(
        index_t surface_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        builder->delete_polygons( to_delete, remove_isolated_vertices );
    }
    void GeoModelBuilderGeometry::delete_region_cells(
        index_t region_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        builder->delete_cells( to_delete, remove_isolated_vertices );
    }

    void GeoModelBuilderGeometry::set_surface_element_adjacency(
        index_t surface_id,
        index_t polygon_id,
        const std::vector< index_t >& adjacents )
    {
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        for( index_t polygon_edge = 0; polygon_edge < adjacents.size();
            polygon_edge++ ) {
            builder->set_polygon_adjacent( polygon_id, polygon_edge,
                adjacents[polygon_edge] );
        }
    }

    void GeoModelBuilderGeometry::compute_surface_adjacencies(
        index_t surface_id,
        bool recompute_adjacency )
    {
        const Surface& surface = geomodel_.surface( surface_id );
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );

        if( recompute_adjacency ) {
            for( index_t p = 0; p < surface.nb_mesh_elements(); p++ ) {
                for( index_t v = 0; v < surface.nb_mesh_element_vertices( p );
                    v++ ) {
                    builder->set_polygon_adjacent( p, v, NO_ID );
                }
            }
        }
        builder->connect_polygons();
    }

    void GeoModelBuilderGeometry::compute_region_adjacencies(
        index_t region_id,
        bool recompute_adjacency )
    {
        const Region& region = geomodel_.region( region_id );
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );

        if( recompute_adjacency ) {
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                for( index_t f = 0; f < region.nb_cell_facets( c ); f++ ) {
                    builder->set_cell_adjacent( c, f, NO_ID );
                }
            }
        }
        builder->connect_cells();
    }

    void GeoModelBuilderGeometry::cut_surfaces_by_internal_lines()
    {
        for( index_t s = 0; s < geomodel_.nb_surfaces(); s++ ) {
            Surface& surface =
                dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                    gmme_id( Surface::type_name_static(), s ) ) );
            std::set< index_t > cutting_lines;
            get_internal_borders( surface, cutting_lines );
            for( index_t line_id : cutting_lines ) {
                cut_surface_by_line( s, line_id );
            }
            if( !cutting_lines.empty() ) {
                std::unique_ptr< SurfaceMeshBuilder > surface_mesh_builder =
                    create_surface_builder( s );
                surface_mesh_builder->remove_isolated_vertices();
            }
        }
    }

    void GeoModelBuilderGeometry::cut_regions_by_internal_surfaces()
    {
        for( index_t r = 0; r < geomodel_.nb_regions(); r++ ) {
            Region& region =
                dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                    gmme_id( Region::type_name_static(), r ) ) );
            if( region.nb_mesh_elements() == 0 ) continue;
            std::set< index_t > cutting_surfaces;
            get_internal_borders( region, cutting_surfaces );
            for( index_t surface_id : cutting_surfaces ) {
                cut_region_by_surface( r, surface_id );
            }
            if( !cutting_surfaces.empty() ) {
                std::unique_ptr< VolumeMeshBuilder > region_mesh_builder =
                    create_region_builder( r );
                region_mesh_builder->remove_isolated_vertices();
            }
        }
    }

    void GeoModelBuilderGeometry::cut_surface_by_line(
        index_t surface_id,
        index_t line_id )
    {
        index_t nb_disconnected_edges = disconnect_surface_polygons_along_line_edges(
            surface_id, line_id );
        if( nb_disconnected_edges > 0 ) {
            duplicate_surface_vertices_along_line( surface_id, line_id );
        }
    }

    void GeoModelBuilderGeometry::cut_region_by_surface(
        index_t region_id,
        index_t surface_id )
    {
        index_t nb_disconnected_polygons =
            disconnect_region_cells_along_surface_polygons( region_id, surface_id );
        if( nb_disconnected_polygons > 0 ) {
            duplicate_region_vertices_along_surface( region_id, surface_id );
        }
    }

    struct ElementVertex {
        index_t element_;
        index_t vertex_;
    };

    void GeoModelBuilderGeometry::duplicate_surface_vertices_along_line(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );
        ringmesh_assert( line_id < geomodel_.nb_lines() );

        gmme_id surface_gme( Surface::type_name_static(), surface_id );
        const Surface& surface = geomodel_.surface( surface_id );
        const Line& line = geomodel_.line( line_id );

        std::vector< ElementVertex > polygon_vertices( line.nb_vertices() );
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v );

            index_t& polygon_vertex = polygon_vertices[v].vertex_;
            index_t& polygon = polygon_vertices[v].element_;
            bool found = find_polygon_from_vertex( surface, p, polygon,
                polygon_vertex );
            ringmesh_unused( found );
            ringmesh_assert( found && polygon != NO_ID && polygon_vertex != NO_ID );
        }

        index_t vertex_id = create_mesh_entity_vertices( surface_gme,
            line.nb_vertices() );
        std::unique_ptr< SurfaceMeshBuilder > surface_mesh_builder =
            create_surface_builder( surface_id );
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v );
            const index_t& polygon_vertex = polygon_vertices[v].vertex_;
            const index_t& polygon = polygon_vertices[v].element_;

            std::vector< index_t > polygons = surface.polygons_around_vertex(
                polygon_vertex, false, polygon );
            update_polygon_vertex( surface_id, polygons, polygon_vertex, vertex_id );
            surface_mesh_builder->set_vertex( vertex_id, p );
            vertex_id++;
        }
    }

    void GeoModelBuilderGeometry::duplicate_region_vertices_along_surface(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() );
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );

        gmme_id region_gme( Region::type_name_static(), region_id );
        const Region& region = geomodel_.region( region_id );
        const Surface& surface = geomodel_.surface( surface_id );

        std::vector< ElementVertex > cell_vertices( surface.nb_vertices() );
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v );

            index_t& cell = cell_vertices[v].element_;
            index_t& cell_vertex = cell_vertices[v].vertex_;
            bool found = find_cell_from_vertex( region, p, cell, cell_vertex );
            ringmesh_unused( found );
            ringmesh_assert( found && cell != NO_ID && cell_vertex != NO_ID );

        }

        index_t vertex_id = create_mesh_entity_vertices( region_gme,
            surface.nb_vertices() );
        std::unique_ptr< VolumeMeshBuilder > region_mesh_builder = create_region_builder(
            region_id );
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v );
            const index_t& cell = cell_vertices[v].element_;
            const index_t& cell_vertex = cell_vertices[v].vertex_;

            std::vector< index_t > cells = region.cells_around_vertex( cell_vertex,
                cell );
            update_cell_vertex( region_id, cells, cell_vertex, vertex_id );
            region_mesh_builder->set_vertex( vertex_id, p );
            vertex_id++;
        }
    }

    index_t GeoModelBuilderGeometry::disconnect_surface_polygons_along_line_edges(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );
        ringmesh_assert( line_id < geomodel_.nb_lines() );

        const Surface& surface = geomodel_.surface( surface_id );
        const Line& line = geomodel_.line( line_id );
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        index_t nb_disconnected_edges = 0;
        for( index_t i = 0; i < line.nb_mesh_elements(); ++i ) {
            const vec3& p0 = line.vertex( i );
            const vec3& p1 = line.vertex( i + 1 );

            index_t p = NO_ID;
            index_t e = NO_ID;
            bool found = find_polygon_from_edge_vertices( surface, p0, p1, p, e );
            ringmesh_unused( found );
            ringmesh_assert( found && p != NO_ID && e != NO_ID );

            index_t adj_f = surface.polygon_adjacent_index( p, e );
            if( adj_f != NO_ID ) {
                index_t adj_e = edge_index_from_polygon_and_edge_vertex_indices(
                    surface, adj_f, p0, p1 );
                ringmesh_assert( adj_e != NO_ID );
                builder->set_polygon_adjacent( p, e, NO_ID );
                builder->set_polygon_adjacent( adj_f, adj_e, NO_ID );
                nb_disconnected_edges++;
            }
        }
        return nb_disconnected_edges;
    }
    index_t GeoModelBuilderGeometry::disconnect_region_cells_along_surface_polygons(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() );
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );

        const Region& region = geomodel_.region( region_id );
        const Surface& surface = geomodel_.surface( surface_id );
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        index_t nb_disconnected_polygons = 0;
        for( index_t polygon = 0; polygon < surface.nb_mesh_elements(); ++polygon ) {
            index_t cell = NO_ID;
            index_t cell_facet = NO_ID;
            bool found = find_cell_facet_from_polygon( region, surface, polygon,
                cell, cell_facet );
            ringmesh_unused( found );
            ringmesh_assert( found && cell != NO_ID && cell_facet != NO_ID );

            index_t adj_cell = region.cell_adjacent_index( cell, cell_facet );
            if( adj_cell != NO_ID ) {
                index_t adj_cell_facet = cell_facet_index_from_cell_and_polygon(
                    region, adj_cell, surface, polygon );
                ringmesh_assert( adj_cell_facet != NO_ID );
                builder->set_cell_adjacent( cell, cell_facet, NO_ID );
                builder->set_cell_adjacent( adj_cell, adj_cell_facet, NO_ID );
                nb_disconnected_polygons++;
            }
        }
        return nb_disconnected_polygons;
    }

    void GeoModelBuilderGeometry::update_polygon_vertex(
        index_t surface_id,
        const std::vector< index_t >& polygons,
        index_t old_vertex,
        index_t new_vertex )
    {
        const Surface& surface = geomodel_.surface( surface_id );
        std::unique_ptr< SurfaceMeshBuilder > builder = create_surface_builder(
            surface_id );
        for( index_t cur_p : polygons ) {
            for( index_t cur_v = 0;
                cur_v < surface.nb_mesh_element_vertices( cur_p ); cur_v++ ) {
                if( surface.mesh_element_vertex_index( cur_p, cur_v )
                    == old_vertex ) {
                    builder->set_polygon_vertex( cur_p, cur_v, new_vertex );
                }
            }
        }
    }

    void GeoModelBuilderGeometry::update_cell_vertex(
        index_t region_id,
        const std::vector< index_t >& cells,
        index_t old_vertex,
        index_t new_vertex )
    {
        const Region& region = geomodel_.region( region_id );
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        for( index_t cur_c : cells ) {
            for( index_t cur_v = 0; cur_v < region.nb_mesh_element_vertices( cur_c );
                cur_v++ ) {
                if( region.mesh_element_vertex_index( cur_c, cur_v )
                    == old_vertex ) {
                    builder->set_cell_vertex( cur_c, cur_v, new_vertex );
                }
            }
        }
    }

    void GeoModelBuilderGeometry::copy_meshes(
        const GeoModel& from,
        const MeshEntityType& entity_type )
    {
        for( index_t i = 0; i < geomodel_.nb_mesh_entities( entity_type ); ++i ) {
            copy_mesh( from, gmme_id( entity_type, i ) );
        }
    }

    void GeoModelBuilderGeometry::copy_mesh(
        const GeoModel& from,
        const gmme_id& mesh_entity )
    {
        const GeoModelMeshEntityConstAccess from_E_const_access(
            from.mesh_entity( mesh_entity ) );
        assign_mesh_to_entity( *from_E_const_access.mesh(), mesh_entity );
    }

    void GeoModelBuilderGeometry::assign_mesh_to_entity(
        const MeshBase& mesh,
        const gmme_id& to )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( to );
        GeoModelMeshEntityAccess gmme_access( E );
        std::unique_ptr< MeshBaseBuilder > builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->copy( mesh, true );
    }

    void GeoModelBuilderGeometry::assign_region_tet_mesh(
        index_t region_id,
        const std::vector< index_t >& tet_vertices )
    {
        std::unique_ptr< VolumeMeshBuilder > builder = create_region_builder(
            region_id );
        builder->assign_cell_tet_mesh( tet_vertices );
        builder->connect_cells();
    }

}
