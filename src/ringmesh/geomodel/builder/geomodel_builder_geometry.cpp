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

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/builder/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>
namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    std::set< index_t > get_internal_borders(
        const GeoModelMeshEntity< DIMENSION >& entity )
    {
        std::set< index_t > internal_borders;
        for( auto i : range( entity.nb_boundaries() ) )
        {
            const auto& border = entity.boundary( i );
            if( border.is_inside_border( entity ) )
            {
                internal_borders.insert( border.index() );
            }
        }
        return internal_borders;
    }

    /*!
     * Finds a polygon and its edge index that are colocalised with an edge
     * defined by its two geomodel vertex indices
     * @param[in] surface the surface where to find the polygon
     * @param[in] v0 the first vertex of the edge
     * @param[in] v1 the second vertex of the edge
     * @param[out] polygon the found polygon index
     * @param[out] edge the found edge index
     * @return True if the polygon and the edge indices are found
     */
    template < index_t DIMENSION >
    std::tuple< bool, index_t, index_t > find_polygon_from_edge_vertices(
        const Surface< DIMENSION >& surface,
        const vecn< DIMENSION >& v0,
        const vecn< DIMENSION >& v1 )
    {
        index_t polygon = NO_ID;
        index_t edge = NO_ID;
        auto v_bary = 0.5 * ( v0 + v1 );
        bool result{ false };
        surface.polygon_nn_search().get_neighbors( v_bary,
            [&surface, &v0, &v1, &result, &edge, &polygon]( index_t i ) {
                for( auto j : range( surface.nb_mesh_element_vertices( i ) ) )
                {
                    if( inexact_equal( surface.mesh_element_vertex( { i, j } ),
                            v0, surface.geomodel().epsilon() ) )
                    {
                        index_t j_next = surface.mesh()
                                             .next_polygon_vertex( { i, j } )
                                             .local_vertex_id;
                        if( inexact_equal(
                                surface.mesh_element_vertex( { i, j_next } ),
                                v1, surface.geomodel().epsilon() ) )
                        {
                            edge = j;
                            polygon = i;
                            result = true;
                            break;
                        }
                    }
                }
                return result;
            } );
        return std::make_tuple( result, polygon, edge );
    }

    template < index_t DIMENSION >
    bool are_cell_facet_and_polygon_equal( const Region< DIMENSION >& region,
        index_t cell,
        index_t cell_facet,
        const Surface< DIMENSION >& surface,
        index_t polygon )
    {
        auto nb_cell_facet_vertices =
            region.nb_cell_facet_vertices( cell, cell_facet );
        auto nb_polygon_vertices = surface.nb_mesh_element_vertices( polygon );
        if( nb_cell_facet_vertices != nb_polygon_vertices )
        {
            return false;
        }
        auto cell_facet_barycenter =
            region.mesh().cell_facet_barycenter( { cell, cell_facet } );
        auto polygon_barycenter = surface.mesh_element_barycenter( polygon );
        return inexact_equal( cell_facet_barycenter, polygon_barycenter,
            region.geomodel().epsilon() );
    }

    template < index_t DIMENSION >
    std::tuple< bool, index_t, index_t > find_cell_facet_from_polygon(
        const Region< DIMENSION >& region,
        const Surface< DIMENSION >& surface,
        index_t polygon )
    {
        index_t cell{ NO_ID };
        index_t cell_facet{ NO_ID };
        auto v_bary = surface.mesh_element_barycenter( polygon );
        bool result{ false };
        region.cell_nn_search().get_neighbors(
            v_bary, [&region, &surface, polygon, &result, &cell_facet, &cell](
                        index_t i ) {
                for( auto cell_facet_i : range( region.nb_cell_facets( i ) ) )
                {
                    if( are_cell_facet_and_polygon_equal(
                            region, i, cell_facet_i, surface, polygon ) )
                    {
                        cell_facet = cell_facet_i;
                        cell = i;
                        result = true;
                        break;
                    }
                }
                return result;
            } );
        return std::make_tuple( result, cell, cell_facet );
    }

    template < index_t DIMENSION >
    bool find_polygon_from_vertex( const Surface< DIMENSION >& surface,
        const vecn< DIMENSION >& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        bool result{ false };
        surface.polygon_nn_search().get_neighbors( v, [&surface, &v, &result,
                                                          &vertex_id,
                                                          &element_id](
                                                          index_t i ) {
            for( auto j : range( surface.nb_mesh_element_vertices( i ) ) )
            {
                if( inexact_equal( surface.mesh_element_vertex( { i, j } ), v,
                        surface.geomodel().epsilon() ) )
                {
                    vertex_id = surface.mesh_element_vertex_index( { i, j } );
                    element_id = i;
                    result = true;
                    break;
                }
            }
            return result;
        } );
        return result;
    }

    template < index_t DIMENSION >
    index_t edge_index_from_polygon_and_edge_vertex_indices(
        const Surface< DIMENSION >& surface,
        index_t p,
        const vecn< DIMENSION >& v0,
        const vecn< DIMENSION >& v1 )
    {
        const auto& mesh = surface.mesh();
        for( auto v : range( surface.nb_mesh_element_vertices( p ) ) )
        {
            if( !inexact_equal( surface.mesh_element_vertex( { p, v } ), v0,
                    surface.geomodel().epsilon() ) )
            {
                continue;
            }
            auto prev_v = mesh.prev_polygon_vertex( { p, v } ).local_vertex_id;
            auto next_v = mesh.next_polygon_vertex( { p, v } ).local_vertex_id;
            if( inexact_equal( surface.mesh_element_vertex( { p, prev_v } ), v1,
                    surface.geomodel().epsilon() ) )
            {
                return prev_v;
            }
            if( inexact_equal( surface.mesh_element_vertex( { p, next_v } ), v1,
                    surface.geomodel().epsilon() ) )
            {
                return v;
            }
        }
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t cell_facet_index_from_cell_and_polygon(
        const Region< DIMENSION >& region,
        index_t cell,
        const Surface< DIMENSION >& surface,
        index_t polygon )
    {
        const auto& mesh = region.mesh();
        auto polygon_barycenter = surface.mesh_element_barycenter( polygon );
        for( auto f : range( region.nb_cell_facets( cell ) ) )
        {
            auto cell_facet_barycenter =
                mesh.cell_facet_barycenter( { cell, f } );
            if( inexact_equal( cell_facet_barycenter, polygon_barycenter,
                    surface.geomodel().epsilon() ) )
            {
                return f;
            }
        }
        return NO_ID;
    }

    template < index_t DIMENSION >
    void check_and_initialize_corner_vertex(
        GeoModel< DIMENSION >& geomodel, index_t corner_id )
    {
        if( geomodel.corner( corner_id ).nb_vertices() == 0 )
        {
            GeoModelBuilder< DIMENSION > builder( geomodel );
            builder.geometry.create_mesh_entity_vertices(
                { corner_type_name_static(), corner_id }, 1 );
        }
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBuilderGeometryBase< DIMENSION >::GeoModelBuilderGeometryBase(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : builder_( builder ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::clear_geomodel_mesh()
    {
        geomodel_.mesh.vertices.clear();
    }

    template < index_t DIMENSION >
    std::unique_ptr< PointSetMeshBuilder< DIMENSION > >
        GeoModelBuilderGeometryBase< DIMENSION >::create_corner_builder(
            index_t corner_id )
    {
        gmme_id id{ corner_type_name_static(), corner_id };
        auto& corner = geomodel_access_.modifiable_mesh_entity( id );
        GeoModelMeshEntityAccess< DIMENSION > corner_access( corner );
        auto& corner_mesh = dynamic_cast< PointSetMesh< DIMENSION >& >(
            *corner_access.modifiable_mesh() );
        return PointSetMeshBuilder< DIMENSION >::create_builder( corner_mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< LineMeshBuilder< DIMENSION > >
        GeoModelBuilderGeometryBase< DIMENSION >::create_line_builder(
            index_t line_id )
    {
        gmme_id id{ line_type_name_static(), line_id };
        auto& line = geomodel_access_.modifiable_mesh_entity( id );
        GeoModelMeshEntityAccess< DIMENSION > line_access( line );
        auto& line_mesh = dynamic_cast< LineMesh< DIMENSION >& >(
            *line_access.modifiable_mesh() );
        return LineMeshBuilder< DIMENSION >::create_builder( line_mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >
        GeoModelBuilderGeometryBase< DIMENSION >::create_surface_builder(
            index_t surface_id )
    {
        gmme_id id{ surface_type_name_static(), surface_id };
        auto& surface = geomodel_access_.modifiable_mesh_entity( id );
        GeoModelMeshEntityAccess< DIMENSION > surface_access( surface );
        auto& surface_mesh = dynamic_cast< SurfaceMesh< DIMENSION >& >(
            *surface_access.modifiable_mesh() );
        return SurfaceMeshBuilder< DIMENSION >::create_builder( surface_mesh );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::copy_meshes(
        const GeoModel< DIMENSION >& from )
    {
        copy_meshes( from, corner_type_name_static() );
        copy_meshes( from, line_type_name_static() );
        copy_meshes( from, surface_type_name_static() );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_mesh_entity_vertex(
        const gmme_id& entity_id,
        index_t v,
        const vecn< DIMENSION >& point,
        bool update )
    {
        auto& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        ringmesh_assert( v < E.nb_vertices() );
        if( update )
        {
            auto& geomodel_vertices = geomodel_.mesh.vertices;
            set_mesh_entity_vertex(
                geomodel_vertices.geomodel_vertex_id( E.gmme(), v ), point );
        }
        else
        {
            GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
            auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
                *gmme_access.modifiable_mesh() );
            builder->set_vertex( v, point );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_mesh_entity_vertex(
        index_t geomodel_vertex_id, const vecn< DIMENSION >& point )
    {
        auto& geomodel_vertices = geomodel_.mesh.vertices;
        geomodel_vertices.set_point( geomodel_vertex_id, point );

        const auto& gme_v =
            geomodel_vertices.gme_vertices( geomodel_vertex_id );
        for( const auto& info : gme_v )
        {
            set_mesh_entity_vertex( info.gmme, info.v_index, point, false );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_mesh_entity_vertex(
        const gmme_id& entity_id, index_t v, index_t geomodel_vertex )
    {
        auto& geomodel_vertices = geomodel_.mesh.vertices;
        set_mesh_entity_vertex(
            entity_id, v, geomodel_vertices.vertex( geomodel_vertex ), false );

        ringmesh_assert( v < geomodel_.mesh_entity( entity_id ).nb_vertices() );
        geomodel_vertices.update_vertex_mapping(
            entity_id, v, geomodel_vertex );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_mesh_entity_vertices(
        const gmme_id& entity_id,
        const std::vector< vecn< DIMENSION > >& points,
        bool clear )
    {
        auto& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        // Clear the mesh, but keep the attributes and the space
        if( clear )
        {
            builder->clear( true, true );
        }
        if( !points.empty() )
        {
            auto nb_points = static_cast< index_t >( points.size() );
            auto start = builder->create_vertices( nb_points );
            for( auto v : range( nb_points ) )
            {
                builder->set_vertex( start + v, points[v] );
            }
        }
    }

    template < index_t DIMENSION >
    index_t
        GeoModelBuilderGeometryBase< DIMENSION >::create_mesh_entity_vertices(
            const gmme_id& entity_id, index_t nb_vertices )
    {
        auto& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        return builder->create_vertices( nb_vertices );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_mesh_entity_vertices(
        const gmme_id& entity_id,
        const std::vector< index_t >& geomodel_vertices,
        bool clear )
    {
        auto& E = geomodel_access_.modifiable_mesh_entity( entity_id );
        GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        // Clear the mesh, but keep the attributes and the space
        if( clear )
        {
            builder->clear( true, true );
        }
        auto nb_model_vertices =
            static_cast< index_t >( geomodel_vertices.size() );
        auto start = builder->create_vertices( nb_model_vertices );
        for( auto v : range( nb_model_vertices ) )
        {
            set_mesh_entity_vertex(
                entity_id, start + v, geomodel_vertices[v] );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_corner(
        index_t corner_id, const vecn< DIMENSION >& point )
    {
        check_and_initialize_corner_vertex( geomodel_, corner_id );
        set_mesh_entity_vertex(
            { corner_type_name_static(), corner_id }, 0, point, false );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_line(
        index_t line_id, const std::vector< vecn< DIMENSION > >& vertices )
    {
        set_mesh_entity_vertices(
            { line_type_name_static(), line_id }, vertices, true );

        auto& line = dynamic_cast< Line< DIMENSION >& >(
            geomodel_access_.modifiable_mesh_entity(
                { line_type_name_static(), line_id } ) );
        auto builder = create_line_builder( line_id );
        for( auto e : range( 1, line.nb_vertices() ) )
        {
            builder->create_edge( e - 1, e );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_surface_geometry(
        index_t surface_id,
        const std::vector< vecn< DIMENSION > >& surface_vertices,
        const std::vector< index_t >& surface_polygons,
        const std::vector< index_t >& surface_polygon_ptr )
    {
        set_mesh_entity_vertices( { surface_type_name_static(), surface_id },
            surface_vertices, true );
        set_surface_geometry(
            surface_id, surface_polygons, surface_polygon_ptr );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& polygons,
        const std::vector< index_t >& polygon_ptr )
    {
        auto builder = create_surface_builder( surface_id );
        builder->create_polygons( polygons, polygon_ptr );
        compute_surface_adjacencies( surface_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_corner(
        index_t corner_id, index_t geomodel_vertex_id )
    {
        check_and_initialize_corner_vertex( geomodel_, corner_id );
        set_mesh_entity_vertex(
            { corner_type_name_static(), corner_id }, 0, geomodel_vertex_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_line(
        index_t line_id, const std::vector< index_t >& unique_vertices )
    {
        bool clear_vertices{ false };
        auto& E = geomodel_access_.modifiable_mesh_entity(
            { line_type_name_static(), line_id } );

        ringmesh_assert( E.nb_vertices() == 0 );
        // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_mesh_entity_vertices( E.gmme(), unique_vertices, clear_vertices );

        auto builder = create_line_builder( line_id );
        for( auto e : range( 1, E.nb_vertices() ) )
        {
            builder->create_edge( e - 1, e );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::set_surface_element_geometry(
        index_t surface_id,
        index_t polygon_id,
        const std::vector< index_t >& corners )
    {
        auto builder = create_surface_builder( surface_id );
        for( auto polygon_vertex : range( corners.size() ) )
        {
            builder->set_polygon_vertex(
                { polygon_id, polygon_vertex }, corners[polygon_vertex] );
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderGeometryBase< DIMENSION >::create_surface_polygon(
        index_t surface_id, const std::vector< index_t >& vertex_indices )
    {
        auto builder = create_surface_builder( surface_id );
        return builder->create_polygon( vertex_indices );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::delete_mesh_entity_mesh(
        const gmme_id& E_id )
    {
        GeoModelMeshEntityAccess< DIMENSION > gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->clear( true, false );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase<
        DIMENSION >::delete_mesh_entity_isolated_vertices( const gmme_id& E_id )
    {
        if( geomodel_.entity_type_manager().mesh_entity_manager.is_line(
                E_id.type() ) )
        {
            auto builder = create_line_builder( E_id.index() );
            builder->remove_isolated_vertices();
        }
        else if( geomodel_.entity_type_manager().mesh_entity_manager.is_surface(
                     E_id.type() ) )
        {
            auto builder = create_surface_builder( E_id.index() );
            builder->remove_isolated_vertices();
        }
        else if( geomodel_.entity_type_manager().mesh_entity_manager.is_corner(
                     E_id.type() ) )
        {
            auto builder = create_corner_builder( E_id.index() );
            builder->remove_isolated_vertices();
        }
        else
        {
            ringmesh_assert_not_reached;
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::delete_mesh_entity_vertices(
        const gmme_id& E_id, const std::vector< bool >& to_delete )
    {
        GeoModelMeshEntityAccess< DIMENSION > gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->delete_vertices( to_delete );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::delete_corner_vertex(
        index_t corner_id )
    {
        gmme_id corner{ corner_type_name_static(), corner_id };
        std::vector< bool > to_delete;
        to_delete.reserve( 1 );
        to_delete.push_back( true );
        delete_mesh_entity_vertices( corner, to_delete );
    }
    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::delete_line_edges(
        index_t line_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        auto builder = create_line_builder( line_id );
        builder->delete_edges( to_delete, remove_isolated_vertices );
    }
    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::delete_surface_polygons(
        index_t surface_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        auto builder = create_surface_builder( surface_id );
        builder->delete_polygons( to_delete, remove_isolated_vertices );
    }

    template < index_t DIMENSION >
    void
        GeoModelBuilderGeometryBase< DIMENSION >::set_surface_element_adjacency(
            index_t surface_id,
            index_t polygon_id,
            const std::vector< index_t >& adjacents )
    {
        auto builder = create_surface_builder( surface_id );
        for( auto polygon_edge : range( adjacents.size() ) )
        {
            builder->set_polygon_adjacent(
                { polygon_id, polygon_edge }, adjacents[polygon_edge] );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::compute_surface_adjacencies(
        index_t surface_id, bool recompute_adjacency )
    {
        const auto& surface = geomodel_.surface( surface_id );
        auto builder = create_surface_builder( surface_id );

        if( recompute_adjacency )
        {
            for( auto p : range( surface.nb_mesh_elements() ) )
            {
                for( auto v : range( surface.nb_mesh_element_vertices( p ) ) )
                {
                    builder->set_polygon_adjacent( { p, v }, NO_ID );
                }
            }
        }
        builder->connect_polygons();
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase<
        DIMENSION >::cut_surfaces_by_internal_lines()
    {
        for( const auto& surface : geomodel_.surfaces() )
        {
            auto cutting_lines = get_internal_borders( surface );
            for( auto line_id : cutting_lines )
            {
                cut_surface_by_line( surface.index(), line_id );
            }
            if( !cutting_lines.empty() )
            {
                auto surface_mesh_builder =
                    create_surface_builder( surface.index() );
                surface_mesh_builder->remove_isolated_vertices();
            }
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::cut_surface_by_line(
        index_t surface_id, index_t line_id )
    {
        auto nb_disconnected_edges =
            disconnect_surface_polygons_along_line_edges( surface_id, line_id );
        if( nb_disconnected_edges > 0 )
        {
            duplicate_surface_vertices_along_line( surface_id, line_id );
        }
    }

    struct ElementVertex
    {
        index_t element_;
        index_t vertex_;
    };

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase<
        DIMENSION >::duplicate_surface_vertices_along_line( index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );
        ringmesh_assert( line_id < geomodel_.nb_lines() );

        gmme_id surface_gme{ surface_type_name_static(), surface_id };
        const auto& surface = geomodel_.surface( surface_id );
        const auto& line = geomodel_.line( line_id );

        std::vector< ElementVertex > polygon_vertices( line.nb_vertices() );
        for( auto v : range( line.nb_vertices() ) )
        {
            const auto& p = line.vertex( v );

            auto& polygon_vertex = polygon_vertices[v].vertex_;
            auto& polygon = polygon_vertices[v].element_;
            bool found{ find_polygon_from_vertex(
                surface, p, polygon, polygon_vertex ) };
            ringmesh_unused( found );
            ringmesh_assert(
                found && polygon != NO_ID && polygon_vertex != NO_ID );
        }

        auto vertex_id =
            create_mesh_entity_vertices( surface_gme, line.nb_vertices() );
        auto surface_mesh_builder = create_surface_builder( surface_id );
        const auto& mesh = surface.mesh();
        for( auto v : range( line.nb_vertices() ) )
        {
            const auto& p = line.vertex( v );
            const auto polygon_vertex = polygon_vertices[v].vertex_;
            const auto polygon = polygon_vertices[v].element_;

            auto polygons =
                mesh.polygons_around_vertex( polygon_vertex, false, polygon );
            update_polygon_vertex(
                surface_id, polygons, polygon_vertex, vertex_id );
            surface_mesh_builder->set_vertex( vertex_id, p );
            vertex_id++;
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderGeometryBase< DIMENSION >::
        disconnect_surface_polygons_along_line_edges(
            index_t surface_id, index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );
        ringmesh_assert( line_id < geomodel_.nb_lines() );

        const auto& surface = geomodel_.surface( surface_id );
        const auto& line = geomodel_.line( line_id );
        auto builder = create_surface_builder( surface_id );
        index_t nb_disconnected_edges{ 0 };
        for( auto i : range( line.nb_mesh_elements() ) )
        {
            const auto& p0 = line.vertex( i );
            const auto& p1 = line.vertex( i + 1 );

            bool found{ false };
            index_t p{ NO_ID };
            index_t e{ NO_ID };
            std::tie( found, p, e ) =
                find_polygon_from_edge_vertices( surface, p0, p1 );
            ringmesh_unused( found );
            ringmesh_assert( found && p != NO_ID && e != NO_ID );

            auto adj_f = surface.polygon_adjacent_index( { p, e } );
            if( adj_f != NO_ID )
            {
                auto adj_e = edge_index_from_polygon_and_edge_vertex_indices(
                    surface, adj_f, p0, p1 );
                ringmesh_assert( adj_e != NO_ID );
                builder->set_polygon_adjacent( { p, e }, NO_ID );
                builder->set_polygon_adjacent( { adj_f, adj_e }, NO_ID );
                nb_disconnected_edges++;
            }
        }
        return nb_disconnected_edges;
    }
    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::update_polygon_vertex(
        index_t surface_id,
        const std::vector< index_t >& polygons,
        index_t old_vertex,
        index_t new_vertex )
    {
        const auto& surface = geomodel_.surface( surface_id );
        auto builder = create_surface_builder( surface_id );
        for( auto cur_p : polygons )
        {
            for( auto cur_v :
                range( surface.nb_mesh_element_vertices( cur_p ) ) )
            {
                if( surface.mesh_element_vertex_index( { cur_p, cur_v } )
                    == old_vertex )
                {
                    builder->set_polygon_vertex( { cur_p, cur_v }, new_vertex );
                }
            }
        }
    }
    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::copy_meshes(
        const GeoModel< DIMENSION >& from, const MeshEntityType& entity_type )
    {
        for( auto i : range( geomodel_.nb_mesh_entities( entity_type ) ) )
        {
            copy_mesh( from, { entity_type, i } );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::copy_mesh(
        const GeoModel< DIMENSION >& from, const gmme_id& mesh_entity )
    {
        const GeoModelMeshEntityConstAccess< DIMENSION > from_E_const_access(
            from.mesh_entity( mesh_entity ) );
        assign_mesh_to_entity( *from_E_const_access.mesh(), mesh_entity );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeometryBase< DIMENSION >::assign_mesh_to_entity(
        const MeshBase< DIMENSION >& mesh, const gmme_id& to )
    {
        auto& E = geomodel_access_.modifiable_mesh_entity( to );
        GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
        auto builder = MeshBaseBuilder< DIMENSION >::create_builder(
            *gmme_access.modifiable_mesh() );
        builder->copy( mesh, true );
    }

    void GeoModelBuilderGeometry< 3 >::update_cell_vertex( index_t region_id,
        const std::vector< index_t >& cells,
        index_t old_vertex,
        index_t new_vertex )
    {
        const auto& region = geomodel_.region( region_id );
        auto builder = create_region_builder( region_id );
        for( auto cur_c : cells )
        {
            for( auto cur_v :
                range( region.nb_mesh_element_vertices( cur_c ) ) )
            {
                if( region.mesh_element_vertex_index( { cur_c, cur_v } )
                    == old_vertex )
                {
                    builder->set_cell_vertex( { cur_c, cur_v }, new_vertex );
                }
            }
        }
    }

    void GeoModelBuilderGeometry< 3 >::assign_region_tet_mesh(
        index_t region_id, const std::vector< index_t >& tet_vertices )
    {
        auto builder = create_region_builder( region_id );
        builder->assign_cell_tet_mesh( tet_vertices );
        builder->connect_cells();
    }

    std::unique_ptr< VolumeMeshBuilder3D >
        GeoModelBuilderGeometry< 3 >::create_region_builder( index_t region_id )
    {
        gmme_id id{ region_type_name_static(), region_id };
        auto& region = geomodel_access_.modifiable_mesh_entity( id );
        GeoModelMeshEntityAccess3D region_access( region );
        auto& region_mesh =
            dynamic_cast< VolumeMesh3D& >( *region_access.modifiable_mesh() );
        return VolumeMeshBuilder3D::create_builder( region_mesh );
    }

    index_t GeoModelBuilderGeometry<
        3 >::disconnect_region_cells_along_surface_polygons( index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() );
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );

        const auto& region = geomodel_.region( region_id );
        const auto& surface = geomodel_.surface( surface_id );
        auto builder = create_region_builder( region_id );
        index_t nb_disconnected_polygons{ 0 };
        for( auto polygon : range( surface.nb_mesh_elements() ) )
        {
            bool found{ false };
            index_t cell{ NO_ID };
            index_t cell_facet{ NO_ID };
            std::tie( found, cell, cell_facet ) =
                find_cell_facet_from_polygon( region, surface, polygon );
            ringmesh_unused( found );
            ringmesh_assert( found && cell != NO_ID && cell_facet != NO_ID );

            auto adj_cell = region.cell_adjacent_index( cell, cell_facet );
            if( adj_cell != NO_ID )
            {
                auto adj_cell_facet = cell_facet_index_from_cell_and_polygon(
                    region, adj_cell, surface, polygon );
                ringmesh_assert( adj_cell_facet != NO_ID );
                builder->set_cell_adjacent( { cell, cell_facet }, NO_ID );
                builder->set_cell_adjacent(
                    { adj_cell, adj_cell_facet }, NO_ID );
                nb_disconnected_polygons++;
            }
        }
        return nb_disconnected_polygons;
    }

    void GeoModelBuilderGeometry< 3 >::duplicate_region_vertices_along_surface(
        index_t region_id, index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() );
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() );

        gmme_id region_gme{ region_type_name_static(), region_id };
        const auto& region = geomodel_.region( region_id );
        const auto& surface = geomodel_.surface( surface_id );

        std::vector< ElementVertex > cell_vertices( surface.nb_vertices() );
        for( auto v : range( surface.nb_vertices() ) )
        {
            const auto& p = surface.vertex( v );

            auto& cell = cell_vertices[v].element_;
            auto& cell_vertex = cell_vertices[v].vertex_;
            auto element_local_vertex =
                region.find_cell_from_colocated_vertex_if_any( p );
            cell = element_local_vertex.element_id;
            cell_vertex = element_local_vertex.local_vertex_id;

            ringmesh_assert( cell != NO_ID && cell_vertex != NO_ID );
        }

        GEO::vector< std::string > names;
        geomodel_.region( region_id )
            .vertex_attribute_manager()
            .list_attribute_names( names );

        auto& E = geomodel_access_.modifiable_mesh_entity( region_gme );
        auto vertices_nb = E.nb_vertices();

        auto vertex_id =
            create_mesh_entity_vertices( region_gme, surface.nb_vertices() );

        auto region_mesh_builder = create_region_builder( region_id );
        const auto& mesh = region.mesh();
        for( auto v : range( surface.nb_vertices() ) )
        {
            const auto& p = surface.vertex( v );
            const auto cell = cell_vertices[v].element_;
            const auto cell_vertex = cell_vertices[v].vertex_;

            auto cells = mesh.cells_around_vertex( cell_vertex, cell );
            update_cell_vertex( region_id, cells, cell_vertex, vertex_id );
            region_mesh_builder->set_vertex( vertex_id, p );
            vertex_id++;

            for( const auto& name : names )
            {
                if( name != "model_vertex_map" && name != "point" )
                {
                    GEO::Attribute< double > attr(
                        geomodel_.region( region_id )
                            .vertex_attribute_manager(),
                        name );
                    auto dim_nb = attr.dimension();
                    for( auto dim : range( dim_nb ) )
                    {
                        attr[( vertices_nb + v ) * dim_nb + dim] =
                            attr[cell_vertex * dim_nb + dim];
                    }
                }
            }
        }
    }

    void GeoModelBuilderGeometry< 3 >::cut_region_by_surface(
        index_t region_id, index_t surface_id )
    {
        auto nb_disconnected_polygons =
            disconnect_region_cells_along_surface_polygons(
                region_id, surface_id );
        if( nb_disconnected_polygons > 0 )
        {
            duplicate_region_vertices_along_surface( region_id, surface_id );
        }
    }

    void GeoModelBuilderGeometry< 3 >::cut_regions_by_internal_surfaces()
    {
        for( const auto& region : geomodel_.regions() )
        {
            if( region.nb_mesh_elements() == 0 )
            {
                continue;
            }
            auto cutting_surfaces = get_internal_borders( region );
            for( auto surface_id : cutting_surfaces )
            {
                cut_region_by_surface( region.index(), surface_id );
            }
            if( !cutting_surfaces.empty() )
            {
                auto region_mesh_builder =
                    create_region_builder( region.index() );
                region_mesh_builder->remove_isolated_vertices();
            }
        }
    }

    void GeoModelBuilderGeometry< 3 >::set_region_geometry( index_t region_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& tetras )
    {
        set_mesh_entity_vertices(
            { region_type_name_static(), region_id }, points, true );
        assign_region_tet_mesh( region_id, tetras );
    }

    void GeoModelBuilderGeometry< 3 >::compute_region_adjacencies(
        index_t region_id, bool recompute_adjacency )
    {
        const auto& region = geomodel_.region( region_id );
        auto builder = create_region_builder( region_id );
        if( recompute_adjacency )
        {
            for( auto c : range( region.nb_mesh_elements() ) )
            {
                for( auto f : range( region.nb_cell_facets( c ) ) )
                {
                    builder->set_cell_adjacent( { c, f }, NO_ID );
                }
            }
        }
        builder->connect_cells();
    }

    void GeoModelBuilderGeometry< 3 >::delete_region_cells( index_t region_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        auto builder = create_region_builder( region_id );
        builder->delete_cells( to_delete, remove_isolated_vertices );
    }

    index_t GeoModelBuilderGeometry< 3 >::create_region_cells(
        index_t region_id, CellType type, index_t nb_cells )
    {
        auto builder = create_region_builder( region_id );
        return builder->create_cells( nb_cells, type );
    }

    void GeoModelBuilderGeometry< 3 >::set_region_element_geometry(
        index_t region_id,
        index_t cell_id,
        const std::vector< index_t >& corners )
    {
        auto builder = create_region_builder( region_id );
        for( auto cell_vertex : range( corners.size() ) )
        {
            builder->set_cell_vertex(
                { cell_id, cell_vertex }, corners[cell_vertex] );
        }
    }
    void GeoModelBuilderGeometry< 3 >::delete_mesh_entity_isolated_vertices(
        const gmme_id& E_id )
    {
        if( geomodel_.entity_type_manager().mesh_entity_manager.is_region(
                E_id.type() ) )
        {
            auto builder = create_region_builder( E_id.index() );
            builder->remove_isolated_vertices();
        }
        else
        {
            GeoModelBuilderGeometryBase<
                3 >::delete_mesh_entity_isolated_vertices( E_id );
        }
    }

    index_t GeoModelBuilderGeometry< 3 >::create_region_cell( index_t region_id,
        CellType type,
        const std::vector< index_t >& vertex_indices )
    {
        auto cell_id = create_region_cells( region_id, type, 1 );
        set_region_element_geometry( region_id, cell_id, vertex_indices );
        return cell_id;
    }

    void GeoModelBuilderGeometry< 3 >::copy_meshes( const GeoModel3D& geomodel )
    {
        GeoModelBuilderGeometryBase3D::copy_meshes( geomodel );
        GeoModelBuilderGeometryBase3D::copy_meshes(
            geomodel, region_type_name_static() );
    }

    template class geomodel_builder_api GeoModelBuilderGeometry< 2 >;
    template class geomodel_builder_api GeoModelBuilderGeometryBase< 2 >;

    template class geomodel_builder_api GeoModelBuilderGeometryBase< 3 >;
} // namespace RINGMesh
