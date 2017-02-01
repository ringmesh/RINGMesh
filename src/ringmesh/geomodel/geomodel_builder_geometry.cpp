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
    using namespace RINGMesh ;

    void get_entity_vertices_and_update_corners(
        std::vector< index_t >& corners,
        std::vector< index_t >& vertices )
    {
        std::vector< index_t > corners_to_vertices ;
        get_unique_input_values_and_mapping< index_t >( corners, vertices,
            corners_to_vertices ) ;
        corners = corners_to_vertices ;
    }

    void get_internal_borders(
        const GeoModelMeshEntity& entity,
        std::set< index_t >& internal_borders )
    {
        for( index_t i = 0; i < entity.nb_boundaries(); ++i ) {
            const GeoModelMeshEntity& border = entity.boundary( i ) ;
            if( border.is_inside_border( entity ) ) {
                internal_borders.insert( border.index() ) ;
            }
        }
    }

    bool inexact_equal( const vec3& v1, const vec3& v2, double epsilon )
    {
        return length( v2 - v1 ) < epsilon ;
    }

    /*!
     * Finds a facet and its edge index that are colocalised with an edge
     * defined by its two geomodel vertex indices
     * @param[in] surface the surface where to find the facet
     * @param[in] v0 the first vertex of the edge
     * @param[in] v1 the second vertex of the edge
     * @param[out] f the found facet index
     * @param[out] e the found edge index
     * @return True if the facet and the edge indices are found
     * @todo RENAME these parameters and split in smaller functions !! [JP]
     */
    bool find_facet_from_edge_vertices(
        const Surface& surface,
        const vec3& v0,
        const vec3& v1,
        index_t& f,
        index_t& e )
    {
        vec3 v_bary = 0.5 * ( v0 + v1 ) ;
        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = surface.facet_nn_search().get_neighbors( v_bary,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                f = neighbors[i] ;
                for( index_t j = 0; j < surface.nb_mesh_element_vertices( f );
                    j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( f, j ), v0,
                        surface.geomodel().epsilon() ) ) {
                        index_t j_next = surface.next_facet_vertex_index( f, j ) ;
                        if( inexact_equal( surface.mesh_element_vertex( f, j_next ),
                            v1, surface.geomodel().epsilon() ) ) {
                            e = j ;
                            return true ;
                        }
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor ) ;

        f = NO_ID ;
        e = NO_ID ;
        return false ;
    }

    bool are_cell_facet_and_facet_equal(
        const Region& region,
        index_t cell,
        index_t cell_facet,
        const Surface& surface,
        index_t facet )
    {
        index_t nb_cell_facet_vertices = region.nb_cell_facet_vertices( cell,
            cell_facet ) ;
        index_t nb_facet_vertices = surface.nb_mesh_element_vertices( facet ) ;
        if( nb_cell_facet_vertices != nb_facet_vertices ) {
            return false ;
        }
        vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell,
            cell_facet ) ;
        vec3 facet_barycenter = surface.mesh_element_barycenter( facet ) ;
        return inexact_equal( cell_facet_barycenter, facet_barycenter,
            region.geomodel().epsilon() ) ;
    }

    bool find_cell_facet_from_facet(
        const Region& region,
        const Surface& surface,
        index_t facet,
        index_t& cell,
        index_t& cell_facet )
    {
        vec3 v_bary = surface.mesh_element_barycenter( facet ) ;
        index_t nb_neighbors = std::min( index_t( 5 ), region.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, region.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = region.cell_nn_search().get_neighbors( v_bary,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                cell = neighbors[i] ;
                for( cell_facet = 0 ; cell_facet < region.nb_cell_facets( cell );
                    cell_facet++ ) {
                    if( are_cell_facet_and_facet_equal( region, cell, cell_facet,
                        surface, facet ) ) {
                        return true ;
                    }
                }
            }
        } while( region.nb_mesh_elements() != cur_neighbor ) ;

        cell = NO_ID ;
        cell_facet = NO_ID ;
        return false ;
    }

    bool find_facet_from_vertex(
        const Surface& surface,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = surface.facet_nn_search().get_neighbors( v, cur_neighbor,
                neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i] ;
                for( index_t j = 0;
                    j < surface.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( element_id, j ),
                        v, surface.geomodel().epsilon() ) ) {
                        vertex_id = surface.mesh_element_vertex_index( element_id,
                            j ) ;
                        return true ;
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor ) ;

        element_id = NO_ID ;
        vertex_id = NO_ID ;
        return false ;
    }

    bool find_cell_from_vertex(
        const Region& entity,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        index_t nb_neighbors = std::min( index_t( 5 ), entity.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, entity.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = entity.cell_nn_search().get_neighbors( v, cur_neighbor,
                neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i] ;
                for( index_t j = 0;
                    j < entity.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( entity.mesh_element_vertex( element_id, j ),
                        v, entity.geomodel().epsilon() ) ) {
                        vertex_id = entity.mesh_element_vertex_index( element_id,
                            j ) ;
                        return true ;
                    }
                }
            }
        } while( entity.nb_mesh_elements() != cur_neighbor ) ;

        element_id = NO_ID ;
        vertex_id = NO_ID ;
        return false ;
    }

    index_t edge_index_from_facet_and_edge_vertex_indices(
        const Surface& surface,
        index_t f,
        const vec3& v0,
        const vec3& v1 )
    {
        for( index_t v = 0; v < surface.nb_mesh_element_vertices( f ); v++ ) {
            if( !inexact_equal( surface.mesh_element_vertex( f, v ), v0,
                surface.geomodel().epsilon() ) ) continue ;
            index_t prev_v = surface.prev_facet_vertex_index( f, v ) ;
            index_t next_v = surface.next_facet_vertex_index( f, v ) ;
            if( inexact_equal( surface.mesh_element_vertex( f, prev_v ), v1,
                surface.geomodel().epsilon() ) ) {
                return prev_v ;
            } else if( inexact_equal( surface.mesh_element_vertex( f, next_v ), v1,
                surface.geomodel().epsilon() ) ) {
                return v ;
            }
        }
        return NO_ID ;
    }

    index_t cell_facet_index_from_cell_and_facet(
        const Region& region,
        index_t cell,
        const Surface& surface,
        index_t facet )
    {
        vec3 facet_barycenter = surface.mesh_element_barycenter( facet ) ;
        for( index_t f = 0; f < region.nb_cell_facets( cell ); f++ ) {
            vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell, f ) ;
            if( inexact_equal( cell_facet_barycenter, facet_barycenter,
                surface.geomodel().epsilon() ) ) {
                return f ;
            }
        }
        return NO_ID ;
    }
}

namespace RINGMesh {
    GeoModelBuilderGeometry::GeoModelBuilderGeometry(
        GeoModelBuilder2& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderGeometry::recompute_geomodel_mesh()
    {
        geomodel_.mesh.vertices.clear() ;
        geomodel_.mesh.vertices.test_and_initialize() ;
    }

    void GeoModelBuilderGeometry::copy_meshes( const GeoModel& geomodel )
    {
        copy_meshes( geomodel, Corner::type_name_static() ) ;
        copy_meshes( geomodel, Line::type_name_static() ) ;
        copy_meshes( geomodel, Surface::type_name_static() ) ;
        copy_meshes( geomodel, Region::type_name_static() ) ;
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertex(
        const gme_t& t,
        index_t v,
        const vec3& point,
        bool update )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( t ) ;
        GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        if( update ) {
            geomodel_vertices.update_point(
                geomodel_vertices.geomodel_vertex_id( E.gme_id(), v ), point ) ;
        } else {
            GeoModelMeshEntityAccess gmme_access( E ) ;
            MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder(
                *gmme_access.modifiable_mesh() ) ;
            builder->set_vertex( v, point ) ;
        }
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertex(
        const gme_t& entity_id,
        index_t v,
        index_t geomodel_vertex )
    {
        GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices ;
        set_mesh_entity_vertex( entity_id, v,
            geomodel_vertices.vertex( geomodel_vertex ), false ) ;

        ringmesh_assert( v < geomodel_.mesh_entity( entity_id ).nb_vertices() ) ;
        geomodel_vertices.update_vertex_mapping( entity_id, v, geomodel_vertex ) ;
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertices(
        const gme_t& id,
        const std::vector< vec3 >& points,
        bool clear )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( id ) ;
        GeoModelMeshEntityAccess gmme_access( E ) ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true ) ;
        }
        if( !points.empty() ) {
            index_t nb_points = static_cast< index_t >( points.size() ) ;
            index_t start = builder->create_vertices( nb_points ) ;
            for( index_t v = 0; v < nb_points; v++ ) {
                builder->set_vertex( start + v, points[v] ) ;
            }
        }
    }

    index_t GeoModelBuilderGeometry::create_mesh_entity_vertices(
        const gme_t& entity_id,
        index_t nb_vertices )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            entity_id ) ;
        GeoModelMeshEntityAccess gmme_access( E ) ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() ) ;
        return builder->create_vertices( nb_vertices ) ;
    }

    void GeoModelBuilderGeometry::set_mesh_entity_vertices(
        const gme_t& entity_id,
        const std::vector< index_t >& geomodel_vertices,
        bool clear )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            entity_id ) ;
        GeoModelMeshEntityAccess gmme_access( E ) ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true ) ;
        }
        index_t nb_model_vertices =
            static_cast< index_t >( geomodel_vertices.size() ) ;
        index_t start = builder->create_vertices( nb_model_vertices ) ;
        for( index_t v = 0; v < nb_model_vertices; v++ ) {
            set_mesh_entity_vertex( entity_id, start + v, geomodel_vertices[v] ) ;
        }
    }

    void GeoModelBuilderGeometry::set_corner( index_t corner_id, const vec3& point )
    {
        set_mesh_entity_vertex( gme_t( Corner::type_name_static(), corner_id ), 0,
            point, false ) ;
    }

    void GeoModelBuilderGeometry::set_line(
        index_t line_id,
        const std::vector< vec3 >& vertices )
    {
        set_mesh_entity_vertices( gme_t( Line::type_name_static(), line_id ),
            vertices, true ) ;

        Line& line = dynamic_cast< Line& >( geomodel_access_.modifiable_mesh_entity(
            gme_t( Line::type_name_static(), line_id ) ) ) ;
        Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
            line.low_level_mesh_storage() ) ;
        for( index_t e = 1; e < line.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e ) ;
        }
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< vec3 >& surface_vertices,
        const std::vector< index_t >& surface_facets,
        const std::vector< index_t >& surface_facet_ptr )
    {
        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            surface_vertices, true ) ;
        assign_surface_mesh_facets( surface_id, surface_facets, surface_facet_ptr ) ;
    }

    void GeoModelBuilderGeometry::set_region_geometry(
        index_t region_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& tetras )
    {
        set_mesh_entity_vertices( gme_t( Region::type_name_static(), region_id ),
            points, true ) ;
        assign_region_tet_mesh( region_id, tetras ) ;
    }

    void GeoModelBuilderGeometry::set_corner(
        index_t corner_id,
        index_t geomodel_vertex_id )
    {
        set_mesh_entity_vertex( gme_t( Corner::type_name_static(), corner_id ), 0,
            geomodel_vertex_id ) ;
    }

    void GeoModelBuilderGeometry::set_line(
        index_t line_id,
        const std::vector< index_t >& unique_vertices )
    {
        bool clear_vertices = false ;
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Line::type_name_static(), line_id ) ) ;

        ringmesh_assert( E.nb_vertices() == 0 ) ;
        // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_mesh_entity_vertices( E.gme_id(), unique_vertices, clear_vertices ) ;

        Line& line = dynamic_cast< Line& >( E ) ;
        Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
            line.low_level_mesh_storage() ) ;
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e ) ;
        }
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& geomodel_vertex_ids,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            geomodel_vertex_ids, false ) ;
        assign_surface_mesh_facets( surface_id, facets, facet_ptr ) ;
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_facets( facets ) ;
        get_entity_vertices_and_update_corners( new_facets, vertices ) ;
        set_surface_geometry( surface_id, vertices, new_facets, facet_ptr ) ;
    }

    void GeoModelBuilderGeometry::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& triangle_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_triangle_corners( triangle_corners ) ;
        get_entity_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            vertices, false ) ;
        assign_surface_triangle_mesh( surface_id, new_triangle_corners ) ;
    }

    void GeoModelBuilderGeometry::set_surface_geometry_with_adjacencies(
        index_t surface_id,
        const std::vector< index_t >& triangle_corners,
        const std::vector< index_t >& adjacent_triangles )
    {
        /// @todo Reorganize to remove duplicated code in the class
        std::vector< index_t > vertices ;
        std::vector< index_t > new_triangle_corners( triangle_corners ) ;
        get_entity_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            vertices, false ) ;

        assign_surface_triangle_mesh( surface_id, new_triangle_corners,
            adjacent_triangles ) ;
    }

    void GeoModelBuilderGeometry::set_surface_element_geometry(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Surface::type_name_static(), surface_id ) ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;

        for( index_t facet_vertex = 0; facet_vertex < corners.size();
            facet_vertex++ ) {
            builder->set_facet_vertex( facet_id, facet_vertex,
                corners[facet_vertex] ) ;
        }
    }

    void GeoModelBuilderGeometry::set_surface_element_adjacency(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& adjacents )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Surface::type_name_static(), surface_id ) ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;

        for( index_t facet_edge = 0; facet_edge < adjacents.size(); facet_edge++ ) {
            builder->set_facet_adjacent( facet_id, facet_edge,
                adjacents[facet_edge] ) ;
        }
    }

    void GeoModelBuilderGeometry::set_region_geometry(
        index_t region_id,
        const std::vector< index_t >& tet_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_tet_corners( tet_corners ) ;
        get_entity_vertices_and_update_corners( new_tet_corners, vertices ) ;

        set_mesh_entity_vertices( gme_t( Region::type_name_static(), region_id ),
            vertices, false ) ;
        assign_region_tet_mesh( region_id, new_tet_corners ) ;
    }

    void GeoModelBuilderGeometry::set_region_element_geometry(
        index_t region_id,
        index_t cell_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Region::type_name_static(), region_id ) ) ;

        Region& region = dynamic_cast< Region& >( E ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;

        for( index_t cell_vertex = 0; cell_vertex < corners.size(); cell_vertex++ ) {
            builder->set_cell_vertex( cell_id, cell_vertex, corners[cell_vertex] ) ;
        }
    }

    index_t GeoModelBuilderGeometry::create_surface_facet(
        index_t surface_id,
        const std::vector< index_t >& vertex_indices )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Surface::type_name_static(), surface_id ) ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        return builder->create_facet_polygon( vertex_indices ) ;
    }

    index_t GeoModelBuilderGeometry::create_region_cells(
        index_t region_id,
        GEO::MeshCellType type,
        index_t nb_cells )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
            gme_t( Region::type_name_static(), region_id ) ) ;
        Region& region = dynamic_cast< Region& >( E ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        return builder->create_cells( nb_cells, type ) ;
    }

    index_t GeoModelBuilderGeometry::create_region_cell(
        index_t region_id,
        GEO::MeshCellType type,
        const std::vector< index_t >& vertex_indices )
    {
        index_t cell_id = create_region_cells( region_id, type, 1 ) ;
        set_region_element_geometry( region_id, cell_id, vertex_indices ) ;
        return cell_id ;
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_mesh( const gme_t& E_id )
    {
        GeoModelMeshEntityAccess gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) ) ;
        MeshBase& M = *gmme_access.modifiable_mesh() ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder( M ) ;
        builder->clear( true, false ) ;
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_isolated_vertices(
        const gme_t& E_id )
    {
        if( geomodel_.entity_type_manager().is_line( E_id.type ) ) {
            Line& line =
                dynamic_cast< Line& >( geomodel_access_.modifiable_mesh_entity(
                    E_id ) ) ;
            Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
                line.low_level_mesh_storage() ) ;
            builder->remove_isolated_vertices() ;
        } else if( geomodel_.entity_type_manager().is_surface( E_id.type ) ) {
            Surface& surface =
                dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                    E_id ) ) ;
            Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
                surface.low_level_mesh_storage() ) ;
            builder->remove_isolated_vertices() ;
        } else if( geomodel_.entity_type_manager().is_region( E_id.type ) ) {
            Region& region =
                dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                    E_id ) ) ;
            Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
                region.low_level_mesh_storage() ) ;
            builder->remove_isolated_vertices() ;
        } else if( geomodel_.entity_type_manager().is_corner( E_id.type ) ) {
            Corner& corner =
                dynamic_cast< Corner& >( geomodel_access_.modifiable_mesh_entity(
                    E_id ) ) ;
            Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder(
                corner.low_level_mesh_storage() ) ;
            builder->remove_isolated_vertices() ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }

    void GeoModelBuilderGeometry::delete_mesh_entity_vertices(
        const gme_t& E_id,
        const std::vector< bool >& to_delete )
    {
        GeoModelMeshEntityAccess gmme_access(
            geomodel_access_.modifiable_mesh_entity( E_id ) ) ;
        MeshBase& M = *gmme_access.modifiable_mesh() ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder( M ) ;
        builder->delete_vertices( to_delete ) ;
    }

    void GeoModelBuilderGeometry::delete_corner_vertex( index_t corner_id )
    {
        gme_t corner( Corner::type_name_static(), corner_id ) ;
        std::vector< bool > to_delete ;
        to_delete.push_back( true ) ;
        delete_mesh_entity_vertices( corner, to_delete ) ;
    }
    void GeoModelBuilderGeometry::delete_line_edges(
        index_t line_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        Line& line = dynamic_cast< Line& >( geomodel_access_.modifiable_mesh_entity(
            gme_t( Line::type_name_static(), line_id ) ) ) ;
        Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
            line.low_level_mesh_storage() ) ;
        builder->delete_edges( to_delete, remove_isolated_vertices ) ;
    }
    void GeoModelBuilderGeometry::delete_surface_facets(
        index_t surface_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        builder->delete_facets( to_delete, remove_isolated_vertices ) ;
    }
    void GeoModelBuilderGeometry::delete_region_cells(
        index_t region_id,
        const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        Region& region =
            dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Region::type_name_static(), region_id ) ) ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        builder->delete_cells( to_delete, remove_isolated_vertices ) ;
    }

    void GeoModelBuilderGeometry::set_surface_facet_adjacencies(
        index_t surface_id,
        const std::vector< index_t >& facets_id,
        const std::vector< index_t >& edges_id,
        const std::vector< index_t >& adjacent_triangles )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        ringmesh_assert( facets_id.size() == edges_id.size() &&
            facets_id.size() == adjacent_triangles.size() ) ;
        for( index_t i = 0; i < facets_id.size(); ++i ) {
            builder->set_facet_adjacent( facets_id[i], edges_id[i],
                adjacent_triangles[i] ) ;
        }
    }

    void GeoModelBuilderGeometry::compute_surface_adjacencies(
        index_t surface_id,
        bool recompute_adjacency )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;

        if( recompute_adjacency ) {
            for( index_t f = 0; f < surface.nb_mesh_elements(); f++ ) {
                for( index_t v = 0; v < surface.nb_mesh_element_vertices( f );
                    v++ ) {
                    builder->set_facet_adjacent( f, v, NO_ID ) ;
                }
            }
        }
        builder->connect_facets() ;
    }

    void GeoModelBuilderGeometry::compute_region_adjacencies(
        index_t region_id,
        bool recompute_adjacency )
    {
        Region& region =
            dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Region::type_name_static(), region_id ) ) ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;

        if( recompute_adjacency ) {
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                for( index_t f = 0; f < region.nb_cell_facets( c ); f++ ) {
                    builder->set_cell_adjacent( c, f, NO_ID ) ;
                }
            }
        }
        builder->connect_cells() ;
    }

    void GeoModelBuilderGeometry::cut_surfaces_by_internal_lines()
    {
        for( index_t s = 0; s < geomodel_.nb_surfaces(); s++ ) {
            Surface& surface =
                dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                    gme_t( Surface::type_name_static(), s ) ) ) ;
            std::set< index_t > cutting_lines ;
            get_internal_borders( surface, cutting_lines ) ;
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                cut_surface_by_line( s, *it ) ;
            }
            if( !cutting_lines.empty() ) {
                Mesh2DBuilder_var surface_mesh_builder =
                    Mesh2DBuilder::create_builder(
                        surface.low_level_mesh_storage() ) ;
                surface_mesh_builder->remove_isolated_vertices() ;
            }
        }
    }

    void GeoModelBuilderGeometry::cut_regions_by_internal_surfaces()
    {
        for( index_t r = 0; r < geomodel_.nb_regions(); r++ ) {
            Region& region =
                dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                    gme_t( Region::type_name_static(), r ) ) ) ;
            if( region.nb_mesh_elements() == 0 ) continue ;
            std::set< index_t > cutting_surfaces ;
            get_internal_borders( region, cutting_surfaces ) ;
            for( std::set< index_t >::iterator it = cutting_surfaces.begin();
                it != cutting_surfaces.end(); ++it ) {
                cut_region_by_surface( r, *it ) ;
            }
            if( !cutting_surfaces.empty() ) {
                Mesh3DBuilder_var region_mesh_builder =
                    Mesh3DBuilder::create_builder(
                        region.low_level_mesh_storage() ) ;
                region_mesh_builder->remove_isolated_vertices() ;
            }
        }
    }

    void GeoModelBuilderGeometry::cut_surface_by_line(
        index_t surface_id,
        index_t line_id )
    {
        index_t nb_disconnected_edges = disconnect_surface_facets_along_line_edges(
            surface_id, line_id ) ;
        if( nb_disconnected_edges > 0 ) {
            duplicate_surface_vertices_along_line( surface_id, line_id ) ;
        }
    }

    void GeoModelBuilderGeometry::cut_region_by_surface(
        index_t region_id,
        index_t surface_id )
    {
        index_t nb_disconnected_facets =
            disconnect_region_cells_along_surface_facets( region_id, surface_id ) ;
        if( nb_disconnected_facets > 0 ) {
            duplicate_region_vertices_along_surface( region_id, surface_id ) ;
        }
    }

    struct ElementVertex {
        index_t element_ ;
        index_t vertex_ ;
    } ;

    void GeoModelBuilderGeometry::duplicate_surface_vertices_along_line(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() ) ;
        ringmesh_assert( line_id < geomodel_.nb_lines() ) ;

        gme_t surface_gme( Surface::type_name_static(), surface_id ) ;
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                surface_gme ) ) ;
        const Line& line = geomodel_.line( line_id ) ;

        std::vector< ElementVertex > facet_vertices( line.nb_vertices() ) ;
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v ) ;

            index_t& facet_vertex = facet_vertices[v].vertex_ ;
            index_t& facet = facet_vertices[v].element_ ;
            bool found = find_facet_from_vertex( surface, p, facet, facet_vertex ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && facet != NO_ID && facet_vertex != NO_ID ) ;
        }

        index_t vertex_id = create_mesh_entity_vertices( surface_gme,
            line.nb_vertices() ) ;
        Mesh2DBuilder_var surface_mesh_builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v ) ;
            const index_t& facet_vertex = facet_vertices[v].vertex_ ;
            const index_t& facet = facet_vertices[v].element_ ;

            std::vector< index_t > facets ;
            surface.facets_around_vertex( facet_vertex, facets, false, facet ) ;
            update_facet_vertex( surface, facets, facet_vertex, vertex_id ) ;
            surface_mesh_builder->set_vertex( vertex_id, p ) ;
            vertex_id++ ;
        }
    }

    void GeoModelBuilderGeometry::duplicate_region_vertices_along_surface(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() ) ;
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() ) ;

        gme_t region_gme( Region::type_name_static(), region_id ) ;
        Region& region =
            dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                region_gme ) ) ;
        const Surface& surface = geomodel_.surface( surface_id ) ;

        std::vector< ElementVertex > cell_vertices( surface.nb_vertices() ) ;
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v ) ;

            index_t& cell = cell_vertices[v].element_ ;
            index_t& cell_vertex = cell_vertices[v].vertex_ ;
            bool found = find_cell_from_vertex( region, p, cell, cell_vertex ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && cell != NO_ID && cell_vertex != NO_ID ) ;

        }

        index_t vertex_id = create_mesh_entity_vertices( region_gme,
            surface.nb_vertices() ) ;
        Mesh3DBuilder_var region_mesh_builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v ) ;
            const index_t& cell = cell_vertices[v].element_ ;
            const index_t& cell_vertex = cell_vertices[v].vertex_ ;

            std::vector< index_t > cells ;
            region.cells_around_vertex( cell_vertex, cells, cell ) ;
            update_cell_vertex( region, cells, cell_vertex, vertex_id ) ;
            region_mesh_builder->set_vertex( vertex_id, p ) ;
            vertex_id++ ;
        }
    }

    index_t GeoModelBuilderGeometry::disconnect_surface_facets_along_line_edges(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() ) ;
        ringmesh_assert( line_id < geomodel_.nb_lines() ) ;

        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        const Line& line = geomodel_.line( line_id ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        index_t nb_disconnected_edges = 0 ;
        for( index_t i = 0; i < line.nb_mesh_elements(); ++i ) {
            const vec3& p0 = line.vertex( i ) ;
            const vec3& p1 = line.vertex( i + 1 ) ;

            index_t f = NO_ID ;
            index_t e = NO_ID ;
            bool found = find_facet_from_edge_vertices( surface, p0, p1, f, e ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && f != NO_ID && e != NO_ID ) ;

            index_t adj_f = surface.facet_adjacent_index( f, e ) ;
            if( adj_f != NO_ID ) {
                index_t adj_e = edge_index_from_facet_and_edge_vertex_indices(
                    surface, adj_f, p0, p1 ) ;
                ringmesh_assert( adj_e != NO_ID ) ;
                builder->set_facet_adjacent( f, e, NO_ID ) ;
                builder->set_facet_adjacent( adj_f, adj_e, NO_ID ) ;
                nb_disconnected_edges++ ;
            }
        }
        return nb_disconnected_edges ;
    }
    index_t GeoModelBuilderGeometry::disconnect_region_cells_along_surface_facets(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < geomodel_.nb_regions() ) ;
        ringmesh_assert( surface_id < geomodel_.nb_surfaces() ) ;

        Region& region =
            dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Region::type_name_static(), region_id ) ) ) ;
        const Surface& surface = geomodel_.surface( surface_id ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        index_t nb_disconnected_facets = 0 ;
        for( index_t facet = 0; facet < surface.nb_mesh_elements(); ++facet ) {
            index_t cell = NO_ID ;
            index_t cell_facet = NO_ID ;
            bool found = find_cell_facet_from_facet( region, surface, facet, cell,
                cell_facet ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && cell != NO_ID && cell_facet != NO_ID ) ;

            index_t adj_cell = region.cell_adjacent_index( cell, cell_facet ) ;
            if( adj_cell != NO_ID ) {
                index_t adj_cell_facet = cell_facet_index_from_cell_and_facet(
                    region, adj_cell, surface, facet ) ;
                ringmesh_assert( adj_cell_facet != NO_ID ) ;
                builder->set_cell_adjacent( cell, cell_facet, NO_ID ) ;
                builder->set_cell_adjacent( adj_cell, adj_cell_facet, NO_ID ) ;
                nb_disconnected_facets++ ;
            }
        }
        return nb_disconnected_facets ;
    }

    void GeoModelBuilderGeometry::update_facet_vertex(
        Surface& surface,
        const std::vector< index_t >& facets,
        index_t old_vertex,
        index_t new_vertex )
    {
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t cur_f = facets[i] ;
            for( index_t cur_v = 0;
                cur_v < surface.nb_mesh_element_vertices( cur_f ); cur_v++ ) {
                if( surface.mesh_element_vertex_index( cur_f, cur_v )
                    == old_vertex ) {
                    builder->set_facet_vertex( cur_f, cur_v, new_vertex ) ;
                }
            }
        }
    }

    void GeoModelBuilderGeometry::update_cell_vertex(
        Region& region,
        const std::vector< index_t >& cells,
        index_t old_vertex,
        index_t new_vertex )
    {
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        for( index_t i = 0; i < cells.size(); ++i ) {
            index_t cur_c = cells[i] ;
            for( index_t cur_v = 0; cur_v < region.nb_mesh_element_vertices( cur_c );
                cur_v++ ) {
                if( region.mesh_element_vertex_index( cur_c, cur_v )
                    == old_vertex ) {
                    builder->set_cell_vertex( cur_c, cur_v, new_vertex ) ;
                }
            }
        }
    }

    void GeoModelBuilderGeometry::copy_meshes(
        const GeoModel& from,
        const std::string& entity_type )
    {
        for( index_t i = 0; i < geomodel_.nb_mesh_entities( entity_type ); ++i ) {
            copy_mesh( from, gme_t( entity_type, i ) ) ;
        }
    }

    void GeoModelBuilderGeometry::copy_mesh(
        const GeoModel& from,
        const gme_t& mesh_entity )
    {
        const GeoModelMeshEntityConstAccess from_E_const_access(
            from.mesh_entity( mesh_entity ) ) ;
        assign_mesh_to_entity( *from_E_const_access.mesh(), mesh_entity ) ;
    }

    void GeoModelBuilderGeometry::assign_mesh_to_entity(
        const MeshBase& mesh,
        const gme_t& to )
    {
        GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity( to ) ;
        GeoModelMeshEntityAccess gmme_access( E ) ;
        MeshBaseBuilder_var builder = MeshBaseBuilder::create_builder(
            *gmme_access.modifiable_mesh() ) ;
        builder->copy( mesh, true ) ;
    }

    void GeoModelBuilderGeometry::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        builder->assign_facet_triangle_mesh( triangle_vertices, true ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilderGeometry::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices,
        const std::vector< index_t >& adjacent_triangles )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        builder->assign_facet_triangle_mesh( triangle_vertices, true ) ;

        ringmesh_assert( adjacent_triangles.size() == surface.nb_mesh_elements() * 3 ) ;
        for( index_t f = 0; f < surface.nb_mesh_elements(); f++ ) {
            for( index_t v = 0; v < 3; v++ ) {
                builder->set_facet_adjacent( f, v, adjacent_triangles[3 * f + v] ) ;
            }
        }
    }

    void GeoModelBuilderGeometry::assign_surface_mesh_facets(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        Surface& surface =
            dynamic_cast< Surface& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
            surface.low_level_mesh_storage() ) ;
        builder->create_facet_polygons( facets, facet_ptr ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilderGeometry::assign_region_tet_mesh(
        index_t region_id,
        const std::vector< index_t >& tet_vertices )
    {
        Region& region =
            dynamic_cast< Region& >( geomodel_access_.modifiable_mesh_entity(
                gme_t( Region::type_name_static(), region_id ) ) ) ;
        ringmesh_assert( region.nb_vertices() > 0 ) ;
        Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
            region.low_level_mesh_storage() ) ;
        builder->assign_cell_tet_mesh( tet_vertices, true ) ;
        builder->connect_cells() ;
    }

    GeoModelBuilderGeology::GeoModelBuilderGeology(
        GeoModelBuilder2& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderGeology::fill_mesh_entities_parent( const EntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const std::vector< EntityType > parent_types(
            geomodel_.entity_type_manager().parent_types( type ) ) ;
        for( index_t i = 0; i < parent_types.size(); ++i ) {
            const EntityType& parent_type = parent_types[i] ;
            if( EntityTypeManager::is_defined_type( parent_type ) ) {
                for( index_t j = 0;
                    j < geomodel_.nb_geological_entities( parent_type ); ++j ) {
                    const GeoModelGeologicalEntity& parent =
                        geomodel_.geological_entity( parent_type, j ) ;
                    for( index_t k = 0; k < parent.nb_children(); ++k ) {
                        add_mesh_entity_parent( parent.child_gme( k ),
                            parent.gme_id() ) ;
                    }
                }
            }
        }
    }

    void GeoModelBuilderGeology::fill_geological_entities_children(
        const EntityType& type )
    {
        if( geomodel_.nb_geological_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& c_type =
            geomodel_.geological_entity( type, 0 ).child_type_name() ;
        if( EntityTypeManager::is_defined_type( c_type ) ) {
            for( index_t i = 0; i < geomodel_.nb_mesh_entities( c_type ); ++i ) {
                const GeoModelMeshEntity& p = geomodel_.mesh_entity( c_type, i ) ;
                for( index_t j = 0; j < p.nb_parents(); j++ ) {
                    add_geological_entity_child( p.parent_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelBuilderGeology::complete_mesh_entities_geol_feature_from_first_parent(
        const EntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const std::vector< EntityType > parents =
            geomodel_.entity_type_manager().parent_types( type ) ;
        if( parents.size() == 0 ) {
            return ;
        } else {
            for( index_t i = 0; i < geomodel_.nb_mesh_entities( type ); ++i ) {
                GeoModelMeshEntity& E = geomodel_access_.modifiable_mesh_entity(
                    gme_t( type, i ) ) ;
                if( !E.has_geological_feature() ) {
                    if( E.nb_parents() > 0
                        && E.parent( 0 ).has_geological_feature() ) {
                        GeoModelMeshEntityAccess gmme_access( E ) ;
                        gmme_access.modifiable_geol_feature() =
                            E.parent( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }

    void GeoModelBuilderGeology::complete_geological_entities_geol_feature_from_first_child(
        const EntityType& type )
    {
        if( geomodel_.nb_geological_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& child_type = geomodel_.entity_type_manager().child_type(
            type ) ;
        if( EntityTypeManager::is_defined_type( child_type ) ) {
            for( index_t i = 0; i < geomodel_.nb_geological_entities( type ); ++i ) {
                GeoModelGeologicalEntity& parent =
                    geomodel_access_.modifiable_geological_entity(
                        gme_t( type, i ) ) ;
                if( !parent.has_geological_feature() ) {
                    if( parent.nb_children() > 0
                        && parent.child( 0 ).has_geological_feature() ) {
                        GeoModelGeologicalEntityAccess gmge_access( parent ) ;
                        gmge_access.modifiable_geol_feature() =
                            parent.child( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }
}
