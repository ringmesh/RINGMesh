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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh.h>

#include <stack>
#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/basic/algorithm.h>

namespace RINGMesh {

    std::unique_ptr< PointSetMesh > PointSetMesh::create_mesh( const MeshType type )
    {
        MeshType new_type = type;
        if( new_type.empty() ) {
            new_type = GeogramPointSetMesh::type_name_static();
        }
        PointSetMesh* mesh = PointSetMeshFactory::create_object( new_type );
        if( !mesh ) {
            Logger::warn( "PointSetMesh", "Could not create mesh data structure: ",
                new_type );
            Logger::warn( "PointSetMesh", "Falling back to GeogramPointSetMesh data structure" );

            mesh = new GeogramPointSetMesh;
        }
        return std::unique_ptr< PointSetMesh >( mesh );
    }

    std::unique_ptr< LineMesh > LineMesh::create_mesh( const MeshType type )
    {
        MeshType new_type = type;
        if( new_type.empty() ) {
            new_type = GeogramLineMesh::type_name_static();
        }
        LineMesh* mesh = LineMeshFactory::create_object( new_type );
        if( !mesh ) {
            Logger::warn( "LineMesh", "Could not create mesh data structure: ",
                new_type );
            Logger::warn( "LineMesh", "Falling back to GeogramLineMesh data structure" );

            mesh = new GeogramLineMesh;
        }
        return std::unique_ptr< LineMesh >( mesh );
    }

    std::unique_ptr< SurfaceMesh > SurfaceMesh::create_mesh( const MeshType type )
    {
        MeshType new_type = type;
        if( new_type.empty() ) {
            new_type = GeogramSurfaceMesh::type_name_static();
        }
        SurfaceMesh* mesh = SurfaceMeshFactory::create_object( new_type );
        if( !mesh ) {
            Logger::warn( "SurfaceMesh", "Could not create mesh data structure: ",
                new_type );
            Logger::warn( "SurfaceMesh", "Falling back to GeogramSurfaceMesh data structure" );

            mesh = new GeogramSurfaceMesh;
        }
        return std::unique_ptr< SurfaceMesh >( mesh );
    }

    void SurfaceMesh::next_on_border(
        index_t p,
        index_t e,
        index_t& next_p,
        index_t& next_e ) const
    {
        ringmesh_assert( e < nb_polygon_vertices( p ) );
        ringmesh_assert( is_edge_on_border( p, e ) );

        // Global indices in the surfaces
        index_t next_v_id = polygon_vertex( p, next_polygon_vertex( p, e ) );

        // Get the polygons around the shared vertex (next_v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > polygons_around_next_v_id = polygons_around_vertex(
            next_v_id, true, p );
        index_t nb_around = static_cast< index_t >( polygons_around_next_v_id.size() );
        ringmesh_assert( nb_around == 1 || nb_around == 2 );

        next_p = polygons_around_next_v_id[0];

        if( nb_around == 2 ) {
            if( next_p == p ) {
                next_p = polygons_around_next_v_id[1];
            }
            ringmesh_assert( next_p != NO_ID );
            ringmesh_assert( is_polygon_on_border( next_p ) );

            // Local index of next vertex in the next polygon
            next_e = vertex_index_in_polygon( next_p, next_v_id );
            ringmesh_assert( is_edge_on_border( next_p, next_e ) );
        } else if( nb_around == 1 ) {
            // next_v_id must be in two border edges of polygon p
            next_e = vertex_index_in_polygon( next_p, next_v_id );
            ringmesh_assert( is_edge_on_border( next_p, next_e ) );
        }
    }

    void SurfaceMesh::prev_on_border(
        index_t p,
        index_t e,
        index_t& prev_p,
        index_t& prev_e ) const
    {
        ringmesh_assert( e < nb_polygon_vertices( p ) );
        ringmesh_assert( is_edge_on_border( p, e ) );

        // Global indices in the surfaces
        index_t v_id = polygon_vertex( p, e );

        // Get the polygons around the shared vertex (v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > polygons_around_v_id = polygons_around_vertex( v_id, true,
            p );
        index_t nb_around = static_cast< index_t >( polygons_around_v_id.size() );
        ringmesh_assert( nb_around == 1 || nb_around == 2 );

        prev_p = polygons_around_v_id[0];

        if( nb_around == 2 ) {
            if( prev_p == p ) {
                prev_p = polygons_around_v_id[1];
            }
            ringmesh_assert( prev_p != NO_ID );
            ringmesh_assert( is_polygon_on_border( prev_p ) );

            // Local index of given vertex in the prev polygon
            index_t v_in_prev_f = vertex_index_in_polygon( prev_p, v_id );
            // Local index of previous vertex in the prev polygon
            prev_e = prev_polygon_vertex( prev_p, v_in_prev_f );
            ringmesh_assert( is_edge_on_border( prev_p, prev_e ) );
        } else if( nb_around == 1 ) {
            // v_id must be in two border edges of polygon p
            index_t v_in_next_polygon = vertex_index_in_polygon( prev_p, v_id );
            prev_e = prev_polygon_vertex( prev_p, v_in_next_polygon );
            ringmesh_assert( is_edge_on_border( prev_p, prev_e ) );
        }
    }

    index_t SurfaceMesh::polygon_from_vertex_ids( index_t in0, index_t in1 ) const
    {
        ringmesh_assert( in0 < nb_vertices() && in1 < nb_vertices() );

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in polygons_ and check that they are in the same polygon

        // Check if the edge is in one of the polygon
        for( index_t poly = 0; poly < nb_polygons(); ++poly ) {
            bool found = false;
            index_t prev = polygon_vertex( poly, nb_polygon_vertices( poly ) - 1 );
            for( index_t v = 0; v < nb_polygon_vertices( poly ); ++v ) {
                index_t p = polygon_vertex( poly, v );
                if( ( prev == in0 && p == in1 ) || ( prev == in1 && p == in0 ) ) {
                    found = true;
                    break;
                }
                prev = p;
            }
            if( found ) {
                return poly;
            }
        }
        return NO_ID;
    }

    index_t SurfaceMesh::vertex_index_in_polygon(
        index_t polygon_index,
        index_t vertex_id ) const
    {
        ringmesh_assert( polygon_index < nb_polygons() );
        for( index_t v = 0; v < nb_polygon_vertices( polygon_index ); v++ ) {
            if( polygon_vertex( polygon_index, v ) == vertex_id ) {
                return v;
            }
        }
        return NO_ID;
    }

    index_t SurfaceMesh::closest_vertex_in_polygon( index_t p, const vec3& v ) const
    {
        index_t result = 0;
        double dist = DBL_MAX;
        for( index_t v_id = 0; v_id < nb_polygon_vertices( p ); v_id++ ) {
            double distance = length2( v - vertex( polygon_vertex( p, v_id ) ) );
            if( dist > distance ) {
                dist = distance;
                result = v_id;
            }
        }
        return result;
    }

    std::vector< index_t > SurfaceMesh::polygons_around_vertex(
        index_t surf_vertex_id,
        bool border_only,
        index_t p0 ) const
    {
        index_t cur_p = 0;
        while( p0 == NO_ID && cur_p < nb_polygons() ) {
            for( index_t lv = 0; lv < nb_polygon_vertices( cur_p ); lv++ ) {
                if( polygon_vertex( cur_p, lv ) == surf_vertex_id ) {
                    p0 = cur_p;
                    break;
                }
            }
            cur_p++;
        }
        ringmesh_assert( p0 != NO_ID );

        // Flag the visited polygons
        std::vector< index_t > visited;
        visited.reserve( 10 );

        // Stack of the adjacent polygons
        std::stack< index_t > S;
        S.push( p0 );
        visited.push_back( p0 );

        std::vector< index_t > result;
        result.reserve( 10 );
        do {
            index_t p = S.top();
            S.pop();

            for( index_t v = 0; v < nb_polygon_vertices( p ); ++v ) {
                if( polygon_vertex( p, v ) == surf_vertex_id ) {
                    index_t adj_P = polygon_adjacent( p, v );
                    index_t prev = prev_polygon_vertex( p, v );
                    index_t adj_prev = polygon_adjacent( p, prev );

                    if( adj_P != NO_ID ) {
                        // The edge starting at P is not on the boundary
                        if( !contains( visited, adj_P ) ) {
                            S.push( adj_P );
                            visited.push_back( adj_P );
                        }
                    }
                    if( adj_prev != NO_ID ) {
                        // The edge ending at P is not on the boundary
                        if( !contains( visited, adj_prev ) ) {
                            S.push( adj_prev );
                            visited.push_back( adj_prev );
                        }
                    }

                    if( border_only ) {
                        if( adj_P == NO_ID || adj_prev == NO_ID ) {
                            result.push_back( p );
                        }
                    } else {
                        result.push_back( p );
                    }

                    // We are done with this polygon
                    break;
                }
            }
        } while( !S.empty() );

        return result;
    }

    std::unique_ptr< VolumeMesh > VolumeMesh::create_mesh( const MeshType type )
    {
        MeshType new_type = type;
        if( new_type.empty() ) {
            new_type = GeogramVolumeMesh::type_name_static();
        }
        VolumeMesh* mesh = VolumeMeshFactory::create_object( new_type );
        if( !mesh ) {
            Logger::warn( "VolumeMesh", "Could not create mesh data structure: ",
                new_type );
            Logger::warn( "VolumeMesh", "Falling back to GeogramVolumeMesh data structure" );

            mesh = new GeogramVolumeMesh;
        }
        return std::unique_ptr< VolumeMesh >( mesh );
    }

    /// TODO border_only parameter [BC] and add a test
    std::vector< index_t > VolumeMesh::cells_around_vertex(
        index_t vertex_id,
        index_t cell_hint ) const
    {
        if( cell_hint == NO_ID ) {
            cell_hint = find_first_cell_owing_vertex( vertex_id );
            if( cell_hint == NO_ID ) {
                return std::vector< index_t >();
            }
        }

        // Flag the visited cells
        std::vector< index_t > visited;
        visited.reserve( 10 );

        // Stack of the adjacent cells
        std::stack< index_t > S;
        S.push( cell_hint );
        visited.push_back( cell_hint );

        std::vector< index_t > result;
        result.reserve( 10 );

        do {
            index_t c = S.top();
            S.pop();

            bool cell_includes_vertex = false;
            for( index_t v = 0; v < nb_cell_vertices( c ); v++ ) {
                if( cell_vertex( c, v ) == vertex_id ) {
                    result.push_back( c );
                    cell_includes_vertex = true;
                    break;
                }
            }
            if( !cell_includes_vertex ) {
                continue;
            }

            for( index_t f = 0; f < nb_cell_facets( c ); f++ ) {
                for( index_t v = 0; v < nb_cell_facet_vertices( c, f ); v++ ) {
                    index_t vertex = cell_facet_vertex( c, f, v );
                    if( vertex == vertex_id ) {
                        index_t adj_P = cell_adjacent( c, f );

                        if( adj_P != NO_ID ) {
                            if( !contains( visited, adj_P ) ) {
                                S.push( adj_P );
                                visited.push_back( adj_P );
                            }
                        }
                        break;
                    }
                }
            }
        } while( !S.empty() );

        return result;
    }

    index_t VolumeMesh::find_first_cell_owing_vertex(
        index_t vertex_id_in_region ) const
    {
        const NNSearch& ann_cells = cell_nn_search();
        const vec3& vertex_pos = vertex( vertex_id_in_region );

        index_t nb_neighbors = std::min( index_t( 5 ), nb_cells() );
        std::vector< index_t > neighbors;
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, nb_cells() );
            neighbors.resize( cur_neighbor );
            neighbors = ann_cells.get_neighbors( vertex_pos, cur_neighbor );
            // nb_neighbors can be less than cur_neighbor.
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                index_t c = neighbors[i];
                for( index_t j = 0; j < nb_cell_vertices( c ); j++ ) {
                    if( cell_vertex( c, j ) == vertex_id_in_region ) {
                        return c;
                    }
                }
            }
        } while( nb_cells() != cur_neighbor );

        return NO_ID;
    }

    MeshSet::MeshSet()
    {
        create_point_set_mesh( "" );
        create_line_mesh( "" );
        create_surface_mesh( "" );
        create_volume_mesh( "" );
    }

    void MeshSet::create_point_set_mesh( const MeshType type )
    {
        point_set_mesh = PointSetMesh::create_mesh( type );
    }

    void MeshSet::create_line_mesh( const MeshType type )
    {
        line_mesh = LineMesh::create_mesh( type );
    }

    void MeshSet::create_surface_mesh( const MeshType type )
    {
        surface_mesh = SurfaceMesh::create_mesh( type );
    }

    void MeshSet::create_volume_mesh( const MeshType type )
    {
        volume_mesh = VolumeMesh::create_mesh( type );
    }

} // namespace
