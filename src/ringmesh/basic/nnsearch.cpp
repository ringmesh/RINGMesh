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

#include <ringmesh/basic/nnsearch.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/geogram_extension/geogram_extension.h>

/*!
 * @file Nearest Neighbor requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace RINGMesh {

    NNSearch::NNSearch(
        const GEO::Mesh& mesh,
        const MeshLocation& location,
        bool copy )
        : nn_points_( nullptr ), delete_points_( true )
    {
        nn_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" );
        switch( location ) {
            case VERTICES: {
                build_nn_search_vertices( mesh, copy );
                break;
            }
            case EDGES: {
                build_nn_search_edges( mesh );
                break;
            }
            case FACETS: {
                build_nn_search_polygons( mesh );
                break;
            }
            case CELLS: {
                build_nn_search_cells( mesh );
                break;
            }
            case CELL_FACETS: {
                build_nn_search_cell_facets( mesh );
                break;
            }
            default:
                ringmesh_assert_not_reached;
                break;
        }
    }

    NNSearch::NNSearch( const std::vector< vec3 >& vertices, bool copy )
    {
        index_t nb_vertices = static_cast< index_t >( vertices.size() );
        nn_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" );
        if( copy ) {
            nn_points_ = new double[nb_vertices * 3];
            delete_points_ = true;
            GEO::Memory::copy( nn_points_, vertices.data()->data(),
                3 * nb_vertices * sizeof(double) );
        } else {
            nn_points_ = const_cast< double* >( vertices.data()->data() );
            delete_points_ = false;
        }
        nn_tree_->set_points( nb_vertices, nn_points_ );
    }

    index_t NNSearch::get_colocated_index_mapping(
        double epsilon,
        std::vector< index_t >& index_map ) const
    {
        index_map.resize( nn_tree_->nb_points() );
        for( index_t i = 0; i < index_map.size(); i++ ) {
            index_map[i] = i;
        }
        index_t nb_threads = static_cast< index_t >( omp_get_max_threads() );
        std::vector< index_t > nb_colocalised_per_thread( nb_threads, 0 );
        for( index_t i = 0; i < index_map.size(); i++ ) {
            std::vector< index_t > results = get_neighbors( point( i ), epsilon );
            index_t id = *std::min_element( results.begin(), results.end() );
            if( id < i ) {
                index_map[i] = id;
                index_t thread_id = static_cast< index_t >( omp_get_thread_num() );
                nb_colocalised_per_thread[thread_id]++;
            }
        }

        index_t nb_colocalised_vertices = 0;
        for( index_t nb_colocalised : nb_colocalised_per_thread ) {
            nb_colocalised_vertices += nb_colocalised;
        }
        return nb_colocalised_vertices;
    }

    index_t NNSearch::get_colocated_index_mapping(
        double epsilon,
        std::vector< index_t >& index_map,
        std::vector< vec3 >& unique_points ) const
    {
        index_t nb_colocalised_vertices = get_colocated_index_mapping( epsilon,
            index_map );
        unique_points.reserve( nb_points() - nb_colocalised_vertices );
        index_t offset = 0;
        for( index_t p = 0; p < index_map.size(); p++ ) {
            if( index_map[p] == p ) {
                unique_points.push_back( point( p ) );
                index_map[p] = p - offset;
            } else {
                offset++;
                index_map[p] = index_map[index_map[p]];
            }
        }
        ringmesh_assert( offset == nb_colocalised_vertices );
        return offset;
    }

    std::vector< index_t > NNSearch::get_neighbors(
        const vec3& v,
        double threshold_distance ) const
    {
        std::vector< index_t > result;
        index_t nb_points = nn_tree_->nb_points();
        if( nb_points != 0 ) {
            double threshold_distance_sq = threshold_distance * threshold_distance;
            index_t nb_neighbors = std::min( index_t( 5 ), nb_points );
            index_t cur_neighbor = 0;
            index_t prev_neighbor = 0;
            do {
                prev_neighbor = cur_neighbor;
                cur_neighbor += nb_neighbors;
                result.reserve( cur_neighbor );
                std::vector< index_t > neighbors = get_neighbors( v, cur_neighbor );
                nb_neighbors = static_cast< index_t >( neighbors.size() );
                for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                    if( length2( v - point( neighbors[i] ) )
                        > threshold_distance_sq ) {
                        break;
                    }
                    result.push_back( neighbors[i] );
                }
            } while( result.size() == cur_neighbor && result.size() < nb_points );
        }
        return result;

    }

    std::vector< index_t > NNSearch::get_neighbors(
        const vec3& v,
        index_t nb_neighbors ) const
    {
        std::vector< index_t > result;
        if( nn_tree_->nb_points() != 0 ) {
            nb_neighbors = std::min( nb_neighbors, nn_tree_->nb_points() );
            std::vector< double > distances( nb_neighbors );
            result.resize( nb_neighbors );
            nn_tree_->get_nearest_neighbors( nb_neighbors, v.data(), &result[0],
                &distances[0] );
        }
        return result;
    }

    void NNSearch::build_nn_search_vertices( const GEO::Mesh& mesh, bool copy )
    {
        const GEO::MeshVertices& mesh_vertices = mesh.vertices;
        index_t nb_vertices = mesh_vertices.nb();
        if( nb_vertices == 0 ) {
            return;
        }
        if( !copy ) {
            nn_points_ = const_cast< double* >( mesh_vertices.point_ptr( 0 ) );
            delete_points_ = false;
        } else {
            nn_points_ = new double[nb_vertices * 3];
            GEO::Memory::copy( nn_points_, mesh_vertices.point_ptr( 0 ),
                nb_vertices * 3 * sizeof(double) );
        }
        nn_tree_->set_points( nb_vertices, nn_points_ );
    }

    void NNSearch::build_nn_search_edges( const GEO::Mesh& mesh )
    {
        const GEO::MeshEdges& mesh_edges = mesh.edges;
        index_t nb_edges = mesh_edges.nb();
        if( nb_edges == 0 ) {
            return;
        }
        nn_points_ = new double[nb_edges * 3];
        for( index_t i = 0; i < nb_edges; i++ ) {
            index_t first_vertex_id = mesh_edges.vertex( i, 0 );
            const vec3& first_vertex_vec = mesh.vertices.point( first_vertex_id );
            index_t second_vertex_id = mesh.edges.vertex( i, 1 );
            const vec3& second_vertex_vec = mesh.vertices.point( second_vertex_id );

            vec3 center = ( first_vertex_vec + second_vertex_vec ) / 2.;
            index_t index_in_nn_search = 3 * i;
            fill_nn_search_points( index_in_nn_search, center );
        }
        nn_tree_->set_points( nb_edges, nn_points_ );
    }

    void NNSearch::build_nn_search_polygons( const GEO::Mesh& mesh )
    {
        index_t nb_polygons = mesh.facets.nb();
        if( nb_polygons == 0 ) {
            return;
        }
        nn_points_ = new double[nb_polygons * 3];
        for( index_t i = 0; i < nb_polygons; i++ ) {
            vec3 center = GEO::Geom::mesh_facet_center( mesh, i );
            index_t index_in_nn_search = 3 * i;
            fill_nn_search_points( index_in_nn_search, center );
        }
        nn_tree_->set_points( nb_polygons, nn_points_ );
    }

    void NNSearch::build_nn_search_cell_facets( const GEO::Mesh& mesh )
    {
        index_t nb_cell_facets = mesh.cell_facets.nb();
        nn_points_ = new double[nb_cell_facets * 3];
        index_t index_in_nn_search = 0;
        for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
            for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                vec3 center = mesh_cell_facet_barycenter( mesh, c, f );
                fill_nn_search_points( index_in_nn_search, center );
                index_in_nn_search += 3;
            }
        }
        nn_tree_->set_points( nb_cell_facets, nn_points_ );
    }

    void NNSearch::build_nn_search_cells( const GEO::Mesh& mesh )
    {
        index_t nb_cells = mesh.cells.nb();
        if( nb_cells == 0 ) {
            return;
        }
        nn_points_ = new double[nb_cells * 3];
        for( index_t i = 0; i < nb_cells; i++ ) {
            vec3 center = mesh_cell_barycenter( mesh, i );
            index_t index_in_nn_search = 3 * i;
            fill_nn_search_points( index_in_nn_search, center );
        }
        nn_tree_->set_points( nb_cells, nn_points_ );
    }

    void NNSearch::fill_nn_search_points(
        index_t index_in_nn_search,
        const vec3& center )
    {
        nn_points_[index_in_nn_search] = center.x;
        nn_points_[index_in_nn_search + 1] = center.y;
        nn_points_[index_in_nn_search + 2] = center.z;
    }
}
