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

/*!
 * @file Nearest Neighbor requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

#include <ringmesh/basic/nn_search.h>

#include <numeric>

namespace RINGMesh {

    template< index_t DIMENSION >
    NNSearch< DIMENSION >::NNSearch(
        const std::vector< vecn< DIMENSION > >& vertices,
        bool copy )
    {
        index_t nb_vertices = static_cast< index_t >( vertices.size() );
        nn_tree_ = GEO::NearestNeighborSearch::create( DIMENSION, "BNN" );
        if( copy ) {
            nn_points_ = new double[nb_vertices * DIMENSION];
            delete_points_ = true;
            GEO::Memory::copy( nn_points_, vertices.data()->data(),
                DIMENSION * nb_vertices * sizeof(double) );
        } else {
            nn_points_ = const_cast< double* >( vertices.data()->data() );
            delete_points_ = false;
        }
        nn_tree_->set_points( nb_vertices, nn_points_ );
    }

    template< index_t DIMENSION >
    index_t NNSearch< DIMENSION >::get_colocated_index_mapping(
        double epsilon,
        std::vector< index_t >& index_map ) const
    {
        index_map.resize( nn_tree_->nb_points() );
        std::iota( index_map.begin(), index_map.end(), 0 );
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

    template< index_t DIMENSION >
    index_t NNSearch< DIMENSION >::get_colocated_index_mapping(
        double epsilon,
        std::vector< index_t >& index_map,
        std::vector< vecn< DIMENSION > >& unique_points ) const
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

    template< index_t DIMENSION >
    std::vector< index_t > NNSearch< DIMENSION >::get_neighbors(
        const vecn< DIMENSION >& v,
        double threshold_distance ) const
    {
        double threshold_distance_sq = threshold_distance * threshold_distance;
        return get_neighbors( v, [this, &v, threshold_distance_sq]( index_t i ) {
            return length2( v - point( i ) ) > threshold_distance_sq;} );
    }

    template< index_t DIMENSION >
    std::vector< index_t > NNSearch< DIMENSION >::get_neighbors(
        const vecn< DIMENSION >& v,
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

    template< index_t DIMENSION >
    void NNSearch< DIMENSION >::fill_nn_search_points(
        index_t index_in_nn_search,
        const vecn< DIMENSION >& center )
    {
        for( index_t i = 0; i < DIMENSION; i++ ) {
            nn_points_[index_in_nn_search + i] = center[i];
        }
    }

    template class RINGMESH_API NNSearch< 2 > ;
    template class RINGMESH_API NNSearch< 3 > ;
}
