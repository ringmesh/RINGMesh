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

#pragma once

#include <ringmesh/basic/common.h>

#include <geogram/points/kd_tree.h>

/*!
 * @file Nearest Neighbor requests
 * @author Arnaud Botella
 */

namespace RINGMesh {

    template< index_t DIMENSION >
    class NNSearch {
    ringmesh_disable_copy( NNSearch );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        NNSearch(
            const std::vector< vecn< DIMENSION > >& vertices,
            bool copy = true );

        ~NNSearch()
        {
            if( delete_points_ ) delete[] nn_points_;
        }

        /*!
         * @brief Gets the \p index_map that link all the duplicated points
         * to their first occurancy
         * @return the number of colocated vertices
         * Example:
         *     vertices = [P1, P2, P1, P3, P2, P4]
         *     index_map = [0, 1, 0, 3, 1, 5]
         *     return 2
         */
        std::tuple< index_t, std::vector< index_t > > get_colocated_index_mapping(
            double epsilon ) const;
        /*!
         * @brief Gets the \p index_map that link all the points
         * to a no duplicated list of index in the list of \p unique_points.
         * @return the number of colocated vertices
         * Example:
         *     vertices = [P1, P2, P1, P3, P2, P4]
         *     unique_points = [P1, P2, P3, P4]
         *     index_map = [0, 1, 0, 2, 1, 3]
         *     return 2
         */
        std::tuple< index_t, std::vector< index_t >, std::vector< vecn< DIMENSION > > >
        get_colocated_index_mapping_and_unique_points( double epsilon ) const;
        /*!
         * Gets the closest neighbor point
         * @param[in] v the point to test
         * return returns the index of the closest point
         */
        index_t get_closest_neighbor( const vecn< DIMENSION >& v ) const
        {
            index_t nb_neighbors { 1 };
            return get_neighbors( v, nb_neighbors ).front();
        }

        /*!
         * Compute the neighbors of a given point, point closer than \param threshold_distance
         * @param[in] v the point to test
         * @param[in] threshold_distance distance defining the neighborhood
         * @return the point indices
         */
        std::vector< index_t > get_neighbors(
            const vecn< DIMENSION >& v,
            double threshold_distance ) const;

        /*!
         * Compute the neighbors of a given point according the \param test
         * @param[in] v the point to test
         * @tparam[in] test This functor takes an index of a neighbor point as parameter
         * and returns true to stop the search, false to continue.
         * @return the point indices
         */
        template< typename TEST >
        std::vector< index_t > get_neighbors(
            const vecn< DIMENSION >& v,
            const TEST& test ) const
        {
            std::vector< index_t > result;
            index_t nb_points = nn_tree_->nb_points();
            if( nb_points != 0 ) {
                index_t nb_neighbors { std::min( index_t( 5 ), nb_points ) };
                index_t cur_neighbor { 0 };
                index_t prev_neighbor { 0 };
                do {
                    prev_neighbor = cur_neighbor;
                    cur_neighbor += nb_neighbors;
                    result.reserve( cur_neighbor );
                    std::vector< index_t > neighbors = get_neighbors( v,
                        cur_neighbor );
                    nb_neighbors = static_cast< index_t >( neighbors.size() );
                    for( index_t i : range( prev_neighbor, cur_neighbor ) ) {
                        if( test( neighbors[i] ) ) {
                            break;
                        }
                        result.push_back( neighbors[i] );
                    }
                } while( result.size() == cur_neighbor && result.size() < nb_points );
            }
            return result;
        }

        /*!
         * Gets the neighboring points of a given one sorted by increasing distance
         * @param[in] v the point to test
         * @param[in] nb_neighbors the number of neighbors to return
         * @return the point indices (can be less than \p nb_neighbors
         * if there is not enough points)
         */
        std::vector< index_t > get_neighbors(
            const vecn< DIMENSION >& v,
            index_t nb_neighbors ) const;

        vecn< DIMENSION > point( index_t v ) const
        {
            vecn< DIMENSION > result;
            for( index_t i : range( DIMENSION ) ) {
                result[i] = nn_points_[DIMENSION * v + i];
            }
            return result;
        }

        index_t nb_points() const
        {
            return nn_tree_->nb_points();
        }

    private:
        void fill_nn_search_points(
            index_t index_in_nn,
            const vecn< DIMENSION >& center );

    private:
        /// KdTree to compute the nearest neighbor search
        GEO::NearestNeighborSearch_var nn_tree_;
        /// Array of the points (size of DIMENSIONxnumber of points)
        double* nn_points_;
        /*!
         * @brief Indicates if ann_points_ should be deleted.
         * @details No need to delete nn_points_ if it is a simple pointer
         * to the mesh vertex array.
         */
        bool delete_points_;
    };
    ALIAS_2D_AND_3D( NNSearch );

}
