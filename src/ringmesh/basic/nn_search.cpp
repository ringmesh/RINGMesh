/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

/*!
 * @file Nearest Neighbor requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

#include <ringmesh/basic/nn_search.h>
#include <ringmesh/basic/pimpl_impl.h>
#include <ringmesh/basic/task_handler.h>

#include <geogram/points/kd_tree.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    class NNSearch< DIMENSION >::Impl
    {
    public:
        Impl( const std::vector< vecn< DIMENSION > >& vertices, bool copy )
            : nn_tree_( GEO::NearestNeighborSearch::create( DIMENSION, "BNN" ) )
        {
            auto nb_vertices = static_cast< index_t >( vertices.size() );
            if( copy )
            {
                nn_points_ = new double[nb_vertices * DIMENSION];
                delete_points_ = true;
                GEO::Memory::copy( nn_points_, vertices.data()->data(),
                    DIMENSION * nb_vertices * sizeof( double ) );
            }
            else
            {
                nn_points_ = const_cast< double* >( vertices.data()->data() );
                delete_points_ = false;
            }
            nn_tree_->set_points( nb_vertices, nn_points_ );
        }

        ~Impl()
        {
            if( delete_points_ )
            {
                delete[] nn_points_;
            }
        }

        void fill_nn_search_points(
            index_t index_in_nn_search, const vecn< DIMENSION >& center )
        {
            for( auto i : range( DIMENSION ) )
            {
                nn_points_[index_in_nn_search + i] = center[i];
            }
        }

        vecn< DIMENSION > point( index_t v ) const
        {
            vecn< DIMENSION > result;
            for( auto i : range( DIMENSION ) )
            {
                result[i] = nn_points_[DIMENSION * v + i];
            }
            return result;
        }

        index_t nb_points() const
        {
            return nn_tree_->nb_points();
        }

        std::vector< index_t > get_neighbors(
            const vecn< DIMENSION >& v, index_t nb_neighbors ) const
        {
            std::vector< index_t > result;
            if( nb_points() != 0 )
            {
                nb_neighbors = std::min( nb_neighbors, nb_points() );
                std::vector< double > distances( nb_neighbors );
                result.resize( nb_neighbors );
                nn_tree_->get_nearest_neighbors(
                    nb_neighbors, v.data(), &result[0], &distances[0] );
            }
            return result;
        }

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

    template < index_t DIMENSION >
    NNSearch< DIMENSION >::NNSearch(
        const std::vector< vecn< DIMENSION > >& vertices, bool copy )
        : impl_( vertices, copy )
    {
    }

    template < index_t DIMENSION >
    NNSearch< DIMENSION >::~NNSearch()
    {
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > NNSearch< DIMENSION >::point( index_t v ) const
    {
        return impl_->point( v );
    }

    template < index_t DIMENSION >
    index_t NNSearch< DIMENSION >::nb_points() const
    {
        return impl_->nb_points();
    }

    template < index_t DIMENSION >
    std::tuple< index_t, std::vector< index_t > >
        NNSearch< DIMENSION >::get_colocated_index_mapping(
            double epsilon ) const
    {
        std::vector< index_t > index_map( nb_points() );
        std::atomic< index_t > nb_colocalised_vertices{ 0 };
        parallel_for( nb_points(), [this, &index_map, &nb_colocalised_vertices,
                                       &epsilon]( index_t i ) {
            auto results = get_neighbors( point( i ), epsilon );
            index_t id{ *std::min_element( results.begin(), results.end() ) };
            index_map[i] = id;
            if( id < i )
            {
                nb_colocalised_vertices++;
            }
        } );
        return std::make_tuple( nb_colocalised_vertices.load(), index_map );
    }

    template < index_t DIMENSION >
    std::tuple< index_t,
        std::vector< index_t >,
        std::vector< vecn< DIMENSION > > >
        NNSearch< DIMENSION >::get_colocated_index_mapping_and_unique_points(
            double epsilon ) const
    {
        index_t nb_colocalised_vertices;
        std::vector< index_t > index_map;
        std::tie( nb_colocalised_vertices, index_map ) =
            get_colocated_index_mapping( epsilon );
        std::vector< vecn< DIMENSION > > unique_points;
        unique_points.reserve( nb_points() - nb_colocalised_vertices );
        index_t offset{ 0 };
        for( auto p : range( index_map.size() ) )
        {
            if( index_map[p] == p )
            {
                unique_points.emplace_back( point( p ) );
                index_map[p] = p - offset;
            }
            else
            {
                offset++;
                index_map[p] = index_map[index_map[p]];
            }
        }
        ringmesh_assert( offset == nb_colocalised_vertices );
        return std::make_tuple( offset, index_map, unique_points );
    }

    template < index_t DIMENSION >
    std::vector< index_t > NNSearch< DIMENSION >::get_neighbors(
        const vecn< DIMENSION >& v, double threshold_distance ) const
    {
        double threshold_distance_sq{ threshold_distance * threshold_distance };
        return get_neighbors(
            v, [this, &v, threshold_distance_sq]( index_t i ) {
                return length2( v - point( i ) ) > threshold_distance_sq;
            } );
    }

    template < index_t DIMENSION >
    std::vector< index_t > NNSearch< DIMENSION >::get_neighbors(
        const vecn< DIMENSION >& v, index_t nb_neighbors ) const
    {
        return impl_->get_neighbors( v, nb_neighbors );
    }

    template class basic_api NNSearch< 2 >;
    template class basic_api NNSearch< 3 >;
} // namespace RINGMesh
