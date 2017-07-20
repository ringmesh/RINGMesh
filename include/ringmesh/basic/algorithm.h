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

#include <algorithm>
#include <vector>

/*!
 * @file ringmesh/algorithm.h
 * @brief Template function for basic operations on std container
 * @author Jeanne Pellerin and Arnaud Botella
 * @todo Rename these functions
 */

namespace RINGMesh {

    /*!
     * @brief Returns the position of the first entity matching @param value
     * in the container, NO_ID if not found.
     */
    template< typename T, typename container >
    inline index_t find( const container& in, const T& value )
    {
        auto it = std::find( in.begin(), in.end(), value );
        if( it == in.end() ) {
            return NO_ID;
        } else {
            return static_cast< index_t >( it - in.begin() );
        }
    }

    /*!
     * @brief Returns the position of the first entity matching @param value
     * in a sorted container, NO_ID if not found.
     */
    template< typename T, typename container >
    inline index_t find_sorted( const container& in, const T& value )
    {
        auto low = std::lower_bound( in.begin(), in.end(), value );
        if( low == in.end() || value < *low ) {
            return NO_ID;
        } else {
            return static_cast< index_t >( low - in.begin() );
        }
    }

    template< typename T, typename container >
    inline bool contains( const container& in, const T& value, bool sorted = false )
    {
        if( sorted ) {
            return find_sorted( in, value ) != NO_ID;
        } else {
            return find( in, value ) != NO_ID;
        }
    }

    /*!
     * @brief Bubble sorting of input and output vectors according to values of input.
     * @note Not efficient.
     */
    template< typename T1, typename T2 >
    inline void indirect_sort( std::vector< T1 >& input, std::vector< T2 >& output )
    {
        if( input.size() < 2 ) {
            return;
        }
        for( index_t it1 : range( input.size() - 1) ) {
            index_t ref_index = it1;
            T1 ref_value = input[it1];
            for( index_t it2 : range( it1 + 1, input.size() ) ) {
                index_t new_index = it2;
                T1 new_value = input[it2];
                if( ref_value > new_value ) {
                    ref_value = new_value;
                    ref_index = new_index;
                }
            }
            std::iter_swap( input.begin() + it1, input.begin() + ref_index );
            std::iter_swap( output.begin() + it1, output.begin() + ref_index );
        }
    }

    /*!
     * @brief Comparator of indices relying on values token in a vector
     * @note To be used in unique_values function
     */
    template< typename T >
    class CompareIndexFromValue {
    public:
        CompareIndexFromValue( const std::vector< T >& values )
            : values_( values )
        {
        }
        /*! @brief Compare two indices based on the stored values at these indices
         */
        inline bool operator()( index_t i, index_t j ) const
        {
            ringmesh_assert( i < values_.size() );
            ringmesh_assert( j < values_.size() );
            if( equal_values( i, j ) ) {
                return i < j;
            } else {
                return values_[i] < values_[j];
            }
        }
        inline bool equal_values( index_t i, index_t j ) const
        {
            return values_[i] == values_[j];
        }
    private:
        const std::vector< T >& values_;
    };

    /*!
     * @brief Determines unique occurences of values in input vector
     * and fill a mapping to the first occurence of each unique value
     * @param[in] input values  example: 1 4 6 0 4 1 5 6
     * @param[out] unique_value_indices example: 0 1 2 3 1 0 6 2
     * Maps each index to the index of the first occurence of the value in the vector.
     * @return Number of unique values in input
     * @note Tricky algorithm, used over and over in Geogram
     */
    template< typename T >
    inline index_t determine_unique_values_indices(
        const std::vector< T >& input_values,
        std::vector< index_t >& unique_value_indices )
    {
        index_t nb_values = static_cast< index_t >( input_values.size() );
        std::vector< index_t > sorted_indices( nb_values );
        std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
        CompareIndexFromValue< T > comparator( input_values );
        // Sort the indices according to the values token in the input vector
        std::sort( sorted_indices.begin(), sorted_indices.end(), comparator );

        unique_value_indices.resize( nb_values, NO_ID );
        index_t nb_unique_values = 0;
        index_t i = 0;
        while( i != nb_values ) {
            nb_unique_values++;
            unique_value_indices[sorted_indices[i]] = sorted_indices[i];
            index_t j = i + 1;
            while( j < nb_values
                && comparator.equal_values( sorted_indices[i], sorted_indices[j] ) ) {
                unique_value_indices[sorted_indices[j]] = sorted_indices[i];
                j++;
            }
            i = j;
        }
        return nb_unique_values;
    }

    /*!
     * @brief Determines unique occurences of input values plus a mapping.
     * Example:
     * Input  : input = 1 3 25 8 3 8
     * Output : unique_values = 1 3 25 8
     *          input2unique_values = 0 1 2 3 1 3
     */
    template< typename T >
    inline void get_unique_input_values_and_mapping(
        const std::vector< T >& input_values,
        std::vector< T >& unique_values,
        std::vector< index_t >& input2unique_values )
    {
        unique_values.resize( 0 );
        index_t nb_values = static_cast< index_t >( input_values.size() );
        input2unique_values.resize( nb_values, NO_ID );

        std::vector< index_t > unique_value_indices;
        index_t nb_unique_values = determine_unique_values_indices( input_values,
            unique_value_indices );
        unique_values.reserve( nb_unique_values );

        for( index_t i : range( nb_values ) ) {
            if( unique_value_indices[i] == i ) {
                input2unique_values[i] =
                    static_cast< index_t >( unique_values.size() );
                unique_values.push_back( input_values[i] );
            } else {
                input2unique_values[i] =
                    input2unique_values[unique_value_indices[i]];
            }
        }
    }

    /*!
     * @brief Sorts a container and suppresses all duplicated entities.
     * @param[in,out] container the container to sort
     * @param[in] cmp a comparator function
     */
    template< typename CONTAINER, typename CMP >
    inline void sort_unique( CONTAINER& container, const CMP& cmp )
    {
        std::sort( container.begin(), container.end(), cmp );
        container.erase( std::unique( container.begin(), container.end(), cmp ),
            container.end() );
    }

    /*!
     * @brief Sorts a container and suppresses all duplicated entities.
     * @param[in,out] container the container to sort
     */
    template< typename CONTAINER >
    inline void sort_unique( CONTAINER& container )
    {
        std::sort( container.begin(), container.end() );
        container.erase( std::unique( container.begin(), container.end() ),
            container.end() );
    }
}
