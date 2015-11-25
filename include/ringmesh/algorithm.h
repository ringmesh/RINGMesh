/*
* Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*  Contacts:
*     Arnaud.Botella@univ-lorraine.fr
*     Antoine.Mazuyer@univ-lorraine.fr
*     Jeanne.Pellerin@wias-berlin.de
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/

#ifndef __RINGMESH_ALGORITHM__
#define __RINGMESH_ALGORITHM__

#include <vector>
#include <algorithm>

/*!
 * @file ringmesh/algorithm.h
 * @brief Template function for basic operations on vectors
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {

    /*!
     * @brief Returns the position of the first element matching @param t 
     * in the vector, NO_ID if not found. 
     */
    template< typename T, typename container >
    index_t find( const container& v, const T& t )
    {
        typename container::const_iterator it = std::find( v.begin(), v.end(), t ) ;
        if( it == v.end() ) {
            return NO_ID ;
        } else {
            return static_cast<index_t>(it - v.begin()) ;
        }
    }

    /*!
     * @brief Returns the position of the first element matching t 
     * in a sorted vector, NO_ID if not found. 
     */
    template< typename T, typename container >
    index_t find_sorted( const container& v, const T& t )
    {
        typename container::const_iterator low = std::lower_bound( v.begin(), v.end(), t ) ;
        if( low == v.end() || t < *low ) {
            return NO_ID ;
        } else {
            return static_cast<index_t>(low - v.begin()) ;
        }
    }


    template< typename T, typename container >
    bool contains( const container& v, const T& t, bool sorted = false )
    {
        if( sorted ) {
            return find_sorted( v, t ) != NO_ID ;
        } else {
            return find( v, t ) != NO_ID ;
        }
    }


    /*!
     * \brief Indirect sorting of two vectors.
     * @todo Comment what is indirect sorting.
     */
    template< class T1, class T2 >
    void indirect_sort( std::vector< T1 >& input, std::vector< T2 >& output )
    {
        if( input.size() < 2 ) {
            return ;
        }
        for( index_t it1 = 0; it1+1 < input.size(); it1++ ) {
            index_t ref_index = it1 ;
            T1 ref_value = input[ it1 ] ;
            for( index_t it2 = it1 + 1; it2 < input.size(); it2++ ) {
                index_t new_index = it2 ;
                T1 new_value = input[ it2 ] ;
                if( ref_value > new_value ) {
                    ref_value = new_value ;
                    ref_index = new_index ;
                }
            }
            std::iter_swap( input.begin() + it1, input.begin() + ref_index ) ;
            std::iter_swap( output.begin() + it1, output.begin() + ref_index ) ;
        }
    }

    /*! 
     * @brief Comparator of indices relying on values token in a vector
     * @note To be used in unique_values function
     */
    template< class T >
    class CompareIndexFromValue
    {
    public:
        CompareIndexFromValue( const std::vector<T>& values ) :
            values_( values )
        {}
        /*! @brief Compare two indices based on the stored values at these indices
         */
        bool operator()( index_t i, index_t j )
        {
            if( are_values_equal( i, j ) ) {
                return i < j ;
            } else {
                return values_[ i ] < values_[ j ] ;
            }
        }
        bool are_values_equal( index_t i, index_t j )
        {
            return values_[ i ] == values_[ j ] ;
        }
    private:
        const std::vector<T>& values_ ;
    };

    /*!
     * @brief Determine unique occurences of values in input vector
     * and fill a mapping to the first occurence of each unique value
     * @param[in] input values  example: 1 4 6 0 4 1 5 6
     * @param[out] input2unique example: 0 1 2 3 1 0 6 2
     * Maps each index to the index of the first occurence of the value in the vector. 
     * @return Number of unique values in input
     * @note Tricky algorithm, used over and over in Geogram
     */
    template< class T >
    index_t unique_values(
        const std::vector< T >& input,
        std::vector< index_t >& unique_value_index )
    {
        index_t nb_values = input.size() ;
        std::vector< index_t > sorted_indices( nb_values ) ;
        for( index_t i = 0; i < nb_values; ++i ) {
            sorted_indices[ i ] = i ;
        }
        CompareIndexFromValue<T> comparator( input ) ;
        // Sort the indices according to the values token in the input vector
        std::sort( sorted_indices.begin(), sorted_indices.end(), comparator ) ;

        unique_value_index.resize( nb_values, NO_ID ) ;
        index_t nb_unique_values = 0 ;
        index_t i = 0 ;
        while( i != nb_values ) {
            nb_unique_values++ ;
            unique_value_index[ sorted_indices[ i ] ] = sorted_indices[ i ] ;
            index_t j = i + 1 ;
            while( j < input.size() &&
                   comparator.are_values_equal( sorted_indices[ i ], sorted_indices[ j ] )
            ) {
                unique_value_index[ sorted_indices[ j ] ] = sorted_indices[ i ] ;
                j++ ;
            }
            i = j ;
        }
        return nb_unique_values ;
    }
}

#endif
