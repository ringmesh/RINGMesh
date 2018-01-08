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

#pragma once

#include <ringmesh/basic/common.h>

#include <algorithm>
#include <vector>

/*!
 * @file ringmesh/algorithm.h
 * @brief Template function for basic operations on std container
 * @author Jeanne Pellerin and Arnaud Botella
 * @todo Rename these functions
 */

namespace RINGMesh
{
    /*!
     * @brief Returns the position of the first entity matching @param value
     * in the container, NO_ID if not found.
     */
    template < typename T, typename container >
    index_t find( const container& in, const T& value )
    {
        auto it = std::find( in.begin(), in.end(), value );
        if( it == in.end() )
        {
            return NO_ID;
        }
        else
        {
            return static_cast< index_t >( it - in.begin() );
        }
    }

    /*!
     * @brief Returns the position of the first entity matching @param value
     * in a sorted container, NO_ID if not found.
     */
    template < typename T, typename container >
    index_t find_sorted( const container& in, const T& value )
    {
        auto low = std::lower_bound( in.begin(), in.end(), value );
        if( low == in.end() || value < *low )
        {
            return NO_ID;
        }
        else
        {
            return static_cast< index_t >( low - in.begin() );
        }
    }

    template < typename T, typename container >
    bool contains( const container& in, const T& value )
    {
        return find( in, value ) != NO_ID;
    }

    template < typename T, typename container >
    bool contains_sorted( const container& in, const T& value )
    {
        return find_sorted( in, value ) != NO_ID;
    }

    /*!
     * @brief Bubble sorting of input and output vectors according to values of
     * input.
     * @note Not efficient.
     */
    template < typename T1, typename T2 >
    void indirect_sort( std::vector< T1 >& input, std::vector< T2 >& output )
    {
        if( input.size() < 2 )
        {
            return;
        }
        for( auto it1 : range( input.size() - 1 ) )
        {
            index_t ref_index = it1;
            T1 ref_value = input[it1];
            for( auto it2 : range( it1 + 1, input.size() ) )
            {
                index_t new_index = it2;
                T1 new_value = input[it2];
                if( ref_value > new_value )
                {
                    ref_value = new_value;
                    ref_index = new_index;
                }
            }
            std::iter_swap( input.begin() + it1, input.begin() + ref_index );
            std::iter_swap( output.begin() + it1, output.begin() + ref_index );
        }
    }

    /*!
     * @brief Sorts a container and suppresses all duplicated entities.
     * @param[in,out] container the container to sort
     * @param[in] cmp a comparator function
     */
    template < typename CONTAINER, typename CMP >
    void sort_unique( CONTAINER& container, const CMP& cmp )
    {
        std::sort( container.begin(), container.end(), cmp );
        container.erase( std::unique( container.begin(), container.end(), cmp ),
            container.end() );
    }

    /*!
     * @brief Sorts a container and suppresses all duplicated entities.
     * @param[in,out] container the container to sort
     */
    template < typename CONTAINER >
    void sort_unique( CONTAINER& container )
    {
        std::sort( container.begin(), container.end() );
        container.erase( std::unique( container.begin(), container.end() ),
            container.end() );
    }
} // namespace RINGMesh
