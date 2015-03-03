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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_PERMUTATION__
#define __RINGMESH_PERMUTATION__

#include <ringmesh/common.h>

namespace RINGMesh {
    namespace Permutation {
        /**
         * @brief checks whether a specified vector encodes
         *  a valid permutation.
         */
        inline bool is_valid( const std::vector< signed_index_t >& permutation )
        {
            std::vector< bool > visited( permutation.size(), false ) ;

            for( index_t i = 0; i < permutation.size(); i++ ) {
                if( permutation[ i ] < 0
                    || permutation[ i ] >= signed_index_t( permutation.size() ) )
                {
                    return false ;
                }

                if( visited[ i ] ) {
                    return false ;
                }

                visited[ i ] = true ;
            }

            return true ;
        }


        /**
         * @brief used internally by apply_permutation()
         */
        inline bool is_marked(
            const std::vector< signed_index_t >& permutation,
            signed_index_t i )
        {
            ringmesh_debug_assert( i >= 0 && i <
                signed_index_t( permutation.size() ) ) ;
            return ( permutation[ i ] < 0 ) ;
        }


        /**
         * @brief used internally by apply_permutation()
         */
        inline void mark(
            std::vector< signed_index_t >& permutation,
            signed_index_t i )
        {
            ringmesh_debug_assert( i >= 0 && i <
                signed_index_t( permutation.size() ) ) ;
            ringmesh_debug_assert( !is_marked( permutation, i ) ) ;
            permutation[ i ] = - permutation[ i ] - 1 ;
        }


        /**
         * @brief used internally by apply_permutation()
         */
        inline void unmark(
            std::vector< signed_index_t >& permutation,
            signed_index_t i )
        {
            ringmesh_debug_assert( i >= 0 && i <
                signed_index_t( permutation.size() ) ) ;
            ringmesh_debug_assert( is_marked( permutation, i ) ) ;
            permutation[ i ] = - permutation[ i ] - 1 ;
        }


        /**
         * @brief applies a permutation in-place.
         * It is equivalent to:
         * @code
         * for(i=0; i<permutation.size(); i++) {
         *    data2[i] = data[permutation[i]]
         * }
         * data = data2 ;
         * @endcode
         * @param data the vector to permute
         * @param [in] permutation_in the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         */
        inline void apply(
            void* data_in,
            std::vector< signed_index_t >& permutation_in,
            index_t elemsize )
        {
            pointer data = (pointer) ( data_in ) ;
            std::vector< signed_index_t >& permutation =
                const_cast< std::vector< signed_index_t >& >( permutation_in ) ;
            ringmesh_debug_assert( is_valid( permutation ) ) ;
            pointer temp = static_cast< pointer >( alloca( elemsize ) ) ;

            for( index_t k = 0; k < permutation.size(); k++ ) {
                if( is_marked( permutation, k ) ) {
                    continue ;
                }

                signed_index_t i = k ;
                ::memcpy( temp, data + i * elemsize, elemsize ) ;
                signed_index_t j = permutation[ k ] ;
                mark( permutation, k ) ;

                while( j != signed_index_t( k ) ) {
                    ::memcpy( data + i * elemsize, data + j * elemsize,
                        elemsize ) ;
                    signed_index_t nj = permutation[ j ] ;
                    mark( permutation, j ) ;
                    i = j ;
                    j = nj ;
                }

                ::memcpy( data + i * elemsize, temp, elemsize ) ;
            }

            for( index_t k = 0; k < permutation.size(); k++ ) {
                unmark( permutation, k ) ;
            }
        }


        /**
         * @brief applies a permutation in-place.
         * It is equivalent to:
         * @code
         * for(i=0; i<permutation.size(); i++) {
         *    data2[i] = data[permutation[i]]
         * }
         * data = data2 ;
         * @endcode
         * @param data the vector to permute
         * @param [in] permutation_in the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         */
        template< class T > inline void apply(
            std::vector< T >& data,
            std::vector< signed_index_t >& permutation_in )
        {
            std::vector< signed_index_t >& permutation =
                const_cast< std::vector< signed_index_t >& >( permutation_in ) ;
            ringmesh_debug_assert( is_valid( permutation ) ) ;
            T temp ;

            for( index_t k = 0; k < permutation.size(); k++ ) {
                if( is_marked( permutation, k ) ) {
                    continue ;
                }

                signed_index_t i = k ;
                temp = data[ i ] ;
                signed_index_t j = permutation[ k ] ;
                mark( permutation, k ) ;

                while( j != signed_index_t( k ) ) {
                    data[ i ] = data[ j ] ;
                    signed_index_t nj = permutation[ j ] ;
                    mark( permutation, j ) ;
                    i = j ;
                    j = nj ;
                }

                data[ i ] = temp ;
            }

            for( index_t k = 0; k < permutation.size(); k++ ) {
                unmark( permutation, k ) ;
            }
        }


        /**
         * @brief inverts a permutation in-place.
         * It is equivalent to:
         * @code
         * for(index_t i=0; i<permutation.size()(); i++) {
         *      perminv[permutation[i]] = i ;
         * }
         * permutation = perminv ;
         * @endcode
         */
        inline void invert( std::vector< signed_index_t >& permutation )
        {
            ringmesh_debug_assert( is_valid( permutation ) ) ;

            for( index_t k = 0; k < permutation.size(); k++ ) {
                if( is_marked( permutation, k ) ) {
                    continue ;
                }

                signed_index_t i = k ;
                signed_index_t j = permutation[ i ] ;

                while( j != signed_index_t( k ) ) {
                    signed_index_t temp = permutation[ j ] ;
                    permutation[ j ] = i ;
                    mark( permutation, j ) ;
                    i = j ;
                    j = temp ;
                }

                permutation[ j ] = i ;
                mark( permutation, j ) ;
            }

            for( index_t k = 0; k < permutation.size(); k++ ) {
                unmark( permutation, k ) ;
            }
        }
    }
}

#endif
