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

#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>

/*!
 * @file Re-definitions of basic types similar to those of Geogram
 */
namespace RINGMesh
{
    // Basic types used in RINGMesh
    // Using definitions of Geogram/basic/numeric.h
    using GEO::Numeric::float32;
    using GEO::Numeric::float64;

    using GEO::Numeric::max_float32;
    using GEO::Numeric::max_float64;
    using GEO::Numeric::min_float32;
    using GEO::Numeric::min_float64;

    static const double global_epsilon = 1E-8;
    static const double global_epsilon_sq = global_epsilon * global_epsilon;
    static const double global_epsilon_3 = global_epsilon_sq * global_epsilon;

    // This is an unsigned int
    using GEO::index_t;
    // This is an int
    using GEO::signed_index_t;
    // This is an array template of doubles
    template < index_t DIMENSION >
    using vecn = GEO::vecng< DIMENSION, double >;
    // This is an array of 3 doubles
    using vec3 = vecn< 3 >;
    // This is an array of 2 doubles
    using vec2 = vecn< 2 >;
    // This is an array template of unsigned integers
    template < index_t DIMENSION >
    using ivecn = GEO::vecng< DIMENSION, index_t >;
    // This is an array of 3 unsigned integers
    using ivec3 = ivecn< 3 >;
    // This is an array of 2 unsigned integers
    using ivec2 = ivecn< 2 >;
    // This is an array template of signed integers
    template < index_t DIMENSION >
    using sivecn = GEO::vecng< DIMENSION, signed_index_t >;
    // This is an array of 3 signed integers
    using sivec3 = sivecn< 3 >;
    // This is an array of 2 signed integers
    using sivec2 = sivecn< 2 >;

    template < index_t DIMENSION >
    bool operator==( const vecn< DIMENSION >& u, const vecn< DIMENSION >& v )
    {
        for( index_t i = 0; i < DIMENSION; i++ )
        {
            if( u[i] != v[i] )
            {
                return false;
            }
        }
        return true;
    }

    template < index_t DIMENSION >
    bool operator==(
        const sivecn< DIMENSION >& u, const sivecn< DIMENSION >& v )
    {
        for( index_t i = 0; i < DIMENSION; i++ )
        {
            if( u[i] != v[i] )
            {
                return false;
            }
        }
        return true;
    }

    template < index_t DIMENSION >
    bool operator!=( const vecn< DIMENSION >& u, const vecn< DIMENSION >& v )
    {
        return !( u == v );
    }

    template < index_t DIMENSION >
    bool operator!=(
        const sivecn< DIMENSION >& u, const sivecn< DIMENSION >& v )
    {
        return !( u == v );
    }

    template < index_t DIMENSION >
    bool inexact_equal( const vecn< DIMENSION >& v1,
        const vecn< DIMENSION >& v2,
        double epsilon )
    {
        double square_length{ 0. };
        for( index_t i = 0; i < DIMENSION; i++ )
        {
            square_length += ( v2[i] - v1[i] ) * ( v2[i] - v1[i] );
        }
        return square_length < epsilon * epsilon;
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > basic_api initialize_vecn_coordinates( double value );

    // This is the value used in RINGMesh for a invalid index
    static const index_t NO_ID = index_t( -1 );

    using MeshType = std::string;

    /*! enum defining the type of cell in region
     *  * CellType::UNCLASSIFIED may be either a connector or more complex cell
     * that is not specified.
     *  * CellType::UNDEFINED means that the cell is not defined and cannot be
     * used.
     */
    enum struct CellType : index_t
    {
        TETRAHEDRON = 0,
        HEXAHEDRON = 1,
        PRISM = 2,
        PYRAMID = 3,
        UNCLASSIFIED = 4,
        UNDEFINED = 5
    };

    /*! enum defining the type of polygon in surface.
     *  * UNCLASSIFIED_POLYGON may be either a connector or more complex polygon
     * that is not specified.
     *  * UNDEFINED_POLYGON means that the polygon is not defined and cannot be
     * used.
     */
    enum struct PolygonType : index_t
    {
        TRIANGLE = 0,
        QUAD = 1,
        UNCLASSIFIED = 2,
        UNDEFINED = 3
    };

    template < typename Enum >
    auto to_underlying_type( Enum e ) ->
        typename std::underlying_type< Enum >::type
    {
        return static_cast< typename std::underlying_type< Enum >::type >( e );
    }

    template < typename Enum >
    struct EnableBitMaskOperators
    {
        static const bool enable = false;
    };

    template < typename Enum >
    typename std::enable_if< EnableBitMaskOperators< Enum >::enable,
        Enum >::type
        operator|( Enum lhs, Enum rhs )
    {
        using underlying = typename std::underlying_type< Enum >::type;
        return static_cast< Enum >( static_cast< underlying >( lhs )
                                    | static_cast< underlying >( rhs ) );
    }

    template < typename Enum >
    typename std::enable_if< EnableBitMaskOperators< Enum >::enable,
        Enum >::type
        operator&( Enum lhs, Enum rhs )
    {
        using underlying = typename std::underlying_type< Enum >::type;
        return static_cast< Enum >( static_cast< underlying >( lhs )
                                    & static_cast< underlying >( rhs ) );
    }

    template < typename Enum >
    typename std::enable_if< EnableBitMaskOperators< Enum >::enable,
        Enum >::type
        operator^( Enum lhs, Enum rhs )
    {
        using underlying = typename std::underlying_type< Enum >::type;
        return static_cast< Enum >( static_cast< underlying >( lhs )
                                    ^ static_cast< underlying >( rhs ) );
    }

    template < typename Enum >
    typename std::enable_if< EnableBitMaskOperators< Enum >::enable,
        Enum >::type
        operator~( Enum lhs )
    {
        using underlying = typename std::underlying_type< Enum >::type;
        return static_cast< Enum >( ~static_cast< underlying >( lhs ) );
    }

    template < typename Enum >
    bool enum_contains( Enum lhs, Enum rhs )
    {
        return ( lhs & rhs ) != Enum::EMPTY;
    }

#define ENABLE_BITMASK_OPERATORS( Enum )                                       \
    template <>                                                                \
    struct EnableBitMaskOperators< Enum >                                      \
    {                                                                          \
        static const bool enable = true;                                       \
    }
} // namespace RINGMesh
