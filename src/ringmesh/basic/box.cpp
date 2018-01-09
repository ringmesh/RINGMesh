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

#include <ringmesh/basic/box.h>

#include <algorithm>

/*!
 * @file Implementation of multi-dimensional Box class
 * @author Arnaud Botella
 */

namespace
{
    template < typename T >
    inline T sqr( T x )
    {
        return x * x;
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    bool Box< DIMENSION >::initialized() const
    {
        return initialized_;
    }

    template < index_t DIMENSION >
    void Box< DIMENSION >::clear()
    {
        initialized_ = false;
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& Box< DIMENSION >::min() const
    {
        return min_;
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& Box< DIMENSION >::max() const
    {
        return max_;
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > Box< DIMENSION >::center() const
    {
        return 0.5 * ( min() + max() );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > Box< DIMENSION >::diagonal() const
    {
        return max() - min();
    }

    template < index_t DIMENSION >
    void Box< DIMENSION >::add_point( const vecn< DIMENSION >& p )
    {
        if( !initialized_ )
        {
            min_ = p;
            max_ = p;
            initialized_ = true;
        }
        else
        {
            for( auto i : range( DIMENSION ) )
            {
                min_[i] = std::min( min_[i], p[i] );
                max_[i] = std::max( max_[i], p[i] );
            }
        }
    }

    template < index_t DIMENSION >
    void Box< DIMENSION >::add_box( const Box< DIMENSION >& b )
    {
        if( b.initialized() )
        {
            add_point( b.min() );
            add_point( b.max() );
        }
    }

    template < index_t DIMENSION >
    bool Box< DIMENSION >::bboxes_overlap( const Box< DIMENSION >& B ) const
    {
        for( auto c : range( DIMENSION ) )
        {
            if( max()[c] < B.min()[c] )
            {
                return false;
            }
            if( min()[c] > B.max()[c] )
            {
                return false;
            }
        }
        return true;
    }

    template < index_t DIMENSION >
    Box< DIMENSION > Box< DIMENSION >::bbox_union(
        const Box< DIMENSION >& B ) const
    {
        Box< DIMENSION > result{ *this };
        result.add_box( B );
        return result;
    }

    template < index_t DIMENSION >
    std::tuple< bool, Box< DIMENSION > > Box< DIMENSION >::bbox_intersection(
        const Box< DIMENSION >& B ) const
    {
        if( !bboxes_overlap( B ) )
        {
            return std::make_tuple( false, Box() );
        }

        Box< DIMENSION > result;
        vecn< DIMENSION > minimal_max;
        vecn< DIMENSION > maximal_min;
        for( auto c : range( DIMENSION ) )
        {
            minimal_max[c] = std::min( this->max()[c], B.max()[c] );
            maximal_min[c] = std::max( this->min()[c], B.min()[c] );
        }
        result.add_point( maximal_min );
        result.add_point( minimal_max );
        return std::make_tuple( true, result );
    }

    template < index_t DIMENSION >
    bool Box< DIMENSION >::contains( const vecn< DIMENSION >& b ) const
    {
        for( auto c : range( DIMENSION ) )
        {
            if( b[c] < min()[c] )
            {
                return false;
            }
            if( b[c] > max()[c] )
            {
                return false;
            }
        }
        return true;
    }

    template < index_t DIMENSION >
    double Box< DIMENSION >::signed_distance( const vecn< DIMENSION >& p ) const
    {
        bool inside{ true };
        double result{ 0.0 };
        for( auto c : range( DIMENSION ) )
        {
            if( p[c] < min()[c] )
            {
                inside = false;
                result += sqr( p[c] - min()[c] );
            }
            else if( p[c] > max()[c] )
            {
                inside = false;
                result += sqr( p[c] - max()[c] );
            }
        }
        if( inside )
        {
            result = sqr( p[0] - min()[0] );
            result = std::min( result, sqr( p[0] - max()[0] ) );
            for( auto c : range( 1, DIMENSION ) )
            {
                result = std::min( result, sqr( p[c] - min()[c] ) );
                result = std::min( result, sqr( p[c] - max()[c] ) );
            }
            result = -result;
        }
        return result;
    }

    template < index_t DIMENSION >
    double Box< DIMENSION >::distance_to_center(
        const vecn< DIMENSION >& p ) const
    {
        double result{ 0.0 };
        for( auto c : range( DIMENSION ) )
        {
            double d = p[c] - 0.5 * ( min()[c] + max()[c] );
            result += sqr( d );
        }
        return result;
    }

    template class basic_api Box< 2 >;
    template class basic_api Box< 3 >;

} // namespace RINGMesh
