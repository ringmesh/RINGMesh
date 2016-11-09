/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/basic/box3d.h>

#include <algorithm> 

#include <geogram/basic/string.h>

/*!
 * @file Implementation of Box3D class
 * @author Arnaud Botella
 * @todo Rename this file.
 */

namespace RINGMesh {

    template< typename T >
    inline T sqr( T x )
    {
        return x * x ;
    }

    void Box3d::add_point( const vec3& p )
    {
        if( !initialized_ ) {
            min_ = p ;
            max_ = p ;
            initialized_ = true ;
        } else {
            for( index_t i = 0; i < 3; i++ ) {
                min_[i] = std::min( min_[i], p[i] ) ;
                max_[i] = std::max( max_[i], p[i] ) ;
            }
        }
    }

    double Box3d::signed_distance( const vec3& p ) const
    {
        bool inside = true ;
        double result = 0.0 ;
        vec3 minimum = min() ;
        vec3 maximum = max() ;
        for( index_t c = 0; c < 3; c++ ) {
            if( p[c] < minimum[c] ) {
                inside = false ;
                result += sqr( p[c] - minimum[c] ) ;
            } else if( p[c] > maximum[c] ) {
                inside = false ;
                result += sqr( p[c] - maximum[c] ) ;
            }
        }
        if( inside ) {
            result = sqr( p[0] - minimum[0] ) ;
            result = std::min( result, sqr( p[0] - maximum[0] ) ) ;
            for( index_t c = 1; c < 3; ++c ) {
                result = std::min( result, sqr( p[c] - minimum[c] ) ) ;
                result = std::min( result, sqr( p[c] - maximum[c] ) ) ;
            }
            result = -result ;
        }
        return result ;
    }

    double Box3d::distance_to_center( const vec3& p ) const
    {
        double result = 0.0 ;
        vec3 minimum = min() ;
        vec3 maximum = max() ;
        for( index_t c = 0; c < 3; ++c ) {
            double d = p[c] - 0.5 * ( minimum[c] + maximum[c] ) ;
            result += sqr( d ) ;
        }
        return result ;
    }

}

