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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_UTILS__
#define __RINGMESH_UTILS__

#include <ringmesh/common.h>

#include <geogram/basic/geometry.h>

namespace RINGMesh {     


    /*! @brief A safer narrow casting function of type S to type T
    *  \return static_cast< T >( in )
    *  \post Check that the result can be cast back to in, if not throws an assertion.
    *  \note cf. The C++ programming language. 4th edition. p299
    */
    template< typename T, typename S >
    T narrow_cast( S in )
    {
        T r = static_cast< T >( in ) ;
        if( static_cast< S >( r ) != in ) {
            ringmesh_assert_not_reached;
        }
        return r ;
    }


    /*!
     * @todo Move this in our geometry file ? [JP]
     */
    class RINGMESH_API Box3d : public GEO::Box {
    public:
        Box3d()
            : initialized_( false )
        {
        }

        bool initialized() const
        {
            return initialized_ ;
        }

        void clear()
        {
            initialized_ = false ;
        }

        float64 width() const
        {
            return xyz_max[0] - xyz_min[0] ;
        }

        float64 height() const
        {
            return xyz_max[1] - xyz_min[1] ;
        }

        float64 depth() const
        {
            return xyz_max[2] - xyz_min[2] ;
        }

        vec3 min() const
        {
            return vec3( xyz_min[0], xyz_min[1], xyz_min[2] ) ;
        }

        vec3 max() const
        {
            return vec3( xyz_max[0], xyz_max[1], xyz_max[2] ) ;
        }

        vec3 center() const
        {
            return 0.5 * ( min() + max() ) ;
        }

        void add_point( const float64* p )
        {
            add_point( vec3( p ) ) ;
        }

        void add_point( const vec3& p ) ;

        void add_box( const Box3d& b )
        {
            if( b.initialized() ) {
                add_point( b.min() ) ;
                add_point( b.max() ) ;
            }
        }

        inline bool bboxes_overlap( const Box3d& B ) const
        {
            vec3 minimum = min() ;
            vec3 b_minimum = B.min() ;
            vec3 maximum = max() ;
            vec3 b_maximum = B.max() ;
            for( index_t c = 0; c < 3; ++c ) {
                if( maximum[c] < b_minimum[c] ) {
                    return false ;
                }
                if( minimum[c] > b_maximum[c] ) {
                    return false ;
                }
            }
            return true ;
        }

        inline Box3d bbox_union( const Box3d& B ) const
        {
            Box3d result = *this ;
            result.add_box( B ) ;
            return result ;
        }

        bool contains( const vec3& b ) const
        {
            vec3 minimum = min() ;
            vec3 maximum = max() ;
            for( index_t c = 0; c < 3; ++c ) {
                if( b[c] < minimum[c] ) {
                    return false ;
                }
                if( b[c] > maximum[c] ) {
                    return false ;
                }
            }
            return true ;
        }

        float64 distance_to_center( const vec3& p ) const ;

        float64 signed_distance( const vec3& p ) const ;

    private:
        bool initialized_ ;
    } ;  

}

#endif
