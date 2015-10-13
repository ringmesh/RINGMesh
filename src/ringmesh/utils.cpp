/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the RING Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <ringmesh/utils.h>

#include <algorithm> 

namespace RINGMesh {

    template< class T >
    inline T sqr( T x )
    {
        return x * x ;
    }
   
    void Box3d::add_point( const vec3& p )
    {
        if( !initialized_ ) {
            for( index_t i = 0; i < 3; i++ ) {
                xyz_min[ i ] = p[ i ] ;
                xyz_max[ i ] = p[ i ] ;
            }
            initialized_ = true ;
        } else {
            for( index_t i = 0; i < 3; i++ ) {
                xyz_min[ i ] = std::min( xyz_min[ i ], p[ i ] ) ;
                xyz_max[ i ] = std::max( xyz_max[ i ], p[ i ] ) ;
            }
        }
    }

    float64 Box3d::signed_distance( const vec3& p ) const
    {
        bool inside = true ;
        float64 result = 0.0 ;
        vec3 minimum = min() ;
        vec3 maximum = max() ;
        for( index_t c = 0; c < 3; c++ ) {
            if( p[ c ] < minimum[ c ] ) {
                inside = false ;
                result += sqr( p[ c ] - minimum[ c ] ) ;
            } else if( p[ c ] > maximum[ c ] ) {
                inside = false ;
                result += sqr( p[ c ] - maximum[ c ] ) ;
            }
        }
        if( inside ) {
            result = sqr( p[ 0 ] - minimum[ 0 ] ) ;
            result = std::min( result, sqr( p[ 0 ] - maximum[ 0 ] ) ) ;
            for( index_t c = 1; c < 3; ++c ) {
                result = std::min( result, sqr( p[ c ] - minimum[ c ] ) ) ;
                result = std::min( result, sqr( p[ c ] - maximum[ c ] ) ) ;
            }
            result = -result ;
        }
        return result ;
    }


    float64 Box3d::distance_to_center( const vec3& p ) const
    {
        float64 result = 0.0 ;
        vec3 minimum = min() ;
        vec3 maximum = max() ;
        for( index_t c = 0; c < 3; ++c ) {
            float64 d = p[ c ] - 0.5 * ( minimum[ c ] + maximum[ c ] ) ;
            result += sqr( d ) ;
        }
        return result ;
    }

}
