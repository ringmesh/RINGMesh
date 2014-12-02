/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * float64his program is a float64rade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#ifndef __GRGMESH_VECN__
#define __GRGMESH_VECN__

#include <grgmesh/common.h>
#include <cstring>
#include <cmath>
#include <iostream>

namespace GRGMesh {

    template< int DIM=3 > class vecn {
    public:
        typedef vecn< DIM > thisclass ;

        vecn()
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] = float64( 0 ) ;
            }
        }

        template< int DIM2 > explicit vecn( const vecn< DIM2 >& rhs )
        {
            grgmesh_debug_assert( DIM == DIM2 ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] = float64( rhs[i] ) ;
            }
        }
        explicit vecn( const float64* rhs )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] = rhs[i] ;
            }
        }
        vecn( const vecn< DIM > &rhs )
        {
            memcpy( data_, rhs.data(), DIM * sizeof(float64) ) ;
        }

        thisclass& operator=( const thisclass& rhs )
        {
            memcpy( data_, rhs.data(), DIM * sizeof(float64) ) ;
            return *this ;
        }

        unsigned int dimension() const
        {
            return (unsigned int) DIM ;
        }

        float64* data()
        {
            return data_ ;
        }

        const float64* data() const
        {
            return data_ ;
        }

        inline float64& operator[]( unsigned int idx )
        {
            grgmesh_debug_assert( idx < DIM ) ;
            return data()[idx] ;
        }

        inline const float64& operator[]( unsigned int idx ) const
        {
            grgmesh_debug_assert( idx < DIM ) ;
            return data()[idx] ;
        }

        inline float64 length2() const
        {
            float64 result = float64( 0 ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result += data_[i] * data_[i] ;
            }
            return result ;
        }

        inline float64 length() const
        {
            return sqrt( length2() ) ;
        }

        inline float64 distance2( const thisclass &rhs ) const
        {
            float64 result = float64( 0 ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                float64 val = rhs.data_[i] - data_[i] ;
                result += val * val ;
            }
            return result ;
        }

        // operators
        inline thisclass& operator+=( const thisclass& v )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] += v.data_[i] ;
            }
            return *this ;
        }

        inline bool operator==( const thisclass& p ) const {
            for( unsigned int i = 0; i < DIM; i++ ) {
                if( data_[i] != p.data_[i] ) return false ;
            }
            return true ;
        }

        inline bool operator!=( const thisclass& p ) const {
            for( unsigned int i = 0; i < DIM; i++ ) {
                if( data_[i] != p.data_[i] ) return true ;
            }
            return false ;
        }

        inline thisclass& operator-=( const thisclass& v )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] -= v.data_[i] ;
            }
            return *this ;
        }

        inline thisclass& operator*=( const thisclass& v )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] *= v.data_[i] ;
            }
            return *this ;
        }

        inline thisclass& operator/=( const thisclass& v )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] /= v.data_[i] ;
            }
            return *this ;
        }

        template< class float642 > inline thisclass& operator*=( float642 s )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] *= float64( s ) ;
            }
            return *this ;
        }

        template< class float642 > inline thisclass& operator/=( float642 s )
        {
            for( unsigned int i = 0; i < DIM; i++ ) {
                data_[i] /= float64( s ) ;
            }
            return *this ;
        }

        inline thisclass operator+( const thisclass& v ) const
        {
            thisclass result( *this ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result.data_[i] += v.data_[i] ;
            }
            return result ;
        }

        inline thisclass operator-( const thisclass& v ) const
        {
            thisclass result( *this ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result.data_[i] -= v.data_[i] ;
            }
            return result ;
        }

        template< class float642 > inline thisclass operator*( float642 s ) const
        {
            thisclass result( *this ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result.data_[i] *= float64( s ) ;
            }
            return result ;
        }

        template< class float642 > inline thisclass operator/( float642 s ) const
        {
            thisclass result( *this ) ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result.data_[i] /= float64( s ) ;
            }
            return result ;
        }

        inline thisclass operator-() const
        {
            thisclass result ;
            for( unsigned int i = 0; i < DIM; i++ ) {
                result.data_[i] = -data_[i] ;
            }
            return result ;
        }

    private:
        float64 data_[DIM] ;
    } ;

    template< int DIM > inline float64 dot(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2 )
    {
        float64 result = 0 ;
        for( unsigned int i = 0; i < DIM; i++ ) {
            result += v1[i] * v2[i] ;
        }
        return result ;
    }

    template< int DIM > inline vecn< DIM > operator-(
        const vecn< DIM >& v1 )
    {
        vecn< DIM > result ;
        for( unsigned int i = 0; i < DIM; i++ ) {
            result[i] = -v1[i] ;
        }
        return result ;
    }

    template< int DIM > inline vecn< DIM > operator*(
        float64 s,
        const vecn< DIM >& v )
    {
        vecn< DIM > result ;
        for( unsigned int i = 0; i < DIM; i++ ) {
            result[i] = float64( s ) * v[i] ;
        }
        return result ;
    }

    template< int DIM > inline vecn< DIM > operator+(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2 )
    {
        vecn< DIM > result ;
        for( unsigned int i = 0; i < DIM; i++ ) {
            result[i] = v1[i] + v2[i] ;
        }
        return result ;
    }

    template< int DIM > inline vecn< DIM > operator-(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2 )
    {
        vecn< DIM > result ;
        for( unsigned int i = 0; i < DIM; i++ ) {
            result[i] = v1[i] - v2[i] ;
        }
        return result ;
    }

    // Compatibility with GLSL
    template< int DIM > inline float64 length( const vecn< DIM >& v )
    {
        return v.length() ;
    }
    template< int DIM > inline float64 length2( const vecn< DIM >& v )
    {
        return v.length2() ;
    }
    template< int DIM > inline float64 distance2(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2 )
    {
        return v2.distance2( v1 ) ;
    }
    template< int DIM > inline float64 distance(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2 )
    {
        return length( v2 - v1 ) ;
    }
    template< int DIM > inline vecn< DIM > normalize(
        const vecn< DIM >& v )
    {
        float64 s = length( v ) ;
        if( s > 1e-30 ) {
            s = float64( 1 ) / s ;
        }
        return s * v ;
    }
    template< int DIM > inline vecn< DIM > mix(
        const vecn< DIM >& v1,
        const vecn< DIM >& v2,
        float64 s )
    {
        return ( float64( 1 ) - s ) * v1 + s * v2 ;
    }

    typedef vecn< 3 > vec3 ;

    template< > class vecn< 3 > {
    public:
        typedef vecn< 3 > thisclass ;

        vecn()
            : x( 0 ), y( 0 ), z( 0 )
        {
        }
        vecn( float64 x_in, float64 y_in, float64 z_in )
            : x( x_in ), y( y_in ), z( z_in )
        {
        }
        vecn( const thisclass& v )
            : x( v.x ), y( v.y ), z( v.z )
        {
        }
        explicit vecn( const float64* v )
            : x( v[0] ), y( v[1] ), z( v[2] )
        {
        }

        inline float64 length2() const
        {
            return x * x + y * y + z * z ;
        }
        inline float64 length() const
        {
            return ::sqrt( x * x + y * y + z * z ) ;
        }
        inline float64 distance2( const thisclass& rhs ) const
        {
            float64 dx = rhs.x - x ;
            float64 dy = rhs.y - y ;
            float64 dz = rhs.z - z ;
            return dx * dx + dy * dy + dz * dz ;
        }

        // operators
        inline thisclass& operator+=( const thisclass& v )
        {
            x += v.x ;
            y += v.y ;
            z += v.z ;
            return *this ;
        }
        inline thisclass& operator-=( const thisclass& v )
        {
            x -= v.x ;
            y -= v.y ;
            z -= v.z ;
            return *this ;
        }
        inline thisclass& operator*=( const thisclass& v )
        {
            x *= v.x ;
            y *= v.y ;
            z *= v.z ;
            return *this ;
        }
        inline thisclass& operator/=( const thisclass& v )
        {
            x /= v.x ;
            y /= v.y ;
            z /= v.z ;
            return *this ;
        }
        inline thisclass& operator*=( float64 s )
        {
            x *= float64( s ) ;
            y *= float64( s ) ;
            z *= float64( s ) ;
            return *this ;
        }
        inline thisclass& operator/=( float64 s )
        {
            x /= float64( s ) ;
            y /= float64( s ) ;
            z /= float64( s ) ;
            return *this ;
        }
        inline bool operator==( const thisclass& p ) const
        {
            for( unsigned int i = 0; i < 3; i++ ) {
                if( (*this)[i] != p[i] ) return false ;
            }
            return true ;
        }

        inline bool operator!=( const thisclass& p ) const
        {
            for( unsigned int i = 0; i < 3; i++ ) {
                if( (*this)[i] != p[i] ) return true ;
            }
            return false ;
        }
        inline thisclass operator+( const thisclass& v ) const
        {
            return thisclass( x + v.x, y + v.y, z + v.z ) ;
        }
        inline thisclass operator-( const thisclass& v ) const
        {
            return thisclass( x - v.x, y - v.y, z - v.z ) ;
        }
        inline thisclass operator*( float64 s ) const
        {
            return thisclass( x * float64( s ), y * float64( s ), z * float64( s ) ) ;
        }
        inline thisclass operator/( float64 s ) const
        {
            return thisclass( x / float64( s ), y / float64( s ), z / float64( s ) ) ;
        }

        inline thisclass operator-() const
        {
            return thisclass( -x, -y, -z ) ;
        }

        unsigned int dimension() const
        {
            return (unsigned int) 3 ;
        }

        float64* data()
        {
            return &x ;
        }
        const float64* data() const
        {
            return &x ;
        }

        inline float64& operator[]( unsigned int idx )
        {
            grgmesh_debug_assert( idx < 3 ) ;
            return data()[idx] ;
        }

        inline const float64& operator[]( unsigned int idx ) const
        {
            grgmesh_debug_assert( idx < 3 ) ;
            return data()[idx] ;
        }

        float64 x ;
        float64 y ;
        float64 z ;
    } ;

    inline float64 dot( const vec3& v1, const vec3& v2 )
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z ;
    }
    inline float64 length( const vec3& v1, const vec3& v2 )
    {
        return (v2 - v1).length() ;
    }

    inline vec3 cross( const vec3& v1, const vec3& v2 )
    {
        return vec3( v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x ) ;
    }
    inline std::ostream& operator<<( std::ostream& out, const vec3& v )
    {
        return out << v.x << "  " << v.y << "  " << v.z ;
    }
    inline std::istream& operator>>( std::istream& in, vec3& v )
    {
        for( uint8 i = 0; i < 3; i++ ) {
            in >> v[i] ;
        }
        return in ;
    }


    static vec3 dummy_vec3 ;
}

#endif

