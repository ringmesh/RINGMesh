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

#include <common.h>

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

}

#endif

