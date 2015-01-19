/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#ifndef __GRGMESH_UTILS__
#define __GRGMESH_UTILS__

#include <grgmesh/common.h>
#include <grgmesh/types.h>

#include <geogram/points/nn_search.h>
#include <geogram/points/kd_tree.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

namespace GEO {
    class Mesh ;
}

namespace GRGMesh {
    class BoundaryModelElement ;
    class Edge ;
    class Surface ;
    class Line ;
}

namespace GRGMesh {

    inline vec3 mesh_cell_center(const GEO::Mesh& M, index_t cell) ;

    static bool operator==( const vec3& u, const vec3& v ) {
        return u.x == v.x && u.y == v.y && u.z == v.z ;
    }

    static bool operator!=( const vec3& u, const vec3& v ) {
        return u.x != v.x || u.y != v.y || u.z != v.z ;
    }

    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    } ;

    template< class T > inline T sqr( T x )
    {
        return x * x ;
    }
    template< class T > inline Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO ) ;
    }
    class Box3d {
    public:
        Box3d()
            :
                initialized_( false ),
                x_min_( 1e30 ),
                y_min_( 1e30 ),
                z_min_( 1e30 ),
                x_max_( -1e30 ),
                y_max_( -1e30 ),
                z_max_( -1e30 )
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
        float64 x_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return x_min_ ;
        }
        float64 y_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return y_min_ ;
        }
        float64 z_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return z_min_ ;
        }
        float64 x_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return x_max_ ;
        }
        float64 y_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return y_max_ ;
        }
        float64 z_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return z_max_ ;
        }
        float64 min( unsigned axis ) const
        {
            return ( axis == 0 ) ? x_min_ : ( ( axis == 1 ) ? y_min_ : z_min_ ) ;
        }
        float64 max( unsigned axis ) const
        {
            return ( axis == 0 ) ? x_max_ : ( ( axis == 1 ) ? y_max_ : z_max_ ) ;
        }
        float64 width() const
        {
            return x_max() - x_min() ;
        }
        float64 height() const
        {
            return y_max() - y_min() ;
        }
        float64 depth() const
        {
            return z_max() - z_min() ;
        }
        vec3 min() const
        {
            return vec3( x_min(), y_min(), z_min() ) ;
        }
        vec3 max() const
        {
            return vec3( x_max(), y_max(), z_max() ) ;
        }
        vec3 center() const
        {
            return vec3( 0.5 * ( x_max() + x_min() ), 0.5 * ( y_max() + y_min() ),
                0.5 * ( z_max() + z_min() ) ) ;
        }
        float64 radius() const
        {
            return 0.5
                * ::sqrt(
                    ::sqrt( x_max() - x_min() ) + ::sqrt( y_max() - y_min() )
                        + ::sqrt( z_max() - z_min() ) ) ;
        }
        void add_point( const float64* p ) {
            add_point( vec3( p ) ) ;
        }
        void add_point( const vec3& p )
        {
            if( !initialized_ ) {
                x_min_ = p[0] ;
                y_min_ = p[1] ;
                z_min_ = p[2] ;
                x_max_ = p[0] ;
                y_max_ = p[1] ;
                z_max_ = p[2] ;
                initialized_ = true ;
            } else {
                x_min_ = std::min( x_min_, p[0] ) ;
                y_min_ = std::min( y_min_, p[1] ) ;
                z_min_ = std::min( z_min_, p[2] ) ;
                x_max_ = std::max( x_max_, p[0] ) ;
                y_max_ = std::max( y_max_, p[1] ) ;
                z_max_ = std::max( z_max_, p[2] ) ;
            }
        }
        void add_box( const Box3d& b )
        {
            if( b.initialized() ) {
                add_point( vec3( b.x_min(), b.y_min(), b.z_min() ) ) ;
                add_point( vec3( b.x_max(), b.y_max(), b.z_max() ) ) ;
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
        inline Box3d bbox_union( const Box3d& B )
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
        float64 distance_to_center( const vec3& p ) const
        {
            float64 result = 0.0 ;
            vec3 minimum = min() ;
            vec3 maximum = max() ;
            for( index_t c = 0; c < 3; ++c ) {
                float64 d = p[c] - 0.5 * ( minimum[c] + maximum[c] ) ;
                result += sqr( d ) ;
            }
            return result ;
        }
        float64 signed_distance( const vec3& p ) const
        {
            bool inside = true ;
            float64 result = 0.0 ;
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

    private:
        bool initialized_ ;
        float64 x_min_ ;
        float64 y_min_ ;
        float64 z_min_ ;
        float64 x_max_ ;
        float64 y_max_ ;
        float64 z_max_ ;
    } ;

    class GRGMESH_API Utils {
    public:
        static bool has_edge(
            GEO::Mesh& mesh,
            index_t t,
            index_t p0,
            index_t p1,
            index_t& edge ) ;
        static signed_index_t next_arround_edge(
            GEO::Mesh& mesh,
            index_t t,
            index_t prev,
            index_t p0,
            index_t p1 ) ;
        static void edges_arround_edge(
            GEO::Mesh& mesh,
            index_t t,
            index_t p0,
            index_t p1,
            std::vector< index_t >& result ) ;
        static index_t get_nearest_vertex_index(
            const GEO::Mesh& mesh,
            const vec3& p,
            signed_index_t t ) ;
        static bool facets_have_same_orientation(
            const GEO::Mesh& mesh,
            index_t f1,
            index_t c11,
            index_t f2 ) ;
        static void check_and_repair_mesh_consistency(
            const BoundaryModelElement& region,
            GEO::Mesh& mesh,
            bool check_duplicated_facet = false ) ;
        static float64 triangle_area( const vec3& p1, const vec3& p2, const vec3& p3 )
        {
            return 0.5 * length( cross( p2 - p1, p3 - p1 ) ) ;
        }
        template< class VEC > static  VEC random_point_in_triangle(
            const VEC& p1,
            const VEC& p2,
            const VEC& p3 )
        {
            float64 l1 = std::rand() ;
            float64 l2 = std::rand() ;
            if( l1 + l2 > 1.0 ) {
                l1 = 1.0 - l1 ;
                l2 = 1.0 - l2 ;
            }
            float64 l3 = 1.0 - l1 - l2 ;
            return l1 * p1 + l2 * p2 + l3 * p3 ;
        }
        template< class T >static bool check_order( std::vector< T > v1, T p1, T p2 ) {
        	v1.resize( v1.size()+ 1 ) ;
        	v1[ v1.size() -1 ] = v1[0] ;
        	for( index_t i = 0 ; i < v1.size()-1 ; i++ ){
        		if( v1[ i ] == p1 && v1[ i + 1 ] == p2 )
        		{
        			return true ;
        		}
        		return false ;

        	}
        return false ;
        } ;
        template< class T > static std::vector< T > intersect(
            const std::vector< T >& v1,
            const  std::vector< T >& v2 )
        {
			std::vector< T > intersect ( v1.size() + v2.size() ) ;
			std::vector< index_t >::iterator it ;
			it = std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), intersect.begin()) ;
			intersect.resize( it - intersect.begin() ) ;
			return intersect ;
        }

        template< class T > static bool contains(
            const std::vector< T >& v,
            const T& t )
        {
            return find( v, t ) != -1 ;
        }      
        template< class T > static signed_index_t find( const std::vector< T >& v, const T& t )
        {
            for( index_t i = 0; i < v.size(); i++ ) {
                if( v[i] == t ) return i ;
            }
            return -1 ;
        }
        template< class T1, class T2 > static bool inexact_equal(
            const T1& v1,
            const T2& v2 )
        {
            for( index_t i = 0; i < 3; i++ ) {
                float64 diff( v1[i] - v2[i] ) ;
                if( diff > epsilon || diff < -epsilon ) {
                    return false ;
                }
            }
            return true ;
        }
        template< class T1, class T2 > static bool triple_equal(
            const T1& rhs1,
            const T1& rhs2,
            const T1& rhs3,
            const T2& lhs1,
            const T2& lhs2,
            const T2& lhs3 )
        {
            if( rhs1 == lhs1 ) {
                if( rhs2 == lhs2 && rhs3 == lhs3 ) {
                    return true ;
                } else if( rhs2 == lhs3 && rhs3 == lhs2 ) {
                    return true ;
                }
            } else if( rhs1 == lhs2 ) {
                if( rhs2 == lhs1 && rhs3 == lhs3 ) {
                    return true ;
                } else if( rhs2 == lhs3 && rhs3 == lhs1 ) {
                    return true ;
                }
            } else if( rhs1 == lhs3 ) {
                if( rhs2 == lhs1 && rhs3 == lhs2 ) {
                    return true ;
                } else if( rhs2 == lhs2 && rhs3 == lhs1 ) {
                    return true ;
                }
            }
            return false ;
        }
        // See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
        template< class VEC > static float64 point_triangle_distance(
            const VEC& point,
            const VEC& V0,
            const VEC& V1,
            const VEC& V2,
            VEC& closest_point,
            float64& lambda0 = dummy_float64,
            float64& lambda1 = dummy_float64,
            float64& lambda2 = dummy_float64 )
        {
            VEC diff = V0 - point ;
            VEC edge0 = V1 - V0 ;
            VEC edge1 = V2 - V0 ;
            float64 a00 = length2( edge0 ) ;
            float64 a01 = dot( edge0, edge1 ) ;
            float64 a11 = length2( edge1 ) ;
            float64 b0 = dot( diff, edge0 ) ;
            float64 b1 = dot( diff, edge1 ) ;
            float64 c = length2( diff ) ;
            float64 det = ::fabs( a00 * a11 - a01 * a01 ) ;
            float64 s = a01 * b1 - a11 * b0 ;
            float64 t = a01 * b0 - a00 * b1 ;
            float64 sqrDistance ;

            if( s + t <= det ) {
                if( s < 0.0 ) {
                    if( t < 0.0 ) { // region 4
                        if( b0 < 0.0 ) {
                            t = 0.0 ;
                            if( -b0 >= a00 ) {
                                s = 1.0 ;
                                sqrDistance = a00 + 2.0 * b0 + c ;
                            } else {
                                s = -b0 / a00 ;
                                sqrDistance = b0 * s + c ;
                            }
                        } else {
                            s = 0.0 ;
                            if( b1 >= 0.0 ) {
                                t = 0.0 ;
                                sqrDistance = c ;
                            } else if( -b1 >= a11 ) {
                                t = 1.0 ;
                                sqrDistance = a11 + 2.0 * b1 + c ;
                            } else {
                                t = -b1 / a11 ;
                                sqrDistance = b1 * t + c ;
                            }
                        }
                    } else { // region 3
                        s = 0.0 ;
                        if( b1 >= 0.0 ) {
                            t = 0.0 ;
                            sqrDistance = c ;
                        } else if( -b1 >= a11 ) {
                            t = 1.0 ;
                            sqrDistance = a11 + 2.0 * b1 + c ;
                        } else {
                            t = -b1 / a11 ;
                            sqrDistance = b1 * t + c ;
                        }
                    }
                } else if( t < 0.0 ) { // region 5
                    t = 0.0 ;
                    if( b0 >= 0.0 ) {
                        s = 0.0 ;
                        sqrDistance = c ;
                    } else if( -b0 >= a00 ) {
                        s = 1.0 ;
                        sqrDistance = a00 + 2.0 * b0 + c ;
                    } else {
                        s = -b0 / a00 ;
                        sqrDistance = b0 * s + c ;
                    }
                } else { // region 0
                    // minimum at interior point
                    float64 invDet = float64( 1.0 ) / det ;
                    s *= invDet ;
                    t *= invDet ;
                    sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                        + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c ;
                }
            } else {
                float64 tmp0, tmp1, numer, denom ;

                if( s < 0.0 ) { // region 2
                    tmp0 = a01 + b0 ;
                    tmp1 = a11 + b1 ;
                    if( tmp1 > tmp0 ) {
                        numer = tmp1 - tmp0 ;
                        denom = a00 - 2.0 * a01 + a11 ;
                        if( numer >= denom ) {
                            s = 1.0 ;
                            t = 0.0 ;
                            sqrDistance = a00 + 2.0 * b0 + c ;
                        } else {
                            s = numer / denom ;
                            t = 1.0 - s ;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c ;
                        }
                    } else {
                        s = 0.0 ;
                        if( tmp1 <= 0.0 ) {
                            t = 1.0 ;
                            sqrDistance = a11 + 2.0 * b1 + c ;
                        } else if( b1 >= 0.0 ) {
                            t = 0.0 ;
                            sqrDistance = c ;
                        } else {
                            t = -b1 / a11 ;
                            sqrDistance = b1 * t + c ;
                        }
                    }
                } else if( t < 0.0 ) { // region 6
                    tmp0 = a01 + b1 ;
                    tmp1 = a00 + b0 ;
                    if( tmp1 > tmp0 ) {
                        numer = tmp1 - tmp0 ;
                        denom = a00 - 2.0 * a01 + a11 ;
                        if( numer >= denom ) {
                            t = 1.0 ;
                            s = 0.0 ;
                            sqrDistance = a11 + 2.0 * b1 + c ;
                        } else {
                            t = numer / denom ;
                            s = 1.0 - t ;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c ;
                        }
                    } else {
                        t = 0.0 ;
                        if( tmp1 <= 0.0 ) {
                            s = 1.0 ;
                            sqrDistance = a00 + 2.0 * b0 + c ;
                        } else if( b0 >= 0.0 ) {
                            s = 0.0 ;
                            sqrDistance = c ;
                        } else {
                            s = -b0 / a00 ;
                            sqrDistance = b0 * s + c ;
                        }
                    }
                } else { // region 1
                    numer = a11 + b1 - a01 - b0 ;
                    if( numer <= 0.0 ) {
                        s = 0.0 ;
                        t = 1.0 ;
                        sqrDistance = a11 + 2.0 * b1 + c ;
                    } else {
                        denom = a00 - 2.0 * a01 + a11 ;
                        if( numer >= denom ) {
                            s = 1.0 ;
                            t = 0.0 ;
                            sqrDistance = a00 + 2.0 * b0 + c ;
                        } else {
                            s = numer / denom ;
                            t = 1.0 - s ;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c ;
                        }
                    }
                }
            }

            // Account for numerical round-off error.
            if( sqrDistance < 0.0 ) {
                sqrDistance = 0.0 ;
            }

            closest_point = V0 + s * edge0 + t * edge1 ;
            lambda0 = 1.0 - s - t ;
            lambda1 = s ;
            lambda2 = t ;
            return sqrt( sqrDistance ) ;
        }
        static float64 point_quad_distance(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            vec3& nearest_p ) ;
        static float64 point_tetra_distance(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            vec3& nearest_p ) ;
        static float64 point_pyramid_distance(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4,
            vec3& nearest_p ) ;
        static float64 point_prism_distance(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4,
            const vec3& p5,
            vec3& nearest_p ) ;
        static float64 point_hexa_distance(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4,
            const vec3& p5,
            const vec3& p6,
            const vec3& p7,
            vec3& nearest_p ) ;
        static float64 nearest_point_segment(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            vec3& nearest_p ) ;
        static bool circle_plane_intersection(
            const vec3& O_plane,
            const vec3& N_plane,
            const vec3& O_circle,
            const vec3& N_circle,
            float64 r,
            std::vector< vec3 >& result ) ;
        static bool circle_triangle_intersection(
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& O_circle,
            const vec3& N_circle,
            float64 r,
            std::vector< vec3 >& result ) ;
        static bool plan_plane_intersection(
            const vec3& O_P0,
            const vec3& N_P0,
            const vec3& O_P1,
            const vec3& N_P1,
            vec3& O_inter,
            vec3& N_inter ) ;
        static bool point_inside_triangle(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;
        static bool point_inside_quad(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3 ) ;
        static bool point_inside_tetra(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3 ) ;
        static bool point_inside_pyramid(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4 ) ;
        static bool point_inside_prism(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4,
            const vec3& p5 ) ;
        static bool point_inside_hexa(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& p3,
            const vec3& p4,
            const vec3& p5,
            const vec3& p6,
            const vec3& p7 ) ;
        static bool point_segment_projection(
            const vec3& p,
            const vec3& p0,
            const vec3& p1,
            vec3& new_p ) ;
        static bool segment_triangle_intersection(
            const vec3& seg0,
            const vec3& seg1,
            const vec3& trgl0,
            const vec3& trgl1,
            const vec3& trgl2,
            vec3& result ) ;
    } ;

    class InputStream {
    public:
        InputStream( std::istream& in )
            : in_( in ), line_in_( nil )
        {
        }
        ~InputStream()
        {
            delete line_in_ ;
            line_in_ = nil ;
        }
        bool eof() const
        {
            return in_.eof() ;
        }
        bool eol() const
        {
            return line_in_ == nil || line_in_->eof() ;
        }
        bool ok() const
        {
            return in_ != 0 ;
        }

        void get_line()
        {
            in_.getline( buffer_, 65536 ) ;
            bool check_multiline = true ;
            signed_index_t total_length = 65536 ;
            char* ptr = buffer_ ;

            // If the line ends with a backslash, append
            // the next line to the current line.
            while( check_multiline ) {
                signed_index_t L = (int) strlen( ptr ) ;
                total_length -= L ;
                ptr = ptr + L - 2 ;
                if( *ptr == '\\' && total_length > 0 ) {
                    *ptr = ' ' ;
                    ptr++ ;
                    in_.getline( ptr, total_length ) ;
                } else {
                    check_multiline = false ;
                }
            }

            if( total_length < 0 ) {
                std::cerr << "MultiLine longer than 65536 bytes" << std::endl ;
            }

            delete line_in_ ;
            line_in_ = new std::istringstream( buffer_ ) ;
        }

        std::istream& line()
        {
            grgmesh_assert( line_in_ != nil ) ;
            return *line_in_ ;
        }

        const char *current_line() const
        {
            return buffer_ ;
        }

        template< class T > InputStream& operator>>( T& param )
        {
            param = T() ; // reset do default value, in case *line_in_ is EOF.
            *line_in_ >> param ;
            return *this ;
        }

    private:
        std::istream& in_ ;
        std::istringstream* line_in_ ;
        char buffer_[65536] ;
    } ;

    class GRGMESH_API MakeUnique {
    public:
        MakeUnique( const std::vector< vec3 >& data ) ;
        template< class T > MakeUnique( const std::vector< T >& data )
        {
            signed_index_t nb_points = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                nb_points += data[i].points().size() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            signed_index_t cur_id = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                for( index_t p = 0; p < data[i].points().size();
                    p++, cur_id++ ) {
                    points_[cur_id] = data[i].points()[p] ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }
        template< class T > MakeUnique(
            const std::vector< T >& data,
            bool T_is_a_pointer )
        {
            signed_index_t nb_points = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                nb_points += data[i]->nb_vertices() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            signed_index_t cur_id = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                for( index_t p = 0; p < data[i]->nb_vertices(); p++, cur_id++ ) {
                    points_[cur_id] = data[i]->vertex( p ) ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }

        void add_edges( const std::vector< Edge >& points ) ;
        void add_points( const std::vector< vec3 >& points ) ;

        void unique( index_t nb_neighbors = 5 ) ;

        const std::vector< vec3 >& points() const
        {
            return points_ ;
        }
        void unique_points( std::vector< vec3 >& results ) const ;
        const std::vector< index_t >& indices() const
        {
            return indices_ ;
        }

    private:
        std::vector< vec3 > points_ ;
        std::vector< index_t > indices_ ;
    } ;



    class GRGMESH_API ColocaterANN {
    public:
        enum MeshLocation {
            VERTICES, FACETS, CELLS
        };

        ColocaterANN( const Surface& mesh ) ;
        ColocaterANN( const Line& mesh ) ;
        ColocaterANN( const GEO::Mesh& mesh, const MeshLocation& location ) ;
        ColocaterANN( const std::vector< vec3 >& vertices ) ;
        ColocaterANN( float64* vertices, index_t nb_vertices ) ;
        ColocaterANN( const std::vector< Edge >& edges ) ;

        ~ColocaterANN()
        {
        }

        bool get_colocated(
            const vec3& v,
            index_t nb_neighbors,
            std::vector< index_t >& result ) const ;

        index_t get_neighbors(
            const vec3& v,
            index_t nb_neighbors,
            std::vector< index_t >& result,
            double * dist = nil ) const ;

        vec3 point( signed_index_t i )
        {
            vec3 p ;
            p.x = ann_tree_->point_ptr(i)[0] ;
            p.y = ann_tree_->point_ptr(i)[1] ;
            p.z = ann_tree_->point_ptr(i)[2] ;
            return p ;
        }

    private:
        GEO::NearestNeighborSearch_var ann_tree_ ;
    } ;

    template< class T, index_t n >
    class Array {
    public:
        void assign( const std::vector< T >& values )
        {
            grgmesh_debug_assert( values.size() < n + 1 ) ;
            for( index_t i = 0; i < values.size(); i++ ) {
                values_[i] = values[i] ;
            }
        }

        T value( signed_index_t i ) const
        {
            return values_[i] ;
        }
        T& value( signed_index_t i )
        {
            return values_[i] ;
        }

        float64 normalized_value( signed_index_t i, float64 max, float64 min, float64 scale ) const
        {
            float64 s = values_[i] ;
            s = ( s - min ) / ( max - min ) ;
            s = std::min( s, 1.0 ) ;
            s = std::max( s, 0.0 ) ;
            s *= scale ;
            return s ;
        }

    protected:
        T values_[n] ;
    } ;

    template< index_t n >
    class intArrayTmpl: public Array< int, n > {
    public:
        intArrayTmpl()
        {
            for( index_t i = 0; i < n; i++ ) {
                this->values_[i] = -1 ;
            }
        }
    } ;
    typedef intArrayTmpl< 6 > intArray ;
    typedef intArrayTmpl< 12 > edgeArray ;

    template< index_t n >
    class boolArrayTmpl: public Array< bool, n > {
    public:
        boolArrayTmpl()
        {
            for( index_t i = 0; i < n; i++ ) {
                this->values_[i] = false ;
            }
        }
    } ;
    typedef boolArrayTmpl< 6 > boolArray ;

    class Edge: public Array< vec3, 2 > {
    public:
        Edge( const vec3& v0, const vec3& v1 )
        {
            values_[0] = v0 ;
            values_[1] = v1 ;
        }
        vec3 barycenter() const
        {
            return ( values_[0] + values_[1] ) / static_cast< float64 >( 2 ) ;
        }
    } ;

}

#endif
