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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_UTILS__
#define __RINGMESH_UTILS__

#include <ringmesh/common.h>

#include <geogram/points/kd_tree.h>
#include <geogram/mesh/mesh_AABB.h>

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    class BoundaryModelElement ;
    class Edge ;
    class Surface ;
    class Line ;
}

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

    static bool operator==( const vec3& u, const vec3& v )
    {
        return u.x == v.x && u.y == v.y && u.z == v.z ;
    }

    static bool operator<( const vec3& u, const vec3& v )
    {
        return u.x < v.x && u.y < v.y && u.z < v.z ;
    }

    static bool operator!=( const vec3& u, const vec3& v )
    {
        return u.x != v.x || u.y != v.y || u.z != v.z ;
    }

    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    } ;

    template< class T >
    inline T sqr( T x )
    {
        return x * x ;
    }

    template< class T >
    inline Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO ) ;
    }

    class Box3d: public GEO::Box {
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

        void add_point( const vec3& p )
        {
            if( !initialized_ ) {
                for( index_t i = 0; i < 3; i++ ) {
                    xyz_min[i] = p[i] ;
                    xyz_max[i] = p[i] ;
                }
                initialized_ = true ;
            } else {
                for( index_t i = 0; i < 3; i++ ) {
                    xyz_min[i] = std::min( xyz_min[i], p[i] ) ;
                    xyz_max[i] = std::max( xyz_max[i], p[i] ) ;
                }
            }
        }

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
    } ;

    class RINGMESH_API Geom {
    public:
        static double mesh_cell_volume( const GEO::Mesh& M, index_t c ) ;

        static vec3 mesh_cell_facet_normal(
            const GEO::Mesh& M,
            index_t c,
            index_t f ) ;

        static vec3 mesh_cell_facet_center(
            const GEO::Mesh& M,
            index_t cell,
            index_t f ) ;

        static vec3 mesh_cell_center( const GEO::Mesh& M, index_t cell ) ;

        static bool has_edge(
            const GEO::Mesh& mesh,
            index_t t,
            index_t p0,
            index_t p1,
            index_t& edge ) ;

        static index_t next_arround_edge(
            const GEO::Mesh& mesh,
            index_t t,
            index_t prev,
            index_t p0,
            index_t p1 ) ;

        static void edges_arround_edge(
            const GEO::Mesh& mesh,
            index_t t,
            index_t p0,
            index_t p1,
            std::vector< index_t >& result ) ;

        static void divide_edge_in_parts(
            const GEO::Mesh& mesh,
            index_t edge,
            index_t nb_parts,
            std::vector< vec3 >& points ) ;

        static void divide_edge_in_parts(
            vec3& node0,
            vec3& node1,
            index_t nb_parts,
            std::vector< vec3 >& points ) ;

    } ;
    class RINGMESH_API Utils {
    public:
        static void print_bounded_attributes( const GEO::Mesh& M ) ;
        static bool compare_file( const std::string& f1, const std::string& f2 ) ;

        static index_t get_nearest_vertex_index(
            const GEO::Mesh& mesh,
            const vec3& p,
            index_t t ) ;
        static bool facets_have_same_orientation(
            const GEO::Mesh& mesh,
            index_t f1,
            index_t c11,
            index_t f2 ) ;

        static void mesh_facet_connect( GEO::Mesh& mesh ) ;
        static void check_and_repair_mesh_consistency(
            const BoundaryModelElement& region,
            GEO::Mesh& mesh,
            bool check_duplicated_facet = false ) ;
        static float64 triangle_area(
            const vec3& p1,
            const vec3& p2,
            const vec3& p3 )
        {
            return 0.5 * length( cross( p2 - p1, p3 - p1 ) ) ;
        }

        template< class VEC >
        static VEC random_point_in_triangle(
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

        template< class T >
        static bool check_order( std::vector< T > v1, T p1, T p2 )
        {
            v1.resize( v1.size() + 1 ) ;
            v1[v1.size() - 1] = v1[0] ;
            for( index_t i = 0; i < v1.size() - 1; i++ ) {
                if( v1[i] == p1 && v1[i + 1] == p2 ) {
                    return true ;
                }
                return false ;
            }
            return false ;
        }

        template< class T >
        static std::vector< T > intersect(
            const std::vector< T >& v1,
            const std::vector< T >& v2 )
        {
            std::vector< T > intersect( v1.size() + v2.size() ) ;
            std::vector< index_t >::iterator it ;
            it = std::set_intersection( v1.begin(), v1.end(), v2.begin(), v2.end(),
                intersect.begin() ) ;
            intersect.resize( it - intersect.begin() ) ;
            return intersect ;
        }

        template< typename T, typename container >
        static bool contains( const container& v, const T& t, bool sorted = false )
        {
            if( sorted )
                return find_sorted( v, t ) != NO_ID ;
            else
                return find( v, t ) != NO_ID ;
        }

        template< typename T, typename container >
        static index_t find( const container& v, const T& t )
        {
            typename container::const_iterator it = std::find(
                v.begin(), v.end(), t ) ;
            if( it == v.end() )
                return NO_ID ;
            else
                return it - v.begin() ;
        }

        template< typename T, typename container >
        static index_t find_sorted( const container& v, const T& t )
        {
            typename container::const_iterator low = std::lower_bound(
                v.begin(), v.end(), t ) ;
            if( low == v.end() || t < *low )
                return NO_ID ;
            else
                return low - v.begin() ;
        }

        template< class T1, class T2 >
        static bool inexact_equal( const T1& v1, const T2& v2 )
        {
            for( index_t i = 0; i < 3; i++ ) {
                float64 diff( v1[i] - v2[i] ) ;
                if( diff > epsilon || diff < -epsilon ) {
                    return false ;
                }
            }
            return true ;
        }

        template< class T1, class T2 >
        static bool triple_equal(
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
    } ;

    class RINGMESH_API Math {
    public:
        // See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
        template< class VEC >
        static float64 point_triangle_distance(
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

        static void rotation_matrix_about_arbitrary_axis(
            const vec3& origin,
            const vec3& axis,
            float64 theta,
            bool degrees,
            GEO::Matrix< float64, 4 >& rot_mat ) ;
        static void rotate_mesh(
            GEO::Mesh& mesh,
            const GEO::Matrix< float64, 4 >& rot_mat ) ;
    } ;

    /*!
     * Given an array of vec3, this class computes the colocated points
     * and a database to identify which colocated point corresponds to
     */
    class RINGMESH_API MakeUnique {
    ringmesh_disable_copy( MakeUnique ) ;
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
                for( index_t p = 0; p < data[i].points().size(); p++, cur_id++ ) {
                    points_[cur_id] = data[i].points()[p] ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }

        template< class T > MakeUnique(
            const std::vector< T >& data,
            bool T_is_a_pointer )
        {
            index_t nb_points = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                nb_points += data[i]->nb_vertices() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            index_t cur_id = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                for( index_t p = 0; p < data[i]->nb_vertices(); p++, cur_id++ ) {
                    points_[cur_id] = data[i]->vertex( p ) ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }

        void add_edges( const std::vector< Edge >& points ) ;

        void add_points( const std::vector< vec3 >& points ) ;

        void unique() ;

        /*!
         * Gets the input vector of vec3
         */
        const std::vector< vec3 >& points() const
        {
            return points_ ;
        }

        /*!
         * Gets the number of points in the database
         * @return returns the corresponding number
         */
        index_t nb_points() const
        {
            return points_.size() ;
        }

        void unique_points( std::vector< vec3 >& results ) const ;

        /*!
         * Gets the computed database that maps
         * the colocated point to the unique one
         */
        const std::vector< index_t >& indices() const
        {
            return indices_ ;
        }

    private:
        /// Input vector of vec3
        std::vector< vec3 > points_ ;
        /// computed database that maps the colocated point to the unique one
        std::vector< index_t > indices_ ;
    } ;

    class RINGMESH_API ColocaterANN {
    ringmesh_disable_copy( ColocaterANN ) ;
    public:
        enum MeshLocation {
            VERTICES, FACETS, CELLS
        } ;
        ColocaterANN() ;
        ColocaterANN(
            const Surface& mesh,
            const MeshLocation& location = VERTICES ) ;
        ColocaterANN( const Line& mesh ) ;
        ColocaterANN( const GEO::Mesh& mesh, const MeshLocation& location ) ;
        ColocaterANN( const std::vector< vec3 >& vertices ) ;
        ColocaterANN( float64* vertices, index_t nb_vertices ) ;
        ColocaterANN( const std::vector< Edge >& edges ) ;

        ~ColocaterANN()
        {
            delete[] ann_points_ ;
        }

        void set_points( const std::vector< vec3 >& vertices ) ;

        bool get_colocated( const vec3& v, std::vector< index_t >& result ) const ;

        /*!
         * Gets the closest neighbor point
         * @param[in] v the point to test
         * @param[out] dist the square distance to the closest point
         * return returns the index of the closest point
         */
        index_t get_closest_neighbor(
            const vec3& v,
            double& dist = dummy_float64 ) const
        {
            std::vector< index_t > result ;
            get_neighbors( v, 1, result, &dist ) ;
            return result[0] ;
        }

        index_t get_neighbors(
            const vec3& v,
            index_t nb_neighbors,
            std::vector< index_t >& result,
            double* dist = nil ) const ;

        /*!
         * Gets the point corresponding to the given index
         * @param[in] i the point index
         * return the corresponding point
         */
        vec3 point( index_t i ) const
        {
            return vec3( ann_tree_->point_ptr( i ) ) ;
        }

    private:
        /// KdTree to compute the nearest neighbor search
        GEO::NearestNeighborSearch_var ann_tree_ ;
        /// Array of the points (size of 3xnumber of points)
        double* ann_points_ ;
    } ;

    template< class T, index_t n >
    class Array {
    public:
        void assign( const std::vector< T >& values )
        {
            ringmesh_debug_assert( values.size() < n + 1 ) ;
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

        float64 normalized_value(
            signed_index_t i,
            float64 max,
            float64 min,
            float64 scale ) const
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

    /*!
     * Class to sort two vectors using indirect sorting
     */
    template< class T1, class T2 >
    class IndirectSort {
    public:
        IndirectSort( std::vector< T1 >& input, std::vector< T2 >& output )
            : input_( input ), output_( output )
        {

        }
        void sort()
        {
            if( input_.size() < 2 ) return ;
            for( index_t it1 = 0; it1 < input_.size() - 1; it1++ ) {
                index_t ref_index = it1 ;
                T1& ref_value = input_[it1] ;
                for( index_t it2 = it1 + 1; it2 < input_.size(); it2++ ) {
                    index_t new_index = it2 ;
                    T1& new_value = input_[it2] ;
                    if( ref_value > new_value ) {
                        ref_value = new_value ;
                        ref_index = new_index ;
                    }
                }
                std::iter_swap( input_.begin() + it1, input_.begin() + ref_index ) ;
                std::iter_swap( output_.begin() + it1,
                    output_.begin() + ref_index ) ;

            }
        }

    private:
        std::vector< T1 >& input_ ;
        std::vector< T2 >& output_ ;
    } ;

/******************************************************************/

}

#endif
