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


#ifndef __RINGMESH_GEOMETRY__
#define __RINGMESH_GEOMETRY__

#include <ringmesh/common.h>

#include <geogram/points/kd_tree.h>
#include <geogram/mesh/mesh_AABB.h>

namespace RINGMesh {

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

    // See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
    template< class VEC >
    float64 point_triangle_distance(
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

    float64 point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p ) ;

    float64 point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p ) ;

    float64 point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p ) ;

    float64 point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p ) ;

    float64 point_hexa_distance(
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

    bool circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result ) ;

    bool circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result ) ;

    bool plan_plane_intersection(
        const vec3& O_P0,
        const vec3& N_P0,
        const vec3& O_P1,
        const vec3& N_P1,
        vec3& O_inter,
        vec3& N_inter ) ;

    bool point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 ) ;

    bool point_inside_quad(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 ) ;

    bool point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 ) ;

    bool point_inside_pyramid(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4 ) ;

    bool point_inside_prism(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5 ) ;

    bool point_inside_hexa(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7 ) ;

    bool point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p ) ;

    bool segment_triangle_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& trgl0,
        const vec3& trgl1,
        const vec3& trgl2,
        vec3& result ) ;

    void RINGMESH_API rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        float64 theta,
        bool degrees,
        GEO::Matrix< float64, 4 >& rot_mat ) ;
        

    template< class VEC >
    VEC random_point_in_triangle(
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




    /*!
    * Given an array of vec3, this class computes the colocated points
    * and a database to identify which colocated point corresponds to
    *
    * @todo Do we really need this class ? [JP]
    * Move it in another file it depends on geogram... We need to compartimentalize.
    *
    */
    class RINGMESH_API MakeUnique {
        ringmesh_disable_copy( MakeUnique ) ;
    public:
        MakeUnique( const std::vector< vec3 >& data ) ;
        template< class T > MakeUnique( const std::vector< T >& data )
        {
            signed_index_t nb_points = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                nb_points += data[ i ].points().size() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            signed_index_t cur_id = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                for( index_t p = 0; p < data[ i ].points().size(); p++, cur_id++ ) {
                    points_[ cur_id ] = data[ i ].points()[ p ] ;
                    indices_[ cur_id ] = cur_id ;
                }
            }
        }

        template< class T > MakeUnique(
            const std::vector< T >& data,
            bool T_is_a_pointer )
        {
            index_t nb_points = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                nb_points += data[ i ]->nb_vertices() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            index_t cur_id = 0 ;
            for( index_t i = 0; i < data.size(); i++ ) {
                for( index_t p = 0; p < data[ i ]->nb_vertices(); p++, cur_id++ ) {
                    points_[ cur_id ] = data[ i ]->vertex( p ) ;
                    indices_[ cur_id ] = cur_id ;
                }
            }
        }

        void add_edges( const std::vector< std::pair< vec3, vec3 > > & points ) ;

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


    /*
     * @todo Do we really need this class ? [JP]
     */
    class RINGMESH_API ColocaterANN {
        ringmesh_disable_copy( ColocaterANN ) ;
    public:
        enum MeshLocation {
            VERTICES, FACETS, CELLS
        } ;
        ColocaterANN() ;
        ColocaterANN(
            const GEO::Mesh& mesh,
            const MeshLocation& location = VERTICES,
            bool copy = false ) ;
        ColocaterANN( const std::vector< vec3 >& vertices, bool copy = true ) ;

        ~ColocaterANN()
        {
            if( ann_points_ ) delete[] ann_points_ ;
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
            return result[ 0 ] ;
        }

        index_t get_neighbors(
            const vec3& v,
            index_t nb_neighbors,
            std::vector< index_t >& result,
            double* dist = nil ) const ;

    private:
        /// KdTree to compute the nearest neighbor search
        GEO::NearestNeighborSearch_var ann_tree_ ;
        /// Array of the points (size of 3xnumber of points), possibly nil
        double* ann_points_ ;
    } ;

}

#endif