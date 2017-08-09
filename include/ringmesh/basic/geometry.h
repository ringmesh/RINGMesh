/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#pragma once

#include <ringmesh/basic/common.h>
#include <array>

/*!
 * @file Basic geometrical requests
 * @author Arnaud Botella
 */

namespace RINGMesh {

    template< index_t DIMENSION >
    bool operator==( const vecn< DIMENSION >& u, const vecn< DIMENSION >& v )
    {
        for( index_t i : range( DIMENSION ) ) {
            if( u[i] != v[i] ) return false;
        }
        return true;
    }

    template< index_t DIMENSION >
    bool operator!=( const vecn< DIMENSION >& u, const vecn< DIMENSION >& v )
    {
        return !( u == v );
    }

    template< index_t DIMENSION >
    bool inexact_equal(
        const vecn< DIMENSION >& v1,
        const vecn< DIMENSION >& v2,
        double epsilon )
    {
        return length( v2 - v1 ) < epsilon;
    }

    double RINGMESH_API dot_perp( const vec2& v0, const vec2& v1 );

    /* @warning Duplicate from Geogram/basic/numeric.h */
    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    };
    /* @warning Duplicate from Geogram/basic/numeric.h */
    template< typename T >
    Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO );
    }

    namespace Geometry {

        template< index_t DIMENSION >
        using Point = vecn< DIMENSION >;
        CLASS_DIMENSION_ALIASES( Point );

        template< index_t DIMENSION >
        struct Segment {
            Segment() = default;
            Segment( vecn< DIMENSION > p0, vecn< DIMENSION > p1 )
                : p0_( std::move( p0 ) ), p1_( std::move( p1 ) )
            {
            }
            vecn< DIMENSION > direction() const
            {
                return normalize( p1_ - p0_ );
            }
            vecn< DIMENSION > barycenter() const
            {
                return ( p1_ + p0_ ) / 2.;
            }
            double length() const
            {
                return ( p1_ - p0_ ).length();
            }
            vecn< DIMENSION > p0_;
            vecn< DIMENSION > p1_;
        };
        CLASS_DIMENSION_ALIASES( Segment );

        template< index_t DIMENSION >
        struct Line {
            Line() = default;
            Line( const vecn< DIMENSION >& direction, vecn< DIMENSION > origin )
                :
                    origin_( std::move( origin ) ),
                    direction_( normalize( direction ) )
            {
            }
            Line( Segment< DIMENSION > segment )
                : Line( segment.p1_ - segment.p0_, std::move( segment.p0_ ) )
            {
            }
            vecn< DIMENSION > origin_;
            vecn< DIMENSION > direction_;
        };
        CLASS_DIMENSION_ALIASES( Line );

        struct Plane {
            Plane() = default;
            Plane( const vec3& normal, vec3 origin )
                : normal_( normalize( normal ) ), origin_( std::move( origin ) )
            {
            }
            double plane_constant() const
            {
                double plane_constant { 0.0 };
                for( index_t i : range( 3 ) ) {
                    plane_constant += origin_[i] * normal_[i];
                }
                return plane_constant;
            }
            vec3 normal_;
            vec3 origin_;
        };

        template< index_t DIMENSION >
        struct Triangle {
            Triangle() = default;
            Triangle(
                vecn< DIMENSION > p0,
                vecn< DIMENSION > p1,
                vecn< DIMENSION > p2 )
                :
                    p0_( std::move( p0 ) ),
                    p1_( std::move( p1 ) ),
                    p2_( std::move( p2 ) )
            {
            }
            vecn< DIMENSION > p0_;
            vecn< DIMENSION > p1_;
            vecn< DIMENSION > p2_;
        };
        template< >
        struct Triangle< 3 > {
            Triangle() = default;
            Triangle( vec3 p0, vec3 p1, vec3 p2 )
                :
                    p0_( std::move( p0 ) ),
                    p1_( std::move( p1 ) ),
                    p2_( std::move( p2 ) )
            {
            }
            Plane plane() const
            {
                return {cross( p1_ - p0_, p2_ - p0_ ), p0_};
            }
            vec3 p0_;
            vec3 p1_;
            vec3 p2_;
        };
        CLASS_DIMENSION_ALIASES( Triangle );

        struct Tetra {
            Tetra() = default;
            Tetra( vec3 p0, vec3 p1, vec3 p2, vec3 p3 )
                :
                    p0_( std::move( p0 ) ),
                    p1_( std::move( p1 ) ),
                    p2_( std::move( p2 ) ),
                    p3_( std::move( p3 ) )
            {
            }
            vec3 p0_;
            vec3 p1_;
            vec3 p2_;
            vec3 p3_;
        };

        struct Sphere {
            Sphere() = default;
            Sphere( vec3 origin, double radius )
                : origin_( std::move( origin ) ), radius_( std::move( radius ) )
            {
            }
            vec3 origin_;
            double radius_ { 0 };
        };

        struct Circle {
            Circle() = default;
            Circle( Plane plane, double radius )
                : plane_( std::move( plane ) ), radius_( std::move( radius ) )
            {
            }
            Plane plane_;
            double radius_ { 0 };
        };

        using Disk = Circle;
    }

    namespace Distance {
        /*!
         * See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
         */
        template< index_t DIMENSION >
        std::tuple< double, vecn< DIMENSION > > point_to_segment(
            const Geometry::Point< DIMENSION >& point,
            const Geometry::Segment< DIMENSION >& segment );

        /*!
         * Computes the smallest distance between a point and a triangle
         * @return a tuple containing the following elements (in this order):
         * - the smallest distance
         * - the closest point on the triangle
         */
        template< index_t DIMENSION >
        std::tuple< double, vecn< DIMENSION > > point_to_triangle(
            const Geometry::Point< DIMENSION >& point,
            const Geometry::Triangle< DIMENSION >& triangle );

        /*!
         * Computes the distance between a point and a tetrahedron
         * @return a tuple containing:
         * - the distance between the point and the tetrahedron facets.
         * - the nearest point on the tetrahedron.
         */
        std::tuple< double, vec3 > RINGMESH_API point_to_tetra(
            const Geometry::Point3D& point,
            const Geometry::Tetra& tetra );
        /*!
         * Computes the distance between a point and a plane
         * @return a tuple containing:
         * - the distance between the point and the plane.
         * - the nearest point on the plane.
         */
        std::tuple< double, vec3 > RINGMESH_API point_to_plane(
            const Geometry::Point3D& point,
            const Geometry::Plane& plane );
    }

    namespace Intersection {
        /*!
         * Computes the intersection between a plane and a line
         * @return returns a tuple containing a boolean (true if there is an intersection)
         * and the intersected point if any.
         */
        std::tuple< bool, vec3 > RINGMESH_API line_plane(
            const Geometry::Line3D& line,
            const Geometry::Plane& plane );

        /*!
         * Computes the intersection(s) between a sphere and a line
         * @return returns a tuple containing a boolean (true if there is at least one intersection)
         * and the intersected points.
         */
        std::tuple< bool, std::vector< vec3 > > RINGMESH_API line_sphere(
            const Geometry::Line3D& line,
            const Geometry::Sphere& sphere );

        /*!
         * Computes the intersection(s) between a sphere and a segment
         * @return returns a tuple containing a boolean (true if there is at least one intersection)
         * and the intersected points.
         */
        std::tuple< bool, std::vector< vec3 > > RINGMESH_API segment_sphere(
            const Geometry::Segment3D& segment,
            const Geometry::Sphere& sphere );

        /*!
         * Computes the intersection between a plane and a segment
         * @return returns a tuple containing a boolean (true if there is an intersection)
         * and the intersected point if any.
         */
        std::tuple< bool, vec3 > RINGMESH_API segment_plane(
            const Geometry::Segment3D& segment,
            const Geometry::Plane& plane );

        /*!
         * Computes the intersection of a segment and a triangle
         * @return a tuple containing a boolean (true is there is an intersection)
         * and the intersected point if any.
         */
        std::tuple< bool, vec3 > RINGMESH_API segment_triangle(
            const Geometry::Segment3D& segment,
            const Geometry::Triangle3D& triangle );

        /*!
         * Computes the intersection(s) between a circle and a plane
         * @return returns a tuple containing a boolean (true if there is at least one intersection)
         * and the intersected points if any.
         */
        std::tuple< bool, std::vector< vec3 > > RINGMESH_API circle_plane(
            const Geometry::Circle& circle,
            const Geometry::Plane& plane );

        /*!
         * Computes the intersection between a disk and a segment
         * @return returns a tuple containing a boolean (true if there is an intersection)
         * and the intersected point if any.
         */
        std::tuple< bool, vec3 > RINGMESH_API segment_disk(
            const Geometry::Segment3D& segment,
            const Geometry::Disk& disk );

        /*!
         * Computes the intersection(s) between a circle and a triangle
         * @return returns a tuple containing a boolean (true if there is at least one intersection)
         * and the intersected points if any.
         */
        std::tuple< bool, std::vector< vec3 > > RINGMESH_API triangle_circle(
            const Geometry::Triangle3D& triangle,
            const Geometry::Circle& circle );

        /*!
         * Computes the intersection between two planes
         * @return a tuple containing:
         * - a boolean: true is there is an intersection between the planes.
         * - the intersected line if any.
         */
        std::tuple< bool, Geometry::Line3D > RINGMESH_API plane_plane(
            const Geometry::Plane& plane0,
            const Geometry::Plane& plane1 );

        /*!
         * Computes the intersection between two lines
         * @return a tuple containing:
         * - a boolean: true if there is an intersection.
         * - the intersected point if any.
         */
        std::tuple< bool, vec2 > RINGMESH_API line_line(
            const Geometry::Line2D& line0,
            const Geometry::Line2D& line1 );

        /*!
         * Computes the intersection between two segments
         * @return a tuple containing:
         * - a boolean: true if there is an intersection.
         * - the intersected point if any.
         */
        std::tuple< bool, vec2 > RINGMESH_API segment_segment(
            const Geometry::Segment2D& segment0,
            const Geometry::Segment2D& segment1 );

        /*!
         * Computes the intersection between a segment and a line
         * @return a tuple containing:
         * - a boolean: true if there is an intersection.
         * - the intersected point if any.
         */
        std::tuple< bool, vec2 > RINGMESH_API segment_line(
            const Geometry::Segment2D& segment,
            const Geometry::Line2D& line );
    }

    namespace Position {

        /*!
         * @brief Tests if a point is on a segment
         * @return returns true if the point is inside
         */
        bool RINGMESH_API point_inside_segment(
            const Geometry::Point3D& point,
            const Geometry::Segment3D& segment );

        /*!
         * @brief Tests if a point is inside a triangle
         * @details if it is inside a prism based on the triangle and its normal
         * @return returns true if the point is inside
         */
        template< index_t DIMENSION >
        bool point_inside_triangle(
            const Geometry::Point< DIMENSION >& point,
            const Geometry::Triangle< DIMENSION >& triangle );

        /*!
         * Tests if a point is inside a tetrahedron
         * @return returns true if the point is inside the tetrahedron
         */
        bool RINGMESH_API point_inside_tetra(
            const Geometry::Point3D& point,
            const Geometry::Tetra& tetra );

        /*!
         * Returns the point side to a segment
         */
        Sign RINGMESH_API point_side_to_segment(
            const Geometry::Point2D& point,
            const Geometry::Segment2D& segment );

        /*!
         * Returns the point side to a plane
         */
        Sign RINGMESH_API point_side_to_plane(
            const Geometry::Point3D& point,
            const Geometry::Plane& plane );

    }

    double RINGMESH_API triangle_signed_area(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& triangle_normal );

    /*!
     * Computes the orthogonal projection of a point on a segment
     * @param[in] p the point to project
     * @param[in] p0 the first vertex of the segment
     * @param[in] p1 the second vertex of the segment
     * @return returns a tuple containing a boolean (true if the projection is possible)
     * and the projected point if any.
     */
    template< index_t DIMENSION >
    std::tuple< bool, vecn< DIMENSION > > point_segment_projection(
        const vecn< DIMENSION >& p,
        const vecn< DIMENSION >& p0,
        const vecn< DIMENSION >& p1 );

    /*!
     * Computes barycentric coordinates of \p p
     * @param[in] p the query point
     * @param[in] p0 the first tetra vertex
     * @param[in] p1 the second tetra vertex
     * @param[in] p2 the third tetra vertex
     * @param[in] p3 the fourth tetra vertex
     * @return a tuple containing:
     * - a boolean (false if the computation failed because of too small tetrahedron volume).
     * - the parametric coordinates corresponding to points
     */
    std::tuple< bool, std::array< double, 4 > > RINGMESH_API tetra_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 );

    /*!
     * Computes barycentric coordinates of \p p
     * @param[in] p the query point
     * @param[in] p0 the first triangle vertex
     * @param[in] p1 the second triangle vertex
     * @param[in] p2 the third triangle vertex
     * @return a tuple containing:
     * - a boolean (false if the computation failed because of too small triangle area).
     * - the parametric coordinates corresponding to points.
     */
    std::tuple< bool, std::array< double, 3 > > RINGMESH_API triangle_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 );

    /*!
     * Computes barycentric coordinates of \p p
     * @param[in] p the query point
     * @param[in] p0 the first triangle vertex
     * @param[in] p1 the second triangle vertex
     * @param[in] p2 the third triangle vertex
     * @return a tuple containing:
     * - a boolean (false if the computation failed because of too small triangle area).
     * - the parametric coordinates corresponding to points.
     */
    std::tuple< bool, std::array< double, 3 > > RINGMESH_API triangle_barycentric_coordinates(
        const vec2& p,
        const vec2& p0,
        const vec2& p1,
        const vec2& p2 );

    /*!
     * @brief Builds a rotational matrix about an arbitrary axis.
     *
     * Mathematical development: http://paulbourke.net/geometry/rotate/.
     * @param[in] origin point in which passes the rotation axis.
     * @param[in] axis vector which defines the rotation axis.
     * @param[in] theta rotation angle (in radians or degrees).
     * @param[in] degrees true is \p theta is in degrees, false
     * if in radians.
     * @return the matrix which defines the rotation
     * of a point around the axis defined by point \p origin
     * and vector \p axis by an angle \p theta.
     * New coordinates of a point (x,y,z) are:
     * (x',y',z') = rot_mat*(x,y,z)
     */
    GEO::Matrix< 4, double > RINGMESH_API rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        double theta,
        bool degrees );

}
