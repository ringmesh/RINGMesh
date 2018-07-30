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

#include <ringmesh/basic/geometry.h>

#include <geogram/mesh/mesh.h>

#include <geogram/numerics/predicates.h>

/*!
 * @file Basic geometrical requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace
{
    using namespace RINGMesh;

    bool is_almost_zero( double value )
    {
        return value < global_epsilon && value > -global_epsilon;
    }

    bool point_inside_segment_exact(
        const Geometry::Point3D& point, const Geometry::Segment3D& segment )
    {
        Sign s1{ Position::point_side_to_plane(
            point, { segment.direction(), segment.p0 } ) };
        Sign s2{ Position::point_side_to_plane(
            point, { segment.direction(), segment.p1 } ) };
        return s1 == ZERO || s2 == ZERO || s1 != s2;
    }

    bool point_inside_segment_exact(
        const Geometry::Point2D& point, const Geometry::Segment2D& segment )
    {
        return sign( GEO::PCK::orient_2d( segment.p0, segment.p1, point ) )
               == ZERO;
    }

    template < index_t DIMENSION >
    bool point_inside_segment_approx( const Geometry::Point< DIMENSION >& point,
        const Geometry::Segment< DIMENSION >& segment )
    {
        double distance;
        std::tie( distance, std::ignore ) =
            Distance::point_to_segment( point, segment );
        if( distance > global_epsilon )
        {
            return false;
        }
        double half_length{ segment.length() / 2. };
        vecn< DIMENSION > segment_center{ segment.barycenter() };
        double point_distance{ length( point - segment_center ) };
        if( point_distance < half_length - global_epsilon )
        {
            return true;
        }
        if( point_distance > half_length + global_epsilon )
        {
            return false;
        }
        return point_inside_segment_exact( point, segment );
    }

    bool point_inside_triangle_exact(
        const Geometry::Point2D& point, const Geometry::Triangle2D& triangle )
    {
        Sign s1{ Position::point_side_to_segment(
            point, { triangle.p0, triangle.p1 } ) };
        Sign s2{ Position::point_side_to_segment(
            point, { triangle.p1, triangle.p2 } ) };
        Sign s3{ Position::point_side_to_segment(
            point, { triangle.p2, triangle.p0 } ) };

        if( s1 == ZERO )
        {
            if( s2 == ZERO || s3 == ZERO )
            {
                // Case where p is exactly equal to one triangle vertex
                return true;
            }
            return s2 == s3;
        }
        if( s2 == ZERO )
        {
            if( s1 == ZERO || s3 == ZERO )
            {
                return true;
            }
            return s1 == s3;
        }
        if( s3 == ZERO )
        {
            if( s1 == ZERO || s2 == ZERO )
            {
                return true;
            }
            return s1 == s2;
        }
        return s1 == s2 && s2 == s3;
    }

    bool point_inside_triangle_approx(
        const Geometry::Point2D& point, const Geometry::Triangle2D& triangle )
    {
        double area1{ GEO::Geom::triangle_signed_area(
            point, triangle.p0, triangle.p1 ) };
        if( is_almost_zero( area1 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s1{ sign( area1 ) };
        double area2{ GEO::Geom::triangle_signed_area(
            point, triangle.p1, triangle.p2 ) };
        if( is_almost_zero( area2 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s2{ sign( area2 ) };
        double area3{ GEO::Geom::triangle_signed_area(
            point, triangle.p2, triangle.p0 ) };
        if( is_almost_zero( area3 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s3{ sign( area3 ) };
        return s1 == s2 && s2 == s3;
    }

    Geometry::Plane plane_from_triangle_normal_and_edge(
        const vec3& normal, const vec3& e0, const vec3& e1 )
    {
        return { cross( normal, normalize( e1 - e0 ) ), e0 };
    }

    bool point_inside_triangle_exact(
        const Geometry::Point3D& point, const Geometry::Triangle3D& triangle )
    {
        vec3 triangle_normal{ triangle.plane().normal };
        Sign s1{ Position::point_side_to_plane(
            point, plane_from_triangle_normal_and_edge(
                       triangle_normal, triangle.p0, triangle.p1 ) ) };
        Sign s2{ Position::point_side_to_plane(
            point, plane_from_triangle_normal_and_edge(
                       triangle_normal, triangle.p1, triangle.p2 ) ) };
        Sign s3{ Position::point_side_to_plane(
            point, plane_from_triangle_normal_and_edge(
                       triangle_normal, triangle.p2, triangle.p0 ) ) };

        if( s1 == ZERO )
        {
            if( s2 == ZERO || s3 == ZERO )
            {
                // Case where p is exactly equal to one triangle vertex
                return true;
            }
            return s2 == s3;
        }
        if( s2 == ZERO )
        {
            if( s1 == ZERO || s3 == ZERO )
            {
                return true;
            }
            return s1 == s3;
        }
        if( s3 == ZERO )
        {
            if( s1 == ZERO || s2 == ZERO )
            {
                return true;
            }
            return s1 == s2;
        }
        return s1 == s2 && s2 == s3;
    }

    bool point_inside_triangle_approx(
        const Geometry::Point3D& point, const Geometry::Triangle3D& triangle )
    {
        // Get another point not in the triangle plane (using its normal)
        vec3 translated_point{ point + triangle.plane().normal };

        double vol1{ GEO::Geom::tetra_signed_volume(
            point, translated_point, triangle.p0, triangle.p1 ) };
        if( is_almost_zero( vol1 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s1{ sign( vol1 ) };
        double vol2{ GEO::Geom::tetra_signed_volume(
            point, translated_point, triangle.p1, triangle.p2 ) };
        if( is_almost_zero( vol2 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s2{ sign( vol2 ) };
        double vol3{ GEO::Geom::tetra_signed_volume(
            point, translated_point, triangle.p2, triangle.p0 ) };
        if( is_almost_zero( vol3 ) )
        {
            return point_inside_triangle_exact( point, triangle );
        }
        Sign s3{ sign( vol3 ) };
        return s1 == s2 && s2 == s3;
    }

    bool point_inside_tetra_exact(
        const vec3& p, std::array< vec3, 4 >& vertices )
    {
        std::array< Sign, 4 > signs;
        for( auto f : range( 4 ) )
        {
            signs[f] = sign( GEO::PCK::orient_3d( p.data(),
                vertices[Geometry::Tetra::tetra_facet_vertex[f][0]].data(),
                vertices[Geometry::Tetra::tetra_facet_vertex[f][1]].data(),
                vertices[Geometry::Tetra::tetra_facet_vertex[f][2]].data() ) );
        }
        return ( signs[0] >= 0 && signs[1] >= 0 && signs[2] >= 0
                   && signs[3] >= 0 )
               || ( signs[0] <= 0 && signs[1] <= 0 && signs[2] <= 0
                      && signs[3] <= 0 );
    }

    bool point_inside_tetra_approx(
        const vec3& p, std::array< vec3, 4 >& vertices )
    {
        std::array< Sign, 4 > signs;
        for( const index_t f : range( 4 ) )
        {
            double volume{ GEO::Geom::tetra_signed_volume( p,
                vertices[Geometry::Tetra::tetra_facet_vertex[f][0]],
                vertices[Geometry::Tetra::tetra_facet_vertex[f][1]],
                vertices[Geometry::Tetra::tetra_facet_vertex[f][2]] ) };
            if( is_almost_zero( volume ) )
            {
                return point_inside_tetra_exact( p, vertices );
            }
            signs[f] = sign( volume );
        }
        return ( signs[0] >= 0 && signs[1] >= 0 && signs[2] >= 0
                   && signs[3] >= 0 )
               || ( signs[0] <= 0 && signs[1] <= 0 && signs[2] <= 0
                      && signs[3] <= 0 );
    }
} // namespace

namespace RINGMesh
{
    namespace Position
    {
        bool point_inside_tetra(
            const Geometry::Point3D& point, const Geometry::Tetra& tetra )
        {
            std::array< vec3, 4 > vertices{ { tetra.p0, tetra.p1, tetra.p2,
                tetra.p3 } };
            return point_inside_tetra_approx( point, vertices );
        }

        template < index_t DIMENSION >
        bool point_inside_triangle( const Geometry::Point< DIMENSION >& point,
            const Geometry::Triangle< DIMENSION >& triangle )
        {
            return point_inside_triangle_approx( point, triangle );
        }

        template < index_t DIMENSION >
        bool point_inside_segment( const Geometry::Point< DIMENSION >& point,
            const Geometry::Segment< DIMENSION >& segment )
        {
            return point_inside_segment_approx( point, segment );
        }

        Sign point_side_to_segment(
            const Geometry::Point2D& point, const Geometry::Segment2D& segment )
        {
            return sign( GEO::PCK::orient_2d( point, segment.p0, segment.p1 ) );
        }

        Sign point_side_to_plane(
            const Geometry::Point3D& point, const Geometry::Plane& plane )
        {
            double distance;
            vec3 projected_point;
            std::tie( distance, projected_point ) =
                Distance::point_to_plane( point, plane );

            vec3 point_on_plane{ projected_point };
            double translation{ std::max( 1.0, distance ) };
            for( auto d : range( 3 ) )
            {
                if( std::fabs( plane.normal[d] ) > global_epsilon )
                {
                    index_t d1{ ( d + 1 ) % 3 };
                    index_t d2{ ( d + 2 ) % 3 };
                    point_on_plane[d1] += translation;
                    point_on_plane[d2] += translation;

                    point_on_plane[d] =
                        -( plane.plane_constant()
                            + plane.normal[d1] * point_on_plane[d1]
                            + plane.normal[d2] * point_on_plane[d2] )
                        / plane.normal[d];
                    break;
                }
            }
            ringmesh_assert( point_on_plane != projected_point );

            vec3 u{ normalize( point_on_plane ) };
            vec3 v{ cross( plane.normal, u ) };

            vec3 p0{ projected_point + distance * u };
            vec3 p1{ projected_point
                     + distance
                           * ( std::cos( 2 * M_PI / 3 ) * u
                                 - std::sin( 2 * M_PI / 3 ) * v ) };
            vec3 p2{ projected_point
                     + distance
                           * ( std::cos( 2 * M_PI / 3 ) * u
                                 + std::sin( 2 * M_PI / 3 ) * v ) };

            return sign( GEO::PCK::orient_3d(
                point.data(), p0.data(), p1.data(), p2.data() ) );
        }

        double basic_api segment_angle( const Geometry::Segment2D& segment1,
            const Geometry::Segment2D& segment2 )
        {
            vec2 seg1{ segment1.direction() };
            vec2 seg2{ segment2.direction() };
            double angle_between_pi_and_minus_pi{
                std::atan2( seg1.y, seg1.x ) - std::atan2( seg2.y, seg2.x )
            };
            if( angle_between_pi_and_minus_pi < 0 )
            {
                return angle_between_pi_and_minus_pi + 2 * M_PI;
            }
            return angle_between_pi_and_minus_pi;
        }

        template bool basic_api point_inside_triangle(
            const Geometry::Point< 2 >&, const Geometry::Triangle< 2 >& );

        template bool basic_api point_inside_segment(
            const Geometry::Point< 2 >&, const Geometry::Segment< 2 >& );

        template bool basic_api point_inside_triangle(
            const Geometry::Point< 3 >&, const Geometry::Triangle< 3 >& );

        template bool basic_api point_inside_segment(
            const Geometry::Point< 3 >&, const Geometry::Segment< 3 >& );
    } // namespace Position
} // namespace RINGMesh
