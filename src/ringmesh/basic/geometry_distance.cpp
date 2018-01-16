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

#include <array>

#include <geogram/mesh/mesh.h>

/*!
 * @file Basic geometrical requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace RINGMesh
{
    namespace Distance
    {
        template <>
        std::tuple< double, vec3 > point_to_triangle(
            const Geometry::Point3D& point,
            const Geometry::Triangle3D& triangle )
        {
            const vec3& V0 = triangle.p0;
            const vec3& V1 = triangle.p1;
            const vec3& V2 = triangle.p2;
            vec3 diff{ V0 - point };
            vec3 edge0{ V1 - V0 };
            vec3 edge1{ V2 - V0 };
            double a00{ length2( edge0 ) };
            double a01{ dot( edge0, edge1 ) };
            double a11{ length2( edge1 ) };
            double b0{ dot( diff, edge0 ) };
            double b1{ dot( diff, edge1 ) };
            double c{ length2( diff ) };
            double det{ std::fabs( a00 * a11 - a01 * a01 ) };
            double s{ a01 * b1 - a11 * b0 };
            double t{ a01 * b0 - a00 * b1 };
            double sqrDistance;

            if( s + t <= det )
            {
                if( s < 0.0 )
                {
                    if( t < 0.0 )
                    { // region 4
                        if( b0 < 0.0 )
                        {
                            t = 0.0;
                            if( -b0 >= a00 )
                            {
                                s = 1.0;
                                sqrDistance = a00 + 2.0 * b0 + c;
                            }
                            else
                            {
                                s = -b0 / a00;
                                sqrDistance = b0 * s + c;
                            }
                        }
                        else
                        {
                            s = 0.0;
                            if( b1 >= 0.0 )
                            {
                                t = 0.0;
                                sqrDistance = c;
                            }
                            else if( -b1 >= a11 )
                            {
                                t = 1.0;
                                sqrDistance = a11 + 2.0 * b1 + c;
                            }
                            else
                            {
                                t = -b1 / a11;
                                sqrDistance = b1 * t + c;
                            }
                        }
                    }
                    else
                    { // region 3
                        s = 0.0;
                        if( b1 >= 0.0 )
                        {
                            t = 0.0;
                            sqrDistance = c;
                        }
                        else if( -b1 >= a11 )
                        {
                            t = 1.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        }
                        else
                        {
                            t = -b1 / a11;
                            sqrDistance = b1 * t + c;
                        }
                    }
                }
                else if( t < 0.0 )
                { // region 5
                    t = 0.0;
                    if( b0 >= 0.0 )
                    {
                        s = 0.0;
                        sqrDistance = c;
                    }
                    else if( -b0 >= a00 )
                    {
                        s = 1.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    }
                    else
                    {
                        s = -b0 / a00;
                        sqrDistance = b0 * s + c;
                    }
                }
                else
                { // region 0
                    // minimum at interior point
                    double invDet = double( 1.0 ) / det;
                    s *= invDet;
                    t *= invDet;
                    sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                  + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
                }
            }
            else
            {
                double tmp0;
                double tmp1;
                double numer;
                double denom;

                if( s < 0.0 )
                { // region 2
                    tmp0 = a01 + b0;
                    tmp1 = a11 + b1;
                    if( tmp1 > tmp0 )
                    {
                        numer = tmp1 - tmp0;
                        denom = a00 - 2.0 * a01 + a11;
                        if( numer >= denom )
                        {
                            s = 1.0;
                            t = 0.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        }
                        else
                        {
                            s = numer / denom;
                            t = 1.0 - s;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                          + t * ( a01 * s + a11 * t + 2.0 * b1 )
                                          + c;
                        }
                    }
                    else
                    {
                        s = 0.0;
                        if( tmp1 <= 0.0 )
                        {
                            t = 1.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        }
                        else if( b1 >= 0.0 )
                        {
                            t = 0.0;
                            sqrDistance = c;
                        }
                        else
                        {
                            t = -b1 / a11;
                            sqrDistance = b1 * t + c;
                        }
                    }
                }
                else if( t < 0.0 )
                { // region 6
                    tmp0 = a01 + b1;
                    tmp1 = a00 + b0;
                    if( tmp1 > tmp0 )
                    {
                        numer = tmp1 - tmp0;
                        denom = a00 - 2.0 * a01 + a11;
                        if( numer >= denom )
                        {
                            t = 1.0;
                            s = 0.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        }
                        else
                        {
                            t = numer / denom;
                            s = 1.0 - t;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                          + t * ( a01 * s + a11 * t + 2.0 * b1 )
                                          + c;
                        }
                    }
                    else
                    {
                        t = 0.0;
                        if( tmp1 <= 0.0 )
                        {
                            s = 1.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        }
                        else if( b0 >= 0.0 )
                        {
                            s = 0.0;
                            sqrDistance = c;
                        }
                        else
                        {
                            s = -b0 / a00;
                            sqrDistance = b0 * s + c;
                        }
                    }
                }
                else
                { // region 1
                    numer = a11 + b1 - a01 - b0;
                    if( numer <= 0.0 )
                    {
                        s = 0.0;
                        t = 1.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    }
                    else
                    {
                        denom = a00 - 2.0 * a01 + a11;
                        if( numer >= denom )
                        {
                            s = 1.0;
                            t = 0.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        }
                        else
                        {
                            s = numer / denom;
                            t = 1.0 - s;
                            sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                                          + t * ( a01 * s + a11 * t + 2.0 * b1 )
                                          + c;
                        }
                    }
                }
            }

            // Account for numerical round-off error.
            if( sqrDistance < 0.0 )
            {
                sqrDistance = 0.0;
            }

            vec3 closest_point{ V0 + s * edge0 + t * edge1 };
            return std::make_tuple( std::sqrt( sqrDistance ), closest_point );
        }

        template <>
        std::tuple< double, vec2 > point_to_segment(
            const Geometry::Point2D& point, const Geometry::Segment2D& segment )
        {
            // The direction vector is not unit length.  The normalization is
            // deferred
            // until it is needed.
            vec2 direction = segment.p1 - segment.p0;
            vec2 diff = point - segment.p1;
            double t = dot( direction, diff );
            vec2 nearest_p;
            if( t >= global_epsilon )
            {
                nearest_p = segment.p1;
            }
            else
            {
                diff = point - segment.p0;
                t = dot( direction, diff );
                if( t <= global_epsilon )
                {
                    nearest_p = segment.p0;
                }
                else
                {
                    double sqrLength = dot( direction, direction );
                    if( sqrLength > global_epsilon )
                    {
                        t /= sqrLength;
                        nearest_p = segment.p0 + t * direction;
                    }
                    else
                    {
                        nearest_p = segment.p0;
                    }
                }
            }

            diff = point - nearest_p;
            return std::make_tuple( std::sqrt( dot( diff, diff ) ), nearest_p );
        }

        template <>
        std::tuple< double, vec2 > point_to_triangle(
            const Geometry::Point2D& point,
            const Geometry::Triangle2D& triangle )
        {
            double result;
            vec2 closest_point;
            if( Position::point_inside_triangle( point, triangle ) )
            {
                closest_point = point;
                result = 0.0;
            }
            else
            {
                std::array< vec2, 3 > closest;
                std::array< double, 3 > distance;
                std::tie( distance[0], closest[0] ) = point_to_segment(
                    point, Geometry::Segment2D{ triangle.p0, triangle.p1 } );
                std::tie( distance[1], closest[1] ) = point_to_segment(
                    point, Geometry::Segment2D{ triangle.p1, triangle.p2 } );
                std::tie( distance[2], closest[2] ) = point_to_segment(
                    point, Geometry::Segment2D{ triangle.p2, triangle.p0 } );
                if( distance[0] < distance[1] )
                {
                    if( distance[0] < distance[2] )
                    {
                        result = distance[0];
                        closest_point = closest[0];
                    }
                    else
                    {
                        result = distance[2];
                        closest_point = closest[2];
                    }
                }
                else
                {
                    if( distance[1] < distance[2] )
                    {
                        result = distance[1];
                        closest_point = closest[1];
                    }
                    else
                    {
                        result = distance[2];
                        closest_point = closest[2];
                    }
                }
            }
            return std::make_tuple( result, closest_point );
        }

        std::tuple< double, vec3 > point_to_tetra(
            const Geometry::Point3D& point, const Geometry::Tetra& tetra )
        {
            std::array< vec3, 4 > vertices{ { tetra.p0, tetra.p1, tetra.p2,
                tetra.p3 } };
            double dist{ max_float64() };
            vec3 nearest_p;
            for( auto f : range( 4 ) )
            {
                double distance{ max_float64() };
                vec3 cur_p;
                std::tie( distance, cur_p ) = point_to_triangle( point,
                    Geometry::Triangle3D{
                        vertices[Geometry::Tetra::tetra_facet_vertex[f][0]],
                        vertices[Geometry::Tetra::tetra_facet_vertex[f][1]],
                        vertices[Geometry::Tetra::tetra_facet_vertex[f][2]] } );
                if( distance < dist )
                {
                    dist = distance;
                    nearest_p = cur_p;
                }
            }
            return std::make_tuple( dist, nearest_p );
        }

        template <>
        std::tuple< double, vec3 > point_to_segment(
            const Geometry::Point3D& point, const Geometry::Segment3D& segment )
        {
            bool can_point_be_projected;
            vec3 nearest_p;
            std::tie( can_point_be_projected, nearest_p ) =
                point_segment_projection( point, segment.p0, segment.p1 );
            if( can_point_be_projected )
            {
                return std::make_tuple(
                    length( nearest_p - point ), nearest_p );
            }
            double p0_sq{ length2( segment.p0 - point ) };
            double p1_sq{ length2( segment.p1 - point ) };
            if( p0_sq < p1_sq )
            {
                return std::make_tuple( std::sqrt( p0_sq ), segment.p0 );
            }
            return std::make_tuple( std::sqrt( p1_sq ), segment.p1 );
        }

        std::tuple< double, vec3 > point_to_plane(
            const Geometry::Point3D& point, const Geometry::Plane& plane )
        {
            vec3 v{ point - plane.origin };
            double distance{ dot( v, plane.normal ) };
            vec3 projected_p{ point - distance * plane.normal };
            return std::make_tuple( distance, projected_p );
        }

        template std::tuple< double, vecn< 2 > > basic_api point_to_segment(
            const Geometry::Point< 2 >&, const Geometry::Segment< 2 >& );
        template std::tuple< double, vecn< 2 > > basic_api point_to_triangle(
            const Geometry::Point< 2 >&, const Geometry::Triangle< 2 >& );

        template std::tuple< double, vecn< 3 > > basic_api point_to_segment(
            const Geometry::Point< 3 >&, const Geometry::Segment< 3 >& );
        template std::tuple< double, vecn< 3 > > basic_api point_to_triangle(
            const Geometry::Point< 3 >&, const Geometry::Triangle< 3 >& );
    } // namespace Distance
} // namespace RINGMesh
