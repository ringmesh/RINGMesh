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

/*!
 * @file Basic geometrical requests
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace RINGMesh
{
    namespace Intersection
    {
        std::tuple< bool, std::vector< vec3 > > circle_plane(
            const Geometry::Circle& circle, const Geometry::Plane& plane )
        {
            bool does_plane_intersect_plane;
            Geometry::Line3D inter;
            std::tie( does_plane_intersect_plane, inter ) =
                plane_plane( plane, circle.plane );
            if( !does_plane_intersect_plane )
            {
                return std::make_tuple( false, std::vector< vec3 >() );
            }

            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // Locate one or two points that are on the circle and line.  If the
            // line is t*D+P, the circle center is C, and the circle radius is
            // r,
            // then r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2.
            // This
            // is a quadratic equation of the form:  a2*t^2 + 2*a1*t + a0 = 0.
            // In our case, a2 = 1 because direction plane is normalized
            vec3 diff{ inter.origin - circle.plane.origin };
            double a1{ dot( diff, inter.direction ) };
            double a0{ diff.length2() - circle.radius * circle.radius };

            double discr{ a1 * a1 - a0 };
            if( discr < 0.0 )
            {
                return std::make_tuple( false, std::vector< vec3 >() );
            }

            std::vector< vec3 > result;
            if( discr < global_epsilon )
            {
                result.emplace_back( inter.origin - a1 * inter.direction );
            }
            else
            {
                double root{ std::sqrt( discr ) };
                result.emplace_back(
                    inter.origin - ( a1 + root ) * inter.direction );
                result.emplace_back(
                    inter.origin - ( a1 - root ) * inter.direction );
            }
            return std::make_tuple( true, result );
        }

        std::tuple< bool, Geometry::Line3D > plane_plane(
            const Geometry::Plane& plane0, const Geometry::Plane& plane1 )
        {
            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // If N0 and N1 are parallel, either the planes are parallel and
            // separated
            // or the same plane.  In both cases, 'false' is returned.
            // Otherwise,
            // the intersection line is
            //   L(t) = t*Cross(N0,N1)/|Cross(N0,N1)| + c0*N0 + c1*N1
            // for some coefficients c0 and c1 and for t any real number (the
            // line
            // parameter).  Taking dot products with the normals,
            //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1) = c0 + c1*d
            //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1) = c0*d + c1
            // where d = Dot(N0,N1).  These are two equations in two unknowns.
            // The
            // solution is
            //   c0 = (d0 - d*d1)/det
            //   c1 = (d1 - d*d0)/det
            // where det = 1 - d^2.

            double norm_d{ dot( plane0.normal, plane1.normal ) };

            // Planes are parallel
            if( std::fabs( std::fabs( norm_d ) - 1 ) < global_epsilon )
            {
                return std::make_tuple( false, Geometry::Line3D{} );
            }

            double invDet{ 1.0 / ( 1.0 - norm_d * norm_d ) };
            double const_P0{ dot( plane0.normal, plane0.origin ) };
            double const_P1{ dot( plane1.normal, plane1.origin ) };
            double c0{ ( const_P0 - norm_d * const_P1 ) * invDet };
            double c1{ ( const_P1 - norm_d * const_P0 ) * invDet };
            vec3 O_inter{ c0 * plane0.normal + c1 * plane1.normal };
            vec3 D_inter{ cross( plane0.normal, plane1.normal ) };
            return std::make_tuple(
                true, Geometry::Line3D{ D_inter, O_inter } );
        }

        std::tuple< bool, vec2 > line_line(
            const Geometry::Line2D& line0, const Geometry::Line2D& line1 )
        {
            // The intersection of two lines is a solution to P0 + s0*D0 = P1 +
            // s1*D1.
            // Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If DotPerp(D0, D1))
            // = 0,
            // the lines are parallel.  Additionally, if DotPerp(Q, D1)) = 0,
            // the
            // lines are the same.  If Dotperp(D0, D1)) is not zero, then
            //   s0 = DotPerp(Q, D1))/DotPerp(D0, D1))
            // produces the point of intersection.  Also,
            //   s1 = DotPerp(Q, D0))/DotPerp(D0, D1))

            vec2 diff{ line1.origin - line0.origin };
            double D0DotPerpD1{ dot_perp( line0.direction, line1.direction ) };
            if( std::fabs( D0DotPerpD1 ) < global_epsilon )
            {
                // The lines are parallel.
                return std::make_tuple( false, vec2() );
            }

            double invD0DotPerpD1{ 1.0 / D0DotPerpD1 };
            double diffDotPerpD1{ dot_perp( diff, line1.direction ) };
            double s0{ diffDotPerpD1 * invD0DotPerpD1 };
            vec2 result{ line0.origin + s0 * line0.direction };
            return std::make_tuple( true, result );
        }

        std::tuple< bool, vec2 > segment_segment(
            const Geometry::Segment2D& segment0,
            const Geometry::Segment2D& segment1 )
        {
            bool does_segment_intersect_segment;
            vec2 line_intersection_result;
            std::tie(
                does_segment_intersect_segment, line_intersection_result ) =
                line_line( Geometry::Line2D{ segment0 },
                    Geometry::Line2D{ segment1 } );
            if( does_segment_intersect_segment )
            {
                // Test whether the line-line intersection is on the segments.
                Sign s0_seg0{ Position::point_side_to_segment(
                    segment0.p0, segment1 ) };
                Sign s1_seg0{ Position::point_side_to_segment(
                    segment0.p1, segment1 ) };
                if( s0_seg0 != ZERO && ( s0_seg0 == s1_seg0 ) )
                {
                    return std::make_tuple( false, vec2() );
                }
                Sign s0_seg1{ Position::point_side_to_segment(
                    segment1.p0, segment0 ) };
                Sign s1_seg1{ Position::point_side_to_segment(
                    segment1.p1, segment0 ) };
                if( s0_seg1 != ZERO && ( s0_seg1 == s1_seg1 ) )
                {
                    return std::make_tuple( false, vec2() );
                }

                if( ( s0_seg0 == ZERO || s1_seg0 == ZERO
                        || ( s0_seg0 != s1_seg0 ) )
                    && ( s0_seg1 == ZERO || s1_seg1 == ZERO
                           || ( s0_seg1 != s1_seg1 ) ) )
                {
                    return std::make_tuple( true, line_intersection_result );
                }
            }
            return std::make_tuple( false, vec2() );
        }

        std::tuple< bool, vec2 > segment_line(
            const Geometry::Segment2D& segment, const Geometry::Line2D& line )
        {
            bool does_segment_intersect_line;
            vec2 line_intersection_result;
            std::tie( does_segment_intersect_line, line_intersection_result ) =
                line_line( Geometry::Line2D{ segment }, line );
            if( does_segment_intersect_line )
            {
                // Test whether the line-line intersection is on the segment.
                Geometry::Segment2D line_segment{ { line_intersection_result
                                                      - line.direction },
                    { line_intersection_result + line.direction } };
                Sign s0{ Position::point_side_to_segment(
                    segment.p0, line_segment ) };
                Sign s1{ Position::point_side_to_segment(
                    segment.p1, line_segment ) };
                if( s0 == ZERO || s1 == ZERO || ( s0 != s1 ) )
                {
                    return std::make_tuple( true, line_intersection_result );
                }
            }
            return std::make_tuple( false, vec2() );
        }

        std::tuple< bool, vec3 > line_plane(
            const Geometry::Line3D& line, const Geometry::Plane& plane )
        {
            double dot_directions{ dot( line.direction, plane.normal ) };
            if( std::fabs( dot_directions ) > global_epsilon )
            {
                double signed_distance{ dot( plane.normal, line.origin )
                                        + plane.plane_constant() };
                vec3 result{ line.origin
                             - signed_distance * line.direction
                                   / dot_directions };
                return std::make_tuple( true, result );
            }
            // line is parallel to the plane
            return std::make_tuple( false, vec3() );
        }

        std::tuple< bool, vec3 > segment_plane(
            const Geometry::Segment3D& segment, const Geometry::Plane& plane )
        {
            bool does_line_intersect_plane;
            vec3 line_plane_result;
            std::tie( does_line_intersect_plane, line_plane_result ) =
                line_plane( Geometry::Line3D{ segment }, plane );
            if( does_line_intersect_plane )
            {
                if( Position::point_inside_segment(
                        line_plane_result, segment ) )
                {
                    // result inside the segment
                    return std::make_tuple( true, line_plane_result );
                }
                return std::make_tuple( false, vec3() );
            }
            return std::make_tuple( false, vec3() );
        }

        std::tuple< bool, vec3 > segment_disk(
            const Geometry::Segment3D& segment, const Geometry::Disk& disk )
        {
            bool does_segment_intersect_plane{ false };
            vec3 segment_plane_result;
            std::tie( does_segment_intersect_plane, segment_plane_result ) =
                segment_plane( segment, disk.plane );
            if( does_segment_intersect_plane
                && ( segment_plane_result - disk.plane.origin ).length()
                       <= disk.radius )
            {
                return std::make_tuple( true, segment_plane_result );
            }
            return std::make_tuple( false, vec3() );
        }

        std::tuple< bool, std::vector< vec3 > > triangle_circle(
            const Geometry::Triangle3D& triangle,
            const Geometry::Circle& circle )
        {
            bool does_circle_intersect_plane{ false };
            std::vector< vec3 > inter_circle_plane;
            std::tie( does_circle_intersect_plane, inter_circle_plane ) =
                circle_plane( circle, triangle.plane() );
            std::vector< vec3 > result;
            if( does_circle_intersect_plane )
            {
                for( const vec3& p : inter_circle_plane )
                {
                    if( Position::point_inside_triangle( p, triangle ) )
                    {
                        result.push_back( p );
                    }
                }
            }
            return std::make_tuple( !result.empty(), result );
        }

        std::tuple< bool, vec3 > segment_triangle(
            const Geometry::Segment3D& segment,
            const Geometry::Triangle3D& triangle )
        {
            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // Compute the offset origin, edges, and normal.
            vec3 seg_center{ segment.barycenter() };
            vec3 diff{ seg_center - triangle.p0 };
            vec3 edge1{ triangle.p1 - triangle.p0 };
            vec3 edge2{ triangle.p2 - triangle.p0 };
            vec3 normal{ cross( edge1, edge2 ) };

            // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
            // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
            //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
            //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
            //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
            vec3 D{ segment.direction() };
            double DdN{ dot( D, normal ) };
            signed_index_t sign;
            if( DdN > global_epsilon )
            {
                sign = 1;
            }
            else if( DdN < -global_epsilon )
            {
                sign = -1;
                DdN = -DdN;
            }
            else
            {
                // Segment and triangle are parallel, call it a "no
                // intersection"
                // even if the segment does intersect.
                return std::make_tuple( false, vec3() );
            }

            double DdQxE2{ sign * dot( D, cross( diff, edge2 ) ) };
            if( DdQxE2 >= 0 )
            {
                double DdE1xQ{ sign * dot( D, cross( edge1, diff ) ) };
                if( DdE1xQ >= 0 && DdQxE2 + DdE1xQ <= DdN )
                {
                    // Line intersects triangle, check if segment does.
                    double QdN{ -sign * dot( diff, normal ) };
                    double extDdN{ segment.length() * DdN / 2. };
                    if( -extDdN <= QdN && QdN <= extDdN )
                    {
                        // Segment intersects triangle.
                        double inv{ 1. / DdN };
                        double seg_parameter{ QdN * inv };

                        vec3 result{ seg_center + seg_parameter * D };
                        return std::make_tuple( true, result );
                    }
                    // else: |t| > extent, no intersection
                }
                // else: b1+b2 > 1, no intersection
                // else: b2 < 0, no intersection
            }
            // else: b1 < 0, no intersection
            return std::make_tuple( false, vec3() );
        }

        template < index_t DIMENSION >
        std::tuple< bool, std::vector< vecn< DIMENSION > > > line_sphere(
            const Geometry::Line< DIMENSION >& line,
            const Geometry::Sphere< DIMENSION >& sphere )
        {
            // The sphere is (X-C)^T*(X-C)-1 = 0 and the line is X = P+t*D.
            // Substitute the line equation into the sphere equation to obtain a
            // quadratic equation Q(t) = t^2 + 2*a1*t + a0 = 0, where a1 =
            // D^T*(P-C),
            // and a0 = (P-C)^T*(P-C)-1.
            vecn< DIMENSION > diff{ line.origin - sphere.origin };
            double a0{ dot( diff, diff ) - sphere.radius * sphere.radius };
            double a1{ dot( line.direction, diff ) };

            // Intersection occurs when Q(t) has real roots.
            double discr{ a1 * a1 - a0 };
            std::vector< vecn< DIMENSION > > results;
            if( discr > global_epsilon )
            {
                double root{ std::sqrt( discr ) };
                results.reserve( 2 );
                results.emplace_back(
                    line.origin + ( -a1 - root ) * line.direction );
                results.emplace_back(
                    line.origin + ( -a1 + root ) * line.direction );
            }
            else if( discr > -global_epsilon )
            {
                results.reserve( 1 );
                results.emplace_back( line.origin + -a1 * line.direction );
            }
            return std::make_tuple( !results.empty(), results );
        }

        template < index_t DIMENSION >
        std::tuple< bool, std::vector< vecn< DIMENSION > > > segment_sphere(
            const Geometry::Segment< DIMENSION >& segment,
            const Geometry::Sphere< DIMENSION >& sphere )
        {
            bool line_intersect;
            std::vector< vecn< DIMENSION > > line_intersections;
            std::tie( line_intersect, line_intersections ) =
                line_sphere( Geometry::Line< DIMENSION >{ segment }, sphere );

            std::vector< vecn< DIMENSION > > segment_intersections;
            if( line_intersect )
            {
                segment_intersections.reserve( line_intersections.size() );
                for( auto& point : line_intersections )
                {
                    if( Position::point_inside_segment( point, segment ) )
                    {
                        segment_intersections.emplace_back( point );
                    }
                }
            }
            return std::make_tuple(
                !segment_intersections.empty(), segment_intersections );
        }

        template std::tuple< bool, std::vector< vec2 > >
            basic_api line_sphere( const Geometry::Line2D& segment,
                const Geometry::Sphere2D& sphere );

        template std::tuple< bool, std::vector< vec2 > >
            basic_api segment_sphere( const Geometry::Segment2D& segment,
                const Geometry::Sphere2D& sphere );

        template std::tuple< bool, std::vector< vec3 > >
            basic_api line_sphere( const Geometry::Line3D& segment,
                const Geometry::Sphere3D& sphere );

        template std::tuple< bool, std::vector< vec3 > >
            basic_api segment_sphere( const Geometry::Segment3D& segment,
                const Geometry::Sphere3D& sphere );
    } // namespace Intersection
} // namespace RINGMesh
