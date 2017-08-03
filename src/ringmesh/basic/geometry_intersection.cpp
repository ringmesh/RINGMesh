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

#include <ringmesh/basic/geometry.h>

#include <geogram/numerics/predicates.h>

/*!
 * @file Basic geometrical requests 
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace RINGMesh {
    namespace Intersection {

        std::tuple< bool, std::vector< vec3 > > circle_plane(
            const vec3& O_plane,
            const vec3& N_plane,
            const vec3& O_circle,
            const vec3& N_circle,
            double r )
        {
            vec3 norm_N_plane { normalize( N_plane ) };
            vec3 norm_N_circle { normalize( N_circle ) };
            bool does_plane_intersect_plane { false };
            vec3 O_inter;
            vec3 D_inter;
            std::tie( does_plane_intersect_plane, O_inter, D_inter ) = plane_plane(
                O_plane, norm_N_plane, O_circle, norm_N_circle );
            if( !does_plane_intersect_plane ) {
                return std::make_tuple( false, std::vector< vec3 >() );
            }

            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // Locate one or two points that are on the circle and line.  If the
            // line is t*D+P, the circle center is C, and the circle radius is r,
            // then r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2.  This
            // is a quadratic equation of the form:  a2*t^2 + 2*a1*t + a0 = 0.
            vec3 diff { O_inter - O_circle };
            double a2 { D_inter.length2() };
            double a1 { dot( diff, D_inter ) };
            double a0 { diff.length2() - r * r };

            double discr { a1 * a1 - a0 * a2 };
            if( discr < 0.0 ) {
                return std::make_tuple( false, std::vector< vec3 >() );
            }

            if( std::fabs( a2 ) < global_epsilon ) {
                return std::make_tuple( false, std::vector< vec3 >() );
            }
            double inv { 1.0 / a2 };
            std::vector< vec3 > result;
            if( discr < global_epsilon ) {
                result.emplace_back( O_inter - ( a1 * inv ) * D_inter );
            } else {
                double root { std::sqrt( discr ) };
                result.emplace_back( O_inter - ( ( a1 + root ) * inv ) * D_inter );
                result.emplace_back( O_inter - ( ( a1 - root ) * inv ) * D_inter );
            }
            return std::make_tuple( true, result );
        }

        std::tuple< bool, vec3, vec3 > plane_plane(
            const vec3& O_P0,
            const vec3& N_P0,
            const vec3& O_P1,
            const vec3& N_P1 )
        {
            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // If N0 and N1 are parallel, either the planes are parallel and separated
            // or the same plane.  In both cases, 'false' is returned.  Otherwise,
            // the intersection line is
            //   L(t) = t*Cross(N0,N1)/|Cross(N0,N1)| + c0*N0 + c1*N1
            // for some coefficients c0 and c1 and for t any real number (the line
            // parameter).  Taking dot products with the normals,
            //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1) = c0 + c1*d
            //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1) = c0*d + c1
            // where d = Dot(N0,N1).  These are two equations in two unknowns.  The
            // solution is
            //   c0 = (d0 - d*d1)/det
            //   c1 = (d1 - d*d0)/det
            // where det = 1 - d^2.

            vec3 norm_N_P0 { normalize( N_P0 ) };
            vec3 norm_N_P1 { normalize( N_P1 ) };
            double norm_d { dot( norm_N_P0, norm_N_P1 ) };

            // Planes are parallel
            if( std::fabs( std::fabs( norm_d ) - 1 ) < global_epsilon ) {
                return std::make_tuple( false, vec3(), vec3() );
            }

            double invDet { 1.0 / ( 1.0 - norm_d * norm_d ) };
            double const_P0 { dot( norm_N_P0, O_P0 ) };
            double const_P1 { dot( norm_N_P1, O_P1 ) };
            double c0 { ( const_P0 - norm_d * const_P1 ) * invDet };
            double c1 { ( const_P1 - norm_d * const_P0 ) * invDet };
            vec3 O_inter { c0 * norm_N_P0 + c1 * norm_N_P1 };
            vec3 D_inter { cross( norm_N_P0, norm_N_P1 ) };
            return std::make_tuple( true, O_inter, D_inter );
        }

        std::tuple< bool, vec2 > line_line(
            const vec2& O_line0,
            const vec2& D_line0,
            const vec2& O_line1,
            const vec2& D_line1 )
        {
            // The intersection of two lines is a solution to P0 + s0*D0 = P1 + s1*D1.
            // Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If DotPerp(D0, D1)) = 0,
            // the lines are parallel.  Additionally, if DotPerp(Q, D1)) = 0, the
            // lines are the same.  If Dotperp(D0, D1)) is not zero, then
            //   s0 = DotPerp(Q, D1))/DotPerp(D0, D1))
            // produces the point of intersection.  Also,
            //   s1 = DotPerp(Q, D0))/DotPerp(D0, D1))

            vec2 norm_D_line0 { normalize( D_line0 ) };
            vec2 norm_D_line1 { normalize( D_line1 ) };
            vec2 diff { O_line1 - O_line0 };
            double D0DotPerpD1 { dot_perp( norm_D_line0, norm_D_line1 ) };
            if( std::fabs( D0DotPerpD1 ) < global_epsilon ) {
                // The lines are parallel.
                return std::make_tuple( false, vec2() );
            }

            double invD0DotPerpD1 { 1.0 / D0DotPerpD1 };
            double diffDotPerpD1 { dot_perp( diff, norm_D_line1 ) };
            double s0 { diffDotPerpD1 * invD0DotPerpD1 };
            vec2 result { O_line0 + s0 * norm_D_line0 };
            return std::make_tuple( true, result );
        }

        std::tuple< bool, vec2 > segment_segment(
            const vec2& p0_seg0,
            const vec2& p1_seg0,
            const vec2& p0_seg1,
            const vec2& p1_seg1 )
        {
            vec2 O_seg0 { ( p0_seg0 + p1_seg0 ) / 2. };
            vec2 D_seg0 { p1_seg0 - p0_seg0 };
            vec2 O_seg1 { ( p0_seg1 + p1_seg1 ) / 2. };
            vec2 D_seg1 { p1_seg1 - p0_seg1 };
            bool does_segment_intersect_segment { false };
            vec2 line_intersection_result;
            std::tie( does_segment_intersect_segment, line_intersection_result ) =
                line_line( O_seg0, D_seg0, O_seg1, D_seg1 );
            if( does_segment_intersect_segment ) {
                // Test whether the line-line intersection is on the segments.
                Sign s0_seg0 = sign(
                    GEO::PCK::orient_2d( p0_seg0, p0_seg1, p1_seg1 ) );
                Sign s1_seg0 = sign(
                    GEO::PCK::orient_2d( p1_seg0, p0_seg1, p1_seg1 ) );
                if( s0_seg0 != ZERO && ( s0_seg0 == s1_seg0 ) ) {
                    return std::make_tuple( false, vec2() );
                }
                Sign s0_seg1 = sign(
                    GEO::PCK::orient_2d( p0_seg1, p0_seg0, p1_seg0 ) );
                Sign s1_seg1 = sign(
                    GEO::PCK::orient_2d( p1_seg1, p0_seg0, p1_seg0 ) );
                if( s0_seg1 != ZERO && ( s0_seg1 == s1_seg1 ) ) {
                    return std::make_tuple( false, vec2() );
                }

                if( s0_seg0 == ZERO || s1_seg0 == ZERO || ( s0_seg0 != s1_seg0 ) ) {
                    if( s0_seg1 == ZERO || s1_seg1 == ZERO
                        || ( s0_seg1 != s1_seg1 ) ) {
                        return std::make_tuple( true, line_intersection_result );
                    }
                }
            }
            return std::make_tuple( false, vec2() );
        }

        std::tuple< bool, vec2 > segment_line(
            const vec2& p0_seg,
            const vec2& p1_seg,
            const vec2& O_line,
            const vec2& D_line )
        {
            vec2 O_seg { ( p0_seg + p1_seg ) / 2. };
            vec2 D_seg { p1_seg - p0_seg };
            bool does_segment_intersect_line { false };
            vec2 line_intersection_result;
            std::tie( does_segment_intersect_line, line_intersection_result ) =
                line_line( O_seg, D_seg, O_line, D_line );
            if( does_segment_intersect_line ) {
                // Test whether the line-line intersection is on the segment.
                vec2 minus_direction { line_intersection_result - D_line };
                vec2 plus_direction { line_intersection_result + D_line };
                Sign s0 = sign(
                    GEO::PCK::orient_2d( p0_seg, minus_direction, plus_direction ) );
                Sign s1 = sign(
                    GEO::PCK::orient_2d( p1_seg, minus_direction, plus_direction ) );
                if( s0 == ZERO || s1 == ZERO || ( s0 != s1 ) ) {
                    return std::make_tuple( true, line_intersection_result );
                }
            }
            return std::make_tuple( false, vec2() );
        }

        std::tuple< bool, vec3 > line_plane(
            const vec3& O_line,
            const vec3& D_line,
            const vec3& O_plane,
            const vec3& N_plane )
        {
            double dot_directions { dot( D_line, N_plane ) };
            if( std::fabs( dot_directions ) > global_epsilon ) {
                double plane_constant { 0.0 };
                for( index_t i : range( 3 ) ) {
                    plane_constant += O_plane[i] * N_plane[i];
                }
                double signed_distance { dot( N_plane, O_line ) - plane_constant };
                vec3 result { O_line - signed_distance * D_line / dot_directions };
                return std::make_tuple( true, result );
            } else {
                // line is parallel to the plane
                return std::make_tuple( false, vec3() );
            }
        }

        std::tuple< bool, vec3 > segment_plane(
            const vec3& seg0,
            const vec3& seg1,
            const vec3& O_plane,
            const vec3& N_plane )
        {
            vec3 segment_direction { normalize( seg1 - seg0 ) };
            vec3 segment_barycenter { 0.5 * ( seg0 + seg1 ) };
            bool does_line_intersect_plane { false };
            vec3 line_plane_result;
            std::tie( does_line_intersect_plane, line_plane_result ) = line_plane(
                segment_barycenter, segment_direction, O_plane, N_plane );
            if( does_line_intersect_plane ) {
                if( ( line_plane_result - segment_barycenter ).length2()
                    > ( seg0 - segment_barycenter ).length2() + global_epsilon ) {
                    // result outside the segment
                    return std::make_tuple( false, vec3() );
                } else {
                    return std::make_tuple( true, line_plane_result );
                }
            } else {
                return std::make_tuple( false, vec3() );
            }
        }

        std::tuple< bool, vec3 > disk_segment(
            const vec3& p0,
            const vec3& p1,
            const vec3& O_circle,
            const vec3& N_circle,
            double r )
        {
            bool does_segment_intersect_plane { false };
            vec3 segment_plane_result;
            std::tie( does_segment_intersect_plane, segment_plane_result ) =
                segment_plane( p0, p1, O_circle, N_circle );
            if( does_segment_intersect_plane ) {
                if( ( segment_plane_result - O_circle ).length() <= r ) {
                    return std::make_tuple( true, segment_plane_result );
                }
            }
            return std::make_tuple( false, vec3() );
        }

        std::tuple< bool, std::vector< vec3 > > circle_triangle(
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            const vec3& O_circle,
            const vec3& N_circle,
            double r )
        {
            vec3 N_triangle { normalize( cross( p1 - p0, p2 - p0 ) ) };
            vec3 barycenter { ( p0 + p1 + p2 ) / 3. };
            bool does_circle_intersect_plane { false };
            std::vector< vec3 > inter_circle_plane;
            std::tie( does_circle_intersect_plane, inter_circle_plane ) =
                circle_plane( barycenter, N_triangle, O_circle, N_circle, r );
            std::vector< vec3 > result;
            if( does_circle_intersect_plane ) {
                for( const vec3& p : inter_circle_plane ) {
                    if( point_inside_triangle( p, p0, p1, p2 ) ) {
                        result.push_back( p );
                    }
                }
            }
            return std::make_tuple( !result.empty(), result );
        }

        std::tuple< bool, vec3 > segment_triangle(
            const vec3& seg0,
            const vec3& seg1,
            const vec3& trgl0,
            const vec3& trgl1,
            const vec3& trgl2 )
        {
            // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
            // Compute the offset origin, edges, and normal.
            vec3 seg_center { ( seg0 + seg1 ) / 2 };
            vec3 diff { seg_center - trgl0 };
            vec3 edge1 { trgl1 - trgl0 };
            vec3 edge2 { trgl2 - trgl0 };
            vec3 normal { cross( edge1, edge2 ) };

            // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
            // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
            //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
            //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
            //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
            vec3 D { normalize( seg1 - seg0 ) };
            double DdN { dot( D, normal ) };
            signed_index_t sign;
            if( DdN > global_epsilon ) {
                sign = 1;
            } else if( DdN < -global_epsilon ) {
                sign = -1;
                DdN = -DdN;
            } else {
                // Segment and triangle are parallel, call it a "no intersection"
                // even if the segment does intersect.
                return std::make_tuple( false, vec3() );
            }

            double DdQxE2 { sign * dot( D, cross( diff, edge2 ) ) };
            if( DdQxE2 >= 0 ) {
                double DdE1xQ { sign * dot( D, cross( edge1, diff ) ) };
                if( DdE1xQ >= 0 ) {
                    if( DdQxE2 + DdE1xQ <= DdN ) {
                        // Line intersects triangle, check if segment does.
                        double QdN { -sign * dot( diff, normal ) };
                        double extDdN { length( seg1 - seg0 ) * DdN / 2. };
                        if( -extDdN <= QdN && QdN <= extDdN ) {
                            // Segment intersects triangle.
                            double inv { 1. / DdN };
                            double seg_parameter { QdN * inv };

                            vec3 result { seg_center + seg_parameter * D };
                            return std::make_tuple( true, result );
                        }
                        // else: |t| > extent, no intersection
                    }
                    // else: b1+b2 > 1, no intersection
                }
                // else: b2 < 0, no intersection
            }
            // else: b1 < 0, no intersection
            return std::make_tuple( false, vec3() );
        }

        std::tuple< bool, std::vector< vec3 > > line_sphere(
            const vec3& O_line,
            const vec3& D_line,
            const vec3& O_sphere,
            double radius )
        {
            vec3 D_line_normalized { normalize( D_line ) };
            // The sphere is (X-C)^T*(X-C)-1 = 0 and the line is X = P+t*D.
            // Substitute the line equation into the sphere equation to obtain a
            // quadratic equation Q(t) = t^2 + 2*a1*t + a0 = 0, where a1 = D^T*(P-C),
            // and a0 = (P-C)^T*(P-C)-1.
            vec3 diff { O_line - O_sphere };
            double a0 { dot( diff, diff ) - radius * radius };
            double a1 { dot( D_line_normalized, diff ) };

            // Intersection occurs when Q(t) has real roots.
            double discr { a1 * a1 - a0 };
            std::vector< vec3 > results;
            if( discr > global_epsilon ) {
                double root { std::sqrt( discr ) };
                results.reserve( 2 );
                results.emplace_back( O_line + ( -a1 - root ) * D_line_normalized );
                results.emplace_back( O_line + ( -a1 + root ) * D_line_normalized );
            } else if( discr > -global_epsilon ) {
                results.reserve( 1 );
                results.emplace_back( O_line + -a1 * D_line_normalized );
            }
            return std::make_tuple( !results.empty(), results );
        }

    }
}
