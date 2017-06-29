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

#include <geogram/mesh/mesh.h>

#include <geogram/numerics/predicates.h>

/*!
 * @file Basic geometrical requests 
 * @author Arnaud Botella
 *
 * @todo Comment on the robustness of the tests
 */

namespace {
    using namespace RINGMesh;

    bool is_almost_zero( const double& value )
    {
        return value < global_epsilon && value > -global_epsilon;
    }

}

namespace RINGMesh {

    double dot_perp( const vec2& v0, const vec2& v1 )
    {
        return dot( v0, vec2( v1.y, -v1.x ) );
    }

    double triangle_signed_area(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& triangle_normal )
    {
        double area = GEO::Geom::triangle_area( p0, p1, p2 );
        vec3 area_normal = cross( p0 - p2, p1 - p2 );
        if( dot( triangle_normal, area_normal ) < 0 ) {
            area = -area;
        }
        return area;
    }

    bool operator==( const vec3& u, const vec3& v )
    {
        return u.x == v.x && u.y == v.y && u.z == v.z;
    }

    bool operator!=( const vec3& u, const vec3& v )
    {
        return u.x != v.x || u.y != v.y || u.z != v.z;
    }

    double point_triangle_distance(
        const vec3& point,
        const vec3& V0,
        const vec3& V1,
        const vec3& V2,
        vec3& closest_point )
    {
        vec3 diff = V0 - point;
        vec3 edge0 = V1 - V0;
        vec3 edge1 = V2 - V0;
        double a00 = length2( edge0 );
        double a01 = dot( edge0, edge1 );
        double a11 = length2( edge1 );
        double b0 = dot( diff, edge0 );
        double b1 = dot( diff, edge1 );
        double c = length2( diff );
        double det = std::fabs( a00 * a11 - a01 * a01 );
        double s = a01 * b1 - a11 * b0;
        double t = a01 * b0 - a00 * b1;
        double sqrDistance;

        if( s + t <= det ) {
            if( s < 0.0 ) {
                if( t < 0.0 ) { // region 4
                    if( b0 < 0.0 ) {
                        t = 0.0;
                        if( -b0 >= a00 ) {
                            s = 1.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        } else {
                            s = -b0 / a00;
                            sqrDistance = b0 * s + c;
                        }
                    } else {
                        s = 0.0;
                        if( b1 >= 0.0 ) {
                            t = 0.0;
                            sqrDistance = c;
                        } else if( -b1 >= a11 ) {
                            t = 1.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        } else {
                            t = -b1 / a11;
                            sqrDistance = b1 * t + c;
                        }
                    }
                } else { // region 3
                    s = 0.0;
                    if( b1 >= 0.0 ) {
                        t = 0.0;
                        sqrDistance = c;
                    } else if( -b1 >= a11 ) {
                        t = 1.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    } else {
                        t = -b1 / a11;
                        sqrDistance = b1 * t + c;
                    }
                }
            } else if( t < 0.0 ) { // region 5
                t = 0.0;
                if( b0 >= 0.0 ) {
                    s = 0.0;
                    sqrDistance = c;
                } else if( -b0 >= a00 ) {
                    s = 1.0;
                    sqrDistance = a00 + 2.0 * b0 + c;
                } else {
                    s = -b0 / a00;
                    sqrDistance = b0 * s + c;
                }
            } else { // region 0
                // minimum at interior point
                double invDet = double( 1.0 ) / det;
                s *= invDet;
                t *= invDet;
                sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                    + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
            }
        } else {
            double tmp0, tmp1, numer, denom;

            if( s < 0.0 ) { // region 2
                tmp0 = a01 + b0;
                tmp1 = a11 + b1;
                if( tmp1 > tmp0 ) {
                    numer = tmp1 - tmp0;
                    denom = a00 - 2.0 * a01 + a11;
                    if( numer >= denom ) {
                        s = 1.0;
                        t = 0.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    } else {
                        s = numer / denom;
                        t = 1.0 - s;
                        sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                            + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
                    }
                } else {
                    s = 0.0;
                    if( tmp1 <= 0.0 ) {
                        t = 1.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    } else if( b1 >= 0.0 ) {
                        t = 0.0;
                        sqrDistance = c;
                    } else {
                        t = -b1 / a11;
                        sqrDistance = b1 * t + c;
                    }
                }
            } else if( t < 0.0 ) { // region 6
                tmp0 = a01 + b1;
                tmp1 = a00 + b0;
                if( tmp1 > tmp0 ) {
                    numer = tmp1 - tmp0;
                    denom = a00 - 2.0 * a01 + a11;
                    if( numer >= denom ) {
                        t = 1.0;
                        s = 0.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    } else {
                        t = numer / denom;
                        s = 1.0 - t;
                        sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                            + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
                    }
                } else {
                    t = 0.0;
                    if( tmp1 <= 0.0 ) {
                        s = 1.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    } else if( b0 >= 0.0 ) {
                        s = 0.0;
                        sqrDistance = c;
                    } else {
                        s = -b0 / a00;
                        sqrDistance = b0 * s + c;
                    }
                }
            } else { // region 1
                numer = a11 + b1 - a01 - b0;
                if( numer <= 0.0 ) {
                    s = 0.0;
                    t = 1.0;
                    sqrDistance = a11 + 2.0 * b1 + c;
                } else {
                    denom = a00 - 2.0 * a01 + a11;
                    if( numer >= denom ) {
                        s = 1.0;
                        t = 0.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    } else {
                        s = numer / denom;
                        t = 1.0 - s;
                        sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 )
                            + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
                    }
                }
            }
        }

        // Account for numerical round-off error.
        if( sqrDistance < 0.0 ) {
            sqrDistance = 0.0;
        }

        closest_point = V0 + s * edge0 + t * edge1;
        return std::sqrt( sqrDistance );
    }

    double point_segment_distance(
        const vec2& p,
        const vec2& p0,
        const vec2& p1,
        vec2& nearest_p )
    {
        // The direction vector is not unit length.  The normalization is deferred
        // until it is needed.
        vec2 direction = p1 - p0;
        vec2 diff = p - p1;
        double t = dot( direction, diff );
        if( t >= global_epsilon ) {
            nearest_p = p1;
        } else {
            diff = p - p0;
            t = dot( direction, diff );
            if( t <= global_epsilon ) {
                nearest_p = p0;
            } else {
                double sqrLength = dot( direction, direction );
                if( sqrLength > global_epsilon ) {
                    t /= sqrLength;
                    nearest_p = p0 + t * direction;
                } else {
                    nearest_p = p0;
                }
            }
        }

        diff = p - nearest_p;
        return std::sqrt( dot( diff, diff ) );
    }

    double point_triangle_distance(
        const vec2& point,
        const vec2& V0,
        const vec2& V1,
        const vec2& V2,
        vec2& closest_point )
    {
        double result = max_float64();
        if( point_inside_triangle( point, V0, V1, V2 ) ) {
            closest_point = point;
            result = 0.0;
        } else {
            vec2 closest[3];
            double distance[3];
            distance[0] = point_segment_distance( point, V0, V1, closest[0] );
            distance[1] = point_segment_distance( point, V1, V2, closest[1] );
            distance[2] = point_segment_distance( point, V2, V0, closest[2] );
            if( distance[0] < distance[1] ) {
                if( distance[0] < distance[2] ) {
                    result = distance[0];
                    closest_point = closest[0];
                } else {
                    result = distance[2];
                    closest_point = closest[2];
                }
            } else {
                if( distance[1] < distance[2] ) {
                    result = distance[1];
                    closest_point = closest[1];
                } else {
                    result = distance[2];
                    closest_point = closest[2];
                }
            }
        }
        return result;
    }

    double point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        vec3 vertices[4] = { p0, p1, p2, p3 };
        double dist = max_float64();
        for( index_t f = 0; f < GEO::MeshCellDescriptors::tet_descriptor.nb_facets;
            f++ ) {
            vec3 cur_p;
            double distance =
                point_triangle_distance( p,
                    vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][0]],
                    vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][1]],
                    vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][2]],
                    cur_p );
            if( distance < dist ) {
                dist = distance;
                nearest_p = cur_p;
            }
        }
        return dist;
    }

    double point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p )
    {
        vec3 vertices[5] = { p0, p1, p2, p3, p4 };
        double dist = max_float64();
        for( index_t f = 0;
            f < GEO::MeshCellDescriptors::pyramid_descriptor.nb_facets; f++ ) {
            vec3 cur_p;
            double distance = max_float64();
            index_t nb_vertices =
                GEO::MeshCellDescriptors::pyramid_descriptor.nb_vertices_in_facet[f];
            if( nb_vertices == 3 ) {
                distance =
                    point_triangle_distance( p,
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][0]],
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][1]],
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][2]],
                        cur_p );
            } else if( nb_vertices == 4 ) {
                distance =
                    point_quad_distance( p,
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][0]],
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][1]],
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][2]],
                        vertices[GEO::MeshCellDescriptors::pyramid_descriptor.facet_vertex[f][3]],
                        cur_p );
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance;
                nearest_p = cur_p;
            }
        }
        return dist;
    }

    double point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p )
    {
        vec3 vertices[6] = { p0, p1, p2, p3, p4, p5 };
        double dist = max_float64();
        for( index_t f = 0; f < GEO::MeshCellDescriptors::prism_descriptor.nb_facets;
            f++ ) {
            vec3 cur_p;
            double distance = max_float64();
            index_t nb_vertices =
                GEO::MeshCellDescriptors::prism_descriptor.nb_vertices_in_facet[f];
            if( nb_vertices == 3 ) {
                distance =
                    point_triangle_distance( p,
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][0]],
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][1]],
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][2]],
                        cur_p );
            } else if( nb_vertices == 4 ) {
                distance =
                    point_quad_distance( p,
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][0]],
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][1]],
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][2]],
                        vertices[GEO::MeshCellDescriptors::prism_descriptor.facet_vertex[f][3]],
                        cur_p );
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance;
                nearest_p = cur_p;
            }
        }
        return dist;
    }

    double point_hexa_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7,
        vec3& nearest_p )
    {
        /// Review: Why not input an array ?
        vec3 vertices[8] = { p0, p1, p2, p3, p4, p5, p6, p7 };
        double dist = max_float64();
        for( index_t f = 0; f < GEO::MeshCellDescriptors::hex_descriptor.nb_facets;
            f++ ) {
            vec3 cur_p;
            double distance =
                point_quad_distance( p,
                    vertices[GEO::MeshCellDescriptors::hex_descriptor.facet_vertex[f][0]],
                    vertices[GEO::MeshCellDescriptors::hex_descriptor.facet_vertex[f][1]],
                    vertices[GEO::MeshCellDescriptors::hex_descriptor.facet_vertex[f][2]],
                    vertices[GEO::MeshCellDescriptors::hex_descriptor.facet_vertex[f][3]],
                    cur_p );
            if( distance < dist ) {
                dist = distance;
                nearest_p = cur_p;
            }
        }
        return dist;
    }

    bool point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        bool exact_predicates )
    {
        vec3 vertices[4] = { p0, p1, p2, p3 };
        Sign signs[4];
        if( !exact_predicates ) {
            for( index_t f = 0;
                f < GEO::MeshCellDescriptors::tet_descriptor.nb_facets; f++ ) {
                double volume =
                    GEO::Geom::tetra_signed_volume( p,
                        vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][0]],
                        vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][1]],
                        vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][2]] );
                if( is_almost_zero( volume ) ) {
                    return point_inside_tetra( p, p0, p1, p2, p3, true );
                }
                signs[f] = sign( volume );
            }
        } else {
            for( index_t f = 0;
                f < GEO::MeshCellDescriptors::tet_descriptor.nb_facets; f++ ) {
                signs[f] =
                    sign(
                        GEO::PCK::orient_3d( p.data(),
                            vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][0]].data(),
                            vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][1]].data(),
                            vertices[GEO::MeshCellDescriptors::tet_descriptor.facet_vertex[f][2]].data() ) );
            }
        }
        return ( signs[0] >= 0 && signs[1] >= 0 && signs[2] >= 0 && signs[3] >= 0 )
            || ( signs[0] <= 0 && signs[1] <= 0 && signs[2] <= 0 && signs[3] <= 0 );
    }

    bool circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result )
    {
        vec3 O_inter, D_inter;
        vec3 norm_N_plane = normalize( N_plane );
        vec3 norm_N_circle = normalize( N_circle );
        if( !plane_plane_intersection( O_plane, norm_N_plane, O_circle,
            norm_N_circle, O_inter, D_inter ) ) {
            return false;
        }

        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Locate one or two points that are on the circle and line.  If the
        // line is t*D+P, the circle center is C, and the circle radius is r,
        // then r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2.  This
        // is a quadratic equation of the form:  a2*t^2 + 2*a1*t + a0 = 0.
        vec3 diff = O_inter - O_circle;
        double a2 = D_inter.length2();
        double a1 = dot( diff, D_inter );
        double a0 = diff.length2() - r * r;

        double discr = a1 * a1 - a0 * a2;
        if( discr < 0.0 ) {
            return false;
        }

        if( std::fabs( a2 ) < global_epsilon ) {
            return false;
        }
        double inv = 1.0 / a2;
        if( discr < global_epsilon ) {
            result.emplace_back( O_inter - ( a1 * inv ) * D_inter );
        } else {
            double root = sqrt( discr );
            result.emplace_back( O_inter - ( ( a1 + root ) * inv ) * D_inter );
            result.emplace_back( O_inter - ( ( a1 - root ) * inv ) * D_inter );
        }
        return true;
    }

    bool plane_plane_intersection(
        const vec3& O_P0,
        const vec3& N_P0,
        const vec3& O_P1,
        const vec3& N_P1,
        vec3& O_inter,
        vec3& D_inter )
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

        vec3 norm_N_P0 = normalize( N_P0 );
        vec3 norm_N_P1 = normalize( N_P1 );
        double norm_d = dot( norm_N_P0, norm_N_P1 );

        // Planes are parallel
        if( std::fabs( std::fabs( norm_d ) - 1 ) < global_epsilon ) {
            return false;
        }

        double invDet = 1.0 / ( 1.0 - norm_d * norm_d );
        double const_P0 = dot( norm_N_P0, O_P0 );
        double const_P1 = dot( norm_N_P1, O_P1 );
        double c0 = ( const_P0 - norm_d * const_P1 ) * invDet;
        double c1 = ( const_P1 - norm_d * const_P0 ) * invDet;
        O_inter = c0 * norm_N_P0 + c1 * norm_N_P1;
        D_inter = cross( norm_N_P0, norm_N_P1 );
        return true;
    }

    bool line_line_intersection(
        const vec2& O_line0,
        const vec2& D_line0,
        const vec2& O_line1,
        const vec2& D_line1,
        vec2& result )
    {
        // The intersection of two lines is a solution to P0 + s0*D0 = P1 + s1*D1.
        // Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If DotPerp(D0, D1)) = 0,
        // the lines are parallel.  Additionally, if DotPerp(Q, D1)) = 0, the
        // lines are the same.  If Dotperp(D0, D1)) is not zero, then
        //   s0 = DotPerp(Q, D1))/DotPerp(D0, D1))
        // produces the point of intersection.  Also,
        //   s1 = DotPerp(Q, D0))/DotPerp(D0, D1))

        vec2 norm_D_line0 = normalize( D_line0 );
        vec2 norm_D_line1 = normalize( D_line1 );
        vec2 diff = O_line1 - O_line0;
        double D0DotPerpD1 = dot_perp( norm_D_line0, norm_D_line1 );
        if( std::fabs( D0DotPerpD1 ) < global_epsilon ) {
            // The lines are parallel.
            return false;
        }

        double invD0DotPerpD1 = 1.0 / D0DotPerpD1;
        double diffDotPerpD1 = dot_perp( diff, norm_D_line1 );
        double s0 = diffDotPerpD1 * invD0DotPerpD1;
        result = O_line0 + s0 * norm_D_line0;
        return true;
    }

    bool segment_segment_intersection(
        const vec2& p0_seg0,
        const vec2& p1_seg0,
        const vec2& p0_seg1,
        const vec2& p1_seg1,
        vec2& result )
    {
        vec2 O_seg0( ( p0_seg0 + p1_seg0 ) / 2. );
        vec2 D_seg0( p1_seg0 - p0_seg0 );
        vec2 O_seg1( ( p0_seg1 + p1_seg1 ) / 2. );
        vec2 D_seg1( p1_seg1 - p0_seg1 );
        vec2 line_intersection_result;
        if( line_line_intersection( O_seg0, D_seg0, O_seg1, D_seg1,
            line_intersection_result ) ) {
            // Test whether the line-line intersection is on the segments.
            if( length( line_intersection_result - O_seg0 )
                <= 0.5 * D_seg0.length() + global_epsilon
                && length( line_intersection_result - O_seg1 )
                    <= 0.5 * D_seg1.length() + global_epsilon ) {
                result = line_intersection_result;
                return true;
            }
        }
        return false;
    }

    bool tetra_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        double lambda[4] )
    {
        double total_volume = GEO::Geom::tetra_signed_volume( p0, p1, p2, p3 );
        if( total_volume < global_epsilon_3 ) {
            for( index_t i = 0; i < 4; i++ ) {
                lambda[i] = 0;
            }
            return false;
        }
        double volume0 = GEO::Geom::tetra_signed_volume( p1, p3, p2, p );
        double volume1 = GEO::Geom::tetra_signed_volume( p0, p2, p3, p );
        double volume2 = GEO::Geom::tetra_signed_volume( p0, p3, p1, p );
        double volume3 = GEO::Geom::tetra_signed_volume( p0, p1, p2, p );

        lambda[0] = volume0 / total_volume;
        lambda[1] = volume1 / total_volume;
        lambda[2] = volume2 / total_volume;
        lambda[3] = volume3 / total_volume;
        return true;
    }

    bool triangle_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        double lambda[3] )
    {
        double total_area = GEO::Geom::triangle_area( p0, p1, p2 );
        if( total_area < global_epsilon_sq ) {
            for( index_t i = 0; i < 3; i++ ) {
                lambda[i] = 0;
            }
            return false;
        }
        vec3 triangle_normal = cross( p2 - p0, p1 - p0 );
        double area0 = triangle_signed_area( p2, p1, p, triangle_normal );
        double area1 = triangle_signed_area( p0, p2, p, triangle_normal );
        double area2 = triangle_signed_area( p1, p0, p, triangle_normal );

        lambda[0] = area0 / total_area;
        lambda[1] = area1 / total_area;
        lambda[2] = area2 / total_area;
        return true;
    }

    bool triangle_barycentric_coordinates(
        const vec2& p,
        const vec2& p0,
        const vec2& p1,
        const vec2& p2,
        double lambda[3] )
    {
        double total_area = GEO::Geom::triangle_signed_area( p2, p1, p0 );
        if( std::fabs( total_area ) < global_epsilon_sq ) {
            for( index_t i = 0; i < 3; i++ ) {
                lambda[i] = 0;
            }
            return false;
        }
        double area0 = GEO::Geom::triangle_signed_area( p2, p1, p );
        double area1 = GEO::Geom::triangle_signed_area( p0, p2, p );
        double area2 = GEO::Geom::triangle_signed_area( p1, p0, p );

        lambda[0] = area0 / total_area;
        lambda[1] = area1 / total_area;
        lambda[2] = area2 / total_area;
        return true;
    }

    bool line_plane_intersection(
        const vec3& O_line,
        const vec3& D_line,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result )
    {
        double dot_directions = dot( D_line, N_plane );
        if( std::fabs( dot_directions ) > global_epsilon ) {
            double plane_constant = 0.0;
            for( index_t i = 0; i < 3; i++ ) {
                plane_constant += O_plane[i] * N_plane[i];
            }
            double signed_distance = dot( N_plane, O_line ) - plane_constant;
            result = O_line - signed_distance * D_line / dot_directions;
            return true;
        } else {
            // line is parallel to the plane
            return false;
        }
    }

    bool segment_plane_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result )
    {
        vec3 segment_direction = normalize( seg1 - seg0 );
        vec3 segment_barycenter = 0.5 * ( seg0 + seg1 );
        vec3 line_plane_result;
        if( line_plane_intersection( segment_barycenter, segment_direction, O_plane,
            N_plane, line_plane_result ) ) {
            if( ( line_plane_result - segment_barycenter ).length2()
                > ( seg0 - segment_barycenter ).length2() + global_epsilon ) {
                // result outside the segment
                return false;
            } else {
                result = line_plane_result;
                return true;
            }
        } else {
            return false;
        }
    }

    bool disk_segment_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        vec3& result )
    {
        vec3 segment_plane_result;
        if( segment_plane_intersection( p0, p1, O_circle, N_circle,
            segment_plane_result ) ) {
            if( ( segment_plane_result - O_circle ).length() <= r ) {
                result = segment_plane_result;
                return true;
            }
        }
        return false;
    }

    bool circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result )
    {
        vec3 N_triangle = normalize( cross( p1 - p0, p2 - p0 ) );
        vec3 barycenter = ( p0 + p1 + p2 ) / 3;
        std::vector< vec3 > inter_circle_plane;
        if( circle_plane_intersection( barycenter, N_triangle, O_circle, N_circle, r,
            inter_circle_plane ) ) {
            for( const vec3& p : inter_circle_plane ) {
                if( point_inside_triangle( p, p0, p1, p2 ) ) {
                    result.push_back( p );
                }
            }
        }
        return !result.empty();
    }

    template< index_t DIMENSION >
    bool point_segment_projection(
        const vecn< DIMENSION >& p,
        const vecn< DIMENSION >& p0,
        const vecn< DIMENSION >& p1,
        vecn< DIMENSION >& new_p )
    {
        vecn< DIMENSION > center = ( p0 + p1 ) * 0.5;
        vecn< DIMENSION > diff = p - center;
        vecn< DIMENSION > edge = p1 - p0;
        double extent = 0.5 * edge.length();
        edge = normalize( edge );
        double d = dot( edge, diff );

        if( std::fabs( d ) <= extent ) {
            new_p = center + d * edge;
            return true;
        }
        return false;
    }

    void point_plane_projection(
        const vec3& p,
        const vec3& N_plane,
        const vec3& O_plane,
        vec3& projected_p )
    {
        vec3 N_unit_plane = normalize( N_plane );
        vec3 v( p - O_plane );
        double distance = dot( v, N_unit_plane );
        projected_p = p - distance * N_unit_plane;
    }

    double point_segment_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& nearest_p )
    {
        if( point_segment_projection( p, p0, p1, nearest_p ) ) {
            return length( nearest_p - p );
        } else {
            double p0_distance_sq = length2( p0 - p );
            double p1_distance_sq = length2( p1 - p );
            if( p0_distance_sq < p1_distance_sq ) {
                nearest_p = p0;
                return std::sqrt( p0_distance_sq );
            } else {
                nearest_p = p1;
                return std::sqrt( p1_distance_sq );
            }
        }
    }

    double point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        const vec3 center( ( p0 + p1 + p2 + p3 ) * 0.25 );
        vec3 edge0( p1 - p0 );
        vec3 edge1( p3 - p0 );
        vec3 axis[2] = { normalize( edge0 ), normalize( edge1 ) };
        double extent[2] = { 0.5 * edge0.length(), 0.5 * edge1.length() };

        vec3 diff = center - p;
        double b0 = dot( diff, axis[0] );
        double b1 = dot( diff, axis[1] );
        double s0 = -b0;
        double s1 = -b1;
        double sqrDistance = dot( diff, diff );

        if( s0 < -extent[0] ) {
            s0 = -extent[0];
        } else if( s0 > extent[0] ) {
            s0 = extent[0];
        }
        sqrDistance += s0 * ( s0 + 2. * b0 );

        if( s1 < -extent[1] ) {
            s1 = -extent[1];
        } else if( s1 > extent[1] ) {
            s1 = extent[1];
        }
        sqrDistance += s1 * ( s1 + 2. * b1 );

        // Account for numerical round-off error.
        if( sqrDistance < 0 ) {
            sqrDistance = 0;
        }

        double distance = std::sqrt( sqrDistance );
        nearest_p = center;
        nearest_p += s0 * axis[0];
        nearest_p += s1 * axis[1];

        return distance;
    }

    bool segment_triangle_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& trgl0,
        const vec3& trgl1,
        const vec3& trgl2,
        vec3& result )
    {
        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Compute the offset origin, edges, and normal.
        vec3 seg_center = ( seg0 + seg1 ) / 2;
        vec3 diff = seg_center - trgl0;
        vec3 edge1 = trgl1 - trgl0;
        vec3 edge2 = trgl2 - trgl0;
        vec3 normal = cross( edge1, edge2 );

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        vec3 D = normalize( seg1 - seg0 );
        double DdN = dot( D, normal );
        signed_index_t sign;
        if( DdN > global_epsilon ) {
            sign = 1;
        } else if( DdN < -global_epsilon ) {
            sign = -1;
            DdN = -DdN;
        } else {
            // Segment and triangle are parallel, call it a "no intersection"
            // even if the segment does intersect.
            return false;
        }

        double DdQxE2 = sign * dot( D, cross( diff, edge2 ) );
        if( DdQxE2 >= 0 ) {
            double DdE1xQ = sign * dot( D, cross( edge1, diff ) );
            if( DdE1xQ >= 0 ) {
                if( DdQxE2 + DdE1xQ <= DdN ) {
                    // Line intersects triangle, check if segment does.
                    double QdN = -sign * dot( diff, normal );
                    double extDdN = length( seg1 - seg0 ) * DdN / 2.;
                    if( -extDdN <= QdN && QdN <= extDdN ) {
                        // Segment intersects triangle.
                        double inv = 1. / DdN;
                        double seg_parameter = QdN * inv;

                        result = seg_center + seg_parameter * D;
                        return true;
                    }
                    // else: |t| > extent, no intersection
                }
                // else: b1+b2 > 1, no intersection
            }
            // else: b2 < 0, no intersection
        }
        // else: b1 < 0, no intersection
        return false;
    }

    void rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        double theta,
        bool degrees,
        GEO::Matrix< 4, double >& rot_mat )
    {
        // Note: Rotation is impossible about an axis with null length.
        ringmesh_assert( axis != vec3() );

        if( degrees ) {
            double pi = 3.141592653589793;
            theta = theta * pi / 180.;
        }

        double axis_length = axis.length();
        ringmesh_assert( axis_length > 0. );
        double x1 = origin[0];
        double y1 = origin[1];
        double z1 = origin[2];
        double a = axis[0] / axis_length;
        double b = axis[1] / axis_length;
        double c = axis[2] / axis_length;
        double d = std::sqrt( b * b + c * c );
        double cos_angle = std::cos( theta );
        double sin_angle = std::sin( theta );

        GEO::Matrix< 4, double > T;
        T( 0, 0 ) = 1;
        T( 0, 1 ) = 0;
        T( 0, 2 ) = 0;
        T( 0, 3 ) = -x1;
        T( 1, 0 ) = 0;
        T( 1, 1 ) = 1;
        T( 1, 2 ) = 0;
        T( 1, 3 ) = -y1;
        T( 2, 0 ) = 0;
        T( 2, 1 ) = 0;
        T( 2, 2 ) = 1;
        T( 2, 3 ) = -z1;
        T( 3, 0 ) = 0;
        T( 3, 1 ) = 0;
        T( 3, 2 ) = 0;
        T( 3, 3 ) = 1;

        GEO::Matrix< 4, double > inv_T;
        inv_T( 0, 0 ) = 1.;
        inv_T( 0, 1 ) = 0.;
        inv_T( 0, 2 ) = 0.;
        inv_T( 0, 3 ) = x1;
        inv_T( 1, 0 ) = 0.;
        inv_T( 1, 1 ) = 1.;
        inv_T( 1, 2 ) = 0.;
        inv_T( 1, 3 ) = y1;
        inv_T( 2, 0 ) = 0.;
        inv_T( 2, 1 ) = 0.;
        inv_T( 2, 2 ) = 1.;
        inv_T( 2, 3 ) = z1;
        inv_T( 3, 0 ) = 0.;
        inv_T( 3, 1 ) = 0.;
        inv_T( 3, 2 ) = 0.;
        inv_T( 3, 3 ) = 1.;

#ifdef RINGMESH_DEBUG
        GEO::Matrix< 4, double > computed_inv_T = T.inverse();
#endif
        ringmesh_assert( inv_T( 0, 0 ) == computed_inv_T( 0, 0 ) );
        ringmesh_assert( inv_T( 0, 1 ) == computed_inv_T( 0, 1 ) );
        ringmesh_assert( inv_T( 0, 2 ) == computed_inv_T( 0, 2 ) );
        ringmesh_assert( inv_T( 0, 3 ) == computed_inv_T( 0, 3 ) );
        ringmesh_assert( inv_T( 1, 0 ) == computed_inv_T( 1, 0 ) );
        ringmesh_assert( inv_T( 1, 1 ) == computed_inv_T( 1, 1 ) );
        ringmesh_assert( inv_T( 1, 2 ) == computed_inv_T( 1, 2 ) );
        ringmesh_assert( inv_T( 1, 3 ) == computed_inv_T( 1, 3 ) );
        ringmesh_assert( inv_T( 2, 0 ) == computed_inv_T( 2, 0 ) );
        ringmesh_assert( inv_T( 2, 1 ) == computed_inv_T( 2, 1 ) );
        ringmesh_assert( inv_T( 2, 2 ) == computed_inv_T( 2, 2 ) );
        ringmesh_assert( inv_T( 2, 3 ) == computed_inv_T( 2, 3 ) );
        ringmesh_assert( inv_T( 3, 0 ) == computed_inv_T( 3, 0 ) );
        ringmesh_assert( inv_T( 3, 1 ) == computed_inv_T( 3, 1 ) );
        ringmesh_assert( inv_T( 3, 2 ) == computed_inv_T( 3, 2 ) );
        ringmesh_assert( inv_T( 3, 3 ) == computed_inv_T( 3, 3 ) );

        // Note: If d = 0, so rotation is along x axis. So Rx = inv_Rx = Id
        GEO::Matrix< 4, double > Rx;
        Rx( 0, 0 ) = 1.;
        Rx( 0, 1 ) = 0.;
        Rx( 0, 2 ) = 0.;
        Rx( 0, 3 ) = 0.;
        Rx( 1, 0 ) = 0.;
        Rx( 1, 3 ) = 0.;
        Rx( 2, 0 ) = 0.;
        Rx( 2, 3 ) = 0.;
        Rx( 3, 0 ) = 0.;
        Rx( 3, 1 ) = 0.;
        Rx( 3, 2 ) = 0.;
        Rx( 3, 3 ) = 1.;
        if( d == 0. ) {
            Rx( 1, 1 ) = 1.;
            Rx( 1, 2 ) = 0.;
            Rx( 2, 1 ) = 0.;
            Rx( 2, 2 ) = 1.;
        } else {
            Rx( 1, 1 ) = c / d;
            Rx( 1, 2 ) = -b / d;
            Rx( 2, 1 ) = b / d;
            Rx( 2, 2 ) = c / d;
        }

        GEO::Matrix< 4, double > inv_Rx;
        inv_Rx( 0, 0 ) = 1.;
        inv_Rx( 0, 1 ) = 0.;
        inv_Rx( 0, 2 ) = 0.;
        inv_Rx( 0, 3 ) = 0.;
        inv_Rx( 1, 0 ) = 0.;
        inv_Rx( 1, 3 ) = 0.;
        inv_Rx( 2, 0 ) = 0.;
        inv_Rx( 2, 3 ) = 0.;
        inv_Rx( 3, 0 ) = 0.;
        inv_Rx( 3, 1 ) = 0.;
        inv_Rx( 3, 2 ) = 0.;
        inv_Rx( 3, 3 ) = 1.;
        if( d == 0. ) {
            inv_Rx( 1, 1 ) = 1.;
            inv_Rx( 1, 2 ) = 0.;
            inv_Rx( 2, 1 ) = 0.;
            inv_Rx( 2, 2 ) = 1.;
        } else {
            inv_Rx( 1, 1 ) = c / d;
            inv_Rx( 1, 2 ) = b / d;
            inv_Rx( 2, 1 ) = -b / d;
            inv_Rx( 2, 2 ) = c / d;
        }

#ifdef RINGMESH_DEBUG
        GEO::Matrix< 4, double > computed_inv_Rx = Rx.inverse();
#endif
        ringmesh_assert( inv_Rx( 0, 0 ) == computed_inv_Rx( 0, 0 ) );
        ringmesh_assert( inv_Rx( 0, 1 ) == computed_inv_Rx( 0, 1 ) );
        ringmesh_assert( inv_Rx( 0, 2 ) == computed_inv_Rx( 0, 2 ) );
        ringmesh_assert( inv_Rx( 0, 3 ) == computed_inv_Rx( 0, 3 ) );
        ringmesh_assert( inv_Rx( 1, 0 ) == computed_inv_Rx( 1, 0 ) );
        ringmesh_assert( inv_Rx( 1, 1 ) == computed_inv_Rx( 1, 1 ) );
        ringmesh_assert( inv_Rx( 1, 2 ) == computed_inv_Rx( 1, 2 ) );
        ringmesh_assert( inv_Rx( 1, 3 ) == computed_inv_Rx( 1, 3 ) );
        ringmesh_assert( inv_Rx( 2, 0 ) == computed_inv_Rx( 2, 0 ) );
        ringmesh_assert( inv_Rx( 2, 1 ) == computed_inv_Rx( 2, 1 ) );
        ringmesh_assert( inv_Rx( 2, 2 ) == computed_inv_Rx( 2, 2 ) );
        ringmesh_assert( inv_Rx( 2, 3 ) == computed_inv_Rx( 2, 3 ) );
        ringmesh_assert( inv_Rx( 3, 0 ) == computed_inv_Rx( 3, 0 ) );
        ringmesh_assert( inv_Rx( 3, 1 ) == computed_inv_Rx( 3, 1 ) );
        ringmesh_assert( inv_Rx( 3, 2 ) == computed_inv_Rx( 3, 2 ) );
        ringmesh_assert( inv_Rx( 3, 3 ) == computed_inv_Rx( 3, 3 ) );

        GEO::Matrix< 4, double > Ry;
        Ry( 0, 0 ) = d;
        Ry( 0, 1 ) = 0.;
        Ry( 0, 2 ) = -a;
        Ry( 0, 3 ) = 0.;
        Ry( 1, 0 ) = 0.;
        Ry( 1, 1 ) = 1.;
        Ry( 1, 2 ) = 0.;
        Ry( 1, 3 ) = 0.;
        Ry( 2, 0 ) = a;
        Ry( 2, 1 ) = 0.;
        Ry( 2, 2 ) = d;
        Ry( 2, 3 ) = 0.;
        Ry( 3, 0 ) = 0.;
        Ry( 3, 1 ) = 0.;
        Ry( 3, 2 ) = 0.;
        Ry( 3, 3 ) = 1.;

        GEO::Matrix< 4, double > inv_Ry;
        inv_Ry( 0, 0 ) = d;
        inv_Ry( 0, 1 ) = 0.;
        inv_Ry( 0, 2 ) = a;
        inv_Ry( 0, 3 ) = 0.;
        inv_Ry( 1, 0 ) = 0.;
        inv_Ry( 1, 1 ) = 1.;
        inv_Ry( 1, 2 ) = 0.;
        inv_Ry( 1, 3 ) = 0.;
        inv_Ry( 2, 0 ) = -a;
        inv_Ry( 2, 1 ) = 0.;
        inv_Ry( 2, 2 ) = d;
        inv_Ry( 2, 3 ) = 0.;
        inv_Ry( 3, 0 ) = 0.;
        inv_Ry( 3, 1 ) = 0.;
        inv_Ry( 3, 2 ) = 0.;
        inv_Ry( 3, 3 ) = 1.;

#ifdef RINGMESH_DEBUG
        GEO::Matrix< 4, double > computed_inv_Ry = Ry.inverse();
#endif
        ringmesh_assert( inv_Ry( 0, 0 ) == computed_inv_Ry( 0, 0 ) );
        ringmesh_assert( inv_Ry( 0, 1 ) == computed_inv_Ry( 0, 1 ) );
        ringmesh_assert( inv_Ry( 0, 2 ) == computed_inv_Ry( 0, 2 ) );
        ringmesh_assert( inv_Ry( 0, 3 ) == computed_inv_Ry( 0, 3 ) );
        ringmesh_assert( inv_Ry( 1, 0 ) == computed_inv_Ry( 1, 0 ) );
        ringmesh_assert( inv_Ry( 1, 1 ) == computed_inv_Ry( 1, 1 ) );
        ringmesh_assert( inv_Ry( 1, 2 ) == computed_inv_Ry( 1, 2 ) );
        ringmesh_assert( inv_Ry( 1, 3 ) == computed_inv_Ry( 1, 3 ) );
        ringmesh_assert( inv_Ry( 2, 0 ) == computed_inv_Ry( 2, 0 ) );
        ringmesh_assert( inv_Ry( 2, 1 ) == computed_inv_Ry( 2, 1 ) );
        ringmesh_assert( inv_Ry( 2, 2 ) == computed_inv_Ry( 2, 2 ) );
        ringmesh_assert( inv_Ry( 2, 3 ) == computed_inv_Ry( 2, 3 ) );
        ringmesh_assert( inv_Ry( 3, 0 ) == computed_inv_Ry( 3, 0 ) );
        ringmesh_assert( inv_Ry( 3, 1 ) == computed_inv_Ry( 3, 1 ) );
        ringmesh_assert( inv_Ry( 3, 2 ) == computed_inv_Ry( 3, 2 ) );
        ringmesh_assert( inv_Ry( 3, 3 ) == computed_inv_Ry( 3, 3 ) );

        GEO::Matrix< 4, double > Rz;
        Rz( 0, 0 ) = cos_angle;
        Rz( 0, 1 ) = -sin_angle;
        Rz( 0, 2 ) = 0.;
        Rz( 0, 3 ) = 0.;
        Rz( 1, 0 ) = sin_angle;
        Rz( 1, 1 ) = cos_angle;
        Rz( 1, 2 ) = 0.;
        Rz( 1, 3 ) = 0.;
        Rz( 2, 0 ) = 0.;
        Rz( 2, 1 ) = 0.;
        Rz( 2, 2 ) = 1.;
        Rz( 2, 3 ) = 0.;
        Rz( 3, 0 ) = 0.;
        Rz( 3, 1 ) = 0.;
        Rz( 3, 2 ) = 0.;
        Rz( 3, 3 ) = 1.;

        rot_mat = inv_T * inv_Rx * inv_Ry * Rz * Ry * Rx * T;
    }

    bool point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool exact_predicates )
    {
        // Get another point not in the triangle plane (using its normal)
        vec3 n = cross( p2 - p0, p1 - p0 );
        vec3 q = p + n;

        Sign s1, s2, s3;
        if( !exact_predicates ) {
            double vol1 = GEO::Geom::tetra_signed_volume( p, q, p0, p1 );
            if( is_almost_zero( vol1 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s1 = sign( vol1 );
            double vol2 = GEO::Geom::tetra_signed_volume( p, q, p1, p2 );
            if( is_almost_zero( vol2 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s2 = sign( vol2 );
            double vol3 = GEO::Geom::tetra_signed_volume( p, q, p2, p0 );
            if( is_almost_zero( vol3 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s3 = sign( vol3 );
        } else {
            s1 = sign(
                GEO::PCK::orient_3d( p.data(), q.data(), p0.data(), p1.data() ) );
            s2 = sign(
                GEO::PCK::orient_3d( p.data(), q.data(), p1.data(), p2.data() ) );
            s3 = sign(
                GEO::PCK::orient_3d( p.data(), q.data(), p2.data(), p0.data() ) );

            if( s1 == ZERO ) {
                if( s2 == ZERO || s3 == ZERO ) {
                    //Case where p is exactly equal to one triangle vertex
                    return true;
                }
                return s2 == s3;
            } else if( s2 == ZERO ) {
                if( s1 == ZERO || s3 == ZERO ) {
                    return true;
                }
                return s1 == s3;
            } else if( s3 == ZERO ) {
                if( s1 == ZERO || s2 == ZERO ) {
                    return true;
                }
                return s1 == s2;
            }
        }

        return s1 == s2 && s2 == s3;
    }

    bool point_inside_triangle(
        const vec2& p,
        const vec2& p0,
        const vec2& p1,
        const vec2& p2,
        bool exact_predicates )
    {
        Sign s1, s2, s3;
        if( !exact_predicates ) {
            double area1 = GEO::Geom::triangle_signed_area( p, p0, p1 );
            if( is_almost_zero( area1 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s1 = sign( area1 );
            double area2 = GEO::Geom::triangle_signed_area( p, p1, p2 );
            if( is_almost_zero( area2 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s2 = sign( area2 );
            double area3 = GEO::Geom::triangle_signed_area( p, p2, p0 );
            if( is_almost_zero( area3 ) ) {
                return point_inside_triangle( p, p0, p1, p2, true );
            }
            s3 = sign( area3 );
        } else {
            s1 = sign( GEO::PCK::orient_2d( p.data(), p0.data(), p1.data() ) );
            s2 = sign( GEO::PCK::orient_2d( p.data(), p1.data(), p2.data() ) );
            s3 = sign( GEO::PCK::orient_2d( p.data(), p2.data(), p0.data() ) );

            if( s1 == ZERO ) {
                if( s2 == ZERO || s3 == ZERO ) {
                    //Case where p is exactly equal to one triangle vertex
                    return true;
                }
                return s2 == s3;
            } else if( s2 == ZERO ) {
                if( s1 == ZERO || s3 == ZERO ) {
                    return true;
                }
                return s1 == s3;
            } else if( s3 == ZERO ) {
                if( s1 == ZERO || s2 == ZERO ) {
                    return true;
                }
                return s1 == s2;
            }
        }

        return s1 == s2 && s2 == s3;
    }

    template bool RINGMESH_API point_segment_projection(
        const vecn< 2 >&,
        const vecn< 2 >&,
        const vecn< 2 >&,
        vecn< 2 >& );

    template bool RINGMESH_API point_segment_projection(
        const vecn< 3 >&,
        const vecn< 3 >&,
        const vecn< 3 >&,
        vecn< 3 >& );
}
