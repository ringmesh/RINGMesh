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

namespace RINGMesh
{
    const std::array< std::array< index_t, 3 >, 4 >
        Geometry::Tetra::tetra_facet_vertex = { { { { 1, 3, 2 } },
            { { 0, 2, 3 } }, { { 3, 1, 0 } }, { { 0, 1, 2 } } } };

    vec2 normalized_perp( const vec2& v )
    {
        return normalize( vec2( v.y, -v.x ) );
    }

    double dot_perp( const vec2& v0, const vec2& v1 )
    {
        return dot( v0, vec2( v1.y, -v1.x ) );
    }

    double triangle_signed_area( const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& triangle_normal )
    {
        double area{ GEO::Geom::triangle_area( p0, p1, p2 ) };
        vec3 area_normal{ cross( p0 - p2, p1 - p2 ) };
        if( dot( triangle_normal, area_normal ) < 0 )
        {
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

    std::tuple< bool, std::array< double, 4 > > tetra_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        double total_volume{ GEO::Geom::tetra_signed_volume( p0, p1, p2, p3 ) };
        if( std::fabs( total_volume ) < global_epsilon_3 )
        {
            std::array< double, 4 > lambdas{ { 0., 0., 0., 0. } };
            return std::make_tuple( false, lambdas );
        }
        double volume0{ GEO::Geom::tetra_signed_volume( p1, p3, p2, p ) };
        double volume1{ GEO::Geom::tetra_signed_volume( p0, p2, p3, p ) };
        double volume2{ GEO::Geom::tetra_signed_volume( p0, p3, p1, p ) };
        double volume3{ GEO::Geom::tetra_signed_volume( p0, p1, p2, p ) };

        double lambda0{ volume0 / total_volume };
        double lambda1{ volume1 / total_volume };
        double lambda2{ volume2 / total_volume };
        double lambda3{ volume3 / total_volume };
        std::array< double, 4 > lambdas{ { lambda0, lambda1, lambda2,
            lambda3 } };
        return std::make_tuple( true, lambdas );
    }

    std::tuple< bool, std::array< double, 3 > >
        triangle_barycentric_coordinates(
            const vec3& p, const vec3& p0, const vec3& p1, const vec3& p2 )
    {
        double total_area{ GEO::Geom::triangle_area( p0, p1, p2 ) };
        if( std::fabs( total_area ) < global_epsilon_sq )
        {
            std::array< double, 3 > lambdas{ { 0., 0., 0. } };
            return std::make_tuple( false, lambdas );
        }
        vec3 triangle_normal{ cross( p2 - p0, p1 - p0 ) };
        double area0{ triangle_signed_area( p2, p1, p, triangle_normal ) };
        double area1{ triangle_signed_area( p0, p2, p, triangle_normal ) };
        double area2{ triangle_signed_area( p1, p0, p, triangle_normal ) };

        double lambda0{ area0 / total_area };
        double lambda1{ area1 / total_area };
        double lambda2{ area2 / total_area };
        std::array< double, 3 > lambdas{ { lambda0, lambda1, lambda2 } };
        return std::make_tuple( true, lambdas );
    }

    std::tuple< bool, std::array< double, 3 > >
        triangle_barycentric_coordinates(
            const vec2& p, const vec2& p0, const vec2& p1, const vec2& p2 )
    {
        double total_area{ GEO::Geom::triangle_signed_area( p2, p1, p0 ) };
        if( std::fabs( total_area ) < global_epsilon_sq )
        {
            std::array< double, 3 > lambdas{ { 0., 0., 0. } };
            return std::make_tuple( false, lambdas );
        }
        double area0{ GEO::Geom::triangle_signed_area( p2, p1, p ) };
        double area1{ GEO::Geom::triangle_signed_area( p0, p2, p ) };
        double area2{ GEO::Geom::triangle_signed_area( p1, p0, p ) };

        double lambda0{ area0 / total_area };
        double lambda1{ area1 / total_area };
        double lambda2{ area2 / total_area };
        std::array< double, 3 > lambdas{ { lambda0, lambda1, lambda2 } };
        return std::make_tuple( true, lambdas );
    }

    template < index_t DIMENSION >
    std::tuple< bool, vecn< DIMENSION > > point_segment_projection(
        const vecn< DIMENSION >& p,
        const vecn< DIMENSION >& p0,
        const vecn< DIMENSION >& p1 )
    {
        vecn< DIMENSION > center{ ( p0 + p1 ) * 0.5 };
        vecn< DIMENSION > diff{ p - center };
        vecn< DIMENSION > edge{ p1 - p0 };
        double extent{ 0.5 * edge.length() };
        edge = normalize( edge );
        double d{ dot( edge, diff ) };

        if( std::fabs( d ) <= extent )
        {
            vecn< DIMENSION > new_p{ center + d * edge };
            return std::make_tuple( true, new_p );
        }
        vecn< DIMENSION > empty_point;
        return std::make_tuple( false, empty_point );
    }

    GEO::Matrix< 4, double > rotation_matrix_about_arbitrary_axis(
        const vec3& origin, const vec3& axis, double theta, bool degrees )
    {
        // Note: Rotation is impossible about an axis with null length.
        ringmesh_assert( axis != vec3() );

        if( degrees )
        {
            theta = theta * M_PI / 180.;
        }

        double axis_length{ axis.length() };
        ringmesh_assert( axis_length > 0. );
        double x1{ origin[0] };
        double y1{ origin[1] };
        double z1{ origin[2] };
        double a{ axis[0] / axis_length };
        double b{ axis[1] / axis_length };
        double c{ axis[2] / axis_length };
        double d{ std::sqrt( b * b + c * c ) };
        double cos_angle{ std::cos( theta ) };
        double sin_angle{ std::sin( theta ) };

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
        if( d == 0. )
        {
            Rx( 1, 1 ) = 1.;
            Rx( 1, 2 ) = 0.;
            Rx( 2, 1 ) = 0.;
            Rx( 2, 2 ) = 1.;
        }
        else
        {
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
        if( d == 0. )
        {
            inv_Rx( 1, 1 ) = 1.;
            inv_Rx( 1, 2 ) = 0.;
            inv_Rx( 2, 1 ) = 0.;
            inv_Rx( 2, 2 ) = 1.;
        }
        else
        {
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

        return inv_T * inv_Rx * inv_Ry * Rz * Ry * Rx * T;
    }

    template std::tuple< bool, vecn< 2 > > basic_api point_segment_projection(
        const vecn< 2 >&, const vecn< 2 >&, const vecn< 2 >& );

    template std::tuple< bool, vecn< 3 > > basic_api point_segment_projection(
        const vecn< 3 >&, const vecn< 3 >&, const vecn< 3 >& );
} // namespace RINGMesh
