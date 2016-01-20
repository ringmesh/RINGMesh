/*
* Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
*
*
*
*
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/

#include <ringmesh/geometry.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geogram_extension.h>

#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh_geometry.h>

namespace {

    using namespace RINGMesh ;

    // Compare the coordinates one by one, returns true
    // if they are all epsilon close.
    bool inexact_equal( const vec3& v1, const vec3& v2 )
    {
        for( index_t i = 0; i < 3; i++ ) {
            if( std::fabs( v1[i] - v2[i] ) > epsilon ) {
                return false ;
            }
        }
        return true ;
    }

}

namespace RINGMesh {

    /*! 
     * @todo  Comment these descriptors [JP]
     */
    GEO::CellDescriptor tetra_descriptor = {
        4,         // nb_vertices
        4,         // nb_facets
        { 3, 3, 3, 3 }, // nb_vertices in facet
        {          // facets
            { 1, 3, 2 },
            { 0, 2, 3 },
            { 3, 1, 0 },
            { 0, 1, 2 }
        },
        6,         // nb_edges
        {          // edges
            { 1, 2 }, { 2, 3 }, { 3, 1 }, { 0, 1 }, { 0, 2 }, { 0, 3 }
        }
    };


    GEO::CellDescriptor hex_descriptor = {
        8,             // nb_vertices
        6,             // nb_facets
        { 4, 4, 4, 4, 4, 4 }, // nb_vertices in facet
        {              // facets
            { 0, 2, 6, 4 },
            { 3, 1, 5, 7 },
            { 1, 0, 4, 5 },
            { 2, 3, 7, 6 },
            { 1, 3, 2, 0 },
            { 4, 6, 7, 5 }
        },
        12,            // nb_edges
        {              // edges
            { 0, 1 }, { 1, 3 }, { 3, 2 }, { 2, 0 }, { 4, 5 }, { 5, 7 },
            { 7, 6 }, { 6, 4 }, { 0, 4 }, { 1, 5 }, { 3, 7 }, { 2, 6 }
        }
    };

    GEO::CellDescriptor prism_descriptor = {
        6,             // nb_vertices
        5,             // nb_facets
        { 3, 3, 4, 4, 4 },   // nb_vertices in facet
        {              // facets
            { 0, 1, 2 },
            { 3, 5, 4 },
            { 0, 3, 4, 1 },
            { 0, 2, 5, 3 },
            { 1, 4, 5, 2 }
        },
        9,             // nb_edges
        {              // edges
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 }, { 5, 3 }, { 0, 3 }, { 1, 4 }, { 2, 5 }
        }
    };


    GEO::CellDescriptor pyramid_descriptor = {
        5,             // nb_vertices
        5,             // nb_facets
        { 4, 3, 3, 3, 3 },   // nb_vertices in facet
        {              // facets
            { 0, 1, 2, 3 },
            { 0, 4, 1 },
            { 0, 3, 4 },
            { 2, 4, 3 },
            { 2, 1, 4 }
        },
        8,             // nb_edges
        {              // edges
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 }
        }
    };


    /*!
    * Computes the distance between a point and a tetrahedron
    * @param[in] p the point
    * @param[in] p0 the first vertex of the tetrahedron
    * @param[in] p1 the second vertex of the tetrahedron
    * @param[in] p2 the third vertex of the tetrahedron
    * @param[in] p3 the fourth vertex of the tetrahedron
    * @param[out] nearest_p the nearest point on the tetrahedron
    * @return the distance between the point and the tetrahedron facets
    */
    float64 point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        vec3 vertices[ 4 ] = { p0, p1, p2, p3 };
        float64 not_used0, not_used1, not_used2 ;
        float64 dist = max_float64() ;
        for( GEO::Numeric::uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_triangle_distance( p,
                vertices[ tetra_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ tetra_descriptor.facet_vertex[ f ][ 1 ] ],
                vertices[ tetra_descriptor.facet_vertex[ f ][ 2 ] ],
                cur_p, not_used0, not_used1, not_used2 
            ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
    * Computes the distance between a point and a pyramid
    * @param[in] p the point
    * @param[in] p0 the first vertex of the pyramid
    * @param[in] p1 the second vertex of the pyramid
    * @param[in] p2 the third vertex of the pyramid
    * @param[in] p3 the fourth vertex of the pyramid
    * @param[in] p4 the fifth vertex of the pyramid
    * @param[out] nearest_p the nearest point on the pyramid
    * @return the distance between the point and the pyramid facets
    */
    float64 point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p )
    {
        vec3 vertices[ 5 ] = { p0, p1, p2, p3, p4 } ;
        float64 not_used0, not_used1, not_used2 ;
        float64 dist = max_float64() ;
        for( GEO::Numeric::uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            GEO::Numeric::uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[ f ] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ],
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 1 ] ],
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 2 ] ], 
                    cur_p, not_used0, not_used1, not_used2
                ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ],
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 1 ] ],
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 2 ] ],
                    vertices[ pyramid_descriptor.facet_vertex[ f ][ 3 ] ], cur_p 
                ) ;
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
    * Computes the distance between a point and a prism
    * @param[in] p the point
    * @param[in] p0 the first vertex of the prism
    * @param[in] p1 the second vertex of the prism
    * @param[in] p2 the third vertex of the prism
    * @param[in] p3 the fourth vertex of the prism
    * @param[in] p4 the fifth vertex of the prism
    * @param[in] p5 the sixth vertex of the prism
    * @param[out] nearest_p the nearest point on the prism
    * @return the distance between the point and the prism facets
    */
    float64 point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p )
    {
        vec3 vertices[ 6 ] = { p0, p1, p2, p3, p4, p5 } ;
        float64 not_used0, not_used1, not_used2 ;

        float64 dist = max_float64() ;
        for( GEO::Numeric::uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            GEO::Numeric::uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[ f ] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ],
                    vertices[ prism_descriptor.facet_vertex[ f ][ 1 ] ],
                    vertices[ prism_descriptor.facet_vertex[ f ][ 2 ] ], 
                    cur_p, not_used0, not_used1, not_used2
                ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ],
                    vertices[ prism_descriptor.facet_vertex[ f ][ 1 ] ],
                    vertices[ prism_descriptor.facet_vertex[ f ][ 2 ] ],
                    vertices[ prism_descriptor.facet_vertex[ f ][ 3 ] ], 
                    cur_p
                ) ;
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }
    /*!
    * Computes the distance between a point and a hexahedron
    * @param[in] p the point
    * @param[in] p0 the first vertex of the hexahedron
    * @param[in] p1 the second vertex of the hexahedron
    * @param[in] p2 the third vertex of the hexahedron
    * @param[in] p3 the fourth vertex of the hexahedron
    * @param[in] p4 the fifth vertex of the hexahedron
    * @param[in] p5 the sixth vertex of the hexahedron
    * @param[in] p6 the seventh vertex of the hexahedron
    * @param[in] p7 the heith vertex of the hexahedron
    * @param[out] nearest_p the nearest point on the hexahedron
    * @return the distance between the point and the hexahedron facets
    */
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
        vec3& nearest_p )
    {
        /// Review: Why not input an array ?
        vec3 vertices[ 8 ] = { p0, p1, p2, p3, p4, p5, p6, p7 } ;
        float64 dist = max_float64() ;
        for( GEO::Numeric::uint8 f = 0; f < hex_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_quad_distance( p,
                vertices[ hex_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ hex_descriptor.facet_vertex[ f ][ 1 ] ],
                vertices[ hex_descriptor.facet_vertex[ f ][ 2 ] ],
                vertices[ hex_descriptor.facet_vertex[ f ][ 3 ] ], 
                cur_p 
            ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
    * Tests if a point is inside a tetrahedron
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the tetrahedron
    * @param[in] p1 the second vertex of the tetrahedron
    * @param[in] p2 the third vertex of the tetrahedron
    * @param[in] p3 the fourth vertex of the tetrahedron
    * @return returns true if the point is inside the tetrahedron
    */
    bool point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3 vertices[ 4 ] = { p0, p1, p2, p3 } ;
        for( GEO::Numeric::uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[ tetra_descriptor.facet_vertex[ f ][ 1 ] ]
                - vertices[ tetra_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ tetra_descriptor.facet_vertex[ f ][ 2 ] ]
                - vertices[ tetra_descriptor.facet_vertex[ f ][ 0 ] ] ) ;
            vec3 n = p
                - ( ( vertices[ tetra_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ tetra_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ tetra_descriptor.facet_vertex[ f ][ 2 ] ] ) / 3. ) ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
    * Tests if a point is inside a pyramid
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the pyramid
    * @param[in] p1 the second vertex of the pyramid
    * @param[in] p2 the third vertex of the pyramid
    * @param[in] p3 the fourth vertex of the pyramid
    * @param[in] p4 the fifth vertex of the pyramid
    * @return returns true if the point is inside the pyramid
    */
    bool point_inside_pyramid(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4 )
    {
        vec3 vertices[ 5 ] = { p0, p1, p2, p3, p4 } ;
        for( GEO::Numeric::uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[ pyramid_descriptor.facet_vertex[ f ][ 1 ] ]
                - vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ pyramid_descriptor.facet_vertex[ f ][ 2 ] ]
                - vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ] ) ;
            GEO::Numeric::uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[ f ] ;
            vec3 barycenter( 0., 0., 0. ) ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ pyramid_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ pyramid_descriptor.facet_vertex[ f ][ 2 ] ] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[ pyramid_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ pyramid_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ pyramid_descriptor.facet_vertex[ f ][ 2 ] ]
                + vertices[ pyramid_descriptor.facet_vertex[ f ][ 3 ] ] ) / 4. ) ;
            else
                ringmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
    * Tests if a point is inside a prism
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the prism
    * @param[in] p1 the second vertex of the prism
    * @param[in] p2 the third vertex of the prism
    * @param[in] p3 the fourth vertex of the prism
    * @param[in] p4 the fifth vertex of the prism
    * @param[in] p5 the sixth vertex of the prism
    * @return returns true if the point is inside the prism
    */
    bool point_inside_prism(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5 )
    {
        vec3 vertices[ 6 ] = { p0, p1, p2, p3, p4, p5 } ;      
        for( GEO::Numeric::uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[ prism_descriptor.facet_vertex[ f ][ 1 ] ]
                - vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ prism_descriptor.facet_vertex[ f ][ 2 ] ]
                - vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ] ) ;
            GEO::Numeric::uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[ f ] ;
            vec3 barycenter(0., 0., 0.) ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ prism_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ prism_descriptor.facet_vertex[ f ][ 2 ] ] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[ prism_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ prism_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ prism_descriptor.facet_vertex[ f ][ 2 ] ]
                + vertices[ prism_descriptor.facet_vertex[ f ][ 3 ] ] ) / 4. ) ;
            else
                ringmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }
    /*!
    * Tests if a point is inside a hexahedron
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the hexahedron
    * @param[in] p1 the second vertex of the hexahedron
    * @param[in] p2 the third vertex of the hexahedron
    * @param[in] p3 the fourth vertex of the hexahedron
    * @param[in] p4 the fifth vertex of the hexahedron
    * @param[in] p5 the sixth vertex of the hexahedron
    * @param[in] p6 the seventh vertex of the hexahedron
    * @param[in] p7 the heigth vertex of the hexahedron
    * @return returns true if the point is inside the hexahedron
    */
    bool point_inside_hexa(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7 )
    {
        vec3 vertices[ 8 ] = { p0, p1, p2, p3, p4, p5, p6, p7 } ;
        for( GEO::Numeric::uint8 f = 0; f < hex_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[ hex_descriptor.facet_vertex[ f ][ 1 ] ]
                - vertices[ hex_descriptor.facet_vertex[ f ][ 0 ] ],
                vertices[ hex_descriptor.facet_vertex[ f ][ 2 ] ]
                - vertices[ hex_descriptor.facet_vertex[ f ][ 0 ] ] ) ;
            vec3 barycenter = ( ( vertices[ hex_descriptor.facet_vertex[ f ][ 0 ] ]
                + vertices[ hex_descriptor.facet_vertex[ f ][ 1 ] ]
                + vertices[ hex_descriptor.facet_vertex[ f ][ 2 ] ]
                + vertices[ hex_descriptor.facet_vertex[ f ][ 3 ] ] ) / 4. ) ;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
    * Computes the intersection(s) between a circle and a plane
    * @param[in] O_plane a point on the plane
    * @param[in] N_plane the normal of the plane
    * @param[in] O_circle the center of the circle
    * @param[in] N_circle the normal of the plane supporting the circle
    * @param[in] r the radius of the circle
    * @param[out] result the intersected points
    * @return returns true if there is at least one intersection
    */
    bool circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result )
    {
        vec3 O_inter, D_inter ;
        if( !plane_plane_intersection( O_plane, N_plane, O_circle, N_circle, O_inter,
            D_inter ) ) {
            return false ;
        }

        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Locate one or two points that are on the circle and line.  If the
        // line is t*D+P, the circle center is C, and the circle radius is r,
        // then r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2.  This
        // is a quadratic equation of the form:  a2*t^2 + 2*a1*t + a0 = 0.
        vec3 diff = O_inter - O_circle ;
        float64 a2 = D_inter.length2() ;
        float64 a1 = dot( diff, D_inter ) ;
        float64 a0 = diff.length2() - r * r ;

        float64 discr = a1 * a1 - a0 * a2 ;
        if( discr < 0.0 ) return false ;

        if( fabs( a2 ) < epsilon ) return false ;
        float64 inv = 1.0 / a2 ;
        if( discr < epsilon ) {
            result.push_back( vec3( O_inter - ( a1 * inv ) * D_inter ) ) ;
        } else {
            float64 root = sqrt( discr ) ;
            result.push_back( vec3( O_inter - ( ( a1 + root ) * inv ) * D_inter ) ) ;
            result.push_back( vec3( O_inter - ( ( a1 - root ) * inv ) * D_inter ) ) ;
        }
        return true ;
    }


    /*!
    * Computes the intersection between two planes
    * @param[in] O_P0 a point on the first plane
    * @param[in] N_P0 the normal of the frst plane
    * @param[in] O_P1 a point on the second plane
    * @param[in] N_P1 the normal of the second plane
    * @param[out] O_inter a point on the intersected line
    * @param[out] D_inter the direction of the intersected line
    * @return true is there is an intersection between the planes
    */
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

        float64 d = dot( N_P0, N_P1 ) ;
        if( fabs( d - 1 ) < epsilon ) return false ;

        float64 invDet = 1.0 / ( 1.0 - d * d ) ;
        float64 const_P0 = dot( N_P0, O_P0 ) ;
        float64 const_P1 = dot( N_P1, O_P1 ) ;
        float64 c0 = ( const_P0 - d * const_P1 ) * invDet ;
        float64 c1 = ( const_P1 - d * const_P0 ) * invDet ;
        O_inter = c0 * N_P0 + c1 * N_P1 ;
        D_inter = cross( N_P0, N_P1 ) ;
        return true ;
    }

    /*!
    * Computes the intersection(s) between a circle and a triangle
    * @param[in] p0 the first vertex of the triangle
    * @param[in] p1 the second vertex of the triangle
    * @param[in] p2 the third vertex of the triangle
    * @param[in] O_circle the center of the circle
    * @param[in] N_circle the normal of the plane supporting the circle
    * @param[in] r the radius of the circle
    * @param[out] result the intersected points
    * @return returns true if there is at least one intersection
    */
    bool circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result )
    {
        vec3 N_triangle = normalize( cross( p1 - p0, p2 - p0 ) ) ;
        vec3 barycenter = ( p0 + p1 + p2 ) / 3 ;
        std::vector< vec3 > inter_circle_plane ;
        if( circle_plane_intersection( barycenter, N_triangle, O_circle, N_circle, r,
            inter_circle_plane ) ) {
            for( index_t i = 0; i < inter_circle_plane.size(); i++ ) {
                const vec3& p = inter_circle_plane[ i ] ;
                if( point_inside_triangle( p, p0, p1, p2 ) ) {
                    result.push_back( p ) ;
                }
            }
        }
        return !result.empty() ;
    }

    /*!
    * Computes the orthogonal projection of a point on a segment
    * @param[in] p the point to project
    * @param[in] p0 the first vertex of the segment
    * @param[in] p1 the second vertex of the segment
    * @param[out] new_p the projected point
    * @return returns true if the projection is possible
    */
    bool point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p )
    {
        vec3 center = ( p0 + p1 ) * 0.5 ;
        vec3 diff = p - center ;
        vec3 edge = p1 - p0 ;
        float64 extent = 0.5 * edge.length() ;
        edge = normalize( edge ) ;
        float64 d = dot( edge, diff ) ;

        if( fabs( d ) <= extent ) {
            new_p = center + d * edge ;
            return true ;
        }
        return false ;
    }

    /*!
    * Computes the smallest distance between a point and a quad
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the quad
    * @param[in] p1 the second vertex of the quad
    * @param[in] p2 the third vertex of the quad
    * @param[in] p3 the fourth vertex of the quad
    * @param[out] nearest_p the closest point on the quad
    * @return the smallest distance
    */
    float64 point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        const vec3 center( ( p0 + p1 + p2 + p3 ) / 4. ) ;
        vec3 edge0( p1 - p0 ) ;
        vec3 edge1( p3 - p0 ) ;
        vec3 axis[ 2 ] = { normalize( edge0 ), normalize( edge1 ) };
        float64 extent[ 2 ] = { 0.5*edge0.length(), 0.5*edge1.length() } ;

        vec3 diff = center - p ;
        float64 b0 = dot( diff, axis[ 0 ] ) ;
        float64 b1 = dot( diff, axis[ 1 ] ) ;
        float64 s0 = -b0 ;
        float64 s1 = -b1 ;
        float64 sqrDistance = dot( diff, diff ) ;

        if( s0 < -extent[ 0 ] ) {
            s0 = -extent[ 0 ] ;
        } else if( s0 > extent[ 0 ] ) {
            s0 = extent[ 0 ] ;
        }
        sqrDistance += s0 * ( s0 + 2. * b0 ) ;

        if( s1 < -extent[ 1 ] ) {
            s1 = -extent[ 1 ] ;
        } else if( s1 > extent[ 1 ] ) {
            s1 = extent[ 1 ] ;
        }
        sqrDistance += s1 * ( s1 + 2. * b1 ) ;

        // Account for numerical round-off error.
        if( sqrDistance < 0 ) {
            sqrDistance = 0 ;
        }

        float64 distance = sqrt( sqrDistance ) ;
        nearest_p = center ;
        nearest_p += s0 * axis[ 0 ] ;
        nearest_p += s1 * axis[ 1 ] ;

        return distance ;
    }

    /*!
    * Computes the intersection of a segment and a triangle
    * @param[in] seg0 the first vertex of the segment
    * @param[in] seg1 the second vertex of the segment
    * @param[in] trgl0 the first vertex of the triangle
    * @param[in] trgl1 the second vertex of the triangle
    * @param[in] trgl2 the third vertex of the triangle
    * @param[out] result the intersected point
    * @return true is there is an intersection
    */
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
        vec3 seg_center = ( seg0 + seg1 ) / 2 ;
        vec3 diff = seg_center - trgl0 ;
        vec3 edge1 = trgl1 - trgl0 ;
        vec3 edge2 = trgl2 - trgl0 ;
        vec3 normal = cross( edge1, edge2 ) ;

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        vec3 D = normalize( seg1 - seg0 ) ;
        float64 DdN = dot( D, normal ) ;
        signed_index_t sign ;
        if( DdN > epsilon ) {
            sign = 1 ;
        } else if( DdN < -epsilon ) {
            sign = -1 ;
            DdN = -DdN ;
        } else {
            // Segment and triangle are parallel, call it a "no intersection"
            // even if the segment does intersect.
            return false ;
        }

        float64 DdQxE2 = sign * dot( D, cross( diff, edge2 ) ) ;
        if( DdQxE2 >= 0 ) {
            float64 DdE1xQ = sign * dot( D, cross( edge1, diff ) ) ;
            if( DdE1xQ >= 0 ) {
                if( DdQxE2 + DdE1xQ <= DdN ) {
                    // Line intersects triangle, check if segment does.
                    float64 QdN = -sign * dot( diff, normal ) ;
                    float64 extDdN = length( seg1 - seg0 ) * DdN / 2. ;
                    if( -extDdN <= QdN && QdN <= extDdN ) {
                        // Segment intersects triangle.
                        float64 inv = 1. / DdN ;
                        float64 seg_parameter = QdN * inv ;

                        result = seg_center + seg_parameter * D ;
                        return true ;
                    }
                    // else: |t| > extent, no intersection
                }
                // else: b1+b2 > 1, no intersection
            }
            // else: b2 < 0, no intersection
        }
        // else: b1 < 0, no intersection
        return false ;
    }


    /*!
    * @brief Builds a rotational matrix about an arbitrary axis.
    *
    * Mathematical development: http://paulbourke.net/geometry/rotate/.
    *
    * @param[in] origin point in which passes the rotation axis.
    *
    * @param[in] axis vector which defines the rotation axis.
    *
    * @param[in] theta rotation angle (in radians or degrees).
    *
    * @param[in] degrees true is \p theta is in degrees, false
    * if in radians.
    *
    * @param[out] rot_mat the matrix which defines the rotation
    * of a point around the axis defined by point \p origin
    * and vector \p axis by an angle \p theta.
    * New coordinates of a point (x,y,z) are:
    * (x',y',z') = rot_mat*(x,y,z)
    */
    void rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        float64 theta,
        bool degrees,
        GEO::Matrix< float64, 4 >& rot_mat )
    {
        // Note: Rotation is impossible about an axis with null length.
        ringmesh_debug_assert( axis != vec3() ) ;

        if( degrees ) {
            float64 pi = 3.141592653589793 ;
            theta = theta * pi / 180. ;
        }

        float64 axis_length = axis.length() ;
        ringmesh_debug_assert( axis_length > 0. ) ;
        float64 x1 = origin[ 0 ] ;
        float64 y1 = origin[ 1 ] ;
        float64 z1 = origin[ 2 ] ;
        float64 a = axis[ 0 ] / axis_length ;
        float64 b = axis[ 1 ] / axis_length ;
        float64 c = axis[ 2 ] / axis_length ;
        float64 d = std::sqrt( b * b + c * c ) ;
        float64 cos_angle = std::cos( theta ) ;
        float64 sin_angle = std::sin( theta ) ;

        GEO::Matrix< float64, 4 > T ;
        T( 0, 0 ) = 1 ;
        T( 0, 1 ) = 0 ;
        T( 0, 2 ) = 0 ;
        T( 0, 3 ) = -x1 ;
        T( 1, 0 ) = 0 ;
        T( 1, 1 ) = 1 ;
        T( 1, 2 ) = 0 ;
        T( 1, 3 ) = -y1 ;
        T( 2, 0 ) = 0 ;
        T( 2, 1 ) = 0 ;
        T( 2, 2 ) = 1 ;
        T( 2, 3 ) = -z1 ;
        T( 3, 0 ) = 0 ;
        T( 3, 1 ) = 0 ;
        T( 3, 2 ) = 0 ;
        T( 3, 3 ) = 1 ;

        GEO::Matrix< float64, 4 > inv_T ;
        inv_T( 0, 0 ) = 1. ;
        inv_T( 0, 1 ) = 0. ;
        inv_T( 0, 2 ) = 0. ;
        inv_T( 0, 3 ) = x1 ;
        inv_T( 1, 0 ) = 0. ;
        inv_T( 1, 1 ) = 1. ;
        inv_T( 1, 2 ) = 0. ;
        inv_T( 1, 3 ) = y1 ;
        inv_T( 2, 0 ) = 0. ;
        inv_T( 2, 1 ) = 0. ;
        inv_T( 2, 2 ) = 1. ;
        inv_T( 2, 3 ) = z1 ;
        inv_T( 3, 0 ) = 0. ;
        inv_T( 3, 1 ) = 0. ;
        inv_T( 3, 2 ) = 0. ;
        inv_T( 3, 3 ) = 1. ;

#ifdef RINGMESH_DEBUG
        GEO::Matrix< float64, 4 > computed_inv_T = T.inverse() ;
#endif
        ringmesh_debug_assert( inv_T( 0, 0 ) == computed_inv_T( 0, 0 ) ) ;
        ringmesh_debug_assert( inv_T( 0, 1 ) == computed_inv_T( 0, 1 ) ) ;
        ringmesh_debug_assert( inv_T( 0, 2 ) == computed_inv_T( 0, 2 ) ) ;
        ringmesh_debug_assert( inv_T( 0, 3 ) == computed_inv_T( 0, 3 ) ) ;
        ringmesh_debug_assert( inv_T( 1, 0 ) == computed_inv_T( 1, 0 ) ) ;
        ringmesh_debug_assert( inv_T( 1, 1 ) == computed_inv_T( 1, 1 ) ) ;
        ringmesh_debug_assert( inv_T( 1, 2 ) == computed_inv_T( 1, 2 ) ) ;
        ringmesh_debug_assert( inv_T( 1, 3 ) == computed_inv_T( 1, 3 ) ) ;
        ringmesh_debug_assert( inv_T( 2, 0 ) == computed_inv_T( 2, 0 ) ) ;
        ringmesh_debug_assert( inv_T( 2, 1 ) == computed_inv_T( 2, 1 ) ) ;
        ringmesh_debug_assert( inv_T( 2, 2 ) == computed_inv_T( 2, 2 ) ) ;
        ringmesh_debug_assert( inv_T( 2, 3 ) == computed_inv_T( 2, 3 ) ) ;
        ringmesh_debug_assert( inv_T( 3, 0 ) == computed_inv_T( 3, 0 ) ) ;
        ringmesh_debug_assert( inv_T( 3, 1 ) == computed_inv_T( 3, 1 ) ) ;
        ringmesh_debug_assert( inv_T( 3, 2 ) == computed_inv_T( 3, 2 ) ) ;
        ringmesh_debug_assert( inv_T( 3, 3 ) == computed_inv_T( 3, 3 ) ) ;

        // Note: If d = 0, so rotation is along x axis. So Rx = inv_Rx = Id
        GEO::Matrix< float64, 4 > Rx ;
        Rx( 0, 0 ) = 1. ;
        Rx( 0, 1 ) = 0. ;
        Rx( 0, 2 ) = 0. ;
        Rx( 0, 3 ) = 0. ;
        Rx( 1, 0 ) = 0. ;
        Rx( 1, 3 ) = 0. ;
        Rx( 2, 0 ) = 0. ;
        Rx( 2, 3 ) = 0. ;
        Rx( 3, 0 ) = 0. ;
        Rx( 3, 1 ) = 0. ;
        Rx( 3, 2 ) = 0. ;
        Rx( 3, 3 ) = 1. ;
        if( d == 0. ) {
            Rx( 1, 1 ) = 1. ;
            Rx( 1, 2 ) = 0. ;
            Rx( 2, 1 ) = 0. ;
            Rx( 2, 2 ) = 1. ;
        } else {
            Rx( 1, 1 ) = c / d ;
            Rx( 1, 2 ) = -b / d ;
            Rx( 2, 1 ) = b / d ;
            Rx( 2, 2 ) = c / d ;
        }

        GEO::Matrix< float64, 4 > inv_Rx ;
        inv_Rx( 0, 0 ) = 1. ;
        inv_Rx( 0, 1 ) = 0. ;
        inv_Rx( 0, 2 ) = 0. ;
        inv_Rx( 0, 3 ) = 0. ;
        inv_Rx( 1, 0 ) = 0. ;
        inv_Rx( 1, 3 ) = 0. ;
        inv_Rx( 2, 0 ) = 0. ;
        inv_Rx( 2, 3 ) = 0. ;
        inv_Rx( 3, 0 ) = 0. ;
        inv_Rx( 3, 1 ) = 0. ;
        inv_Rx( 3, 2 ) = 0. ;
        inv_Rx( 3, 3 ) = 1. ;
        if( d == 0. ) {
            inv_Rx( 1, 1 ) = 1. ;
            inv_Rx( 1, 2 ) = 0. ;
            inv_Rx( 2, 1 ) = 0. ;
            inv_Rx( 2, 2 ) = 1. ;
        } else {
            inv_Rx( 1, 1 ) = c / d ;
            inv_Rx( 1, 2 ) = b / d ;
            inv_Rx( 2, 1 ) = -b / d ;
            inv_Rx( 2, 2 ) = c / d ;
        }

#ifdef RINGMESH_DEBUG
        GEO::Matrix< float64, 4 > computed_inv_Rx = Rx.inverse() ;
#endif
        ringmesh_debug_assert( inv_Rx( 0, 0 ) == computed_inv_Rx( 0, 0 ) ) ;
        ringmesh_debug_assert( inv_Rx( 0, 1 ) == computed_inv_Rx( 0, 1 ) ) ;
        ringmesh_debug_assert( inv_Rx( 0, 2 ) == computed_inv_Rx( 0, 2 ) ) ;
        ringmesh_debug_assert( inv_Rx( 0, 3 ) == computed_inv_Rx( 0, 3 ) ) ;
        ringmesh_debug_assert( inv_Rx( 1, 0 ) == computed_inv_Rx( 1, 0 ) ) ;
        ringmesh_debug_assert( inv_Rx( 1, 1 ) == computed_inv_Rx( 1, 1 ) ) ;
        ringmesh_debug_assert( inv_Rx( 1, 2 ) == computed_inv_Rx( 1, 2 ) ) ;
        ringmesh_debug_assert( inv_Rx( 1, 3 ) == computed_inv_Rx( 1, 3 ) ) ;
        ringmesh_debug_assert( inv_Rx( 2, 0 ) == computed_inv_Rx( 2, 0 ) ) ;
        ringmesh_debug_assert( inv_Rx( 2, 1 ) == computed_inv_Rx( 2, 1 ) ) ;
        ringmesh_debug_assert( inv_Rx( 2, 2 ) == computed_inv_Rx( 2, 2 ) ) ;
        ringmesh_debug_assert( inv_Rx( 2, 3 ) == computed_inv_Rx( 2, 3 ) ) ;
        ringmesh_debug_assert( inv_Rx( 3, 0 ) == computed_inv_Rx( 3, 0 ) ) ;
        ringmesh_debug_assert( inv_Rx( 3, 1 ) == computed_inv_Rx( 3, 1 ) ) ;
        ringmesh_debug_assert( inv_Rx( 3, 2 ) == computed_inv_Rx( 3, 2 ) ) ;
        ringmesh_debug_assert( inv_Rx( 3, 3 ) == computed_inv_Rx( 3, 3 ) ) ;

        GEO::Matrix< float64, 4 > Ry ;
        Ry( 0, 0 ) = d ;
        Ry( 0, 1 ) = 0. ;
        Ry( 0, 2 ) = -a ;
        Ry( 0, 3 ) = 0. ;
        Ry( 1, 0 ) = 0. ;
        Ry( 1, 1 ) = 1. ;
        Ry( 1, 2 ) = 0. ;
        Ry( 1, 3 ) = 0. ;
        Ry( 2, 0 ) = a ;
        Ry( 2, 1 ) = 0. ;
        Ry( 2, 2 ) = d ;
        Ry( 2, 3 ) = 0. ;
        Ry( 3, 0 ) = 0. ;
        Ry( 3, 1 ) = 0. ;
        Ry( 3, 2 ) = 0. ;
        Ry( 3, 3 ) = 1. ;

        GEO::Matrix< float64, 4 > inv_Ry ;
        inv_Ry( 0, 0 ) = d ;
        inv_Ry( 0, 1 ) = 0. ;
        inv_Ry( 0, 2 ) = a ;
        inv_Ry( 0, 3 ) = 0. ;
        inv_Ry( 1, 0 ) = 0. ;
        inv_Ry( 1, 1 ) = 1. ;
        inv_Ry( 1, 2 ) = 0. ;
        inv_Ry( 1, 3 ) = 0. ;
        inv_Ry( 2, 0 ) = -a ;
        inv_Ry( 2, 1 ) = 0. ;
        inv_Ry( 2, 2 ) = d ;
        inv_Ry( 2, 3 ) = 0. ;
        inv_Ry( 3, 0 ) = 0. ;
        inv_Ry( 3, 1 ) = 0. ;
        inv_Ry( 3, 2 ) = 0. ;
        inv_Ry( 3, 3 ) = 1. ;

#ifdef RINGMESH_DEBUG
        GEO::Matrix< float64, 4 > computed_inv_Ry = Ry.inverse() ;
#endif
        ringmesh_debug_assert( inv_Ry( 0, 0 ) == computed_inv_Ry( 0, 0 ) ) ;
        ringmesh_debug_assert( inv_Ry( 0, 1 ) == computed_inv_Ry( 0, 1 ) ) ;
        ringmesh_debug_assert( inv_Ry( 0, 2 ) == computed_inv_Ry( 0, 2 ) ) ;
        ringmesh_debug_assert( inv_Ry( 0, 3 ) == computed_inv_Ry( 0, 3 ) ) ;
        ringmesh_debug_assert( inv_Ry( 1, 0 ) == computed_inv_Ry( 1, 0 ) ) ;
        ringmesh_debug_assert( inv_Ry( 1, 1 ) == computed_inv_Ry( 1, 1 ) ) ;
        ringmesh_debug_assert( inv_Ry( 1, 2 ) == computed_inv_Ry( 1, 2 ) ) ;
        ringmesh_debug_assert( inv_Ry( 1, 3 ) == computed_inv_Ry( 1, 3 ) ) ;
        ringmesh_debug_assert( inv_Ry( 2, 0 ) == computed_inv_Ry( 2, 0 ) ) ;
        ringmesh_debug_assert( inv_Ry( 2, 1 ) == computed_inv_Ry( 2, 1 ) ) ;
        ringmesh_debug_assert( inv_Ry( 2, 2 ) == computed_inv_Ry( 2, 2 ) ) ;
        ringmesh_debug_assert( inv_Ry( 2, 3 ) == computed_inv_Ry( 2, 3 ) ) ;
        ringmesh_debug_assert( inv_Ry( 3, 0 ) == computed_inv_Ry( 3, 0 ) ) ;
        ringmesh_debug_assert( inv_Ry( 3, 1 ) == computed_inv_Ry( 3, 1 ) ) ;
        ringmesh_debug_assert( inv_Ry( 3, 2 ) == computed_inv_Ry( 3, 2 ) ) ;
        ringmesh_debug_assert( inv_Ry( 3, 3 ) == computed_inv_Ry( 3, 3 ) ) ;

        GEO::Matrix< float64, 4 > Rz ;
        Rz( 0, 0 ) = cos_angle ;
        Rz( 0, 1 ) = -sin_angle ;
        Rz( 0, 2 ) = 0. ;
        Rz( 0, 3 ) = 0. ;
        Rz( 1, 0 ) = sin_angle ;
        Rz( 1, 1 ) = cos_angle ;
        Rz( 1, 2 ) = 0. ;
        Rz( 1, 3 ) = 0. ;
        Rz( 2, 0 ) = 0. ;
        Rz( 2, 1 ) = 0. ;
        Rz( 2, 2 ) = 1. ;
        Rz( 2, 3 ) = 0. ;
        Rz( 3, 0 ) = 0. ;
        Rz( 3, 1 ) = 0. ;
        Rz( 3, 2 ) = 0. ;
        Rz( 3, 3 ) = 1. ;

        rot_mat = inv_T * inv_Rx * inv_Ry * Rz * Ry * Rx * T ;
    }



    /*!
    * Tests if a point is inside a triangle, more precisely if it is inside
    * a prism based on the triangle and its normal
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the triangle
    * @param[in] p1 the second vertex of the triangle
    * @param[in] p2 the third vertex of the triangle
    * @return returns true if the point is inside
    */
    bool point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        vec3 n = cross( p2 - p0, p1 - p0 ) ;
        vec3 q = p + n ;

        Sign s1 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p0.data(), p1.data() ) ) ;
        Sign s2 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p1.data(), p2.data() ) ) ;
        Sign s3 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p2.data(), p0.data() ) ) ;
        if( s1 == ZERO || s2 == ZERO || s3 == ZERO ) {
            return true ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 ;
    }

    /*!
    * Tests if a point is inside a quad, more precisely if it is inside the box
    * based on the quad and its normal
    * @param[in] p the point to test
    * @param[in] p0 the first vertex of the quad
    * @param[in] p1 the second vertex of the quad
    * @param[in] p2 the third vertex of the quad
    * @param[in] p3 the fourth vertex of the quad
    * @return returns true if the point is inside
    */
    bool point_inside_quad(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3 n = cross( p2 - p0, p1 - p0 ) ;
        vec3 q = p + n ;

        Sign s1 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p0.data(), p1.data() ) ) ;
        Sign s2 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p1.data(), p2.data() ) ) ;
        Sign s3 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p2.data(), p3.data() ) ) ;
        Sign s4 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p3.data(), p0.data() ) ) ;

        if( s1 == ZERO || s2 == ZERO || s3 == ZERO || s4 == ZERO ) {
            if( inexact_equal( p, p0 ) || inexact_equal( p, p1 )
                || inexact_equal( p, p2 ) || inexact_equal( p, p3 ) ) {
                return true ;
            }
            return false ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 && s3 == s4 ;
    }


    MakeUnique::MakeUnique( const std::vector< vec3 >& points )
        : points_( points )
    {
        index_t nb_points = points_.size() ;
        indices_.resize( nb_points ) ;
        for( index_t i = 0; i < nb_points; i++ ) {
            indices_[ i ] = i ;
        }
    }



    /*!
    * Gets a vector of unique points (initial vector - colocated points)
    * @param[out] results the vector to fill
    */
    void MakeUnique::unique_points( std::vector< vec3 >& results ) const
    {
        results.reserve( indices_.size() ) ;
        signed_index_t offset = 0, cur_id = 0 ;
        for( index_t p = 0; p < indices_.size(); p++ ) {
            if( cur_id == indices_[ p ] ) {
                cur_id++ ;
                results.push_back( points_[ indices_[ p ] + offset ] ) ;
            } else {
                offset++ ;
            }
        }
    }

    /*!
    * Computes the unique database
    */
    void MakeUnique::unique()
    {
        ColocaterANN ann( points_ ) ;
        for( index_t i = 0; i < indices_.size(); i++ ) {
            if( indices_[ i ] != i ) continue ;
            std::vector< index_t > results ;
            ann.get_colocated( points_[ i ], results ) ;
            index_t id = *std::min_element( results.begin(), results.end() ) ;
            for( index_t j = 0; j < results.size(); j++ ) {
                if( id == results[ j ] ) continue ;
                indices_[ results[ j ] ] = id ;
            }
        }
        index_t offset = 0 ;
        for( index_t i = 0; i < indices_.size(); i++ ) {
            if( indices_[ i ] != i ) {
                indices_[ i ] = indices_[ indices_[ i ] ] ;
                offset++ ;
            } else {
                indices_[ i ] -= offset ;
            }
        }
    }
    /*!
    * Add edges to the initial vector
    * @param[in] points the edges to add
    */
    void MakeUnique::add_edges( const std::vector< std::pair< vec3, vec3 > >& points )
    {
        signed_index_t offset = points_.size() ;
        points_.resize( offset + ( points.size() * 2 ) ) ;
        indices_.resize( offset + ( points.size() * 2 ) ) ;
        for( index_t p = 0; p < points.size(); p++ ) {
            points_[ offset ] = points[ p ].first ;
            indices_[ offset ] = offset ;
            offset++ ;
            points_[ offset ] = points[ p ].second ;
            indices_[ offset ] = offset ;
            offset++ ;
        }
    }
    /*!
    * Add points to the initial vector
    * @param[in] points the points to add
    */
    void MakeUnique::add_points( const std::vector< vec3 >& points )
    {
        signed_index_t offset = points_.size() ;
        points_.resize( offset + points.size() ) ;
        indices_.resize( offset + points.size() ) ;
        for( index_t p = 0; p < points.size(); p++, offset++ ) {
            points_[ offset ] = points[ p ] ;
            indices_[ offset ] = offset ;
        }
    }

    ColocaterANN::ColocaterANN()
        : ann_points_( nil )
    {}

    ColocaterANN::ColocaterANN(
        const GEO::Mesh& mesh,
        const MeshLocation& location,
        bool copy )
    {
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        switch( location ) {
            case VERTICES: {
                if( !copy ) {
                    ann_points_ = nil ;
                    ann_tree_->set_points( mesh.vertices.nb(),
                                           mesh.vertices.point_ptr( 0 ) ) ;
                } else {
                    index_t nb_vertices = mesh.vertices.nb() ;
                    ann_points_ = new double[ nb_vertices * 3 ] ;
                    GEO::Memory::copy( ann_points_, mesh.vertices.point_ptr( 0 ),
                                       nb_vertices * 3 * sizeof( double ) ) ;
                    ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                }
                break ;
            }
            case FACETS: {
                index_t nb_vertices = mesh.facets.nb() ;
                ann_points_ = new double[ nb_vertices * 3 ] ;
                for( index_t i = 0; i < mesh.facets.nb(); i++ ) {
                    vec3 center = GEO::Geom::mesh_facet_center( mesh, i ) ;
                    index_t index_in_ann = 3 * i ;
                    ann_points_[ index_in_ann ] = center.x ;
                    ann_points_[ index_in_ann + 1 ] = center.y ;
                    ann_points_[ index_in_ann + 2 ] = center.z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
            case CELLS: {
                index_t nb_vertices = mesh.cells.nb() ;
                ann_points_ = new double[ nb_vertices * 3 ] ;
                for( index_t i = 0; i < mesh.cells.nb(); i++ ) {
                    vec3 center = mesh_cell_center( mesh, i ) ;
                    index_t index_in_ann = 3 * i ;
                    ann_points_[ index_in_ann ] = center.x ;
                    ann_points_[ index_in_ann + 1 ] = center.y ;
                    ann_points_[ index_in_ann + 2 ] = center.z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
        }
    }

    ColocaterANN::ColocaterANN( const std::vector< vec3 >& vertices, bool copy )
    {
        index_t nb_vertices = vertices.size() ;
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        if( copy ) {
            ann_points_ = new double[nb_vertices * 3] ;
            GEO::Memory::copy( ann_points_, vertices.data()->data(),
                3 * nb_vertices * sizeof(double) ) ;
            ann_tree_->set_points( nb_vertices, ann_points_ ) ;
        } else {
            ann_points_ = nil ;
            ann_tree_->set_points( nb_vertices, vertices.data()->data() ) ;
        }
    }

    /*!
    * Compute the colocated point(s) of a given point
    * @param[in] v the point to test
    * @param[out] result the colocated point indices with a precision epsilon.
    * @return return true if there is at least one intersection
    */
    bool ColocaterANN::get_colocated(
        const vec3& v,
        std::vector< index_t >& result ) const
    {
        index_t nb_neighbors = std::min( index_t( 5 ), ann_tree_->nb_points() ) ;
        result.clear() ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*)alloca( sizeof( double ) * cur_neighbor ) ;
            nb_neighbors = get_neighbors( v, cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                if( dist[ i ] > epsilon_sq ) {
                    break ;
                }
                result.push_back( neighbors[ i ] ) ;
            }
        } while( result.size() == cur_neighbor ) ;

        return !result.empty() ;
    }

    /*!
    * Gets the neighboring points of a given one sorted by increasing distance
    * @param[in] v the point to test
    * @param[in] nb_neighbors the number of neighbors to return
    * @param[out] result the neighboring points
    * @param[out] dist the square distance between each neigbhor and the point \p v
    * @return the number of neighbors returned (can be less than \p nb_neighbors
    * if there is not enough points)
    */
    index_t ColocaterANN::get_neighbors(
        const vec3& v,
        index_t nb_neighbors,
        std::vector< index_t >& result,
        double* dist ) const
    {
        if( ann_tree_->nb_points() == 0 ) return 0 ;
        if( !dist ) {
            dist = (double*)alloca( sizeof( double ) * nb_neighbors ) ;
        }
        nb_neighbors = std::min( nb_neighbors, ann_tree_->nb_points() ) ;
        result.resize( nb_neighbors ) ;
        ann_tree_->get_nearest_neighbors( nb_neighbors, v.data(), &result[ 0 ],
                                          dist ) ;
        return nb_neighbors ;
    }
}
