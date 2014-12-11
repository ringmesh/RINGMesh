/*[
* Association Scientifique pour la Geologie et ses Applications (ASGA)
* Copyright (c) 1993-2013 ASGA. All Rights Reserved.
*
* This program is a Trade Secret of the ASGA and it is not to be:
* - reproduced, published, or disclosed to other,
* - distributed or displayed,
* - used for purposes or on Sites other than described
*   in the GOCAD Advancement Agreement,
* without the prior written authorization of the ASGA. Licencee
* agrees to attach or embed this Notice on all copies of the program,
* including partial copies or modified versions thereof.
]*/

#include <grgmesh/utils.h>
#include <grgmesh/boundary_model.h>
#include <grgmesh/boundary_model_element.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_private.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/third_party/shewchuk/shewchuk.h>

#include <iostream>
#include <sstream>
#include <stack>
#include <algorithm>

namespace GRGMesh {

    MakeUnique::MakeUnique( const std::vector< vec3 >& points )
        : points_( points )
    {
        uint32 nb_points = points_.size() ;
        indices_.resize( nb_points ) ;
        for( uint32 i = 0; i < nb_points; i++ ) {
            indices_[i] = i ;
        }
    }

    static bool inexact_equal( const vec3& v1, const vec3& v2 ) {
        for( uint32 i = 0; i < 3; i++ ) {
            double diff( v1[i]-v2[i] ) ;
            if( diff > epsilon || diff < -epsilon ) {
                return false ;
            }
        }
        return true ;
    }

    bool Utils::facets_have_same_orientation(
        const GEO::Mesh& mesh,
        uint32 f1,
        uint32 c11,
        uint32 f2 )
    {
        uint32 c12 = mesh.next_around_facet(f1, c11);
        uint32 v11 = mesh.corner_vertex_index( c11 ) ;
        uint32 v12 = mesh.corner_vertex_index( c12 ) ;
        for( uint32 c21 = mesh.facet_begin( f2 ); c21 < mesh.facet_end( f2 ); c21++ ) {
            uint32 c22 = mesh.next_around_facet( f2, c21 ) ;
            uint32 v21 = mesh.corner_vertex_index( c21 ) ;
            uint32 v22 = mesh.corner_vertex_index( c22 ) ;
            if( v11 == v21 && v12 == v22 ) {
                return false ;
            }
            if( v11 == v22 && v12 == v21 ) {
                return true ;
            }
        }
        return true ;
    }

    void Utils::check_and_repair_mesh_consistency(
        const BoundaryModelElement& region,
        GEO::Mesh& mesh,
        bool check_duplicated_facet )
    {
        /// 1 - Remove duplicated facets (optionnal)
        if( check_duplicated_facet ) {
            std::vector< vec3 > barycenters( mesh.nb_facets(), vec3( 0, 0, 0 ) ) ;
            for( uint32 f = 0; f < mesh.nb_facets(); f++ ) {
                barycenters[f] = GEO::Geom::mesh_facet_center( mesh, f ) ;
            }

            MakeUnique unique( barycenters ) ;
            unique.unique() ;
            const std::vector< int32 > indices = unique.indices() ;
            GEO::vector< uint32 > facet_to_remove( mesh.nb_facets(), 0 ) ;
            int32 cur_id = 0 ;
            for( uint32 f = 0; f < mesh.nb_facets(); f++ ) {
                if( cur_id == indices[f] ) {
                    cur_id++ ;
                } else {
                    facet_to_remove[f] = 1 ;
                }
            }
            mesh.remove_facets( facet_to_remove ) ;

            if( mesh.has_attribute( GEO::MESH_FACET_REGION ) ) {
                GEO::vector< int32 >& attribute = GEO::MeshMutator::facet_regions( mesh ) ;
                int32 offset = 0 ;
                cur_id = 0 ;
                for( uint32 f = 0; f < mesh.nb_facets(); f++ ) {
                    if( cur_id == indices[f] ) {
                        cur_id++ ;
                        attribute[f-offset] = attribute[f] ;
                    } else {
                        offset++ ;
                    }
                }
                attribute.erase( attribute.end()-offset, attribute.end() ) ;
            }
            mesh.update_cached_variables() ;
            mesh.connect_facets() ;
        }

        /// 2 - Reorient in the same direction using propagation
        std::vector< bool > facet_visited( mesh.nb_facets(), false ) ;
        for( uint32 f = 0; f < mesh.nb_facets(); f++ ) {
            if( facet_visited[f] ) continue ;
            uint32 surface_id = mesh.facet_region( f ) ;
            std::stack< uint32 > S ;
            S.push( f ) ;
            do {
                uint32 cur_f = S.top() ;
                S.pop() ;
                for( unsigned int c = mesh.facet_begin( cur_f );
                    c < mesh.facet_end( cur_f ); c++ ) {
                    int32 f_adj = mesh.corner_adjacent_facet( c ) ;
                    if( f_adj == -1 || mesh.facet_region( f_adj ) != surface_id
                        || facet_visited[f_adj] ) continue ;
                    if( !facets_have_same_orientation( mesh, cur_f, c, f_adj ) ) {
                        GEO::MeshMutator::flip_facet( mesh, f_adj ) ;
                    }
                    S.push( f_adj ) ;
                }
            } while( !S.empty() ) ;
        }

        /// 3 - Check for consistent orientation with BoundaryModel
        GEO::MeshFacetsAABB aabb( mesh ) ;
        std::vector< bool > flip_surface( region.model().nb_surfaces(), false ) ;
        bool flip_sthg = false ;
        for( uint32 s = 0; s < region.nb_boundaries(); s++ ) {
            const Surface& surface =
                dynamic_cast< const Surface& >( region.boundary( s ) ) ;
            vec3 barycenter = surface.facet_barycenter( 0 ) ;
            vec3 nearest_point ;
            float64 distance ;
            uint32 f = aabb.nearest_facet( barycenter, nearest_point, distance ) ;
            grgmesh_debug_assert( surface.id() == mesh.facet_region( f ) ) ;

            vec3 ori_normal = surface.facet_normal( 0 ) ;
            vec3 test_normal = GEO::Geom::mesh_facet_normal( mesh, f ) ;
            if( dot( ori_normal, test_normal ) < 0 ) {
                flip_surface[surface.id()] = true ;
                flip_sthg = true ;
            }
        }
        if( flip_sthg ) {
            for( uint32 f = 0; f < mesh.nb_facets(); f++ ) {
                uint32 surface_id = mesh.facet_region( f ) ;
                if( flip_surface[surface_id] ) {
                    GEO::MeshMutator::flip_facet( mesh, f ) ;
                }
            }
        }
    }

    bool Utils::circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result )
    {
        vec3 O_inter, D_inter ;
        if( !plan_plane_intersection( O_plane, N_plane, O_circle, N_circle, O_inter,
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
        float64 a0 = diff.length2() - r*r ;

        float64 discr = a1*a1 - a0*a2 ;
        if( discr < 0.0 ) return false ;

        if( fabs(a2) < epsilon ) return false ;
        float64 inv = 1.0/a2 ;
        if( discr < epsilon ) {
            result.push_back( vec3( O_inter - (a1*inv)*D_inter ) ) ;
        } else {
            float64 root = sqrt( discr ) ;
            result.push_back( vec3( O_inter - ( (a1 + root)*inv )*D_inter ) ) ;
            result.push_back( vec3( O_inter - ( (a1 - root)*inv )*D_inter ) ) ;
        }
        return true ;
    }

    bool Utils::plan_plane_intersection(
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
        if( fabs( d - 1  ) < epsilon ) return false ;

        float64 invDet = 1.0 / ( 1.0 - d*d ) ;
        float64 const_P0 = dot( N_P0, O_P0 ) ;
        float64 const_P1 = dot( N_P1, O_P1 ) ;
        float64 c0 = ( const_P0 - d*const_P1 ) * invDet ;
        float64 c1 = ( const_P1 - d*const_P0 ) * invDet ;
        O_inter = c0*N_P0 + c1*N_P1 ;
        D_inter = cross( N_P0, N_P1 ) ;
        return true ;
    }

    bool Utils::circle_triangle_intersection(
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
            for( uint32 i = 0; i < inter_circle_plane.size(); i++ ) {
                const vec3& p = inter_circle_plane[i] ;
                if( point_inside_triangle( p, p0, p1, p2 ) ) {
                    result.push_back( p ) ;
                }
            }
        }
        return !result.empty() ;
    }

    bool Utils::point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p )
    {
        vec3 center = (p0+p1) * 0.5 ;
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

    float64 Utils::point_quad_distance(
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
        vec3 axis[2] ;
        axis[0] = normalize( edge0 ) ;
        axis[1] = normalize( edge1 ) ;
        float64 extent[2] ;
        extent[0] = edge0.length() / 2. ;
        extent[1] = edge1.length() / 2. ;

        vec3 diff = center - p ;
        float64 b0 = dot( diff, axis[0] ) ;
        float64 b1 = dot( diff, axis[1] ) ;
        float64 s0 = -b0 ;
        float64 s1 = -b1 ;
        float64 sqrDistance = dot( diff, diff ) ;

        if( s0 < -extent[0] ) {
            s0 = -extent[0] ;
        } else if( s0 > extent[0] ) {
            s0 = extent[0] ;
        }
        sqrDistance += s0 * ( s0 + 2. * b0 ) ;

        if( s1 < -extent[1] ) {
            s1 = -extent[1] ;
        } else if( s1 > extent[1] ) {
            s1 = extent[1] ;
        }
        sqrDistance += s1 * ( s1 +  2. * b1 ) ;

        // Account for numerical round-off error.
        if( sqrDistance < 0 ) {
            sqrDistance = 0 ;
        }

        float64 distance = sqrt( sqrDistance ) ;
        nearest_p = center ;
        nearest_p += s0 * axis[0] ;
        nearest_p += s1 * axis[1] ;

        return distance ;
    }

    bool Utils::segment_triangle_intersection(
        const vec3& seg0, const vec3& seg1,
        const vec3& trgl0, const vec3& trgl1, const vec3& trgl2,
        vec3& result )
    {
        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Compute the offset origin, edges, and normal.
        vec3 seg_center = (seg0+seg1)/2 ;
        vec3 diff = seg_center - trgl0 ;
        vec3 edge1 = trgl1 - trgl0 ;
        vec3 edge2 = trgl2 - trgl0 ;
        vec3 normal = cross( edge1, edge2 ) ;

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        vec3 D = normalize( seg1-seg0 ) ;
        float64 DdN = dot( D, normal ) ;
        int32 sign ;
        if( DdN > epsilon ) {
            sign = 1 ;
        } else if( DdN < -epsilon ) {
            sign = - 1 ;
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
                    float64 extDdN = length( seg1-seg0 ) * DdN / 2. ;
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

    bool Utils::point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        vec3& P = const_cast< vec3& >( p ) ;
        vec3& P0 = const_cast< vec3& >( p0 ) ;
        vec3& P1 = const_cast< vec3& >( p1 ) ;
        vec3& P2 = const_cast< vec3& >( p2 ) ;

        // calculer la normale au triangle
        vec3 n = cross( P2 - P0, P1 - P0 ) ;

        // calculer un deuxieme point un peu au dessus du triangle
        vec3 q = P + n ;

        // calculer le signe du volume signé des trois tétraèdres qui
        // s'appuient sur [p,q] et sur les trois aretes du triangle.
        Sign s1 = sign( orient3d( P.data(), q.data(), P0.data(), P1.data() ) ) ;
        Sign s2 = sign( orient3d( P.data(), q.data(), P1.data(), P2.data() ) ) ;
        Sign s3 = sign( orient3d( P.data(), q.data(), P2.data(), P0.data() ) ) ;

        if( s1 == ZERO || s2 == ZERO || s3 == ZERO ) {
            if( inexact_equal( P, P0 ) || inexact_equal( P, P1 )
                || inexact_equal( P, P2 ) ) {
                return true ;
            }
//            std::cerr << "Point on edge... :(" << std::endl ;
            return false ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 ;
    }

    bool Utils::point_inside_quad(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3& P = const_cast< vec3& >( p ) ;
        vec3& P0 = const_cast< vec3& >( p0 ) ;
        vec3& P1 = const_cast< vec3& >( p1 ) ;
        vec3& P2 = const_cast< vec3& >( p2 ) ;
        vec3& P3 = const_cast< vec3& >( p3 ) ;

        // calculer la normale au quad
        vec3 n = cross( P2 - P0, P1 - P0 ) ;

        // calculer un deuxieme point un peu au dessus du quad
        vec3 q = P + n ;

        // calculer le signe du volume signé des quatre tétraèdres qui
        // s'appuient sur [p,q] et sur les quatre aretes du quad.
        Sign s1 = sign( orient3d( P.data(), q.data(), P0.data(), P1.data() ) ) ;
        Sign s2 = sign( orient3d( P.data(), q.data(), P1.data(), P2.data() ) ) ;
        Sign s3 = sign( orient3d( P.data(), q.data(), P2.data(), P3.data() ) ) ;
        Sign s4 = sign( orient3d( P.data(), q.data(), P3.data(), P0.data() ) ) ;

        if( s1 == ZERO || s2 == ZERO || s3 == ZERO || s4 == ZERO ) {
            if( inexact_equal( P, P0 ) || inexact_equal( P, P1 )
                || inexact_equal( P, P2 ) || inexact_equal( P, P3 ) ) {
                return true ;
            }
//            std::cerr << "Point on edge... :(" << std::endl ;
            return false ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 && s3 == s4 ;
    }

    float64 Utils::point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        vec3 vertices[4] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_triangle_distance( p,
                vertices[tetra_descriptor.facet[f][0]],
                vertices[tetra_descriptor.facet[f][1]],
                vertices[tetra_descriptor.facet[f][2]], cur_p ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    float64 Utils::point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p )
    {
        vec3 vertices[5] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[f] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[pyramid_descriptor.facet[f][0]],
                    vertices[pyramid_descriptor.facet[f][1]],
                    vertices[pyramid_descriptor.facet[f][2]], cur_p ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[pyramid_descriptor.facet[f][0]],
                    vertices[pyramid_descriptor.facet[f][1]],
                    vertices[pyramid_descriptor.facet[f][2]],
                    vertices[pyramid_descriptor.facet[f][3]], cur_p ) ;
            } else {
                grgmesh_assert_not_reached ;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    float64 Utils::point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p )
    {
        vec3 vertices[6] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[f] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[prism_descriptor.facet[f][0]],
                    vertices[prism_descriptor.facet[f][1]],
                    vertices[prism_descriptor.facet[f][2]], cur_p ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[prism_descriptor.facet[f][0]],
                    vertices[prism_descriptor.facet[f][1]],
                    vertices[prism_descriptor.facet[f][2]],
                    vertices[prism_descriptor.facet[f][3]], cur_p ) ;
            } else {
                grgmesh_assert_not_reached ;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    float64 Utils::point_hexa_distance(
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
        vec3 vertices[8] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        vertices[6] = p6 ;
        vertices[7] = p7 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < hexa_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_quad_distance( p,
                vertices[hexa_descriptor.facet[f][0]],
                vertices[hexa_descriptor.facet[f][1]],
                vertices[hexa_descriptor.facet[f][2]],
                vertices[hexa_descriptor.facet[f][3]], cur_p ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    bool Utils::point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3 vertices[4] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        for( uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[tetra_descriptor.facet[f][1]]
                    - vertices[tetra_descriptor.facet[f][0]],
                vertices[tetra_descriptor.facet[f][2]]
                    - vertices[tetra_descriptor.facet[f][0]] ) ;
            vec3 n = p
                - ( ( vertices[tetra_descriptor.facet[f][0]]
                    + vertices[tetra_descriptor.facet[f][1]]
                    + vertices[tetra_descriptor.facet[f][2]] ) / 3. ) ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    bool Utils::point_inside_pyramid(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4 )
    {
        vec3 vertices[5] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        for( uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[pyramid_descriptor.facet[f][1]]
                    - vertices[pyramid_descriptor.facet[f][0]],
                vertices[pyramid_descriptor.facet[f][2]]
                    - vertices[pyramid_descriptor.facet[f][0]] ) ;
            uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[f] ;
            vec3 barycenter ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[pyramid_descriptor.facet[f][0]]
                    + vertices[pyramid_descriptor.facet[f][1]]
                    + vertices[pyramid_descriptor.facet[f][2]] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[pyramid_descriptor.facet[f][0]]
                    + vertices[pyramid_descriptor.facet[f][1]]
                    + vertices[pyramid_descriptor.facet[f][2]]
                    + vertices[pyramid_descriptor.facet[f][3]] ) / 4. ) ;
            else
                grgmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    bool Utils::point_inside_prism(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5 )
    {
        vec3 vertices[6] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        for( uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[prism_descriptor.facet[f][1]]
                    - vertices[prism_descriptor.facet[f][0]],
                vertices[prism_descriptor.facet[f][2]]
                    - vertices[prism_descriptor.facet[f][0]] ) ;
            uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[f] ;
            vec3 barycenter ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[prism_descriptor.facet[f][0]]
                    + vertices[prism_descriptor.facet[f][1]]
                    + vertices[prism_descriptor.facet[f][2]] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[prism_descriptor.facet[f][0]]
                    + vertices[prism_descriptor.facet[f][1]]
                    + vertices[prism_descriptor.facet[f][2]]
                    + vertices[prism_descriptor.facet[f][3]] ) / 4. ) ;
            else
                grgmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    bool Utils::point_inside_hexa(
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
        vec3 vertices[7] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        vertices[6] = p6 ;
        vertices[7] = p7 ;
        for( uint8 f = 0; f < hexa_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[hexa_descriptor.facet[f][1]]
                    - vertices[hexa_descriptor.facet[f][0]],
                vertices[hexa_descriptor.facet[f][2]]
                    - vertices[hexa_descriptor.facet[f][0]] ) ;
            vec3 barycenter = ( ( vertices[hexa_descriptor.facet[f][0]]
                + vertices[hexa_descriptor.facet[f][1]]
                + vertices[hexa_descriptor.facet[f][2]]
                + vertices[hexa_descriptor.facet[f][3]] ) / 4. ) ;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    float64 Utils::nearest_point_segment(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& nearest_p )
    {
        bool inside = point_segment_projection( p, p0, p1, nearest_p ) ;
        float64 distance = length( p - nearest_p ) ;
        return inside ? -distance : distance ;
    }

    void MakeUnique::unique_points( std::vector< vec3 >& results ) const
    {
        results.reserve( indices_.size() ) ;
        int32 offset = 0, cur_id = 0 ;
        for( uint32 p = 0; p < indices_.size(); p++ ) {
            if( cur_id == indices_[p] ) {
                cur_id++ ;
                results.push_back( points_[indices_[p] + offset] ) ;
            } else {
                offset++ ;
            }
        }

        std::ofstream ind2( "ind3" ) ;
        for( unsigned int p = 0; p < indices_.size(); p++ ) {
            ind2 << indices_[p] << std::endl ;
        }
    }

    void MakeUnique::unique( int32 nb_neighbors )
    {
        ColocaterANN ann( points_ ) ;
        for( uint32 i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) continue ;
            std::vector< uint32 > results ;
            int32 cur_neighbor = 0 ;
            do {
                cur_neighbor += nb_neighbors ;
                ann.get_colocated( points_[i], results, cur_neighbor ) ;
            } while( results.size() == cur_neighbor ) ;
            uint32 id = *std::min_element( results.begin(), results.end() ) ;
            for( uint32 j = 0; j < results.size(); j++ ) {
                if( id == results[j] ) continue ;
                indices_[results[j]] = id ;
            }
        }
        int32 offset = 0 ;
        for( uint32 i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) {
                indices_[i] = indices_[indices_[i]] ;
                offset++ ;
            } else {
                indices_[i] -= offset ;
            }
        }
    }
    void MakeUnique::add_edges( const std::vector< Edge >& points ) {
        int32 offset = points_.size() ;
        points_.resize( offset+(points.size()*2) ) ;
        indices_.resize( offset+(points.size()*2) ) ;
        for( uint32 p = 0; p < points.size(); p++ ) {
            points_[offset] = points[p].value( 0 ) ;
            indices_[offset] = offset ;
            offset++ ;
            points_[offset] = points[p].value( 1 ) ;
            indices_[offset] = offset ;
            offset++ ;
        }
    }
    void MakeUnique::add_points( const std::vector< vec3 >& points ) {
        int32 offset = points_.size() ;
        points_.resize( offset+points.size() ) ;
        indices_.resize( offset+points.size() ) ;
        for( uint32 p = 0; p < points.size(); p++, offset++ ) {
            points_[offset] = points[p] ;
            indices_[offset] = offset ;
        }
    }

    ColocaterANN::ColocaterANN( const Surface& mesh )
    {
        int32 nb_vertices = mesh.nb_points() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < mesh.nb_points(); i++ ) {
            std::copy( mesh.point( i ).data(), mesh.point( i ).data() + 3,
                &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const Line& mesh )
    {
        int32 nb_vertices = mesh.nb_points() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < mesh.nb_points(); i++ ) {
            std::copy( mesh.point( i ).data(), mesh.point( i ).data() + 3,
                &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const GEO::Mesh& mesh )
    {
        int32 nb_vertices = mesh.nb_vertices() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < nb_vertices; i++ ) {
            std::copy( mesh.vertex_ptr( i ), mesh.vertex_ptr( i ) + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const std::vector< vec3 >& vertices )
    {
        int32 nb_vertices = vertices.size() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < nb_vertices; i++ ) {
            std::copy( vertices[i].data(), vertices[i].data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( float64* vertices, uint32 nb_vertices )
    {
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < nb_vertices; i++ ) {
            ann_points_[i][0] = vertices[3*i  ] ;
            ann_points_[i][1] = vertices[3*i+1] ;
            ann_points_[i][2] = vertices[3*i+2] ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const std::vector< Edge >& edges )
    {
        int32 nb_vertices = edges.size() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( uint32 i = 0; i < nb_vertices; i++ ) {
            vec3 barycenter( ( edges[i].value( 0 ) + edges[i].value( 1 ) ) / 2.0 ) ;
            std::copy( barycenter.data(), barycenter.data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( ann_points_, nb_vertices, 3 ) ;
    }

    void ColocaterANN::get_mapped_colocated(
        vec3& v,
        std::vector< uint32 >& result,
        int32 nb_neighbors )
    {
        grgmesh_debug_assert( !mapped_indices_.empty() ) ;
        get_colocated( v, result, nb_neighbors ) ;
        for( uint32 i = 0; i < result.size(); i++ ) {
            result[i] = mapped_indices_[result[i]] ;
        }
    }

    bool ColocaterANN::get_colocated(
        vec3& v,
        std::vector< uint32 >& result,
        int32 nb_neighbors )
    {
        result.clear() ;
        std::vector< int32 > neighbors(nb_neighbors) ;
        get_neighbors( v, neighbors, nb_neighbors ) ;
        for( int32 i = 0; i < neighbors.size(); ++i ) {
            if( Utils::inexact_equal( v, ann_points_[neighbors[i]] ) ) {
                result.push_back( neighbors[i] ) ;
            }
        }

        return !result.empty() ;
    }

    void ColocaterANN::get_neighbors(
        vec3& v,
        std::vector< int >& result,
        int nb_neighbors ) const
    {
        nb_neighbors = std::min( nb_neighbors, ann_tree_->nPoints() ) ;
        result.resize( nb_neighbors ) ;
        ANNdistArray dist = new ANNdist[nb_neighbors] ;
        ann_tree_->annkSearch( v.data(), nb_neighbors, &result[0], dist ) ;
        delete[] dist ;
    }

}

