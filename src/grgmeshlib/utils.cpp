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

#include <grgmeshlib/utils.h>
#include <grgmeshlib/boundary_model_element.h>
#include <grgmeshlib/mixed_mesh.h>

#include <third_party/shewchuk.h>

#include <iostream>
#include <sstream>
#include <algorithm>

namespace GRGMesh {

    MakeUnique::MakeUnique( const std::vector< vec3 >& points )
        : points_( points )
    {
        int nb_points = points_.size() ;
        indices_.resize( nb_points ) ;
        for( unsigned int i = 0; i < nb_points; i++ ) {
            indices_[i] = i ;
        }
    }

    static bool inexact_equal( const vec3& v1, const vec3& v2 ) {
        for( unsigned int i = 0; i < 3; i++ ) {
            double diff( v1[i]-v2[i] ) ;
            if( diff > epsilon || diff < -epsilon ) {
                return false ;
            }
        }
        return true ;
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
            for( unsigned int i = 0; i < inter_circle_plane.size(); i++ ) {
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
        int sign ;
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



    void MakeUnique::unique_points( std::vector< vec3 >& results ) const
    {
        results.reserve( indices_.size() ) ;
        int offset = 0, cur_id = 0 ;
        for( unsigned int p = 0; p < indices_.size(); p++ ) {
            if( cur_id == indices_[p] ) {
                cur_id++ ;
                results.push_back( points_[indices_[p] + offset] ) ;
            } else {
                offset++ ;
            }
        }
    }

    void MakeUnique::unique( int nb_neighbors )
    {
        ColocaterANN ann( points_ ) ;
        for( unsigned int i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) continue ;
            std::vector< unsigned int > results ;
            if( ann.get_colocated( points_[i], results, nb_neighbors ) ) {
                unsigned int id = *std::min_element( results.begin(),
                    results.end() ) ;
                for( unsigned int j = 0; j < results.size(); j++ ) {
                    if( id == results[j] ) continue ;
                    indices_[results[j]] = id ;
                }
            }
        }
        int offset = 0 ;
        for( unsigned int i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) {
                indices_[i] = indices_[indices_[i]] ;
                offset++ ;
            } else {
                indices_[i] -= offset ;
            }
        }
    }
    void MakeUnique::add_edges( const std::vector< Edge >& points ) {
        int offset = points_.size() ;
        points_.resize( offset+(points.size()*2) ) ;
        indices_.resize( offset+(points.size()*2) ) ;
        for( unsigned int p = 0; p < points.size(); p++ ) {
            points_[offset] = points[p].value( 0 ) ;
            indices_[offset] = offset ;
            offset++ ;
            points_[offset] = points[p].value( 1 ) ;
            indices_[offset] = offset ;
            offset++ ;
        }
    }
    void MakeUnique::add_points( const std::vector< vec3 >& points ) {
        int offset = points_.size() ;
        points_.resize( offset+points.size() ) ;
        indices_.resize( offset+points.size() ) ;
        for( unsigned int p = 0; p < points.size(); p++, offset++ ) {
            points_[offset] = points[p] ;
            indices_[offset] = offset ;
        }
    }

    ColocaterANN::ColocaterANN( const SurfacePart& mesh )
    {
        int nb_vertices = mesh.nb_points() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( unsigned int i = 0; i < mesh.nb_points(); i++ ) {
            std::copy( mesh.point( i ).data(), mesh.point( i ).data() + 3,
                &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const ContactPart& mesh )
    {
        int nb_vertices = mesh.nb_points() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( unsigned int i = 0; i < mesh.nb_points(); i++ ) {
            std::copy( mesh.point( i ).data(), mesh.point( i ).data() + 3,
                &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const MixedMesh& mesh )
    {
        int nb_vertices = mesh.nb_vertices() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( unsigned int i = 0; i < nb_vertices; i++ ) {
            vec3 vertex = mesh.vertex( i ) ;
            std::copy( vertex.data(), vertex.data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const MixedMesh& mesh, bool use_surface_id )
    {
        grgmesh_assert_not_reached ;
        //todo
        /*
        std::vector< vec3 > barycenters ;
        barycenters.reserve( 4*mesh.nb_tetra() ) ;
        for( unsigned int i = 0; i < mesh.nb_tetra(); i++ ) {
            for( unsigned int f = 0; f < 4; f++ ) {
                if( mesh.surface_id( i, f ) != -1 ) {
                    barycenters.push_back(
                        VorteXUtils::barycenter( mesh.vertex( i, f2e[f][0] ),
                            mesh.vertex( i, f2e[f][1] ),
                            mesh.vertex( i, f2e[f][2] ) ) ) ;
                    mapped_indices_.push_back( 4*i + f ) ;
                }
            }
        }
        ann_points_ = annAllocPts( barycenters.size(), 3 ) ;
        for( unsigned int i = 0; i < barycenters.size(); i++ ) {
            std::copy( barycenters[i].data(), barycenters[i].data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], barycenters.size(), 3 ) ;
        */
    }

    ColocaterANN::ColocaterANN( const std::vector< vec3 >& vertices )
    {
        int nb_vertices = vertices.size() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( unsigned int i = 0; i < nb_vertices; i++ ) {
            std::copy( vertices[i].data(), vertices[i].data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
    }

    ColocaterANN::ColocaterANN( const std::vector< Edge >& edges )
    {
        int nb_vertices = edges.size() ;
        ann_points_ = annAllocPts( nb_vertices, 3 ) ;
        for( unsigned int i = 0; i < nb_vertices; i++ ) {
            vec3 barycenter( ( edges[i].value( 0 ) + edges[i].value( 1 ) ) / 2.0 ) ;
            std::copy( barycenter.data(), barycenter.data() + 3, &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
    }

    void ColocaterANN::get_mapped_colocated(
        vec3& v,
        std::vector< unsigned int >& result,
        int nb_neighbors ) const
    {
        grgmesh_debug_assert( !mapped_indices_.empty() ) ;
        get_colocated( v, result, nb_neighbors ) ;
        for( unsigned int i = 0; i < result.size(); i++ ) {
            result[i] = mapped_indices_[result[i]] ;
        }
    }

    bool ColocaterANN::get_colocated(
        vec3& v,
        std::vector< unsigned int >& result,
        int nb_neighbors ) const
    {
        result.clear() ;
        std::vector< int > neighbors(nb_neighbors) ;
        get_neighbors( v, neighbors, nb_neighbors ) ;
        for( int i = 0; i < neighbors.size(); ++i ) {
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

