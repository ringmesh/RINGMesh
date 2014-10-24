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

