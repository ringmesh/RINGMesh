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

#include <iostream>
#include <sstream>
#include <algorithm>

namespace GRGMesh {

    MakeUnique::MakeUnique( const std::vector< vec3 >& points )
        : points_( points )
    {
        int nb_points = points_.size() ;
        indices_.resize( nb_points ) ;
        ann_points_ = annAllocPts( nb_points, 3 ) ;
        for( unsigned int i = 0; i < nb_points; i++ ) {
            indices_[i] = i ;
            std::copy( points_[i].data(), points_[i].data() + 3,
                &ann_points_[i][0] ) ;
        }
        ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_points, 3 ) ;
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

    bool MakeUnique::get_colocated(
        vec3& v,
        std::vector< unsigned int >& result,
        int nb_neighbors ) const
    {
        nb_neighbors = std::min( nb_neighbors, ann_tree_->nPoints() ) ;
        std::vector< int > neighbors( nb_neighbors ) ;
        ANNdistArray dist = new ANNdist[nb_neighbors] ;
        ann_tree_->annkSearch( v.data(), nb_neighbors, &neighbors[0], dist ) ;
        delete[] dist ;

        result.clear() ;
        for( int i = 0; i < neighbors.size(); ++i ) {
            if( inexact_equal( v, vec3( ann_points_[neighbors[i]] ) ) ) {
                result.push_back( neighbors[i] ) ;
            }
        }

        return !result.empty() ;
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
        for( unsigned int i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) continue ;
            std::vector< unsigned int > results ;
            int cur_neighbor = 0 ;
            do {
                cur_neighbor += nb_neighbors ;
                get_colocated( points_[i], results, cur_neighbor ) ;
            } while( results.size() == cur_neighbor ) ;
            unsigned int id = *std::min_element( results.begin(),
                results.end() ) ;
            for( unsigned int j = 0; j < results.size(); j++ ) {
                if( id == results[j] ) continue ;
                indices_[results[j]] = id ;
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

}

