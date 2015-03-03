/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_AABB_TREE__
#define __RINGMESH_AABB_TREE__

#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_element.h>
#include <ringmesh/common.h>
#include <ringmesh/utils.h>

namespace RINGMesh {
    class RINGMESH_API FacetAABBTree {
    public:
        FacetAABBTree( Surface& M ) ;
        /**
         * @brief computes all the pairs of intersecting facets.
         * @param action ACTION::operator(index_t,index_t) is
         *   invoked of all pairs of facets that have overlapping
         *   bounding boxes. triangles_intersection() needs to be
         *   called to detect the actual intersections.
         */
        template< class ACTION > void compute_facet_bbox_intersections(
            ACTION& action ) const
        {
            intersect_recursive( action, 1, 0,
                mesh_.nb_cells(), 1, 0, mesh_.nb_cells() ) ;
        }

        void compute_bbox_intersections(
            const Box3d& box,
            std::vector< index_t >& results ) const
        {
            intersect_recursive( results, box, 1, 0, mesh_.nb_cells() ) ;
        }

        /**
         * @return the index of the facet nearest to point p.
         * @param [in] p query point
         * @param [out] nearest_point nearest point on the surface
         * @param [out] sq_dist squared distance between p and the surface.
         */
        index_t nearest_facet(
            const float64* p,
            float64* nearest_point,
            float64& sq_dist ) const
        {
            vec3 nearest ;
            vec3 query( p ) ;
            index_t nearest_t = nearest_facet( query, nearest, sq_dist ) ;
            std::copy( nearest.data(), nearest.data() + 3, nearest_point ) ;
            return nearest_t ;
        }

        index_t nearest_facet(
            const vec3& p,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            index_t nearest_t ;
            get_nearest_facet_hint( p, nearest_t, nearest_point, sq_dist ) ;
            nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist, 1, 0,
                mesh_.nb_cells() ) ;
            return nearest_t ;
        }

        static float64 get_nearest_point(
            const Surface& M,
            const vec3& p,
            int t,
            vec3& nearest_p ) ;

        index_t nearest_facet( const vec3& p ) const
        {
            vec3 nearest_point ;
            float64 dist ;
            return nearest_facet( p, nearest_point, dist ) ;
        }

        index_t nearest_facet( const float64* p ) const
        {
            vec3 nearest_point ;
            float64 dist ;
            return nearest_facet( vec3( p ), nearest_point, dist ) ;
        }

        /**
         * @brief computes the nearest point and nearest facet from
         * a query point, using user-specified hint.
         *
         * The hint is specified as reasonable initial values of
         * (nearest_f, nearest_p, sq_dist). If multiple queries
         * are done on a set of points that has spatial locality,
         * the hint can be the result of the previous call.
         *
         * @param [in]     p         query point
         * @param [in,out] nearest_f the nearest facet so far,
         * @param [in,out] nearest_p a point in nearest_f
         * @param [in,out] sq_dist   squared distance between p and nearest_p
         *
         */
        void nearest_facet_with_hint(
            const vec3& p,
            index_t& nearest_t,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist, 1, 0,
                mesh_.nb_cells() ) ;
        }

        /**
         * @return the squared distance between p and the surface.
         * @param [in] p query point
         */
        float64 squared_distance( const vec3& p ) const
        {
            vec3 nearest_point ;
            float64 result ;
            nearest_facet( p, nearest_point, result ) ;
            return result ;
        }

    private:
        void reorder_morton() ;

        void init_bboxes_recursive(
            index_t node,
            index_t b,
            index_t e ) ;

        /**
         * @brief computes all the pairs of intersecting facets
         *  for two sub-trees of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * @param action ACTION::operator(index_t,index_t) is
         *   invoked of all pairs of facets that have overlapping
         *   bounding boxes.
         * @param node1 index of the first node of the AABB tree
         * @param b1    index of the first facet in node1
         * @param e1    one position past the index of the last facet in node1
         * @param node2 index of the second node of the AABB tree
         * @param b1    index of the first facet in node2
         * @param e2    one position past the index of the second facet in node2
         */
        template< class ACTION > void intersect_recursive(
            ACTION& action,
            index_t node1,
            index_t b1,
            index_t e1,
            index_t node2,
            index_t b2,
            index_t e2 ) const
        {
            ringmesh_debug_assert( e1 != b1 ) ;
            ringmesh_debug_assert( e2 != b2 ) ;

            // Since we are intersecting the AABBTree with *itself*,
            // we can prune half of the cases by skipping the test
            // whenever node2's facet index interval is greated than
            // node1's facet index interval.
            if( e2 <= b1 ) {
                return ;
            }

            // The acceleration is here:
            if( !bboxes_[ node1 ].bboxes_overlap( bboxes_[ node2 ] ) ) {
                return ;
            }

            // Simple case: leaf - leaf intersection.
            if( b1 + 1 == e1 && b2 + 1 == e2 ) {
                action( b1, b2 ) ;
                return ;
            }

            // If node2 has more facets than node1, then
            //   intersect node2's two children with node1
            // else
            //   intersect node1's two children with node2
            if( e2 - b2 > e1 - b1 ) {
                index_t m2 = b2 + ( e2 - b2 ) / 2 ;
                index_t node2_l = 2 * node2 ;
                index_t node2_r = 2 * node2 + 1 ;
                intersect_recursive( action, node1, b1, e1, node2_l, b2, m2 ) ;
                intersect_recursive( action, node1, b1, e1, node2_r, m2, e2 ) ;
            } else {
                index_t m1 = b1 + ( e1 - b1 ) / 2 ;
                index_t node1_l = 2 * node1 ;
                index_t node1_r = 2 * node1 + 1 ;
                intersect_recursive( action, node1_l, b1, m1, node2, b2, e2 ) ;
                intersect_recursive( action, node1_r, m1, e1, node2, b2, e2 ) ;
            }
        }

        /**
         * @brief computes all the pairs of intersecting facets
         *  for two sub-trees of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * @param action ACTION::operator(index_t) is
         *   invoked of all pairs of facets that have overlapping
         *   bounding boxes.
         * @param box   input box to test
         * @param node  index of the node of the AABB tree
         * @param b     index of the facet in node1
         * @param e     one position past the index of the last facet in node1
         */
        void intersect_recursive(
            std::vector< index_t >& results,
            const Box3d& box,
            index_t node,
            index_t b,
            index_t e ) const
        {
            ringmesh_debug_assert( e != b ) ;

            // The acceleration is here:
            if( !box.bboxes_overlap( bboxes_[ node ] ) ) {
                return ;
            }

            // Simple case: leaf - leaf intersection.
            if( b + 1 == e ) {
                results.push_back( b ) ;
                return ;
            }

            // If node2 has more facets than node1, then
            //   intersect node2's two children with node1
            // else
            //   intersect node1's two children with node2

            index_t m = b + ( e - b ) / 2 ;
            index_t node_l = 2 * node ;
            index_t node_r = 2 * node + 1 ;
            intersect_recursive( results, box, node_l, b, m ) ;
            intersect_recursive( results, box, node_r, m, e ) ;
        }

        /**
         * @brief computes a reasonable initialization for
         *   nearest facet search.
         *
         *  A good initialization makes the algorithm faster,
         *  by allowing early pruning of subtrees that provably
         *  do not contain the nearest neighbor.
         *
         * @param [in]  p         query point
         * @param [out] nearest_f a facet reasonably near p
         * @param [out] nearest_p a point in nearest_f
         * @param [out] sq_dist   squared distance between p and nearest_p
         */
        void get_nearest_facet_hint(
            const vec3& p,
            index_t& nearest_t,
            vec3& nearest_p,
            float64& sq_dist ) const ;

        /**
         * @brief the recursive function used by the implementation
         *   of nearest_facet().
         *
         * The first call may use get_nearest_facet_hint() to initialize
         *  nearest_f, nearest_p and sq_dist, as done in nearest_facet().
         *
         * @param [in]     p         query point
         * @param [in,out] nearest_f the nearest facet so far,
         * @param [in,out] nearest_p a point in nearest_f
         * @param [in,out] sq_dist   squared distance between p and nearest_p
         */
        void nearest_facet_recursive(
            const vec3& p,
            index_t& nearest_t,
            vec3& nearest_point,
            float64& sq_dist,
            index_t n,
            index_t b,
            index_t e ) const ;

    private:
        std::vector< Box3d > bboxes_ ;
        Surface& mesh_ ;
    } ;
}

#endif
