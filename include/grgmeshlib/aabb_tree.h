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

#ifndef __GRGMESH_AABB_TREE__
#define __GRGMESH_AABB_TREE__

#include <grgmeshlib/common.h>
#include <grgmeshlib/mixed_mesh.h>
#include <grgmeshlib/boundary_model.h>
#include <grgmeshlib/boundary_model_element.h>
#include <grgmeshlib/utils.h>

namespace GRGMesh {

    class GRGMESH_API MixedAABBTree {
    public:
        MixedAABBTree( MixedMesh& M ) ;

        void build_tree() ;

        void use_lines( bool b = true ) { use_lines_ = b ; }
        void use_triangles( bool b = true ) { use_trgl_ = b ; }
        void use_quad( bool b = true ) { use_quad_ = b ; }
        void use_tetra( bool b = true ) { use_tetra_ = b ; }
        void use_pyramids( bool b = true ) { use_pyramids_ = b ; }
        void use_prisms( bool b = true ) { use_prisms_ = b ; }
        void use_hexa( bool b = true ) { use_hexa_ = b ; }

        void use_1D_cells( bool b = true ) {
            use_lines( b ) ;
        }
        void use_2D_cells( bool b = true ) {
            use_triangles( b ) ; use_quad( b ) ; 
        }
        void use_3D_cells( bool b = true ) {
            use_tetra( b ) ; use_pyramids( b ) ; use_prisms( b ) ; use_hexa( b ) ;
        }
        void use_all( bool b = true ) {
            use_1D_cells( b ) ; use_2D_cells( b ) ; use_3D_cells( b ) ;
        }

        uint32 nb_vertices() const { return ordered_vertices_.size() ; }
        uint32 vertex( uint32 v ) const { return ordered_vertices_[v] ; }
        uint32 nb_cells() const { return ordered_cells_.size() ; }
        uint32 cell( uint32 e ) const { return ordered_cells_[e] ; }

        /**
         * @brief computes all the pairs of intersecting tets.
         * @param action ACTION::operator(uint32,uint32) is
         *   invoked of all pairs of tets that have overlapping
         *   bounding boxes. triangles_intersection() needs to be
         *   called to detect the actual intersections.
         */
        template< class ACTION > void compute_tet_bbox_intersections(
            ACTION& action ) const
        {
            intersect_recursive( action, 1, 0, nb_cells(), 1, 0, nb_cells() ) ;
        }

        /**
         * @return the index of the tet nearest to point p.
         * @param [in] p query point
         * @param [out] nearest_point nearest point on the surface
         * @param [out] sq_dist squared distance between p and the surface.
         */
        uint32 nearest_cell(
            const vec3& p,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            uint32 nearest_c ;
            get_nearest_cell_hint( p, nearest_c, nearest_point, sq_dist ) ;
            nearest_cell_recursive( p, nearest_c, nearest_point, sq_dist, 1, 0,
                nb_cells() ) ;
            return cell( nearest_c ) ;
        }

        uint32 nearest_cell( const vec3& p ) const
        {
            vec3 nearest_point ;
            float64 dist ;
            return nearest_cell( p, nearest_point, dist ) ;
        }

        /**
         * @brief computes the nearest point and nearest tet from
         * a query point, using user-specified hint.
         *
         * The hint is specified as reasonable initial values of
         * (nearest_f, nearest_p, sq_dist). If multiple queries
         * are done on a set of points that has spatial locality,
         * the hint can be the result of the previous call.
         *
         * @param [in]     p         query point
         * @param [in,out] nearest_f the nearest tet so far,
         * @param [in,out] nearest_p a point in nearest_f
         * @param [in,out] sq_dist   squared distance between p and nearest_p
         *
         */
        void nearest_cell_with_hint(
            const vec3& p,
            uint32& nearest_c,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            nearest_cell_recursive( p, nearest_c, nearest_point, sq_dist, 1, 0,
                nb_cells() ) ;
            nearest_c = cell( nearest_c ) ;
        }

        /**
         * @return the squared distance between p and the surface.
         * @param [in] p query point
         */
        float64 squared_distance( const vec3& p ) const
        {
            vec3 nearest_point ;
            float64 result ;
            nearest_cell( p, nearest_point, result ) ;
            return result ;
        }

    private:
        void reorder_morton() ;

        void init_bboxes_recursive( uint32 node, uint32 b, uint32 e ) ;

        /**
         * @brief computes all the pairs of intersecting tets
         *  for two sub-trees of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) tet indices
         *  sequences that correspond to the tets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * @param action ACTION::operator(uint32,uint32) is
         *   invoked of all pairs of tets that have overlapping
         *   bounding boxes.
         * @param node1 index of the first node of the AABB tree
         * @param b1    index of the first tet in node1
         * @param e1    one position past the index of the last tet in node1
         * @param node2 index of the second node of the AABB tree
         * @param b1    index of the first tet in node2
         * @param e2    one position past the index of the second tet in node2
         */
        template< class ACTION > void intersect_recursive(
            ACTION& action,
            uint32 node1,
            uint32 b1,
            uint32 e1,
            uint32 node2,
            uint32 b2,
            uint32 e2 ) const
        {
            grgmesh_debug_assert( e1 != b1 ) ;
            grgmesh_debug_assert( e2 != b2 ) ;

            // Since we are intersecting the AABBTree with *itself*,
            // we can prune half of the cases by skipping the test
            // whenever node2's tet index interval is greated than
            // node1's tet index interval.
            if( e2 <= b1 ) {
                return ;
            }

            // The acceleration is here:
            if( !bboxes_[node1].bboxes_overlap( bboxes_[node2] ) ) {
                return ;
            }

            // Simple case: leaf - leaf intersection.
            if( b1 + 1 == e1 && b2 + 1 == e2 ) {
                action( b1, b2 ) ;
                return ;
            }

            // If node2 has more tets than node1, then
            //   intersect node2's two children with node1
            // else
            //   intersect node1's two children with node2
            if( e2 - b2 > e1 - b1 ) {
                uint32 m2 = b2 + ( e2 - b2 ) / 2 ;
                uint32 node2_l = 2 * node2 ;
                uint32 node2_r = 2 * node2 + 1 ;
                intersect_recursive( action, node1, b1, e1, node2_l, b2, m2 ) ;
                intersect_recursive( action, node1, b1, e1, node2_r, m2, e2 ) ;
            } else {
                uint32 m1 = b1 + ( e1 - b1 ) / 2 ;
                uint32 node1_l = 2 * node1 ;
                uint32 node1_r = 2 * node1 + 1 ;
                intersect_recursive( action, node1_l, b1, m1, node2, b2, e2 ) ;
                intersect_recursive( action, node1_r, m1, e1, node2, b2, e2 ) ;
            }
        }

        /**
         * @brief computes a reasonable initialization for
         *   nearest tet search.
         *
         *  A good initialization makes the algorithm faster,
         *  by allowing early pruning of subtrees that provably
         *  do not contain the nearest neighbor.
         *
         * @param [in]  p         query point
         * @param [out] nearest_f a tet reasonably near p
         * @param [out] nearest_p a point in nearest_f
         * @param [out] sq_dist   squared distance between p and nearest_p
         */
        void get_nearest_cell_hint(
            const vec3& p,
            uint32& nearest_c,
            vec3& nearest_p,
            float64& sq_dist ) const ;

        /**
         * @brief the recursive function used by the implementation
         *   of nearest_tet().
         *
         * The first call may use get_nearest_tet_hint() to initialize
         *  nearest_f, nearest_p and sq_dist, as done in nearest_tet().
         *
         * @param [in]     p         query point
         * @param [in,out] nearest_f the nearest tet so far,
         * @param [in,out] nearest_p a point in nearest_f
         * @param [in,out] sq_dist   squared distance between p and nearest_p
         */
        void nearest_cell_recursive(
            const vec3& p,
            uint32& nearest_c,
            vec3& nearest_point,
            float64& sq_dist,
            uint32 n,
            uint32 b,
            uint32 e ) const ;

    private:
        std::vector< Box3d > bboxes_ ;
        MixedMesh& mesh_ ;

        std::vector< int32 > ordered_vertices_ ;
        std::vector< int32 > ordered_cells_ ;

        bool use_lines_ ;
        bool use_trgl_ ;
        bool use_quad_ ;
        bool use_tetra_ ;
        bool use_pyramids_ ;
        bool use_prisms_ ;
        bool use_hexa_ ;

    } ;

    class GRGMESH_API FacetAABBTree {
    public:
        FacetAABBTree( SurfacePart& M ) ;
        /**
         * @brief computes all the pairs of intersecting facets.
         * @param action ACTION::operator(uint32,uint32) is
         *   invoked of all pairs of facets that have overlapping
         *   bounding boxes. triangles_intersection() needs to be
         *   called to detect the actual intersections.
         */
        template< class ACTION > void compute_facet_bbox_intersections(
            ACTION& action ) const
        {
            intersect_recursive( action, 1, 0, mesh_.nb_simplices(), 1, 0, mesh_.nb_simplices() ) ;
        }

        void compute_bbox_intersections(
            const Box3d& box, std::vector< uint32 >& results ) const
    {
            intersect_recursive( results, box, 1, 0, mesh_.nb_simplices() ) ;
        }

        /**
         * @return the index of the facet nearest to point p.
         * @param [in] p query point
         * @param [out] nearest_point nearest point on the surface
         * @param [out] sq_dist squared distance between p and the surface.
         */
        uint32 nearest_facet(
            const float64* p,
            float64* nearest_point,
            float64& sq_dist ) const
        {
            vec3 nearest( nearest_point ) ;
            vec3 query( p ) ;
            uint32 nearest_t = nearest_facet( query, nearest, sq_dist ) ;
            nearest_point = nearest.data() ;
            return nearest_t ;

        }
        uint32 nearest_facet(
            const vec3& p,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            uint32 nearest_t ;
            get_nearest_facet_hint( p, nearest_t, nearest_point, sq_dist ) ;
            nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist, 1, 0,
                mesh_.nb_simplices() ) ;
            return nearest_t ;
        }

        static float64 get_nearest_point(
            const SurfacePart& M,
            const vec3& p,
            int t,
            vec3& nearest_p ) ;

        uint32 nearest_facet( const vec3& p ) const
        {
            vec3 nearest_point ;
            float64 dist ;
            return nearest_facet( p, nearest_point, dist ) ;
        }

        uint32 nearest_facet( const float64* p ) const
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
            uint32& nearest_t,
            vec3& nearest_point,
            float64& sq_dist ) const
        {
            nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist, 1, 0,
                mesh_.nb_simplices() ) ;
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

        void init_bboxes_recursive( uint32 node, uint32 b, uint32 e ) ;

        /**
         * @brief computes all the pairs of intersecting facets
         *  for two sub-trees of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * @param action ACTION::operator(uint32,uint32) is
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
            uint32 node1,
            uint32 b1,
            uint32 e1,
            uint32 node2,
            uint32 b2,
            uint32 e2 ) const
        {
            grgmesh_debug_assert( e1 != b1 ) ;
            grgmesh_debug_assert( e2 != b2 ) ;

            // Since we are intersecting the AABBTree with *itself*,
            // we can prune half of the cases by skipping the test
            // whenever node2's facet index interval is greated than
            // node1's facet index interval.
            if( e2 <= b1 ) {
                return ;
            }

            // The acceleration is here:
            if( !bboxes_[node1].bboxes_overlap( bboxes_[node2] ) ) {
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
                uint32 m2 = b2 + ( e2 - b2 ) / 2 ;
                uint32 node2_l = 2 * node2 ;
                uint32 node2_r = 2 * node2 + 1 ;
                intersect_recursive( action, node1, b1, e1, node2_l, b2, m2 ) ;
                intersect_recursive( action, node1, b1, e1, node2_r, m2, e2 ) ;
            } else {
                uint32 m1 = b1 + ( e1 - b1 ) / 2 ;
                uint32 node1_l = 2 * node1 ;
                uint32 node1_r = 2 * node1 + 1 ;
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
         * @param action ACTION::operator(uint32) is
         *   invoked of all pairs of facets that have overlapping
         *   bounding boxes.
         * @param box   input box to test
         * @param node  index of the node of the AABB tree
         * @param b     index of the facet in node1
         * @param e     one position past the index of the last facet in node1
         */
        void intersect_recursive(
            std::vector< uint32 >& results,
            const Box3d& box,
            uint32 node,
            uint32 b,
            uint32 e ) const
        {
            grgmesh_debug_assert( e != b ) ;

            // The acceleration is here:
            if( !box.bboxes_overlap( bboxes_[node] ) ) {
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

            uint32 m = b + ( e - b ) / 2 ;
            uint32 node_l = 2 * node ;
            uint32 node_r = 2 * node + 1 ;
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
            uint32& nearest_t,
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
            uint32& nearest_t,
            vec3& nearest_point,
            float64& sq_dist,
            uint32 n,
            uint32 b,
            uint32 e ) const ;

    private:
        std::vector< Box3d > bboxes_ ;
        SurfacePart& mesh_ ;
    } ;
}

#endif
