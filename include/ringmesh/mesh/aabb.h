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

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/basic/box3d.h>

namespace RINGMesh {
    class MeshBase;
    class LineMesh;
    class SurfaceMesh;
    class VolumeMesh;
}

namespace RINGMesh {

    /*!
     * @brief AABB tree structure
     * @details The tree is store in s single vector following this example:
     *                          ROOT
     *                        /      \
     *                      A1        A2
     *                    /    \     /   \
     *                  B1     B2   B3    B4
     *  where B* are the input bboxes
     *  Storage: |empty|ROOT|A1|A2|B1|B2|B3|B4|
     */
    class RINGMESH_API AABBTree {
    public:
        /// The index where to store the root. It starts to one for algorithm trick.
        static const index_t ROOT_INDEX = 1;

        /*!
         * @brief Saves the tree in a set of files
         * @details Each level of the tree is saved in a .geogram file
         * prefixed by \p name
         * @param[in] name the prefix used for the file naming
         */
        void save_tree( const std::string& name ) const;
        index_t nb_bboxes() const
        {
            return static_cast< index_t >( mapping_morton_.size() );
        }

        /*!
         * @brief Gets the closest element box to a point
         * @param[in] query the point to test
         * @param[in] action the functor to compute the distance between
         * the \p query and the tree element boxes
         * @return a tuple containing:
         * - the index of the closest element box.
         * - the nearest point on the element box.
         * - the distance between the \p query.
         * and \p nearest_point.
         * @tparam EvalDistance this functor should have an operator() defined like this:
         *  void operator()(
         *      const vec3& query,
         *      index_t cur_box,
         *      vec3& nearest_point,
         *      double& distance ) const ;
         * where query is the same than \p query, cur_box is the element box index
         * (e.g. in the case of AABBTree2D, this index is a polygon index) and nearest_point
         * and distance are the value computed using the element in the \p cur_box.
         */
        template< typename EvalDistance >
        std::tuple< index_t, vec3, double > closest_element_box(
            const vec3& query,
            const EvalDistance& action ) const
        {
            index_t nearest_box = NO_ID;
            vec3 nearest_point;
            double distance;
            std::tie( nearest_box, nearest_point, distance ) =
                get_nearest_element_box_hint( query );
            closest_element_box_recursive< EvalDistance >( query, nearest_box,
                nearest_point, distance, ROOT_INDEX, 0, nb_bboxes(), action );
            ringmesh_assert( nearest_box != NO_ID );
            return std::make_tuple( nearest_box, nearest_point, distance );
        }
        /*
         * @brief Computes the intersections between a given
         *  box and the element boxes.
         * @param[in] box the box to test
         * @param[in] action The functor used to compute intersection
         * with the element boxes when they intersect \p box
         * @tparam EvalIntersection this functor should have an operator() defined like this:
         *  void operator()( index_t cur_box ) ;
         * where cur_box is the element box index
         * (e.g. in the case of AABBTree2D, this index is a polygon index)
         */
        template< class EvalIntersection >
        void compute_bbox_element_bbox_intersections(
            const Box3d& box,
            EvalIntersection& action ) const
        {
            bbox_intersect_recursive< EvalIntersection >( box, ROOT_INDEX, 0,
                nb_bboxes(), action );
        }
        /*
         * @brief Computes the self intersections of the element boxes.
         * @param[in] action The functor used to compute intersection
         * with the intersected element boxes
         * @tparam EvalIntersection this functor should have an operator() defined like this:
         *  void operator()( index_t box1, index_t box2 ) ;
         * where box1 and box2 are the element box indices
         * (e.g. in the case of AABBTree2D, this index is a polygon index)
         */
        template< class EvalIntersection >
        void compute_self_element_bbox_intersections(
            EvalIntersection& action ) const
        {
            self_intersect_recursive< EvalIntersection >( ROOT_INDEX, 0, nb_bboxes(),
                ROOT_INDEX, 0, nb_bboxes(), action );
        }
    protected:
        virtual ~AABBTree() = default;

        /*!
         * @brief Builds the tree
         * @details Comptes the morton order and build the tree
         * using the ordered bboxes
         * @param[in] bboxes the set of unordered bboxes
         */
        void initialize_tree( const std::vector< Box3d >& bboxes );

        bool is_leaf( index_t box_begin, index_t box_end ) const
        {
            return box_begin + 1 == box_end;
        }
        void get_recursive_iterators(
            index_t node_index,
            index_t box_begin,
            index_t box_end,
            index_t& middle_box,
            index_t& child_left,
            index_t& child_right ) const
        {
            middle_box = box_begin + ( box_end - box_begin ) / 2;
            child_left = 2 * node_index;
            child_right = 2 * node_index + 1;
        }

    private:
        /*!
         * @brief Gets the number of nodes in the tree subset
         */
        index_t max_node_index(
            index_t node_index,
            index_t box_begin,
            index_t box_end );
        /*!
         * @brief The recursive instruction used in initialize_tree()
         */
        void initialize_tree_recursive(
            const std::vector< Box3d >& bboxes,
            index_t node_index,
            index_t element_begin,
            index_t element_end );

        /*!
         * @brief The recursive instruction used in closest_element_box()
         */
        template< typename ACTION >
        void closest_element_box_recursive(
            const vec3& query,
            index_t& nearest_box,
            vec3& nearest_point,
            double& distance,
            index_t node_index,
            index_t element_begin,
            index_t element_end,
            const ACTION& action ) const;

        template< class ACTION >
        void bbox_intersect_recursive(
            const Box3d& box,
            index_t node_index,
            index_t element_begin,
            index_t element_end,
            ACTION& action ) const;

        template< class ACTION >
        void self_intersect_recursive(
            index_t node_index1,
            index_t element_begin1,
            index_t element_end1,
            index_t node_index2,
            index_t element_begin2,
            index_t element_end2,
            ACTION& action ) const;

        /*!
         * @brief Gets an hint of the result
         * @details Compute the result by approximating each bbox to its barycenter.
         * This result is then used to speed-up the computation by minimizing
         * the distance computation between \p query and the real elements
         * inside the bboxes
         */
        std::tuple< index_t, vec3, double > get_nearest_element_box_hint(
            const vec3& query ) const;
        /*!
         * @brief Gets an element point from its box
         * @details This function is used to get a result from the selected hint box
         */
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const = 0;

    protected:
        std::vector< Box3d > tree_;
        std::vector< index_t > mapping_morton_;
    };

    class RINGMESH_API BoxAABBTree: public AABBTree {
    public:
        BoxAABBTree( const std::vector< Box3d >& boxes );
        virtual ~BoxAABBTree() = default;

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the barycenter of the box
         */
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const override;
    };

    class RINGMESH_API LineAABBTree: public AABBTree {
    public:
        LineAABBTree( const LineMesh& mesh );
        virtual ~LineAABBTree() = default;

        /*!
         * @brief Gets the closest edge to a given point
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest edge index.
         * - nearest_point the nearest point on the closest edge.
         * - distance the distance between \p query and \p nearest_point.
         */
        std::tuple< index_t, vec3, double > closest_edge( const vec3& query ) const;
    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const override;
        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and an edge
         */
        class DistanceToEdge {
        public:
            DistanceToEdge( const LineMesh& mesh )
                : mesh_( mesh )
            {
            }

            void operator()(
                const vec3& query,
                index_t cur_box,
                vec3& nearest_point,
                double& distance ) const;

        private:
            const LineMesh& mesh_;
        };

    private:
        const LineMesh& mesh_;
    };

    class RINGMESH_API SurfaceAABBTree: public AABBTree {
    public:
        SurfaceAABBTree( const SurfaceMesh& mesh );
        virtual ~SurfaceAABBTree() = default;

        /*!
         * @brief Gets the closest triangle to a given point
         * @pre The mesh needs to be triangulated
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest triangle index.
         * - the nearest point on the closest triangle.
         * - the distance between \p query and \p nearest_point.
         */
        std::tuple< index_t, vec3, double > closest_triangle(
            const vec3& query ) const;
    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const override;
        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and a triangle
         */
        class DistanceToTriangle {
        public:
            DistanceToTriangle( const SurfaceMesh& mesh )
                : mesh_( mesh )
            {
            }

            void operator()(
                const vec3& query,
                index_t cur_box,
                vec3& nearest_point,
                double& distance ) const;

        private:
            const SurfaceMesh& mesh_;
        };

    private:
        const SurfaceMesh& mesh_;
    };

    class RINGMESH_API VolumeAABBTree: public AABBTree {
    public:
        VolumeAABBTree( const VolumeMesh& mesh );
        virtual ~VolumeAABBTree() = default;

        /*!
         * @brief Gets the cell contining a point
         * @param[in] query the point to use
         * @return the cell index containing \p query,
         * NO_ID if no cell is corresponding
         */
        index_t containing_cell( const vec3& query ) const;

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const override;
        index_t containing_cell_recursive(
            const vec3& query,
            index_t node_index,
            index_t box_begin,
            index_t box_end ) const;

    private:
        const VolumeMesh& mesh_;
    };

    double inner_point_box_distance( const vec3& p, const Box3d& B );

    double point_box_signed_distance( const vec3& p, const Box3d& B );

    template< typename ACTION >
    void AABBTree::closest_element_box_recursive(
        const vec3& query,
        index_t& nearest_box,
        vec3& nearest_point,
        double& distance,
        index_t node_index,
        index_t box_begin,
        index_t box_end,
        const ACTION& action ) const
    {
        ringmesh_assert( node_index < tree_.size() );
        ringmesh_assert( box_begin != box_end );

        // If node is a leaf: compute point-element distance
        // and replace current if nearer
        if( is_leaf( box_begin, box_end ) ) {
            index_t cur_box = mapping_morton_[box_begin];
            vec3 cur_nearest_point;
            double cur_distance;
            action( query, cur_box, cur_nearest_point, cur_distance );
            if( cur_distance < distance ) {
                nearest_box = cur_box;
                nearest_point = cur_nearest_point;
                distance = cur_distance;
            }
            return;
        }
        index_t box_middle, child_left, child_right;
        get_recursive_iterators( node_index, box_begin, box_end, box_middle,
            child_left, child_right );

        double distance_left = point_box_signed_distance( query, tree_[child_left] );
        double distance_right = point_box_signed_distance( query,
            tree_[child_right] );

        // Traverse the "nearest" child first, so that it has more chances
        // to prune the traversal of the other child.
        if( distance_left < distance_right ) {
            if( distance_left < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_left, box_begin, box_middle,
                    action );
            }
            if( distance_right < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_right, box_middle, box_end,
                    action );
            }
        } else {
            if( distance_right < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_right, box_middle, box_end,
                    action );
            }
            if( distance_left < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_left, box_begin, box_middle,
                    action );
            }
        }
    }

    template< typename ACTION >
    void AABBTree::bbox_intersect_recursive(
        const Box3d& box,
        index_t node_index,
        index_t element_begin,
        index_t element_end,
        ACTION& action ) const
    {
        ringmesh_assert( node_index < tree_.size() );
        ringmesh_assert( element_begin != element_end );

        // Prune sub-tree that does not have intersection
        if( !box.bboxes_overlap( tree_[node_index] ) ) {
            return;
        }

        // Leaf case
        if( is_leaf( element_begin, element_end ) ) {
            index_t cur_box = mapping_morton_[element_begin];
            action( cur_box );
            return;
        }

        index_t box_middle, child_left, child_right;
        get_recursive_iterators( node_index, element_begin, element_end, box_middle,
            child_left, child_right );

        bbox_intersect_recursive< ACTION >( box, child_left, element_begin,
            box_middle, action );
        bbox_intersect_recursive< ACTION >( box, child_right, box_middle,
            element_end, action );
    }

    template< class ACTION >
    void AABBTree::self_intersect_recursive(
        index_t node_index1,
        index_t element_begin1,
        index_t element_end1,
        index_t node_index2,
        index_t element_begin2,
        index_t element_end2,
        ACTION& action ) const
    {
        ringmesh_assert( element_end1 != element_begin1 );
        ringmesh_assert( element_end2 != element_begin2 );

        // Since we are intersecting the AABBTree with *itself*,
        // we can prune half of the cases by skipping the test
        // whenever node2's polygon index interval is greated than
        // node1's polygon index interval.
        if( element_end2 <= element_begin1 ) {
            return;
        }

        // The acceleration is here:
        if( !tree_[node_index1].bboxes_overlap( tree_[node_index2] ) ) {
            return;
        }

        // Simple case: leaf - leaf intersection.
        if( is_leaf( element_begin1, element_end1 )
            && is_leaf( element_begin2, element_end2 ) ) {
            action( mapping_morton_[element_begin1],
                mapping_morton_[element_begin2] );
            return;
        }

        // If node2 has more polygons than node1, then
        //   intersect node2's two children with node1
        // else
        //   intersect node1's two children with node2
        if( element_end2 - element_begin2 > element_end1 - element_begin1 ) {
            index_t middle_box2, child_left2, child_right2;
            get_recursive_iterators( node_index2, element_begin2, element_end2,
                middle_box2, child_left2, child_right2 );
            self_intersect_recursive< ACTION >( node_index1, element_begin1,
                element_end1, child_left2, element_begin2, middle_box2, action );
            self_intersect_recursive< ACTION >( node_index1, element_begin1,
                element_end1, child_right2, middle_box2, element_end2, action );
        } else {
            index_t middle_box1, child_left1, child_right1;
            get_recursive_iterators( node_index1, element_begin1, element_end1,
                middle_box1, child_left1, child_right1 );
            self_intersect_recursive< ACTION >( child_left1, element_begin1,
                middle_box1, node_index2, element_begin2, element_end2, action );
            self_intersect_recursive< ACTION >( child_right1, middle_box1,
                element_end1, node_index2, element_begin2, element_end2, action );
        }
    }
}
