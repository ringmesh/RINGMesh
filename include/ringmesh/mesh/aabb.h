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

#include <ringmesh/basic/box.h>
#include <ringmesh/basic/common.h>

namespace RINGMesh {
    template< index_t DIMENSION > class MeshBase;
    template< index_t DIMENSION > class LineMesh;
    template< index_t DIMENSION > class SurfaceMeshBase;
    template< index_t DIMENSION > class VolumeMesh;
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
    template< index_t DIMENSION >
    class AABBTree {
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        /// The index where to store the root. It starts to one for algorithm trick.
        static const index_t ROOT_INDEX = 1;

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
         *  std::tuple< double, vecn< DIMENSION > > operator()(
         *      const vecn< DIMENSION >& query,
         *      index_t cur_box ) const ;
         * where query is the same than \p query, cur_box is the element box index
         * (e.g. in the case of SurfaceAABBTree, this index is a polygon index).
         * The returned tuple contains the distance of the nearest point and the
         * nearest point computed using the element in the \p cur_box.
         */
        template< typename EvalDistance >
        std::tuple< index_t, vecn< DIMENSION >, double > closest_element_box(
            const vecn< DIMENSION >& query,
            const EvalDistance& action ) const
        {
            index_t nearest_box = NO_ID;
            vecn< DIMENSION > nearest_point;
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
            const Box< DIMENSION > & box,
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
        void initialize_tree( const std::vector< Box< DIMENSION > >& bboxes );

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

        const Box< DIMENSION >& node( index_t i ) const
        {
            ringmesh_assert( i < tree_.size() );
            return tree_[i];
        }

        Box< DIMENSION >& node( index_t i )
        {
            ringmesh_assert( i < tree_.size() );
            return tree_[i];
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
            const std::vector< Box< DIMENSION > >& bboxes,
            index_t node_index,
            index_t element_begin,
            index_t element_end );

        /*!
         * @brief The recursive instruction used in closest_element_box()
         */
        template< typename ACTION >
        void closest_element_box_recursive(
            const vecn< DIMENSION > & query,
            index_t& nearest_box,
            vecn< DIMENSION > & nearest_point,
            double& distance,
            index_t node_index,
            index_t element_begin,
            index_t element_end,
            const ACTION& action ) const;

        template< class ACTION >
        void bbox_intersect_recursive(
            const Box< DIMENSION >& box,
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
        std::tuple< index_t, vecn< DIMENSION >, double > get_nearest_element_box_hint(
            const vecn< DIMENSION >& query ) const;
        /*!
         * @brief Gets an element point from its box
         * @details This function is used to get a result from the selected hint box
         */
        virtual vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box,
            index_t element_id ) const = 0;

    protected:
        std::vector< Box< DIMENSION > > tree_;
        std::vector< index_t > mapping_morton_;
    };

    template< index_t DIMENSION >
    class BoxAABBTree: public AABBTree< DIMENSION > {
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        BoxAABBTree( const std::vector< Box< DIMENSION > >& boxes );
        virtual ~BoxAABBTree() = default;

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the barycenter of the box
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box,
            index_t element_id ) const override;
    };

    CLASS_DIMENSION_ALIASES( BoxAABBTree );

    template< index_t DIMENSION >
    class LineAABBTree: public AABBTree< DIMENSION > {
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        LineAABBTree( const LineMesh< DIMENSION >& mesh );
        virtual ~LineAABBTree() = default;

        /*!
         * @brief Gets the closest edge to a given point
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest edge index.
         * - nearest_point the nearest point on the closest edge.
         * - distance the distance between \p query and \p nearest_point.
         */
        std::tuple< index_t, vecn< DIMENSION >, double > closest_edge( const vecn< DIMENSION >& query ) const;
    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box,
            index_t element_id ) const override;
        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and an edge
         */
        class DistanceToEdge {
        public:
            DistanceToEdge( const LineMesh< DIMENSION >& mesh )
                : mesh_( mesh )
            {
            }

            std::tuple< double, vecn< DIMENSION > > operator()(
                const vecn< DIMENSION >& query,
                index_t cur_box ) const;

        private:
            const LineMesh< DIMENSION >& mesh_;
        };

    private:
        const LineMesh< DIMENSION >& mesh_;
    };

    CLASS_DIMENSION_ALIASES( LineAABBTree );

    template< index_t DIMENSION >
    class SurfaceAABBTree: public AABBTree< DIMENSION > {
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        SurfaceAABBTree( const SurfaceMeshBase< DIMENSION >& mesh );
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
        std::tuple< index_t, vecn< DIMENSION >, double > closest_triangle(
            const vecn< DIMENSION >& query ) const;
    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box,
            index_t element_id ) const override;
        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and a triangle
         */
        class DistanceToTriangle {
        public:
            DistanceToTriangle( const SurfaceMeshBase< DIMENSION >& mesh )
                : mesh_( mesh )
            {
            }

            std::tuple< double, vecn< DIMENSION > > operator()(
                const vecn< DIMENSION >& query,
                index_t cur_box ) const;

        private:
            const SurfaceMeshBase< DIMENSION >& mesh_;
        };

    private:
        const SurfaceMeshBase< DIMENSION >& mesh_;
    };

    CLASS_DIMENSION_ALIASES( SurfaceAABBTree );

    template< index_t DIMENSION >
    class VolumeAABBTree: public AABBTree< DIMENSION > {
        static_assert( DIMENSION == 3, "DIMENSION template should be 3" );
    public:
        VolumeAABBTree( const VolumeMesh< DIMENSION >& mesh );
        virtual ~VolumeAABBTree() = default;

        /*!
         * @brief Gets the cell contining a point
         * @param[in] query the point to use
         * @return the cell index containing \p query,
         * NO_ID if no cell is corresponding
         */
        index_t containing_cell( const vecn< DIMENSION >& query ) const;

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box,
            index_t element_id ) const override;
        index_t containing_cell_recursive(
            const vecn< DIMENSION >& query,
            index_t node_index,
            index_t box_begin,
            index_t box_end ) const;

    private:
        const VolumeMesh< DIMENSION >& mesh_;
    };

    using VolumeAABBTree3D = VolumeAABBTree< 3 >;

    template< index_t DIMENSION >
    double inner_point_box_distance(
        const vecn< DIMENSION >& p,
        const Box< DIMENSION >& B );

    template< index_t DIMENSION >
    double point_box_signed_distance(
        const vecn< DIMENSION >& p,
        const Box< DIMENSION >& B );

    template< index_t DIMENSION >
    template< typename ACTION >
    void AABBTree< DIMENSION >::closest_element_box_recursive(
        const vecn< DIMENSION >& query,
        index_t& nearest_box,
        vecn< DIMENSION >& nearest_point,
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
            vecn< DIMENSION > cur_nearest_point;
            double cur_distance;
            std::tie( cur_distance, cur_nearest_point ) = action( query, cur_box );
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

        double distance_left = point_box_signed_distance( query,
            node( child_left ) );
        double distance_right = point_box_signed_distance( query,
            node( child_right ) );

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

    template< index_t DIMENSION >
    template< typename ACTION >
    void AABBTree< DIMENSION >::bbox_intersect_recursive(
        const Box< DIMENSION >& box,
        index_t node_index,
        index_t element_begin,
        index_t element_end,
        ACTION& action ) const
    {
        ringmesh_assert( node_index < tree_.size() );
        ringmesh_assert( element_begin != element_end );

        // Prune sub-tree that does not have intersection
        if( !box.bboxes_overlap( node( node_index ) ) ) {
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

    template< index_t DIMENSION >
    template< class ACTION >
    void AABBTree< DIMENSION >::self_intersect_recursive(
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
        if( !node( node_index1 ).bboxes_overlap( node( node_index2 ) ) ) {
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
