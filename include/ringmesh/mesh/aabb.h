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

#ifndef RINGMESH_MESH_AABB
#define RINGMESH_MESH_AABB

#include <ringmesh/basic/common.h>

#include <ringmesh/basic/box3d.h>

namespace RINGMesh {
    class MeshBase ;
    class Mesh1D ;
}

namespace RINGMesh {

    class RINGMESH_API AABBTree {
    public:
        static const index_t ROOT_INDEX = 1 ;

        void save_tree( const std::string& name ) const ;
        index_t nb_bboxes() const
        {
            return static_cast< index_t >( mapping_morton_.size() ) ;
        }

        template< typename ACTION >
        index_t closest_element_box(
            const vec3& query,
            vec3& nearest_point,
            double& distance,
            const ACTION& action ) const
        {
            get_nearest_element_box_hint( query, nearest_point, distance ) ;
            index_t nearest_box ;
            closest_element_box_recursive< ACTION >( query, nearest_box,
                nearest_point, distance, ROOT_INDEX, 0, nb_bboxes(), action ) ;
            return nearest_box ;
        }
    protected:
        virtual ~AABBTree()
        {
        }
        void initialize_tree( const std::vector< Box3d >& bboxes ) ;

        bool is_leaf( index_t box_begin, index_t box_end ) const
        {
            return box_begin + 1 == box_end ;
        }
        void get_recursive_iterators(
            index_t node_index,
            index_t box_begin,
            index_t box_end,
            index_t& middle_box,
            index_t& child_left,
            index_t& child_right ) const
        {
            middle_box = box_begin + ( box_end - box_begin ) / 2 ;
            child_left = 2 * node_index ;
            child_right = 2 * node_index + 1 ;
        }

    private:
        index_t max_node_index(
            index_t node_index,
            index_t box_begin,
            index_t box_end ) ;
        void initialize_tree_recursive(
            const std::vector< Box3d >& bboxes,
            index_t node_index,
            index_t element_begin,
            index_t element_end ) ;

        template< typename ACTION >
        void closest_element_box_recursive(
            const vec3& query,
            index_t& nearest_box,
            vec3& nearest_point,
            double& distance,
            index_t node_index,
            index_t element_begin,
            index_t element_end,
            const ACTION& action ) const ;

        void get_nearest_element_box_hint(
            const vec3& query,
            vec3& nearest_point,
            double& distance ) const ;
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const = 0 ;

    private:
        std::vector< Box3d > tree_ ;
        std::vector< index_t > mapping_morton_ ;
    } ;

    class RINGMESH_API AABBTreeBox: public AABBTree {
    public:
        AABBTreeBox( const std::vector< Box3d >& boxes ) ;
        virtual ~AABBTreeBox()
        {
        }

    private:
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const ;
    } ;

    class RINGMESH_API AABBTreeMesh: public AABBTree {
    protected:
        AABBTreeMesh( const MeshBase& mesh ) ;
        virtual ~AABBTreeMesh()
        {
        }
    } ;

    class RINGMESH_API AABBTree1D: public AABBTreeMesh {
    public:
        AABBTree1D( const Mesh1D& mesh ) ;
        virtual ~AABBTree1D()
        {
        }

        index_t closest_edge(
            const vec3& query,
            vec3& nearest_point,
            double& distance ) const ;
    private:
        virtual vec3 get_point_hint_from_box(
            const Box3d& box,
            index_t element_id ) const ;

        class DistanceToEdge {
        public:
            DistanceToEdge( const Mesh1D& mesh )
                : mesh_( mesh )
            {
            }

            void operator()(
                const vec3& query,
                index_t cur_box,
                vec3& nearest_point,
                double& distance ) const ;

        private:
            const Mesh1D& mesh_ ;
        } ;

    private:
        const Mesh1D& mesh_ ;
    } ;

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
        ringmesh_assert( node_index < tree_.size() ) ;
        ringmesh_assert( box_begin != box_end ) ;

        // If node is a leaf: compute point-element distance
        // and replace current if nearer
        if( is_leaf( box_begin, box_end ) ) {
            index_t cur_box = mapping_morton_[box_begin] ;
            vec3 cur_nearest_point ;
            double cur_distance ;
            action( query, cur_box, cur_nearest_point, cur_distance ) ;
            if( cur_distance < distance ) {
                nearest_box = cur_box ;
                nearest_point = cur_nearest_point ;
                distance = cur_distance ;
            }
            return ;
        }
        index_t box_middle, child_left, child_right ;
        get_recursive_iterators( node_index, box_begin, box_end, box_middle,
            child_left, child_right ) ;

        double distance_left = length2( tree_[child_left].center() - query ) ;
        double distance_right = length2( tree_[child_right].center() - query ) ;

        // Traverse the "nearest" child first, so that it has more chances
        // to prune the traversal of the other child.
        if( distance_left < distance_right ) {
            if( distance_left < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_left, box_begin, box_middle,
                    action ) ;
            }
            if( distance_right < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_right, box_middle, box_end,
                    action ) ;
            }
        } else {
            if( distance_right < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_right, box_middle, box_end,
                    action ) ;
            }
            if( distance_left < distance ) {
                closest_element_box_recursive< ACTION >( query, nearest_box,
                    nearest_point, distance, child_left, box_begin, box_middle,
                    action ) ;
            }
        }
    }
}

#endif

