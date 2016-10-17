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

#include <ringmesh/mesh/aabb.h>

#include <ringmesh/mesh/mesh.h>

/// Copied and adapted from Geogram

namespace {

    using namespace RINGMesh ;

    typedef std::vector< index_t >::iterator vector_itr ;

    /**
     * \brief Computes the axis-aligned bounding box of a mesh cell
     * \param[in] mesh the mesh
     * \param[in] element_id the index of the element in mesh \p mesh
     * \param[out] bbox the bounding box of the facet
     */
    void get_element_bbox( const MeshBase& mesh, index_t element_id, Box3d& bbox )
    {
        for( index_t v = 0; v < mesh.nb_mesh_element_vertices( element_id ); v++ ) {
            const vec3& point = mesh.mesh_element_vertex( element_id, v ) ;
            bbox.add_point( point ) ;
        }
    }

    /**
     * \brief Computes the maximum node index in a subtree
     * \param[in] node_index node index of the root of the subtree
     * \param[in] box_begin first box index
     * \param[in] box_end one position past the last box index
     * \return the maximum node index in the subtree rooted at \p node_index
     */
    index_t max_node_index( index_t node_index, index_t box_begin, index_t box_end )
    {
        ringmesh_assert( box_end > box_begin ) ;
        if( box_begin + 1 == box_end ) {
            return node_index ;
        }
        index_t element_middle = box_begin + ( box_end - box_begin ) / 2 ;
        index_t child_left = 2 * node_index ;
        index_t child_right = 2 * node_index + 1 ;
        return std::max( max_node_index( child_left, box_begin, element_middle ),
            max_node_index( child_right, element_middle, box_end ) ) ;
    }

    template< index_t COORD >
    class Morton_cmp {
    public:
        Morton_cmp(const std::vector< Box3d >& bboxes) :
            bboxes_(bboxes) {
        }

        bool operator()( index_t box1, index_t box2 )
        {
            return bboxes_[box1].center()[COORD] < bboxes_[box2].center()[COORD] ;
        }

    private:
        const std::vector< Box3d >& bboxes_;
    };

    /**
     * \brief Splits a sequence into two ordered halves.
     * \details The algorithm shuffles the sequence and
     *  partitions its into two halves with the same number of elements
     *  and such that the elements of the first half are smaller
     *  than the elements of the second half.
     * \param[in] begin an iterator to the first element
     * \param[in] end an iterator one position past the last element
     * \param[in] cmp the comparator object
     * \return an iterator to the middle of the sequence that separates
     *  the two halves
     */
    template< class CMP >
    inline vector_itr split( vector_itr begin, vector_itr end, CMP cmp )
    {
        if( begin >= end ) {
            return begin ;
        }
        vector_itr middle = begin + ( end - begin ) / 2 ;
        std::nth_element( begin, middle, end, cmp ) ;
        return middle ;
    }

    /**
     * \brief Generic class for sorting arbitrary elements in
     *  Hilbert and Morton orders.
     * \details The implementation is inspired by:
     *  - Christophe Delage and Olivier Devillers. Spatial Sorting.
     *   In CGAL User and Reference Manual. CGAL Editorial Board,
     *   3.9 edition, 2011
     * \tparam CMP the comparator class for ordering the elements. CMP
     *  is itself a template parameterized by~:
     *    - COORD the coordinate along which elements should be
     *      sorted
     */
    template< template< index_t COORD > class CMP >
    struct HilbertSort {

        /**
         * \brief Low-level recursive spatial sorting function
         * \details This function is recursive
         * \param[in] M the mesh in which the elements reside
         * \param[in] begin an iterator that points to the
         *  first element of the sequence
         * \param[in] end an interator that points one position past the
         *  last element of the sequence
         */
        template< index_t COORDX >
        static void sort(
            const std::vector< Box3d >& bboxes,
            vector_itr begin,
            vector_itr end )
        {
            if( end - begin <= 1 ) {
                return ;
            }
            const index_t COORDY = (COORDX + 1) % 3, COORDZ = (COORDY + 1) % 3;

            vector_itr m0 = begin, m8 = end ;
            vector_itr m4 = split( m0, m8, CMP< COORDX >( bboxes ) ) ;
            vector_itr m2 = split( m0, m4, CMP< COORDY >( bboxes ) ) ;
            vector_itr m1 = split( m0, m2, CMP< COORDZ >( bboxes ) ) ;
            vector_itr m3 = split( m2, m4, CMP< COORDZ >( bboxes ) ) ;
            vector_itr m6 = split( m4, m8, CMP< COORDY >( bboxes ) ) ;
            vector_itr m5 = split( m4, m6, CMP< COORDZ >( bboxes ) ) ;
            vector_itr m7 = split( m6, m8, CMP< COORDZ >( bboxes ) ) ;
            sort< COORDZ >( bboxes, m0, m1 ) ;
            sort< COORDY >( bboxes, m1, m2 ) ;
            sort< COORDY >( bboxes, m2, m3 ) ;
            sort< COORDX >( bboxes, m3, m4 ) ;
            sort< COORDX >( bboxes, m4, m5 ) ;
            sort< COORDY >( bboxes, m5, m6 ) ;
            sort< COORDY >( bboxes, m6, m7 ) ;
            sort< COORDZ >( bboxes, m7, m8 ) ;
        }

        /**
         * \brief Sorts a sequence of boxes spatially.
         * \details This function does an indirect sort,
         *  in the sense that a sequence
         *  of indices that refer to the elements is sorted.
         * \param[in] bboxes the array of bounding boxes
         * \param[in] box_begin an interator to the first index to be sorted
         * \param[in] box_end an interator one position past the last index
         *  to be sorted
         */
        HilbertSort(
            const std::vector< Box3d >& bboxes,
            vector_itr box_begin,
            vector_itr box_end )
        {
            ringmesh_assert( box_end > box_begin ) ;
            sort< 0 >( bboxes, box_begin, box_end ) ;
        }
    } ;

    void morton_sort(
        const std::vector< Box3d >& bboxes,
        std::vector< index_t >& mapping_morton )
    {
        mapping_morton.resize( bboxes.size() ) ;
        for( index_t i = 0; i < bboxes.size(); i++ ) {
            mapping_morton[i] = i ;
        }
        HilbertSort< Morton_cmp >( bboxes, mapping_morton.begin(), mapping_morton.end() ) ;
    }

}

/****************************************************************************/

namespace RINGMesh {

    AABBTree::AABBTree( const std::vector< Box3d >& bboxes )
    {
        morton_sort( bboxes, mapping_morton_ ) ;
        tree_.resize( max_node_index( 1, 0, bboxes.size() ) + 1 ) ;
        init_bboxes_recursive( bboxes, 1, 0, bboxes.size() ) ;
    }

    /**
     * \brief Computes the hierarchy of bounding boxes recursively.
     * \param[in] bboxes the array of bounding boxes
     * \param[in] node_index the index of the root of the subtree
     * \param[in] box_begin first box index in the vector \p bboxes
     * \param[in] box_end one position past the last box index in the vector \p bboxes
     */
    void AABBTree::init_bboxes_recursive(
        const std::vector< Box3d >& bboxes,
        index_t node_index,
        index_t box_begin,
        index_t box_end )
    {
        ringmesh_assert( node_index < bboxes.size() ) ;
        ringmesh_assert( box_begin != box_end ) ;
        if( box_begin + 1 == box_end ) {
            tree_[node_index] = bboxes[box_begin] ;
            return ;
        }
        index_t element_middle = box_begin + ( box_end - box_begin ) / 2 ;
        index_t child_left = 2 * node_index ;
        index_t child_right = 2 * node_index + 1 ;
        ringmesh_assert( child_left < tree_.size() ) ;
        ringmesh_assert( child_right < tree_.size() ) ;
        init_bboxes_recursive( bboxes, child_left, box_begin, element_middle ) ;
        init_bboxes_recursive( bboxes, child_right, element_middle, box_end ) ;
        ringmesh_assert( child_left < bboxes.size() ) ;
        ringmesh_assert( child_right < bboxes.size() ) ;
        tree_[node_index] = tree_[child_left].bbox_union( tree_[child_right] ) ;
    }

}

