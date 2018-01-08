/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/basic/aabb.h>

#include <algorithm>
#include <numeric>

#include <ringmesh/basic/geometry.h>

/// Copied and adapted from Geogram

namespace
{
    using namespace RINGMesh;

    using const_vector_itr = const std::vector< index_t >::iterator;

    template < index_t DIMENSION >
    class Morton_cmp
    {
    public:
        Morton_cmp(
            const std::vector< Box< DIMENSION > >& bboxes, index_t coord )
            : bboxes_( bboxes ), coord_( coord )
        {
        }

        bool operator()( index_t box1, index_t box2 )
        {
            return bboxes_[box1].center()[coord_]
                   < bboxes_[box2].center()[coord_];
        }

    private:
        const std::vector< Box< DIMENSION > >& bboxes_;
        index_t coord_;
    };

    ALIAS_2D_AND_3D( Morton_cmp );

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
    template < class CMP >
    const_vector_itr split(
        const_vector_itr& begin, const_vector_itr& end, CMP cmp )
    {
        if( begin >= end )
        {
            return begin;
        }
        const_vector_itr middle = begin + ( end - begin ) / 2;
        std::nth_element( begin, middle, end, cmp );
        return middle;
    }

    /**
     * \brief Generic class for sorting arbitrary elements in Morton order.
     * \details The implementation is inspired by:
     *  - Christophe Delage and Olivier Devillers. Spatial Sorting.
     *   In CGAL User and Reference Manual. CGAL Editorial Board,
     *   3.9 edition, 2011
     */
    template < index_t DIMENSION >
    struct MortonSort
    {
        template < index_t COORDX >
        void sort( const std::vector< Box< DIMENSION > >& bboxes,
            const_vector_itr& begin,
            const_vector_itr& end );

        MortonSort( const std::vector< Box< DIMENSION > >& bboxes,
            std::vector< index_t >& mapping_morton )
        {
            sort< 0 >( bboxes, mapping_morton.begin(), mapping_morton.end() );
        }
    };

    template <>
    template < index_t COORDX >
    void MortonSort< 3 >::sort( const std::vector< Box3D >& bboxes,
        const_vector_itr& begin,
        const_vector_itr& end )
    {
        if( end - begin <= 1 )
        {
            return;
        }
        const index_t COORDY{ ( COORDX + 1 ) % 3 };
        const index_t COORDZ{ ( COORDY + 1 ) % 3 };

        auto m0 = begin;
        auto m8 = end;
        auto m4 = split( m0, m8, Morton_cmp3D( bboxes, COORDX ) );
        auto m2 = split( m0, m4, Morton_cmp3D( bboxes, COORDY ) );
        auto m1 = split( m0, m2, Morton_cmp3D( bboxes, COORDZ ) );
        auto m3 = split( m2, m4, Morton_cmp3D( bboxes, COORDZ ) );
        auto m6 = split( m4, m8, Morton_cmp3D( bboxes, COORDY ) );
        auto m5 = split( m4, m6, Morton_cmp3D( bboxes, COORDZ ) );
        auto m7 = split( m6, m8, Morton_cmp3D( bboxes, COORDZ ) );
        sort< COORDZ >( bboxes, m0, m1 );
        sort< COORDY >( bboxes, m1, m2 );
        sort< COORDY >( bboxes, m2, m3 );
        sort< COORDX >( bboxes, m3, m4 );
        sort< COORDX >( bboxes, m4, m5 );
        sort< COORDY >( bboxes, m5, m6 );
        sort< COORDY >( bboxes, m6, m7 );
        sort< COORDZ >( bboxes, m7, m8 );
    }

    template <>
    template < index_t COORDX >
    void MortonSort< 2 >::sort( const std::vector< Box2D >& bboxes,
        const_vector_itr& begin,
        const_vector_itr& end )
    {
        if( end - begin <= 1 )
        {
            return;
        }
        const index_t COORDY = ( COORDX + 1 ) % 2;

        auto m0 = begin;
        auto m4 = end;
        auto m2 = split( m0, m4, Morton_cmp2D( bboxes, COORDX ) );
        auto m1 = split( m0, m2, Morton_cmp2D( bboxes, COORDY ) );
        auto m3 = split( m2, m4, Morton_cmp2D( bboxes, COORDY ) );
        sort< COORDY >( bboxes, m0, m1 );
        sort< COORDX >( bboxes, m1, m2 );
        sort< COORDX >( bboxes, m2, m3 );
        sort< COORDY >( bboxes, m3, m4 );
    }

    template < index_t DIMENSION >
    void morton_sort( const std::vector< Box< DIMENSION > >& bboxes,
        std::vector< index_t >& mapping_morton )
    {
        mapping_morton.resize( bboxes.size() );
        std::iota( mapping_morton.begin(), mapping_morton.end(), 0 );
        MortonSort< DIMENSION >( bboxes, mapping_morton );
    }
} // namespace

/****************************************************************************/

namespace RINGMesh
{
    template < index_t DIMENSION >
    void AABBTree< DIMENSION >::initialize_tree(
        const std::vector< Box< DIMENSION > >& bboxes )
    {
        morton_sort( bboxes, mapping_morton_ );
        auto nb_bboxes = static_cast< index_t >( bboxes.size() );
        tree_.resize( max_node_index( ROOT_INDEX, 0, nb_bboxes ) + ROOT_INDEX );
        initialize_tree_recursive( bboxes, ROOT_INDEX, 0, nb_bboxes );
    }

    template < index_t DIMENSION >
    index_t AABBTree< DIMENSION >::max_node_index(
        index_t node_index, index_t box_begin, index_t box_end )
    {
        ringmesh_assert( box_end > box_begin );
        if( is_leaf( box_begin, box_end ) )
        {
            return node_index;
        }
        index_t element_middle;
        index_t child_left;
        index_t child_right;
        get_recursive_iterators( node_index, box_begin, box_end, element_middle,
            child_left, child_right );
        return std::max(
            max_node_index( child_left, box_begin, element_middle ),
            max_node_index( child_right, element_middle, box_end ) );
    }

    /**
     * \brief Computes the hierarchy of bounding boxes recursively.
     * \param[in] bboxes the array of bounding boxes
     * \param[in] node_index the index of the root of the subtree
     * \param[in] element_begin first box index in the vector \p bboxes
     * \param[in] element_end one position past the last box index in the vector
     * \p
     * bboxes
     */
    template < index_t DIMENSION >
    void AABBTree< DIMENSION >::initialize_tree_recursive(
        const std::vector< Box< DIMENSION > >& bboxes,
        index_t node_index,
        index_t element_begin,
        index_t element_end )
    {
        ringmesh_assert( node_index < tree_.size() );
        ringmesh_assert( element_begin != element_end );
        if( is_leaf( element_begin, element_end ) )
        {
            node( node_index ) = bboxes[mapping_morton_[element_begin]];
            return;
        }
        index_t element_middle;
        index_t child_left;
        index_t child_right;
        get_recursive_iterators( node_index, element_begin, element_end,
            element_middle, child_left, child_right );
        ringmesh_assert( child_left < tree_.size() );
        ringmesh_assert( child_right < tree_.size() );
        initialize_tree_recursive(
            bboxes, child_left, element_begin, element_middle );
        initialize_tree_recursive(
            bboxes, child_right, element_middle, element_end );
        node( node_index ) =
            node( child_left ).bbox_union( node( child_right ) );
    }

    template < index_t DIMENSION >
    std::tuple< index_t, vecn< DIMENSION >, double >
        AABBTree< DIMENSION >::get_nearest_element_box_hint(
            const vecn< DIMENSION >& query ) const
    {
        index_t box_begin = 0;
        index_t box_end = nb_bboxes();
        index_t node_index = ROOT_INDEX;
        while( !is_leaf( box_begin, box_end ) )
        {
            index_t box_middle;
            index_t child_left;
            index_t child_right;
            get_recursive_iterators( node_index, box_begin, box_end, box_middle,
                child_left, child_right );
            if( length2( node( child_left ).center() - query )
                < length2( node( child_right ).center() - query ) )
            {
                box_end = box_middle;
                node_index = child_left;
            }
            else
            {
                box_begin = box_middle;
                node_index = child_right;
            }
        }

        index_t nearest_box = mapping_morton_[box_begin];
        vecn< DIMENSION > nearest_point =
            get_point_hint_from_box( tree_[box_begin], nearest_box );
        double distance = length( query - nearest_point );
        return std::make_tuple( nearest_box, nearest_point, distance );
    }

    template < index_t DIMENSION >
    BoxAABBTree< DIMENSION >::BoxAABBTree(
        const std::vector< Box< DIMENSION > >& bboxes )
    {
        this->initialize_tree( bboxes );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > BoxAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box, index_t element_id ) const
    {
        ringmesh_unused( element_id );
        return box.center();
    }

    template < index_t DIMENSION >
    double inner_point_box_distance(
        const vecn< DIMENSION >& p, const Box< DIMENSION >& B )
    {
        ringmesh_assert( B.contains( p ) );
        double result = std::abs( p[0] - B.min()[0] );
        result = std::min( result, std::abs( p[0] - B.max()[0] ) );
        for( auto c : range( 1, DIMENSION ) )
        {
            result = std::min( result, std::abs( p[c] - B.min()[c] ) );
            result = std::min( result, std::abs( p[c] - B.max()[c] ) );
        }
        return result;
    }

    template < index_t DIMENSION >
    double point_box_signed_distance(
        const vecn< DIMENSION >& p, const Box< DIMENSION >& B )
    {
        bool inside = true;
        vecn< DIMENSION > result;
        for( auto c : range( DIMENSION ) )
        {
            if( p[c] < B.min()[c] )
            {
                inside = false;
                result[c] = p[c] - B.min()[c];
            }
            else if( p[c] > B.max()[c] )
            {
                inside = false;
                result[c] = p[c] - B.max()[c];
            }
        }
        if( inside )
        {
            return -inner_point_box_distance( p, B );
        }
        return result.length();
    }

    template double basic_api inner_point_box_distance(
        const vec2&, const Box2D& );
    template double basic_api point_box_signed_distance(
        const vec2&, const Box2D& );
    template class basic_api AABBTree< 2 >;
    template class basic_api BoxAABBTree< 2 >;

    template double basic_api inner_point_box_distance(
        const vec3&, const Box3D& );
    template double basic_api point_box_signed_distance(
        const vec3&, const Box3D& );
    template class basic_api AABBTree< 3 >;
    template class basic_api BoxAABBTree< 3 >;
} // namespace RINGMesh
