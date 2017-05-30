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

#include <ringmesh/mesh/aabb.h>

#include <geogram/mesh/mesh_io.h>

#include <ringmesh/basic/geometry.h>

#include <ringmesh/mesh/mesh.h>

/// Copied and adapted from Geogram

namespace {

    using namespace RINGMesh;

    typedef const std::vector< index_t >::iterator const_vector_itr;

    template< index_t COORD, index_t DIMENSION >
    class Morton_cmp {
    public:
        Morton_cmp( const std::vector< Box< DIMENSION > >& bboxes )
            : bboxes_( bboxes )
        {
        }

        bool operator()( index_t box1, index_t box2 )
        {
            return bboxes_[box1].center()[COORD] < bboxes_[box2].center()[COORD];
        }

    private:
        const std::vector< Box< DIMENSION > >& bboxes_;
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
    inline const_vector_itr split(
        const_vector_itr& begin,
        const_vector_itr& end,
        CMP cmp )
    {
        if( begin >= end ) {
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
     * \tparam CMP the comparator class for ordering the elements. CMP
     *  is itself a template parameterized by:
     *    - the coordinate along which elements should be sorted
     *    - the dimension of sorted elements
     */
    template< index_t DIMENSION, template< index_t, index_t > class CMP >
    struct MortonSort {

        template< index_t COORDX >
        static void sort(
            const std::vector< Box< DIMENSION > >& bboxes,
            const_vector_itr& begin,
            const_vector_itr& end )
        {
            if( end - begin <= 1 ) {
                return;
            }
            const index_t COORDY = ( COORDX + 1 ) % 3, COORDZ = ( COORDY + 1 ) % 3;

            const_vector_itr m0 = begin, m8 = end;
            const_vector_itr m4 = split( m0, m8,
                CMP< COORDX, DIMENSION >( bboxes ) );
            const_vector_itr m2 = split( m0, m4,
                CMP< COORDY, DIMENSION >( bboxes ) );
            const_vector_itr m1 = split( m0, m2,
                CMP< COORDZ, DIMENSION >( bboxes ) );
            const_vector_itr m3 = split( m2, m4,
                CMP< COORDZ, DIMENSION >( bboxes ) );
            const_vector_itr m6 = split( m4, m8,
                CMP< COORDY, DIMENSION >( bboxes ) );
            const_vector_itr m5 = split( m4, m6,
                CMP< COORDZ, DIMENSION >( bboxes ) );
            const_vector_itr m7 = split( m6, m8,
                CMP< COORDZ, DIMENSION >( bboxes ) );
            sort< COORDZ >( bboxes, m0, m1 );
            sort< COORDY >( bboxes, m1, m2 );
            sort< COORDY >( bboxes, m2, m3 );
            sort< COORDX >( bboxes, m3, m4 );
            sort< COORDX >( bboxes, m4, m5 );
            sort< COORDY >( bboxes, m5, m6 );
            sort< COORDY >( bboxes, m6, m7 );
            sort< COORDZ >( bboxes, m7, m8 );
        }

        MortonSort(
            const std::vector< Box< DIMENSION > >& bboxes,
            std::vector< index_t >& mapping_morton )
        {
            sort< 0 >( bboxes, mapping_morton.begin(), mapping_morton.end() );
        }
    };

    template< index_t DIMENSION >
    void morton_sort(
        const std::vector< Box< DIMENSION > >& bboxes,
        std::vector< index_t >& mapping_morton )
    {
        mapping_morton.resize( bboxes.size() );
        for( index_t i = 0; i < bboxes.size(); i++ ) {
            mapping_morton[i] = i;
        }
        MortonSort< DIMENSION, Morton_cmp >( bboxes, mapping_morton );
    }

    template< index_t DIMENSION >
    void add_cube( GEO::Mesh& M, const Box< DIMENSION >& box, index_t n )
    {
        if( !box.initialized() ) return;
        const vecn< DIMENSION >& min_vertex = box.min();
        const vecn< DIMENSION >& max_vertex = box.max();
        vecn< DIMENSION > width( max_vertex[0] - min_vertex[0], 0, 0 );
        vecn< DIMENSION > height( 0, max_vertex[1] - min_vertex[1], 0 );
        vecn< DIMENSION > depth( 0, 0, max_vertex[2] - min_vertex[2] );
        index_t v0 = M.vertices.create_vertex( min_vertex.data() );
        index_t v1 = M.vertices.create_vertex(
            vecn< DIMENSION >( min_vertex + width ).data() );
        index_t v2 = M.vertices.create_vertex(
            vecn< DIMENSION >( max_vertex - depth ).data() );
        index_t v3 = M.vertices.create_vertex(
            vecn< DIMENSION >( min_vertex + height ).data() );
        index_t v4 = M.vertices.create_vertex(
            vecn< DIMENSION >( min_vertex + depth ).data() );
        index_t v5 = M.vertices.create_vertex(
            vecn< DIMENSION >( max_vertex - height ).data() );
        index_t v6 = M.vertices.create_vertex( max_vertex.data() );
        index_t v7 = M.vertices.create_vertex(
            vecn< DIMENSION >( max_vertex - width ).data() );

        GEO::Attribute< index_t > id( M.edges.attributes(), "id" );
        id[M.edges.create_edge( v0, v1 )] = n;
        id[M.edges.create_edge( v1, v2 )] = n;
        id[M.edges.create_edge( v2, v3 )] = n;
        id[M.edges.create_edge( v3, v0 )] = n;
        id[M.edges.create_edge( v4, v5 )] = n;
        id[M.edges.create_edge( v5, v6 )] = n;
        id[M.edges.create_edge( v6, v7 )] = n;
        id[M.edges.create_edge( v7, v4 )] = n;
        id[M.edges.create_edge( v0, v4 )] = n;
        id[M.edges.create_edge( v1, v5 )] = n;
        id[M.edges.create_edge( v2, v6 )] = n;
        id[M.edges.create_edge( v3, v7 )] = n;
    }

    template< index_t DIMENSION >
    bool mesh_cell_contains_point(
        const VolumeMesh< DIMENSION >& M,
        index_t cell,
        const vecn< DIMENSION >& p )
    {
        switch( M.cell_type( cell ) ) {
            case GEO::MESH_TET: {
                const vecn< DIMENSION >& p0 = M.vertex( M.cell_vertex( cell, 0 ) );
                const vecn< DIMENSION >& p1 = M.vertex( M.cell_vertex( cell, 1 ) );
                const vecn< DIMENSION >& p2 = M.vertex( M.cell_vertex( cell, 2 ) );
                const vecn< DIMENSION >& p3 = M.vertex( M.cell_vertex( cell, 3 ) );
                return point_inside_tetra( p, p0, p1, p2, p3 );
            }
            default:
                ringmesh_assert_not_reached;
                return false;
        }
    }
}

/****************************************************************************/

namespace RINGMesh {

    template< index_t DIMENSION >
    void AABBTree< DIMENSION >::initialize_tree(
        const std::vector< Box< DIMENSION > >& bboxes )
    {
        morton_sort( bboxes, mapping_morton_ );
        index_t nb_bboxes = static_cast< index_t >( bboxes.size() );
        tree_.resize( max_node_index( ROOT_INDEX, 0, nb_bboxes ) + ROOT_INDEX );
        initialize_tree_recursive( bboxes, ROOT_INDEX, 0, nb_bboxes );
    }

    template< index_t DIMENSION >
    index_t AABBTree< DIMENSION >::max_node_index(
        index_t node_index,
        index_t box_begin,
        index_t box_end )
    {
        ringmesh_assert( box_end > box_begin );
        if( is_leaf( box_begin, box_end ) ) {
            return node_index;
        }
        index_t element_middle, child_left, child_right;
        get_recursive_iterators( node_index, box_begin, box_end, element_middle,
            child_left, child_right );
        return std::max( max_node_index( child_left, box_begin, element_middle ),
            max_node_index( child_right, element_middle, box_end ) );
    }

    /**
     * \brief Computes the hierarchy of bounding boxes recursively.
     * \param[in] bboxes the array of bounding boxes
     * \param[in] node_index the index of the root of the subtree
     * \param[in] box_begin first box index in the vector \p bboxes
     * \param[in] box_end one position past the last box index in the vector \p bboxes
     */
    template< index_t DIMENSION >
    void AABBTree< DIMENSION >::initialize_tree_recursive(
        const std::vector< Box< DIMENSION > >& bboxes,
        index_t node_index,
        index_t box_begin,
        index_t box_end )
    {
        ringmesh_assert( node_index < tree_.size() );
        ringmesh_assert( box_begin != box_end );
        if( is_leaf( box_begin, box_end ) ) {
            tree_[node_index] = bboxes[mapping_morton_[box_begin]];
            return;
        }
        index_t element_middle, child_left, child_right;
        get_recursive_iterators( node_index, box_begin, box_end, element_middle,
            child_left, child_right );
        ringmesh_assert( child_left < tree_.size() );
        ringmesh_assert( child_right < tree_.size() );
        initialize_tree_recursive( bboxes, child_left, box_begin, element_middle );
        initialize_tree_recursive( bboxes, child_right, element_middle, box_end );
        tree_[node_index] = tree_[child_left].bbox_union( tree_[child_right] );
    }

    template< index_t DIMENSION >
    void AABBTree< DIMENSION >::save_tree( const std::string& name ) const
    {
        index_t nb_nodes = 0;
        for( double level = 1.; nb_nodes < tree_.size(); level++ ) {
            index_t start_node = static_cast< index_t >( std::pow( 2., level ) );
            nb_nodes = 2 * start_node;
            GEO::Mesh M;
            for( index_t n = start_node; n < nb_nodes; n++ ) {
                add_cube( M, tree_[n], n );
            }
            std::ostringstream oss;
            oss << name << level << ".geogram";
            GEO::mesh_save( M, oss.str() );
        }
    }

    template< index_t DIMENSION >
    void AABBTree< DIMENSION >::get_nearest_element_box_hint(
        const vecn< DIMENSION >& query,
        index_t& nearest_box,
        vecn< DIMENSION >& nearest_point,
        double& distance ) const
    {
        index_t box_begin = 0;
        index_t box_end = nb_bboxes();
        index_t node_index = ROOT_INDEX;
        while( !is_leaf( box_begin, box_end ) ) {
            index_t box_middle, child_left, child_right;
            get_recursive_iterators( node_index, box_begin, box_end, box_middle,
                child_left, child_right );
            if( length2( tree_[child_left].center() - query )
                < length2( tree_[child_right].center() - query ) ) {
                box_end = box_middle;
                node_index = child_left;
            } else {
                box_begin = box_middle;
                node_index = child_right;
            }
        }

        nearest_box = mapping_morton_[box_begin];
        nearest_point = get_point_hint_from_box( tree_[box_begin], nearest_box );
        distance = length( query - nearest_point );
    }

    /****************************************************************************/

    template< index_t DIMENSION >
    BoxAABBTree< DIMENSION >::BoxAABBTree(
        const std::vector< Box< DIMENSION > >& bboxes )
    {
        this->initialize_tree( bboxes );
    }

    template< index_t DIMENSION >
    vecn< DIMENSION > BoxAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box,
        index_t element_id ) const
    {
        ringmesh_unused( element_id );
        return box.center();
    }

    /****************************************************************************/

    template< index_t DIMENSION >
    LineAABBTree< DIMENSION >::LineAABBTree( const LineMesh< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_edges() );
        for( index_t i = 0; i < mesh.nb_edges(); i++ ) {
            for( index_t v = 0; v < 2; v++ ) {
                const vecn< DIMENSION >& point = mesh.vertex(
                    mesh.edge_vertex( i, v ) );
                bboxes[i].add_point( point );
            }
        }
        this->initialize_tree( bboxes );
    }

    template< index_t DIMENSION >
    index_t LineAABBTree< DIMENSION >::closest_edge(
        const vecn< DIMENSION >& query,
        vecn< DIMENSION >& nearest_point,
        double& distance ) const
    {
        DistanceToEdge action( mesh_ );
        return this->closest_element_box( query, nearest_point, distance, action );
    }

    template< index_t DIMENSION >
    void LineAABBTree< DIMENSION >::DistanceToEdge::operator()(
        const vecn< DIMENSION >& query,
        index_t cur_box,
        vecn< DIMENSION >& nearest_point,
        double& distance ) const
    {
        const vecn< DIMENSION >& v0 = mesh_.vertex(
            mesh_.edge_vertex( cur_box, 0 ) );
        const vecn< DIMENSION >& v1 = mesh_.vertex(
            mesh_.edge_vertex( cur_box, 1 ) );
        distance = point_segment_distance( query, v0, v1, nearest_point );
    }

    template< index_t DIMENSION >
    vecn< DIMENSION > LineAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box,
        index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex( mesh_.edge_vertex( element_id, 0 ) );
    }

    /****************************************************************************/

    template< index_t DIMENSION >
    SurfaceAABBTree< DIMENSION >::SurfaceAABBTree(
        const SurfaceMeshBase< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_polygons() );
        for( index_t i = 0; i < mesh.nb_polygons(); i++ ) {
            for( index_t v = 0; v < mesh.nb_polygon_vertices( i ); v++ ) {
                const vecn< DIMENSION >& point = mesh.vertex(
                    mesh.polygon_vertex( i, v ) );
                bboxes[i].add_point( point );
            }
        }
        this->initialize_tree( bboxes );
    }

    template< index_t DIMENSION >
    index_t SurfaceAABBTree< DIMENSION >::closest_triangle(
        const vecn< DIMENSION >& query,
        vecn< DIMENSION >& nearest_point,
        double& distance ) const
    {
        DistanceToTriangle action( mesh_ );
        return this->closest_element_box( query, nearest_point, distance, action );
    }

    template< index_t DIMENSION >
    void SurfaceAABBTree< DIMENSION >::DistanceToTriangle::operator()(
        const vecn< DIMENSION >& query,
        index_t cur_box,
        vecn< DIMENSION >& nearest_point,
        double& distance ) const
    {
        const vecn< DIMENSION >& v0 = mesh_.vertex(
            mesh_.polygon_vertex( cur_box, 0 ) );
        const vecn< DIMENSION >& v1 = mesh_.vertex(
            mesh_.polygon_vertex( cur_box, 1 ) );
        const vecn< DIMENSION >& v2 = mesh_.vertex(
            mesh_.polygon_vertex( cur_box, 2 ) );
        distance = point_triangle_distance( query, v0, v1, v2, nearest_point );
    }

    template< index_t DIMENSION >
    vecn< DIMENSION > SurfaceAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box,
        index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex( mesh_.polygon_vertex( element_id, 0 ) );
    }

    /****************************************************************************/

    template< index_t DIMENSION >
    VolumeAABBTree< DIMENSION >::VolumeAABBTree(
        const VolumeMesh< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_cells() );
        for( index_t i = 0; i < mesh.nb_cells(); i++ ) {
            for( index_t v = 0; v < mesh.nb_cell_vertices( i ); v++ ) {
                const vecn< DIMENSION >& point = mesh.vertex(
                    mesh.cell_vertex( i, v ) );
                bboxes[i].add_point( point );
            }
        }
        this->initialize_tree( bboxes );
    }

    template< index_t DIMENSION >
    vecn< DIMENSION > VolumeAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box,
        index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex( mesh_.cell_vertex( element_id, 0 ) );
    }

    template< index_t DIMENSION >
    index_t VolumeAABBTree< DIMENSION >::containing_cell(
        const vecn< DIMENSION >& query ) const
    {
        return containing_cell_recursive( query, AABBTree< DIMENSION >::ROOT_INDEX,
            0, this->nb_bboxes() );
    }

    template< index_t DIMENSION >
    index_t VolumeAABBTree< DIMENSION >::containing_cell_recursive(
        const vecn< DIMENSION >& query,
        index_t node_index,
        index_t box_begin,
        index_t box_end ) const
    {

        if( !this->tree_[node_index].contains( query ) ) {
            return NO_ID;
        }
        if( box_end == box_begin + 1 ) {
            index_t cell_id = this->mapping_morton_[box_begin];
            if( mesh_cell_contains_point( mesh_, cell_id, query ) ) {
                return cell_id;
            } else {
                return NO_ID;
            }
        }

        index_t box_middle, child_left, child_right;
        this->get_recursive_iterators( node_index, box_begin, box_end, box_middle,
            child_left, child_right );

        index_t result = containing_cell_recursive( query, child_left, box_begin,
            box_middle );
        if( result == NO_ID ) {
            result = containing_cell_recursive( query, child_right, box_middle,
                box_end );
        }
        return result;
    }

    template< index_t DIMENSION >
    double inner_point_box_distance(
        const vecn< DIMENSION >& p,
        const Box< DIMENSION >& B )
    {
        ringmesh_assert( B.contains( p ) );
        double result = std::abs( p[0] - B.min()[0] );
        result = std::min( result, std::abs( p[0] - B.max()[0] ) );
        for( index_t c = 1; c < 3; ++c ) {
            result = std::min( result, std::abs( p[c] - B.min()[c] ) );
            result = std::min( result, std::abs( p[c] - B.max()[c] ) );
        }
        return result;
    }

    template< index_t DIMENSION >
    double point_box_signed_distance(
        const vecn< DIMENSION >& p,
        const Box< DIMENSION >& B )
    {
        bool inside = true;
        vecn< DIMENSION > result;
        for( index_t c = 0; c < DIMENSION; c++ ) {
            if( p[c] < B.min()[c] ) {
                inside = false;
                result[c] = p[c] - B.min()[c];
            } else if( p[c] > B.max()[c] ) {
                inside = false;
                result[c] = p[c] - B.max()[c];
            }
        }
        if( inside ) {
            return -inner_point_box_distance( p, B );
        } else {
            return result.length();
        }
    }

//    template class AABBTree< 2 >;
    template class BoxAABBTree< 2 >;
    template class LineAABBTree< 2 >;
    template class SurfaceAABBTree< 2 >;

    template class AABBTree< 3 >;
    template class BoxAABBTree< 3 >;
    template class LineAABBTree< 3 >;
    template class SurfaceAABBTree< 3 >;
    template class VolumeAABBTree< 3 >;
}

