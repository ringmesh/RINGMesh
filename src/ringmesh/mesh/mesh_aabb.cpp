/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/mesh/mesh_aabb.h>

#include <algorithm>
#include <numeric>

#include <ringmesh/basic/geometry.h>

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_index.h>

/// Copied and adapted from Geogram

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    bool mesh_cell_contains_point( const VolumeMesh< DIMENSION >& M,
        index_t cell,
        const vecn< DIMENSION >& p )
    {
        if( M.cell_type( cell ) != CellType::TETRAHEDRON )
        {
            throw RINGMeshException( "AABB", "VolumeAABBTre containing_cell "
                                             "request only handles "
                                             "tetrahedra." );
        }
        const auto& p0 = M.vertex( M.cell_vertex( { cell, 0 } ) );
        const auto& p1 = M.vertex( M.cell_vertex( { cell, 1 } ) );
        const auto& p2 = M.vertex( M.cell_vertex( { cell, 2 } ) );
        const auto& p3 = M.vertex( M.cell_vertex( { cell, 3 } ) );
        return Position::point_inside_tetra( p, { p0, p1, p2, p3 } );
    }
} // namespace

/****************************************************************************/

namespace RINGMesh
{
    template < index_t DIMENSION >
    LineAABBTree< DIMENSION >::LineAABBTree( const LineMesh< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_edges() );
        for( auto i : range( mesh.nb_edges() ) )
        {
            for( auto v : range( 2 ) )
            {
                bboxes[i].add_point( mesh.vertex(
                    mesh.edge_vertex( ElementLocalVertex( i, v ) ) ) );
            }
        }
        this->initialize_tree( bboxes );
    }

    template < index_t DIMENSION >
    std::tuple< index_t, vecn< DIMENSION >, double >
        LineAABBTree< DIMENSION >::closest_edge(
            const vecn< DIMENSION >& query ) const
    {
        DistanceToEdge action( mesh_ );
        return this->closest_element_box( query, action );
    }

    template < index_t DIMENSION >
    std::tuple< double, vecn< DIMENSION > >
        LineAABBTree< DIMENSION >::DistanceToEdge::operator()(
            const vecn< DIMENSION >& query, index_t cur_box ) const
    {
        const auto& v0 = mesh_.vertex( mesh_.edge_vertex( { cur_box, 0 } ) );
        const auto& v1 = mesh_.vertex( mesh_.edge_vertex( { cur_box, 1 } ) );
        return Distance::point_to_segment(
            query, Geometry::Segment< DIMENSION >{ v0, v1 } );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > LineAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box, index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex(
            mesh_.edge_vertex( ElementLocalVertex( element_id, 0 ) ) );
    }

    /****************************************************************************/

    template < index_t DIMENSION >
    SurfaceAABBTree< DIMENSION >::SurfaceAABBTree(
        const SurfaceMeshBase< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_polygons() );
        for( auto i : range( mesh.nb_polygons() ) )
        {
            for( auto v : range( mesh.nb_polygon_vertices( i ) ) )
            {
                bboxes[i].add_point( mesh.vertex(
                    mesh.polygon_vertex( ElementLocalVertex( i, v ) ) ) );
            }
        }
        this->initialize_tree( bboxes );
    }

    template < index_t DIMENSION >
    std::tuple< index_t, vecn< DIMENSION >, double >
        SurfaceAABBTree< DIMENSION >::closest_triangle(
            const vecn< DIMENSION >& query ) const
    {
        DistanceToTriangle action( mesh_ );
        return this->closest_element_box( query, action );
    }

    template < index_t DIMENSION >
    std::tuple< double, vecn< DIMENSION > >
        SurfaceAABBTree< DIMENSION >::DistanceToTriangle::operator()(
            const vecn< DIMENSION >& query, index_t cur_box ) const
    {
        const auto& v0 = mesh_.vertex( mesh_.polygon_vertex( { cur_box, 0 } ) );
        const auto& v1 = mesh_.vertex( mesh_.polygon_vertex( { cur_box, 1 } ) );
        const auto& v2 = mesh_.vertex( mesh_.polygon_vertex( { cur_box, 2 } ) );
        return Distance::point_to_triangle(
            query, Geometry::Triangle< DIMENSION >{ v0, v1, v2 } );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > SurfaceAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box, index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex(
            mesh_.polygon_vertex( ElementLocalVertex( element_id, 0 ) ) );
    }

    /****************************************************************************/

    template < index_t DIMENSION >
    VolumeAABBTree< DIMENSION >::VolumeAABBTree(
        const VolumeMesh< DIMENSION >& mesh )
        : mesh_( mesh )
    {
        std::vector< Box< DIMENSION > > bboxes;
        bboxes.resize( mesh.nb_cells() );
        for( auto i : range( mesh.nb_cells() ) )
        {
            for( auto v : range( mesh.nb_cell_vertices( i ) ) )
            {
                bboxes[i].add_point( mesh.vertex(
                    mesh.cell_vertex( ElementLocalVertex( i, v ) ) ) );
            }
        }
        this->initialize_tree( bboxes );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > VolumeAABBTree< DIMENSION >::get_point_hint_from_box(
        const Box< DIMENSION >& box, index_t element_id ) const
    {
        ringmesh_unused( box );
        return mesh_.vertex(
            mesh_.cell_vertex( ElementLocalVertex( element_id, 0 ) ) );
    }

    template < index_t DIMENSION >
    index_t VolumeAABBTree< DIMENSION >::containing_cell(
        const vecn< DIMENSION >& query ) const
    {
        return containing_cell_recursive(
            query, AABBTree< DIMENSION >::ROOT_INDEX, 0, this->nb_bboxes() );
    }

    template < index_t DIMENSION >
    index_t VolumeAABBTree< DIMENSION >::containing_cell_recursive(
        const vecn< DIMENSION >& query,
        index_t node_index,
        index_t box_begin,
        index_t box_end ) const
    {
        if( !this->node( node_index ).contains( query ) )
        {
            return NO_ID;
        }
        if( box_end == box_begin + 1 )
        {
            index_t cell_id = this->mapping_morton_[box_begin];
            if( mesh_cell_contains_point( mesh_, cell_id, query ) )
            {
                return cell_id;
            }
            return NO_ID;
        }

        index_t box_middle;
        index_t child_left;
        index_t child_right;
        this->get_recursive_iterators( node_index, box_begin, box_end,
            box_middle, child_left, child_right );

        index_t result = containing_cell_recursive(
            query, child_left, box_begin, box_middle );
        if( result == NO_ID )
        {
            result = containing_cell_recursive(
                query, child_right, box_middle, box_end );
        }
        return result;
    }

    template class mesh_api LineAABBTree< 2 >;
    template class mesh_api SurfaceAABBTree< 2 >;

    template class mesh_api LineAABBTree< 3 >;
    template class mesh_api SurfaceAABBTree< 3 >;
    template class mesh_api VolumeAABBTree< 3 >;
} // namespace RINGMesh
