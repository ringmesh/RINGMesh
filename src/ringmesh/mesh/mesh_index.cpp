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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh_index.h>

namespace
{
    using RINGMesh::index_t;

    bool compare_equal(
        index_t lhs_index1,
        index_t lhs_index2,
        index_t rhs_index1,
        index_t rhs_index2 )
    {
        if( lhs_index1 == rhs_index1 ) {
            return lhs_index2 == rhs_index2;
        }
        return false;
    }
} // namespace

namespace RINGMesh
{
    ElementLocalVertex::ElementLocalVertex( EdgeLocalVertex edge_local_vertex )
        : element_id_( std::move( edge_local_vertex.edge_id_ ) ),
          local_vertex_id_( std::move( edge_local_vertex.local_vertex_id_ ) )
    {
    }

    ElementLocalVertex::ElementLocalVertex(
        PolygonLocalEdge polygon_local_edge )
        : element_id_( std::move( polygon_local_edge.polygon_id_ ) ),
          local_vertex_id_( std::move( polygon_local_edge.local_edge_id_ ) )
    {
    }

    ElementLocalVertex::ElementLocalVertex( CellLocalFacet cell_local_facet )
        : element_id_( std::move( cell_local_facet.cell_id_ ) ),
          local_vertex_id_( std::move( cell_local_facet.local_facet_id_ ) )
    {
    }

    bool ElementLocalVertex::operator==( const ElementLocalVertex& rhs ) const
    {
        return compare_equal( element_id_, local_vertex_id_, rhs.element_id_,
            rhs.local_vertex_id_ );
    }
    bool ElementLocalVertex::operator!=( const ElementLocalVertex& rhs ) const
    {
        return !operator==( rhs );
    }

    bool EdgeLocalVertex::operator==( const EdgeLocalVertex& rhs ) const
    {
        return compare_equal( edge_id_, local_vertex_id_, rhs.edge_id_,
            rhs.local_vertex_id_ );
    }
    bool EdgeLocalVertex::operator!=( const EdgeLocalVertex& rhs ) const
    {
        return !operator==( rhs );
    }

    bool PolygonLocalEdge::operator==( const PolygonLocalEdge& rhs ) const
    {
        return compare_equal( polygon_id_, local_edge_id_, rhs.polygon_id_,
            rhs.local_edge_id_ );
    }
    bool PolygonLocalEdge::operator!=( const PolygonLocalEdge& rhs ) const
    {
        return !operator==( rhs );
    }

    bool CellLocalFacet::operator==( const CellLocalFacet& rhs ) const
    {
        return compare_equal( cell_id_, local_facet_id_, rhs.cell_id_,
            rhs.local_facet_id_ );
    }
    bool CellLocalFacet::operator!=( const CellLocalFacet& rhs ) const
    {
        return !operator==( rhs );
    }
} // namespace RINGMesh
