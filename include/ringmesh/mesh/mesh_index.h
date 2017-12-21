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

#pragma once

#include <ringmesh/mesh/common.h>

namespace RINGMesh
{
    struct EdgeLocalVertex;
    struct PolygonLocalEdge;
    struct CellLocalFacet;
} // namespace RINGMesh

namespace RINGMesh
{
    struct mesh_api ElementLocalVertex
    {
        ElementLocalVertex() = default;
        ElementLocalVertex( index_t element_id, index_t local_vertex_id )
            : element_id( element_id ), local_vertex_id( local_vertex_id )
        {
        }
        ElementLocalVertex( EdgeLocalVertex edge_local_vertex );
        ElementLocalVertex( PolygonLocalEdge polygon_local_edge );
        ElementLocalVertex( CellLocalFacet cell_local_facet );
        bool operator==( const ElementLocalVertex& rhs ) const;
        bool operator!=( const ElementLocalVertex& rhs ) const;
        index_t element_id{ NO_ID };
        index_t local_vertex_id{ NO_ID };
    };

    struct mesh_api EdgeLocalVertex
    {
        EdgeLocalVertex() = default;
        EdgeLocalVertex( index_t edge_id, index_t local_vertex_id )
            : edge_id( edge_id ), local_vertex_id( local_vertex_id )
        {
        }
        EdgeLocalVertex( ElementLocalVertex edge_local_vertex )
            : edge_id( std::move( edge_local_vertex.element_id ) ),
              local_vertex_id( std::move( edge_local_vertex.local_vertex_id ) )
        {
        }
        bool operator==( const EdgeLocalVertex& rhs ) const;
        bool operator!=( const EdgeLocalVertex& rhs ) const;
        index_t edge_id{ NO_ID };
        index_t local_vertex_id{ NO_ID };
    };

    struct mesh_api PolygonLocalEdge
    {
        PolygonLocalEdge() = default;
        PolygonLocalEdge( index_t polygon_id, index_t local_edge_id )
            : polygon_id( polygon_id ), local_edge_id( local_edge_id )
        {
        }
        PolygonLocalEdge( ElementLocalVertex polygon_local_vertex )
            : polygon_id( std::move( polygon_local_vertex.element_id ) ),
              local_edge_id( std::move( polygon_local_vertex.local_vertex_id ) )
        {
        }
        bool operator==( const PolygonLocalEdge& rhs ) const;
        bool operator!=( const PolygonLocalEdge& rhs ) const;
        index_t polygon_id{ NO_ID };
        index_t local_edge_id{ NO_ID };
    };

    struct mesh_api CellLocalFacet
    {
        CellLocalFacet() = default;
        CellLocalFacet( index_t cell_id, index_t local_facet_id )
            : cell_id( cell_id ), local_facet_id( local_facet_id )
        {
        }
        bool operator==( const CellLocalFacet& rhs ) const;
        bool operator!=( const CellLocalFacet& rhs ) const;
        index_t cell_id{ NO_ID };
        index_t local_facet_id{ NO_ID };
    };
} // namespace RINGMesh
