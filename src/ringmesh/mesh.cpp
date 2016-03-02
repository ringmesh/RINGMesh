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

/*! \author Francois Bonneau */

#include <ringmesh/mesh.h>
#include <geogram/mesh/mesh.h>

namespace RINGMesh {

    Mesh::Mesh()
    {
        mesh_ = new GEO::Mesh ;
    }

    Mesh::~Mesh()
    {
        delete mesh_ ;
    }

    const vec3& Mesh::vertex( index_t v_id ) const
    {
        return mesh_->vertices.point( v_id ) ;
    }

    index_t Mesh::nb_vertices() const
    {
        return mesh_->vertices.nb() ;
    }

    index_t Mesh::edge_vertex( index_t edge_id, index_t vertex_id ) const
    {
        return mesh_->edges.vertex( edge_id, vertex_id ) ;
    }

    index_t Mesh::nb_edges() const
    {
        return mesh_->edges.nb() ;
    }

    index_t Mesh::facet_vertex( index_t facet_id, index_t vertex_id ) const
    {
        return mesh_->facets.vertex( facet_id, vertex_id ) ;
    }

    index_t Mesh::facet_adjacent( index_t facet_id, index_t vertex_id ) const
    {
        return mesh_->facets.adjacent( facet_id, vertex_id ) ;
    }

    index_t Mesh::nb_facets() const
    {
        return mesh_->facets.nb() ;
    }

    index_t Mesh::cell_vertex( index_t cell_id, index_t vertex_id ) const
    {
        return mesh_->cells.vertex( cell_id, vertex_id ) ;
    }

    index_t Mesh::cell_facet_vertex(
        index_t cell_id,
        index_t facet_id,
        index_t vertex_id ) const
    {
        return mesh_->cells.facet_vertex( cell_id, facet_id, vertex_id ) ;
    }

    index_t Mesh::cell_adjacent( index_t cell_id, index_t facet_id ) const
    {
        return mesh_->cells.adjacent( cell_id, facet_id ) ;
    }

    index_t Mesh::nb_cell_facets( index_t cell_id ) const
    {
        return mesh_->cells.nb_facets( cell_id ) ;
    }

    index_t Mesh::nb_cell_facet_vertices(
        index_t cell_id,
        index_t facet_id ) const
    {
        return mesh_->cells.facet_nb_vertices( cell_id, facet_id ) ;
    }

    index_t Mesh::nb_cells() const
    {
        return mesh_->cells.nb() ;
    }

} // namespace
