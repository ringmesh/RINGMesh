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

#pragma once

#include <ringmesh/geogram_extension/common.h>

#include <memory>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeogramPointSetMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeogramLineMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeogramSurfaceMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeogramVolumeMeshBuilder );
} // namespace RINGMesh

namespace RINGMesh
{
#define COMMON_GEOGRAM_MESH_IMPLEMENTATION( Class )                            \
    friend class Class##Builder< DIMENSION >;                                  \
    \
public:                                                                        \
    Class() : mesh_( new GEO::Mesh( DIMENSION, false ) ) {}                    \
    static MeshType type_name_static()                                         \
    {                                                                          \
        return #Class;                                                         \
    }                                                                          \
    void save_mesh( const std::string& filename ) const override               \
    {                                                                          \
        GEO::mesh_save( *mesh_, filename, GEO::MeshIOFlags() );                \
    }                                                                          \
    GEO::Mesh& geogram_mesh()                                                  \
    {                                                                          \
        return *mesh_;                                                         \
    }                                                                          \
    const GEO::Mesh& geogram_mesh() const                                      \
    {                                                                          \
        return *mesh_;                                                         \
    }                                                                          \
    GEO::AttributesManager& vertex_attribute_manager() const override          \
    {                                                                          \
        return mesh_->vertices.attributes();                                   \
    }                                                                          \
    MeshType type_name() const override                                        \
    {                                                                          \
        return type_name_static();                                             \
    }                                                                          \
    static std::string default_extension_static()                              \
    {                                                                          \
        return "geogram";                                                      \
    }                                                                          \
    std::string default_extension() const override                             \
    {                                                                          \
        return default_extension_static();                                     \
    }                                                                          \
    const vecn< DIMENSION >& vertex( index_t v_id ) const override             \
    {                                                                          \
        ringmesh_assert( v_id < nb_vertices() );                               \
        double* vertex_ptr = mesh_->vertices.point_ptr( v_id );                \
        return *(vecn< DIMENSION >*) ( vertex_ptr );                           \
    }                                                                          \
    index_t nb_vertices() const override                                       \
    {                                                                          \
        return mesh_->vertices.nb();                                           \
    }                                                                          \
    \
private:                                                                       \
    vecn< DIMENSION >& ref_vertex( index_t v_id )                              \
    {                                                                          \
        ringmesh_assert( v_id < nb_vertices() );                               \
        double* vertex_ptr = mesh_->vertices.point_ptr( v_id );                \
        return *(vecn< DIMENSION >*) ( vertex_ptr );                           \
    }                                                                          \
    \
protected:                                                                     \
    std::unique_ptr< GEO::Mesh > mesh_

    template < index_t DIMENSION >
    class GeogramPointSetMesh : public PointSetMesh< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramPointSetMesh );
    };

    ALIAS_2D_AND_3D( GeogramPointSetMesh );

    template < index_t DIMENSION >
    class GeogramLineMesh : public LineMesh< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramLineMesh );

    public:
        index_t edge_vertex(
            const ElementLocalVertex& edge_local_vertex ) const override
        {
            return mesh_->edges.vertex( edge_local_vertex.element_id,
                edge_local_vertex.local_vertex_id );
        }

        index_t nb_edges() const override
        {
            return mesh_->edges.nb();
        }

        GEO::AttributesManager& edge_attribute_manager() const override
        {
            return mesh_->edges.attributes();
        }
    };

    ALIAS_2D_AND_3D( GeogramLineMesh );

    template < index_t DIMENSION >
    class GeogramSurfaceMesh : public SurfaceMesh< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramSurfaceMesh );

    public:
        index_t polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const override
        {
            return mesh_->facets.vertex( polygon_local_vertex.element_id,
                polygon_local_vertex.local_vertex_id );
        }

        index_t nb_polygons() const override
        {
            return mesh_->facets.nb();
        }

        index_t nb_polygon_vertices( index_t polygon_id ) const override
        {
            return mesh_->facets.nb_vertices( polygon_id );
        }

        index_t polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge ) const override
        {
            return mesh_->facets.adjacent( polygon_local_edge.polygon_id,
                polygon_local_edge.local_edge_id );
        }

        GEO::AttributesManager& polygon_attribute_manager() const override
        {
            return mesh_->facets.attributes();
        }

        bool polygons_are_simplices() const override
        {
            return mesh_->facets.are_simplices();
        }
    };

    ALIAS_2D_AND_3D( GeogramSurfaceMesh );

    template < index_t DIMENSION >
    class GeogramVolumeMesh : public VolumeMesh< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramVolumeMesh );

    public:
        index_t cell_vertex(
            const ElementLocalVertex& cell_local_vertex ) const override
        {
            return mesh_->cells.vertex( cell_local_vertex.element_id,
                cell_local_vertex.local_vertex_id );
        }

        index_t cell_edge_vertex(
            index_t cell_id, index_t edge_id, index_t vertex_id ) const override
        {
            return mesh_->cells.edge_vertex( cell_id, edge_id, vertex_id );
        }

        index_t cell_facet_vertex( const CellLocalFacet& cell_local_facet,
            index_t vertex_id ) const override
        {
            return mesh_->cells.facet_vertex( cell_local_facet.cell_id,
                cell_local_facet.local_facet_id, vertex_id );
        }

        index_t cell_facet(
            const CellLocalFacet& cell_local_facet ) const override
        {
            return mesh_->cells.facet(
                cell_local_facet.cell_id, cell_local_facet.local_facet_id );
        }

        index_t nb_cell_facets( index_t cell_id ) const override
        {
            return mesh_->cells.nb_facets( cell_id );
        }

        index_t nb_cell_facets() const override
        {
            return mesh_->cell_facets.nb();
        }

        index_t nb_cell_edges( index_t cell_id ) const override
        {
            return mesh_->cells.nb_edges( cell_id );
        }

        index_t nb_cell_facet_vertices(
            const CellLocalFacet& cell_local_facet ) const override
        {
            return mesh_->cells.facet_nb_vertices(
                cell_local_facet.cell_id, cell_local_facet.local_facet_id );
        }

        index_t nb_cell_vertices( index_t cell_id ) const override
        {
            return mesh_->cells.nb_vertices( cell_id );
        }

        index_t nb_cells() const override
        {
            return mesh_->cells.nb();
        }

        index_t cell_begin( index_t cell_id ) const override
        {
            return mesh_->cells.corners_begin( cell_id );
        }

        index_t cell_end( index_t cell_id ) const override
        {
            return mesh_->cells.corners_end( cell_id );
        }

        index_t cell_adjacent(
            const CellLocalFacet& cell_local_facet ) const override
        {
            return mesh_->cells.adjacent(
                cell_local_facet.cell_id, cell_local_facet.local_facet_id );
        }

        GEO::AttributesManager& cell_attribute_manager() const override
        {
            return mesh_->cells.attributes();
        }

        GEO::AttributesManager& cell_facet_attribute_manager() const override
        {
            return mesh_->cell_facets.attributes();
        }

        CellType cell_type( index_t cell_id ) const override
        {
            return static_cast< CellType >( mesh_->cells.type( cell_id ) );
        }

        bool cells_are_simplicies() const override
        {
            return mesh_->cells.are_simplices();
        }

        double cell_volume( index_t cell_id ) const override
        {
            return RINGMesh::mesh_cell_volume( *mesh_, cell_id );
        }
    };

    using GeogramVolumeMesh3D = GeogramVolumeMesh< 3 >;
} // namespace RINGMesh
