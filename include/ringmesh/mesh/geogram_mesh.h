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

#include <ringmesh/basic/common.h>

#include <memory>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <ringmesh/mesh/mesh.h>

namespace RINGMesh {
    template< index_t DIMENSION > class GeogramPointSetMeshBuilder;
    template< index_t DIMENSION > class GeogramLineMeshBuilder;
    template< index_t DIMENSION > class GeogramSurfaceMeshBuilder;
    template< index_t DIMENSION > class GeogramVolumeMeshBuilder;
}

namespace RINGMesh {

#define COMMON_GEOGRAM_MESH_IMPLEMENTATION( Class )                                 \
    friend class Class ## Builder< DIMENSION >;                                     \
    public:                                                                         \
        Class()                                                                     \
            : mesh_( new GEO::Mesh( DIMENSION, false ) )                            \
        {                                                                           \
        }                                                                           \
        virtual ~Class() = default;                                                 \
        static MeshType type_name_static()                                          \
        {                                                                           \
            return #Class;                                                          \
        }                                                                           \
        virtual void save_mesh( const std::string& filename ) const override        \
        {                                                                           \
            GEO::mesh_save( *mesh_, filename, GEO::MeshIOFlags() );                 \
        }                                                                           \
        virtual const GEO::Mesh& gfx_mesh() const override                          \
        {                                                                           \
            return *mesh_;                                                          \
        }                                                                           \
        virtual void print_mesh_bounded_attributes() const override                 \
        {                                                                           \
            print_bounded_attributes( *mesh_ );                                     \
        }                                                                           \
        virtual GEO::AttributesManager& vertex_attribute_manager() const override   \
        {                                                                           \
            return mesh_->vertices.attributes();                                    \
        }                                                                           \
        virtual MeshType type_name() const override                                 \
        {                                                                           \
            return type_name_static();                                              \
        }                                                                           \
        static std::string default_extension_static()                               \
        {                                                                           \
            return "geogram";                                                       \
        }                                                                           \
        virtual std::string default_extension() const override                      \
        {                                                                           \
            return default_extension_static();                                      \
        }                                                                           \
        virtual const vecn< DIMENSION >& vertex( index_t v_id ) const override      \
        {                                                                           \
            ringmesh_assert( v_id < nb_vertices() );                                \
            double* vertex_ptr = mesh_->vertices.point_ptr( v_id );                 \
            return *( vecn< DIMENSION >* )( vertex_ptr );                           \
        }                                                                           \
        vecn< DIMENSION >& ref_vertex( index_t v_id )                               \
        {                                                                           \
            ringmesh_assert( v_id < nb_vertices() );                                \
            double* vertex_ptr = mesh_->vertices.point_ptr( v_id );                 \
            return *( vecn< DIMENSION >* )( vertex_ptr );                           \
        }                                                                           \
        virtual index_t nb_vertices() const override                                \
        {                                                                           \
            return mesh_->vertices.nb();                                            \
        }                                                                           \
    protected:                                                                      \
        std::unique_ptr< GEO::Mesh > mesh_

    template< index_t DIMENSION >
    class RINGMESH_API GeogramPointSetMesh: public PointSetMesh< DIMENSION > {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramPointSetMesh );
    };

    using GeogramPointSetMesh2D = GeogramPointSetMesh< 2 >;
    using GeogramPointSetMesh3D = GeogramPointSetMesh< 3 >;

    template< index_t DIMENSION >
    class RINGMESH_API GeogramLineMesh: public LineMesh< DIMENSION > {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramLineMesh );
    public:
        virtual index_t edge_vertex( index_t edge_id, index_t vertex_id ) const override
        {
            return mesh_->edges.vertex( edge_id, vertex_id );
        }

        virtual index_t nb_edges() const override
        {
            return mesh_->edges.nb();
        }

        virtual GEO::AttributesManager& edge_attribute_manager() const override
        {
            return mesh_->edges.attributes();
        }
    };

    using GeogramLineMesh2D = GeogramLineMesh< 2 >;
    using GeogramLineMesh3D = GeogramLineMesh< 3 >;

    template< index_t DIMENSION >
    class RINGMESH_API GeogramSurfaceMesh: public SurfaceMesh< DIMENSION > {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramSurfaceMesh );
    public:
        virtual index_t polygon_vertex( index_t polygon_id, index_t vertex_id ) const override
        {
            return mesh_->facets.vertex( polygon_id, vertex_id );
        }

        virtual index_t nb_polygons() const override
        {
            return mesh_->facets.nb();
        }

        virtual index_t nb_polygon_vertices( index_t polygon_id ) const override
        {
            return mesh_->facets.nb_vertices( polygon_id );
        }

        virtual index_t polygon_adjacent( index_t polygon_id, index_t edge_id ) const override
        {
            return mesh_->facets.adjacent( polygon_id, edge_id );
        }
        virtual GEO::AttributesManager& polygon_attribute_manager() const override
        {
            return mesh_->facets.attributes();
        }

        virtual bool polygons_are_simplicies() const override
        {
            return mesh_->facets.are_simplices();
        }
    };

    using GeogramSurfaceMesh2D = GeogramSurfaceMesh< 2 >;
    using GeogramSurfaceMesh3D = GeogramSurfaceMesh< 3 >;

    template< index_t DIMENSION >
    class RINGMESH_API GeogramVolumeMesh: public VolumeMesh< DIMENSION > {
        COMMON_GEOGRAM_MESH_IMPLEMENTATION( GeogramVolumeMesh );
    public:
        virtual index_t cell_vertex( index_t cell_id, index_t vertex_id ) const override
        {
            return mesh_->cells.vertex( cell_id, vertex_id );
        }

        virtual index_t cell_edge_vertex(
            index_t cell_id,
            index_t edge_id,
            index_t vertex_id ) const override
        {
            return mesh_->cells.edge_vertex( cell_id, edge_id, vertex_id );
        }

        virtual index_t cell_facet_vertex(
            index_t cell_id,
            index_t facet_id,
            index_t vertex_id ) const override
        {
            return mesh_->cells.facet_vertex( cell_id, facet_id, vertex_id );
        }

        virtual index_t cell_facet( index_t cell_id, index_t facet_id ) const override
        {
            return mesh_->cells.facet( cell_id, facet_id );
        }

        virtual index_t nb_cell_facets( index_t cell_id ) const override
        {
            return mesh_->cells.nb_facets( cell_id );
        }

        virtual index_t nb_cell_facets() const override
        {
            return mesh_->cell_facets.nb();
        }

        virtual index_t nb_cell_edges( index_t cell_id ) const override
        {
            return mesh_->cells.nb_edges( cell_id );
        }

        virtual index_t nb_cell_facet_vertices(
            index_t cell_id,
            index_t facet_id ) const override
        {
            return mesh_->cells.facet_nb_vertices( cell_id, facet_id );
        }

        virtual index_t nb_cell_vertices( index_t cell_id ) const override
        {
            return mesh_->cells.nb_vertices( cell_id );
        }

        virtual index_t nb_cells() const override
        {
            return mesh_->cells.nb();
        }

        virtual index_t cell_begin( index_t cell_id ) const override
        {
            return mesh_->cells.corners_begin( cell_id );
        }
        virtual index_t cell_end( index_t cell_id ) const override
        {
            return mesh_->cells.corners_end( cell_id );
        }

        virtual index_t cell_adjacent( index_t cell_id, index_t facet_id ) const override
        {
            return mesh_->cells.adjacent( cell_id, facet_id );
        }
        virtual GEO::AttributesManager& cell_attribute_manager() const override
        {
            return mesh_->cells.attributes();
        }
        virtual GEO::AttributesManager& cell_facet_attribute_manager() const override
        {
            return mesh_->cell_facets.attributes();
        }

        virtual GEO::MeshCellType cell_type( index_t cell_id ) const override
        {
            return mesh_->cells.type( cell_id );
        }

        virtual bool cells_are_simplicies() const override
        {
            return mesh_->cells.are_simplices();
        }

        virtual double cell_volume( index_t cell_id ) const override
        {
            return RINGMesh::mesh_cell_volume( *mesh_, cell_id );
        }
    };

    using GeogramVolumeMesh2D = GeogramVolumeMesh< 2 >;
    using GeogramVolumeMesh3D = GeogramVolumeMesh< 3 >;

    void register_geogram_mesh();
}
