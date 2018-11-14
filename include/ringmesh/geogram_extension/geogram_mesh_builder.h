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

#include <geogram/basic/command_line.h>

#include <ringmesh/geogram_extension/geogram_mesh.h>

#include <ringmesh/mesh/mesh_builder.h>

namespace RINGMesh
{
#define COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( Class )                    \
    \
public:                                                                        \
    void do_copy( const MeshBase< DIMENSION >& rhs, bool copy_attributes )     \
        override                                                               \
    {                                                                          \
        const auto& geogrammesh =                                              \
            dynamic_cast< const Class< DIMENSION >& >( rhs );                  \
        mesh_.mesh_->copy(                                                     \
            *geogrammesh.mesh_, copy_attributes, GEO::MESH_ALL_ELEMENTS );     \
    }                                                                          \
    void load_mesh( const std::string& filename ) override                     \
    {                                                                          \
        GEO::MeshIOFlags ioflags;                                              \
        ioflags.set_attribute( GEO::MESH_ALL_ATTRIBUTES );                     \
        GEO::mesh_load( filename, *mesh_.mesh_, ioflags );                     \
    }                                                                          \
    void do_clear( bool keep_attributes, bool keep_memory ) override           \
    {                                                                          \
        mesh_.mesh_->clear( keep_attributes, keep_memory );                    \
    }                                                                          \
    void do_set_vertex( index_t v_id, const vecn< DIMENSION >& vertex )        \
        override                                                               \
    {                                                                          \
        mesh_.ref_vertex( v_id ) = vertex;                                     \
    }                                                                          \
    index_t do_create_vertex() override                                        \
    {                                                                          \
        return mesh_.mesh_->vertices.create_vertex();                          \
    }                                                                          \
    index_t do_create_vertices( index_t nb ) override                          \
    {                                                                          \
        return mesh_.mesh_->vertices.create_vertices( nb );                    \
    }                                                                          \
    void do_assign_vertices( const std::vector< double >& point_coordinates )  \
        override                                                               \
    {                                                                          \
        GEO::vector< double > point_coordinates_cp =                           \
            copy_std_vector_to_geo_vector( point_coordinates );                \
        mesh_.mesh_->vertices.assign_points(                                   \
            point_coordinates_cp, DIMENSION, false );                          \
    }                                                                          \
    void do_delete_vertices( const std::vector< bool >& to_delete ) override   \
    {                                                                          \
        GEO::vector< index_t > vertices_to_delete =                            \
            copy_std_vector_to_geo_vector< bool, index_t >( to_delete );       \
        mesh_.mesh_->vertices.delete_elements( vertices_to_delete, false );    \
    }                                                                          \
    void do_clear_vertices( bool keep_attributes, bool keep_memory ) override  \
    {                                                                          \
        mesh_.mesh_->vertices.clear( keep_attributes, keep_memory );           \
    }                                                                          \
    void do_permute_vertices( const std::vector< index_t >& permutation )      \
        override                                                               \
    {                                                                          \
        GEO::vector< index_t > geo_vector_permutation =                        \
            copy_std_vector_to_geo_vector( permutation );                      \
        mesh_.mesh_->vertices.permute_elements( geo_vector_permutation );      \
    }                                                                          \
    \
private:                                                                       \
    Class< DIMENSION >& mesh_

    template < index_t DIMENSION >
    class GeogramPointSetMeshBuilder : public PointSetMeshBuilder< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramPointSetMesh );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        explicit GeogramPointSetMeshBuilder( PointSetMesh< DIMENSION >& mesh )
            : PointSetMeshBuilder< DIMENSION >( mesh ),
              mesh_( dynamic_cast< GeogramPointSetMesh< DIMENSION >& >( mesh ) )
        {
        }
    };

    ALIAS_2D_AND_3D( GeogramPointSetMeshBuilder );

    template < index_t DIMENSION >
    class GeogramLineMeshBuilder : public LineMeshBuilder< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramLineMesh );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        explicit GeogramLineMeshBuilder( LineMesh< DIMENSION >& mesh )
            : LineMeshBuilder< DIMENSION >( mesh ),
              mesh_( dynamic_cast< GeogramLineMesh< DIMENSION >& >( mesh ) )
        {
        }

        void do_create_edge( index_t v1_id, index_t v2_id ) override
        {
            mesh_.mesh_->edges.create_edge( v1_id, v2_id );
        }

        index_t do_create_edges( index_t nb_edges ) override
        {
            return mesh_.mesh_->edges.create_edges( nb_edges );
        }

        void do_set_edge_vertex( const EdgeLocalVertex& edge_local_vertex,
            index_t vertex_id ) override
        {
            mesh_.mesh_->edges.set_vertex( edge_local_vertex.edge_id,
                edge_local_vertex.local_vertex_id, vertex_id );
        }

        void do_delete_edges( const std::vector< bool >& to_delete ) override
        {
            GEO::vector< index_t > edges_to_delete =
                copy_std_vector_to_geo_vector< bool, index_t >( to_delete );
            mesh_.mesh_->edges.delete_elements( edges_to_delete, false );
        }

        void do_clear_edges( bool keep_attributes, bool keep_memory ) override
        {
            mesh_.mesh_->edges.clear( keep_attributes, keep_memory );
        }

        void do_permute_edges(
            const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_.mesh_->edges.permute_elements( geo_vector_permutation );
        }
    };

    ALIAS_2D_AND_3D( GeogramLineMeshBuilder );

    template < index_t DIMENSION >
    class GeogramSurfaceMeshBuilder : public SurfaceMeshBuilder< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramSurfaceMesh );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        explicit GeogramSurfaceMeshBuilder( SurfaceMesh< DIMENSION >& mesh )
            : SurfaceMeshBuilder< DIMENSION >( mesh ),
              mesh_( dynamic_cast< GeogramSurfaceMesh< DIMENSION >& >( mesh ) )
        {
        }

		void triangulate_with_geogram_cvt(
		                 const SurfaceMeshBase< DIMENSION >& surface_in )
		{
		   Logger::instance()->set_minimal( true );
		   const auto& geogram_surface_in =
			   dynamic_cast< const RINGMesh::GeogramSurfaceMesh< DIMENSION >& >(surface_in);
		   GEO::CentroidalVoronoiTesselation CVT(geogram_surface_in.mesh_.get(), DIMENSION, GEO::CmdLine::get_arg( "algo:delaunay" ) );
		   CVT.set_points(
		                   mesh_.nb_vertices(), mesh_.mesh_->vertices.point_ptr( 0 ) );
		   CVT.compute_surface( mesh_.mesh_.get(), false );
		   Logger::instance()->set_minimal( false );
		   clear_vertex_linked_objects();
		}

        index_t do_create_polygon(
            const std::vector< index_t >& vertices ) override
        {
            GEO::vector< index_t > polygon_vertices =
                copy_std_vector_to_geo_vector( vertices );
            return mesh_.mesh_->facets.create_polygon( polygon_vertices );
        }

        index_t do_create_triangles( index_t nb_triangles ) override
        {
            return mesh_.mesh_->facets.create_triangles( nb_triangles );
        }

        index_t do_create_quads( index_t nb_quads ) override
        {
            return mesh_.mesh_->facets.create_quads( nb_quads );
        }

        void do_set_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex,
            index_t vertex_id ) override
        {
            mesh_.mesh_->facets.set_vertex( polygon_local_vertex.element_id,
                polygon_local_vertex.local_vertex_id, vertex_id );
        }

        void do_set_polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge,
            index_t specifies ) override
        {
            mesh_.mesh_->facets.set_adjacent( polygon_local_edge.polygon_id,
                polygon_local_edge.local_edge_id, specifies );
        }

        void do_clear_polygons(
            bool keep_attributes, bool keep_memory ) override
        {
            mesh_.mesh_->facets.clear( keep_attributes, keep_memory );
        }

        void do_permute_polygons(
            const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_.mesh_->facets.permute_elements( geo_vector_permutation );
        }

        void do_delete_polygons( const std::vector< bool >& to_delete ) override
        {
            GEO::vector< index_t > polygons_to_delete =
                copy_std_vector_to_geo_vector< bool, index_t >( to_delete );
            mesh_.mesh_->facets.delete_elements( polygons_to_delete, false );
        }
    };

    ALIAS_2D_AND_3D( GeogramSurfaceMeshBuilder );

    template < index_t DIMENSION >
    class GeogramVolumeMeshBuilder : public VolumeMeshBuilder< DIMENSION >
    {
        COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramVolumeMesh );
        ringmesh_template_assert_3d( DIMENSION );

    public:
        explicit GeogramVolumeMeshBuilder( VolumeMesh< DIMENSION >& mesh )
            : VolumeMeshBuilder< DIMENSION >( mesh ),
              mesh_( dynamic_cast< GeogramVolumeMesh< DIMENSION >& >( mesh ) )
        {
        }

        index_t do_create_cells( index_t nb_cells, CellType type ) override
        {
            return mesh_.mesh_->cells.create_cells(
                nb_cells, static_cast< GEO::MeshCellType >( type ) );
        }

        void do_assign_cell_tet_mesh(
            const std::vector< index_t >& tets ) override
        {
            GEO::vector< index_t > copy = copy_std_vector_to_geo_vector( tets );
            mesh_.mesh_->cells.assign_tet_mesh( copy, false );
        }

        void do_set_cell_vertex( const ElementLocalVertex& cell_local_vertex,
            index_t vertex_id ) override
        {
            mesh_.mesh_->cells.set_vertex( cell_local_vertex.element_id,
                cell_local_vertex.local_vertex_id, vertex_id );
        }

        void do_set_cell_corner_vertex_index(
            index_t corner_index, index_t vertex_index ) override
        {
            mesh_.mesh_->cell_corners.set_vertex( corner_index, vertex_index );
        }

        void do_set_cell_adjacent( const CellLocalFacet& cell_local_facet,
            index_t cell_adjacent ) override
        {
            mesh_.mesh_->cells.set_adjacent( cell_local_facet.cell_id,
                cell_local_facet.local_facet_id, cell_adjacent );
        }

        void connect_cells() override
        {
            mesh_.mesh_->cells.connect();
        }

        void do_clear_cells( bool keep_attributes, bool keep_memory ) override
        {
            mesh_.mesh_->cells.clear( keep_attributes, keep_memory );
        }

        void do_permute_cells(
            const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_.mesh_->cells.permute_elements( geo_vector_permutation );
        }

        void do_delete_cells( const std::vector< bool >& to_delete ) override
        {
            GEO::vector< index_t > geo_to_delete =
                copy_std_vector_to_geo_vector< bool, index_t >( to_delete );
            mesh_.mesh_->cells.delete_elements( geo_to_delete, false );
        }
    };

    using GeogramVolumeMeshBuilder3D = GeogramVolumeMeshBuilder< 3 >;

} // namespace RINGMesh
