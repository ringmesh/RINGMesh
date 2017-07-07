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

#include <geogram/basic/command_line.h>

#include <geogram/mesh/mesh_preprocessing.h>

#include <geogram/voronoi/CVT.h>

#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>

namespace RINGMesh {

#define COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( Class )                             \
    public:                                                                             \
        Class ## Builder()                                                              \
            : mesh_( nullptr )                                                          \
        {                                                                               \
        }                                                                               \
        virtual ~Class ## Builder() = default ;                                         \
        void copy( const MeshBase< DIMENSION>& rhs, bool copy_attributes ) override     \
        {                                                                               \
            const Class< DIMENSION >& geogrammesh =                                     \
                dynamic_cast< const Class< DIMENSION >& >( rhs );                       \
            mesh_->mesh_->copy( *geogrammesh.mesh_, copy_attributes,                    \
                GEO::MESH_ALL_ELEMENTS );                                               \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void load_mesh( const std::string& filename ) override                          \
        {                                                                               \
            GEO::MeshIOFlags ioflags;                                                   \
            ioflags.set_attribute( GEO::MESH_ALL_ATTRIBUTES );                          \
            GEO::mesh_load( filename, *mesh_->mesh_, ioflags );                         \
        }                                                                               \
        void clear( bool keep_attributes, bool keep_memory ) override                   \
        {                                                                               \
            mesh_->mesh_->clear( keep_attributes, keep_memory );                        \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void mesh_repair( GEO::MeshRepairMode mode, double colocate_epsilon ) override  \
        {                                                                               \
            GEO::mesh_repair( *mesh_->mesh_, mode, colocate_epsilon );                  \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void set_vertex( index_t v_id, const vecn< DIMENSION >& vertex ) override       \
        {                                                                               \
            mesh_->ref_vertex( v_id ) = vertex;                                         \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        index_t create_vertex() override                                                \
        {                                                                               \
            clear_vertex_linked_objects();                                              \
            return mesh_->mesh_->vertices.create_vertex();                              \
        }                                                                               \
        index_t create_vertices( index_t nb ) override                                  \
        {                                                                               \
            clear_vertex_linked_objects();                                              \
            return mesh_->mesh_->vertices.create_vertices( nb );                        \
        }                                                                               \
        void assign_vertices(                                                           \
            const std::vector< double >& point_coordinates ) override                   \
        {                                                                               \
            GEO::vector< double > point_coordinates_cp =                                \
                copy_std_vector_to_geo_vector( point_coordinates );                     \
            mesh_->mesh_->vertices.assign_points( point_coordinates_cp, DIMENSION,      \
                false );                                                                \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void delete_vertices( const std::vector< bool >& to_delete ) override           \
        {                                                                               \
            GEO::vector< index_t > vertices_to_delete =                                 \
                copy_std_vector_to_geo_vector< bool, index_t >( to_delete );            \
            mesh_->mesh_->vertices.delete_elements( vertices_to_delete, false );        \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void clear_vertices( bool keep_attributes, bool keep_memory ) override          \
        {                                                                               \
            mesh_->mesh_->vertices.clear( keep_attributes, keep_memory );               \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void permute_vertices( const std::vector< index_t >& permutation ) override     \
        {                                                                               \
            GEO::vector< index_t > geo_vector_permutation =                             \
                copy_std_vector_to_geo_vector( permutation );                           \
            mesh_->mesh_->vertices.permute_elements( geo_vector_permutation );          \
            clear_vertex_linked_objects();                                              \
        }                                                                               \
        void set_geogram_mesh( Class< DIMENSION >& mesh )                               \
        {                                                                               \
            mesh_ = &mesh;                                                              \
        }                                                                               \
    protected:                                                                          \
        void delete_vertex_nn_search()                                                  \
        {                                                                               \
            mesh_->vertex_nn_search_.reset();                                           \
        }                                                                               \
    private:                                                                            \
        Class< DIMENSION >* mesh_

    template< index_t DIMENSION >
    class GeogramPointSetMeshBuilder: public PointSetMeshBuilder< DIMENSION > {
    COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramPointSetMesh );ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        void set_mesh( PointSetMesh< DIMENSION >& mesh ) override
        {
            set_geogram_mesh(
                dynamic_cast< GeogramPointSetMesh< DIMENSION > & >( mesh ) );
        }
        void clear_vertex_linked_objects() override
        {
            delete_vertex_nn_search();
        }
    };

    using GeogramPointSetMesh2DBuilder = GeogramPointSetMeshBuilder< 2 >;
    using GeogramPointSetMesh3DBuilder = GeogramPointSetMeshBuilder< 3 >;

    template< index_t DIMENSION >
    class GeogramLineMeshBuilder: public LineMeshBuilder< DIMENSION > {
    COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramLineMesh );ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:

        void set_mesh( LineMesh< DIMENSION >& mesh ) override
        {
            set_geogram_mesh(
                dynamic_cast< GeogramLineMesh< DIMENSION >& >( mesh ) );
        }

        void create_edge( index_t v1_id, index_t v2_id ) override
        {
            mesh_->mesh_->edges.create_edge( v1_id, v2_id );
            clear_edge_linked_objects();
        }

        index_t create_edges( index_t nb_edges ) override
        {
            return mesh_->mesh_->edges.create_edges( nb_edges );
        }

        void set_edge_vertex(
            index_t edge_id,
            index_t local_vertex_id,
            index_t vertex_id ) override
        {
            mesh_->mesh_->edges.set_vertex( edge_id, local_vertex_id, vertex_id );
            clear_edge_linked_objects();
        }

        void delete_edges(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) override
        {
            GEO::vector< index_t > edges_to_delete = copy_std_vector_to_geo_vector<
                bool, index_t >( to_delete );
            mesh_->mesh_->edges.delete_elements( edges_to_delete, false );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_edge_linked_objects();
        }

        void clear_edges( bool keep_attributes, bool keep_memory ) override
        {
            mesh_->mesh_->edges.clear( keep_attributes, keep_memory );
            clear_edge_linked_objects();
        }

        void remove_isolated_vertices() override
        {
            std::vector< bool > to_delete( mesh_->nb_vertices(), true );
            for( index_t e : range( mesh_->nb_edges() ) ) {
                for( index_t v : range( 2 ) ) {
                    index_t vertex_id = mesh_->edge_vertex( e, v );
                    to_delete[vertex_id] = false;
                }
            }
            delete_vertices( to_delete );
        }

        void permute_edges( const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_->mesh_->edges.permute_elements( geo_vector_permutation );
        }

        void clear_vertex_linked_objects() override
        {
            delete_vertex_nn_search();
            clear_edge_linked_objects();
        }

        void clear_edge_linked_objects() override
        {
            delete_edge_nn_search();
        }

    protected:
        /*!
         * @brief Deletes the NNSearch on edges
         */
        void delete_edge_nn_search()
        {
            mesh_->edge_nn_search_.reset();
        }
    };

    using GeogramLineMesh2DBuilder = GeogramLineMeshBuilder< 2 >;
    using GeogramLineMesh3DBuilder = GeogramLineMeshBuilder< 3 >;

    template< index_t DIMENSION >
    class GeogramSurfaceMeshBuilder: public SurfaceMeshBuilder< DIMENSION > {
    COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramSurfaceMesh );ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:

        void set_mesh( SurfaceMeshBase< DIMENSION >& mesh ) override
        {
            set_geogram_mesh(
                dynamic_cast< GeogramSurfaceMesh< DIMENSION >& >( mesh ) );
        }

        void remove_small_connected_components(
            double min_area,
            index_t min_polygons ) override
        {
            GEO::remove_small_connected_components( *mesh_->mesh_, min_area,
                min_polygons );
        }

        void triangulate( const SurfaceMeshBase< DIMENSION >& surface_in ) override
        {
            Logger::instance()->set_minimal( true );
            const GeogramSurfaceMesh< DIMENSION >& geogram_surf_in =
                dynamic_cast< const GeogramSurfaceMesh< DIMENSION >& >( surface_in );
            GEO::CentroidalVoronoiTesselation CVT( geogram_surf_in.mesh_.get(), 3,
                GEO::CmdLine::get_arg( "algo:delaunay" ) );
            CVT.set_points( mesh_->nb_vertices(),
                mesh_->mesh_->vertices.point_ptr( 0 ) );
            CVT.compute_surface( mesh_->mesh_.get(), false );
            Logger::instance()->set_minimal( false );
        }

        void create_polygons(
            const std::vector< index_t >& polygons,
            const std::vector< index_t >& polygon_ptr ) override
        {
            for( index_t p : range( polygon_ptr.size() - 1 ) ) {
                index_t start = polygon_ptr[p];
                index_t end = polygon_ptr[p + 1];
                GEO::vector< index_t > polygon_vertices =
                    copy_std_vector_to_geo_vector( polygons, start, end );
                mesh_->mesh_->facets.create_polygon( polygon_vertices );
            }
            clear_polygon_linked_objects();
        }

        index_t create_polygon( const std::vector< index_t >& vertices ) override
        {
            GEO::vector< index_t > polygon_vertices = copy_std_vector_to_geo_vector(
                vertices );
            index_t index = mesh_->mesh_->facets.create_polygon( polygon_vertices );
            clear_polygon_linked_objects();
            return index;
        }

        index_t create_triangles( index_t nb_triangles ) override
        {
            return mesh_->mesh_->facets.create_triangles( nb_triangles );
        }

        index_t create_quads( index_t nb_quads ) override
        {
            return mesh_->mesh_->facets.create_quads( nb_quads );
        }

        void set_polygon_vertex(
            index_t polygon_id,
            index_t local_vertex_id,
            index_t vertex_id ) override
        {
            mesh_->mesh_->facets.set_vertex( polygon_id, local_vertex_id,
                vertex_id );
            clear_polygon_linked_objects();
        }

        void set_polygon_adjacent(
            index_t polygon_id,
            index_t edge_id,
            index_t specifies ) override
        {
            mesh_->mesh_->facets.set_adjacent( polygon_id, edge_id, specifies );
        }

        void assign_triangle_mesh( const std::vector< index_t >& triangles ) override
        {
            GEO::vector< index_t > geo_triangles = copy_std_vector_to_geo_vector(
                triangles );
            mesh_->mesh_->facets.assign_triangle_mesh( geo_triangles, false );
            clear_polygon_linked_objects();
        }

        void clear_polygons( bool keep_attributes, bool keep_memory ) override
        {
            mesh_->mesh_->facets.clear( keep_attributes, keep_memory );
        }

        void connect_polygons() override
        {
            mesh_->mesh_->facets.connect();
        }

        void permute_polygons( const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_->mesh_->facets.permute_elements( geo_vector_permutation );
        }

        void delete_polygons(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) override
        {
            GEO::vector< index_t > polygons_to_delete =
                copy_std_vector_to_geo_vector< bool, index_t >( to_delete );
            mesh_->mesh_->facets.delete_elements( polygons_to_delete, false );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_polygon_linked_objects();
        }

        void remove_isolated_vertices() override
        {
            std::vector< bool > to_delete( mesh_->nb_vertices(), true );

            for( index_t p : range( mesh_->nb_polygons() ) ) {
                for( index_t v : range( mesh_->nb_polygon_vertices( p ) ) ) {
                    index_t vertex_id = mesh_->polygon_vertex( p, v );
                    to_delete[vertex_id] = false;
                }
            }

            delete_vertices( to_delete );
        }

        void clear_vertex_linked_objects() override
        {
            delete_vertex_nn_search();
            clear_polygon_linked_objects();
        }

        void clear_polygon_linked_objects() override
        {
            delete_polygon_aabb();
            delete_polygon_nn_search();
        }
    protected:
        /*!
         * @brief Deletes the NNSearch on polygons
         */
        void delete_polygon_nn_search()
        {
            mesh_->nn_search_.reset();
        }
        /*!
         * @brief Deletes the AABB on polygons
         */
        void delete_polygon_aabb()
        {
            mesh_->polygon_aabb_.reset();
        }
    };

    using GeogramSurfaceMesh2DBuilder = GeogramSurfaceMeshBuilder< 2 >;
    using GeogramSurfaceMesh3DBuilder = GeogramSurfaceMeshBuilder< 3 >;

    template< index_t DIMENSION >
    class GeogramVolumeMeshBuilder: public VolumeMeshBuilder< DIMENSION > {
    COMMON_GEOGRAM_MESH_BUILDER_IMPLEMENTATION( GeogramVolumeMesh );ringmesh_template_assert_3d( DIMENSION );
    public:

        void set_mesh( VolumeMesh< DIMENSION >& mesh ) override
        {
            set_geogram_mesh(
                dynamic_cast< GeogramVolumeMesh< DIMENSION >& >( mesh ) );
        }

        index_t create_cells( index_t nb_cells, CellType type ) override
        {
            return mesh_->mesh_->cells.create_cells( nb_cells,
                static_cast< GEO::MeshCellType >( type ) );
        }

        void assign_cell_tet_mesh( const std::vector< index_t >& tets ) override
        {
            GEO::vector< index_t > copy = copy_std_vector_to_geo_vector( tets );
            mesh_->mesh_->cells.assign_tet_mesh( copy, false );
            clear_cell_linked_objects();
        }

        void set_cell_vertex(
            index_t cell_id,
            index_t local_vertex_id,
            index_t vertex_id ) override
        {
            mesh_->mesh_->cells.set_vertex( cell_id, local_vertex_id, vertex_id );
            clear_cell_linked_objects();
        }

        void set_cell_corner_vertex_index(
            index_t corner_index,
            index_t vertex_index ) override
        {
            mesh_->mesh_->cell_corners.set_vertex( corner_index, vertex_index );
            clear_cell_linked_objects();
        }

        void set_cell_adjacent(
            index_t cell_index,
            index_t facet_index,
            index_t cell_adjacent ) override
        {
            mesh_->mesh_->cells.set_adjacent( cell_index, facet_index,
                cell_adjacent );
        }

        void connect_cells() override
        {
            mesh_->mesh_->cells.connect();
        }

        void clear_cells( bool keep_attributes, bool keep_memory ) override
        {
            mesh_->mesh_->cells.clear( keep_attributes, keep_memory );
        }

        void permute_cells( const std::vector< index_t >& permutation ) override
        {
            GEO::vector< index_t > geo_vector_permutation =
                copy_std_vector_to_geo_vector( permutation );
            mesh_->mesh_->cells.permute_elements( geo_vector_permutation );
        }

        void delete_cells(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) override
        {
            GEO::vector< index_t > geo_to_delete = copy_std_vector_to_geo_vector<
                bool, index_t >( to_delete );
            mesh_->mesh_->cells.delete_elements( geo_to_delete, false );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_cell_linked_objects();
        }

        void remove_isolated_vertices() override
        {
            std::vector< bool > to_delete( mesh_->nb_vertices(), true );
            for( index_t c : range( mesh_->nb_cells() ) ) {
                for( index_t v : range( mesh_->nb_cell_vertices( c ) ) ) {
                    index_t vertex_id = mesh_->cell_vertex( c, v );
                    to_delete[vertex_id] = false;
                }
            }
            delete_vertices( to_delete );
        }

        void clear_vertex_linked_objects() override
        {
            delete_vertex_nn_search();
            clear_cell_linked_objects();
        }

        void clear_cell_linked_objects() override
        {
            delete_cell_aabb();
            delete_cell_nn_search();
        }

    protected:
        /*!
         * @brief Deletes the NNSearch on cells
         */
        void delete_cell_nn_search()
        {
            mesh_->cell_nn_search_.reset();
            mesh_->cell_facet_nn_search_.reset();
        }
        /*!
         * @brief Deletes the AABB on cells
         */
        void delete_cell_aabb()
        {
            mesh_->cell_aabb_.reset();
        }
    };

    using GeogramVolumeMesh3DBuilder = GeogramVolumeMeshBuilder< 3 >;

}
