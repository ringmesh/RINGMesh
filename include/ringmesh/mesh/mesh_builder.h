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
#include <numeric>

#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/basic/factory.h>

#include <ringmesh/mesh/mesh.h>

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModel;
}

namespace RINGMesh {
    template< index_t DIMENSION >
    class MeshBaseBuilder {
    ringmesh_disable_copy( MeshBaseBuilder );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~MeshBaseBuilder()
        {
        }
        /*!
         * \name general methods
         * @{
         */
        /*!
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @param[in] copy_attributes if true, all attributes are copied.
         * @return a modifiable reference to the point that corresponds to the vertex.
         */
        void copy( const MeshBase< DIMENSION >& rhs, bool copy_attributes )
        {
            do_copy( rhs, copy_attributes );
            clear_vertex_linked_objects();
        }

        virtual void load_mesh( const std::string& filename ) = 0;
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear( bool keep_attributes, bool keep_memory )
        {
            do_clear( keep_attributes, keep_memory );
            clear_vertex_linked_objects();
        }
        /*!
         * \brief Fixes some defaults in a mesh.
         * \param[in] mode a combination of #MeshRepairMode flags.
         *  Combine them with the 'bitwise or' (|) operator.
         * \param[in] colocate_epsilon tolerance used to colocate vertices
         *  (if #MESH_REPAIR_COLOCATE is set in mode).
         */
        void repair( GEO::MeshRepairMode mode, double colocate_epsilon )
        {
            do_repair( mode, colocate_epsilon );
            clear_vertex_linked_objects();
        }
        /*!@}
         * \name Vertex related methods
         * @{
         */
        /*!
         * @brief Sets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @param[in] vertex the vertex coordinates
         * @return reference to the point that corresponds to the vertex.
         */
        void set_vertex( index_t v_id, const vecn< DIMENSION >& vertex )
        {
            do_set_vertex( v_id, vertex );
            clear_vertex_linked_objects();
        }
        /*!
         * @brief Creates a new vertex.
         * @return the index of the created vertex
         */
        index_t create_vertex()
        {
            index_t index = do_create_vertex();
            clear_vertex_linked_objects();
            return index;
        }
        /*!
         * @brief Creates a new vertex.
         * @param[in] coords a pointer to @function dimension() coordinate.
         * @return the index of the created vertex
         */
        index_t create_vertex( const vecn< DIMENSION >& vertex )
        {
            index_t index = create_vertex();
            set_vertex( index, vertex );
            return index;
        }
        /*!
         * @brief Creates a contiguous chunk of vertices.
         * @param[in] nb number of sub-entities to create.
         * @return the index of the first created vertex
         */
        index_t create_vertices( index_t nb )
        {
            index_t index = do_create_vertices( nb );
            clear_vertex_linked_objects();
            return index;
        }

        /*!
         * @brief set vertex coordinates from a std::vector of coordinates
         * @param[in] point_coordinates a set of x, y (, z) coordinates
         */
        void assign_vertices( const std::vector< double >& point_coordinates )
        {
            do_assign_vertices( point_coordinates );
            clear_vertex_linked_objects();
        }

        /*!
         * @brief Deletes a set of vertices.
         * @param[in] to_delete     a vector of size @function nb(). If to_delete[e] is true,
         * then entity e will be destroyed, else it will be kept.
         */
        void delete_vertices( const std::vector< bool >& to_delete )
        {
            do_delete_vertices( to_delete );
            clear_vertex_linked_objects();
        }
        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_vertices( bool keep_attributes, bool keep_memory )
        {
            do_clear_vertices( keep_attributes, keep_memory );
            clear_vertex_linked_objects();
        }
        void permute_vertices( const std::vector< index_t >& permutation )
        {
            do_permute_vertices( permutation );
            clear_vertex_linked_objects();
        }
        /*!@}
         */

        static std::unique_ptr< MeshBaseBuilder< DIMENSION > > create_builder(
            MeshBase< DIMENSION >& mesh );
    protected:
        MeshBaseBuilder( MeshBase< DIMENSION >& mesh )
            : mesh_base_( mesh )
        {
        }

        void delete_vertex_nn_search()
        {
            mesh_base_.vertex_nn_search_.reset();
        }

        /*!
         * @brief Deletes the NNSearch on vertices
         */
        virtual void clear_vertex_linked_objects() = 0;

    private:
        /*!
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @param[in] copy_attributes if true, all attributes are copied.
         * @return a modifiable reference to the point that corresponds to the vertex.
         */
        virtual void do_copy(
            const MeshBase< DIMENSION >& rhs,
            bool copy_attributes ) = 0;
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear( bool keep_attributes, bool keep_memory ) = 0;
        /*!
         * \brief Fixes some defaults in a mesh.
         * \param[in] mode a combination of #MeshRepairMode flags.
         *  Combine them with the 'bitwise or' (|) operator.
         * \param[in] colocate_epsilon tolerance used to colocate vertices
         *  (if #MESH_REPAIR_COLOCATE is set in mode).
         */
        virtual void do_repair(
            GEO::MeshRepairMode mode,
            double colocate_epsilon ) = 0;
        /*!
         * @brief Sets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @param[in] vertex the vertex coordinates
         * @return reference to the point that corresponds to the vertex.
         */
        virtual void do_set_vertex(
            index_t v_id,
            const vecn< DIMENSION >& vertex ) = 0;
        /*!
         * @brief Creates a new vertex.
         * @return the index of the created vertex
         */
        virtual index_t do_create_vertex() = 0;
        /*!
         * @brief Creates a contiguous chunk of vertices.
         * @param[in] nb number of sub-entities to create.
         * @return the index of the first created vertex
         */
        virtual index_t do_create_vertices( index_t nb ) = 0;
        /*!
         * @brief set vertex coordinates from a std::vector of coordinates
         * @param[in] point_coordinates a set of x, y (, z) coordinates
         */
        virtual void do_assign_vertices(
            const std::vector< double >& point_coordinates ) = 0;
        /*!
         * @brief Deletes a set of vertices.
         * @param[in] to_delete     a vector of size @function nb(). If to_delete[e] is true,
         * then entity e will be destroyed, else it will be kept.
         */
        virtual void do_delete_vertices( const std::vector< bool >& to_delete ) = 0;
        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_vertices( bool keep_attributes, bool keep_memory ) = 0;
        virtual void do_permute_vertices(
            const std::vector< index_t >& permutation ) = 0;

    protected:
        MeshBase< DIMENSION >& mesh_base_;
    };

    CLASS_DIMENSION_ALIASES( MeshBaseBuilder );

    template< index_t DIMENSION >
    class PointSetMeshBuilder: public MeshBaseBuilder< DIMENSION > {
    ringmesh_disable_copy( PointSetMeshBuilder );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        virtual ~PointSetMeshBuilder() = default;

        static std::unique_ptr< PointSetMeshBuilder< DIMENSION > > create_builder(
            PointSetMesh< DIMENSION >& mesh );

        virtual void remove_isolated_vertices()
        {
            // All vertices are isolated in a Mesh0D
        }

    protected:
        PointSetMeshBuilder( PointSetMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), pointset_mesh_( mesh )
        {
        }

    private:
        void clear_vertex_linked_objects() final
        {
            this->delete_vertex_nn_search();
        }

    protected:
        PointSetMesh< DIMENSION >& pointset_mesh_;
    };

    CLASS_DIMENSION_ALIASES( PointSetMeshBuilder );

    template< index_t DIMENSION >
    using PointSetMeshBuilderFactory = Factory< MeshType, PointSetMeshBuilder< DIMENSION >, PointSetMesh< DIMENSION >& >;

    CLASS_DIMENSION_ALIASES( PointSetMeshBuilderFactory );

    template< index_t DIMENSION >
    class LineMeshBuilder: public MeshBaseBuilder< DIMENSION > {
    ringmesh_disable_copy( LineMeshBuilder );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        virtual ~LineMeshBuilder() = default;

        static std::unique_ptr< LineMeshBuilder > create_builder(
            LineMesh< DIMENSION >& mesh );

        /*!
         * @brief Create a new edge.
         * @param[in] v1_id index of the starting vertex.
         * @param[in] v2_id index of the ending vertex.
         */
        void create_edge( index_t v1_id, index_t v2_id )
        {
            do_create_edge( v1_id, v2_id );
            clear_edge_linked_objects();
        }
        /*!
         * \brief Creates a contiguous chunk of edges
         * \param[in] nb_edges number of edges to create
         * \return the index of the first edge
         */
        index_t create_edges( index_t nb_edges )
        {
            index_t index = do_create_edges( nb_edges );
            clear_edge_linked_objects();
            return index;
        }
        /*!
         * @brief Sets a vertex of a edge by local vertex index.
         * @param[in] edge_id index of the edge, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the edge.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of edge
         * \param edge_id. Index between 0 and @function nb() - 1.
         */
        void set_edge_vertex(
            index_t edge_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            do_set_edge_vertex( edge_id, local_vertex_id, vertex_id );
            clear_edge_linked_objects();
        }
        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices
         * that are no longer incident to any entity are deleted.
         */
        void delete_edges(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices )
        {
            do_delete_edges( to_delete );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_edge_linked_objects();
        }
        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_edges( bool keep_attributes, bool keep_memory )
        {
            do_clear_edges( keep_attributes, keep_memory );
            clear_edge_linked_objects();
        }
        void permute_edges( const std::vector< index_t >& permutation )
        {
            do_permute_edges( permutation );
            clear_edge_linked_objects();
        }

        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        void remove_isolated_vertices()
        {
            std::vector< bool > to_delete( line_mesh_.nb_vertices(), true );
            for( index_t e : range( line_mesh_.nb_edges() ) ) {
                for( index_t v : range( 2 ) ) {
                    index_t vertex_id = line_mesh_.edge_vertex(
                        ElementLocalVertex( e, v ) );
                    to_delete[vertex_id] = false;
                }
            }
            this->delete_vertices( to_delete );

        }
    protected:
        LineMeshBuilder( LineMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), line_mesh_( mesh )
        {
        }

    private:
        /*!
         * @brief Deletes the NNSearch on edges
         */
        void delete_edge_nn_search()
        {
            line_mesh_.edge_nn_search_.reset();
        }

        void clear_vertex_linked_objects() override
        {
            this->delete_vertex_nn_search();
            clear_edge_linked_objects();
        }

        void clear_edge_linked_objects()
        {
            delete_edge_nn_search();
        }

        /*!
         * @brief Create a new edge.
         * @param[in] v1_id index of the starting vertex.
         * @param[in] v2_id index of the ending vertex.
         */
        virtual void do_create_edge( index_t v1_id, index_t v2_id ) = 0;
        /*!
         * \brief Creates a contiguous chunk of edges
         * \param[in] nb_edges number of edges to create
         * \return the index of the first edge
         */
        virtual index_t do_create_edges( index_t nb_edges ) = 0;
        /*!
         * @brief Sets a vertex of a edge by local vertex index.
         * @param[in] edge_id index of the edge, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the edge.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of edge
         * \param edge_id. Index between 0 and @function nb() - 1.
         */
        virtual void do_set_edge_vertex(
            index_t edge_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0;
        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         */
        virtual void do_delete_edges( const std::vector< bool >& to_delete ) = 0;
        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_edges( bool keep_attributes, bool keep_memory ) = 0;
        virtual void do_permute_edges(
            const std::vector< index_t >& permutation ) = 0;

    protected:
        LineMesh< DIMENSION >& line_mesh_;
    };

    CLASS_DIMENSION_ALIASES( LineMeshBuilder );

    template< index_t DIMENSION >
    using LineMeshBuilderFactory = Factory< MeshType, LineMeshBuilder< DIMENSION >, LineMesh< DIMENSION >& >;

    CLASS_DIMENSION_ALIASES( LineMeshBuilderFactory );

    template< index_t DIMENSION >
    class SurfaceMeshBuilder: public MeshBaseBuilder< DIMENSION > {
    ringmesh_disable_copy( SurfaceMeshBuilder );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        virtual ~SurfaceMeshBuilder() = default;

        static std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > create_builder(
            SurfaceMesh< DIMENSION >& mesh );

        /*!@}
         * \name Polygon related methods
         * @{
         */
        /*!
         * brief create polygons
         * @param[in] polygons is the vector of vertex index for each polygon
         * @param[in] polygon_ptr is the vector addressing the first polygon vertex for each polygon.
         */
        void create_polygons(
            const std::vector< index_t >& polygons,
            const std::vector< index_t >& polygon_ptr )
        {
            do_create_polygons( polygons, polygon_ptr );
            clear_polygon_linked_objects();
        }
        /*!
         * \brief Creates a polygon
         * \param[in] vertices a const reference to a vector that
         *  contains the vertices
         * \return the index of the created polygon
         */
        index_t create_polygon( const std::vector< index_t >& vertices )
        {
            index_t index = do_create_polygon( vertices );
            clear_polygon_linked_objects();
            return index;
        }
        /*!
         * \brief Creates a contiguous chunk of triangles
         * \param[in] nb_triangles number of triangles to create
         * \return the index of the first triangle
         */
        index_t create_triangles( index_t nb_triangles )
        {
            index_t index = do_create_triangles( nb_triangles );
            clear_polygon_linked_objects();
            return index;
        }
        /*!
         * \brief Creates a contiguous chunk of quads
         * \param[in] nb_quads number of quads to create
         * \return the index of the first quad
         */
        index_t create_quads( index_t nb_quads )
        {
            index_t index = do_create_quads( nb_quads );
            clear_polygon_linked_objects();
            return index;
        }
        /*!
         * @brief Sets a vertex of a polygon by local vertex index.
         * @param[in] polygon_id index of the polygon, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the polygon.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of the
         * polygon \param polygon_id. Index between 0 and @function nb() - 1.
         */
        void set_polygon_vertex(
            index_t polygon_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            do_set_polygon_vertex( polygon_id, local_vertex_id, vertex_id );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Sets an adjacent polygon by both its polygon \param polygon_id
         * and its local edge index \param edge_id.
         * @param[in] polygon_id the polygon index
         * @param[in] edge_id the local index of an edge in polygon \p polygon_id
         * @param[in] specifies the polygon adjacent to \param polygon_id along edge
         * \param edge_id or GEO::NO_FACET if the parameter \param edge_id is
         * on the border.
         */
        void set_polygon_adjacent(
            index_t polygon_id,
            index_t edge_id,
            index_t specifies )
        {
            do_set_polygon_adjacent( polygon_id, edge_id, specifies );
        }
        /*!
         * @brief Removes all the polygons and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_polygons( bool keep_attributes, bool keep_memory )
        {
            do_clear_polygons( keep_attributes, keep_memory );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Retrieve the adjacencies of polygons
         */
        void connect_polygons()
        {
            std::vector< index_t > polygons_to_connect(
                surface_mesh_.nb_polygons() );
            std::iota( polygons_to_connect.begin(), polygons_to_connect.end(), 0 );
            connect_polygons( polygons_to_connect );
        }
        void connect_polygons( const std::vector< index_t >& polygons_to_connect )
        {
            struct PolygonLocalVertex {
                PolygonLocalVertex( index_t polygon, index_t local_vertex )
                    : polygon_( polygon ), local_vertex_( local_vertex )
                {
                }
                index_t polygon_ { NO_ID };
                index_t local_vertex_ { NO_ID };
            };

            index_t nb_local_vertices = 0;
            for( index_t polygon : polygons_to_connect ) {
                nb_local_vertices += this->surface_mesh_.nb_polygon_vertices(
                    polygon );
            }

            std::vector< PolygonLocalVertex > polygon_vertices;
            polygon_vertices.reserve( nb_local_vertices );
            for( index_t polygon : polygons_to_connect ) {
                for( index_t v : range(
                    this->surface_mesh_.nb_polygon_vertices( polygon ) ) ) {
                    polygon_vertices.emplace_back( polygon, v );
                }
            }

            std::vector< index_t > next_local_vertex_around_vertex(
                nb_local_vertices, NO_ID );
            std::vector< index_t > vertex2polygon_local_vertex(
                this->surface_mesh_.nb_vertices(), NO_ID );
            index_t local_vertex_count = 0;
            for( index_t polygon : polygons_to_connect ) {
                for( index_t v = 0;
                    v < this->surface_mesh_.nb_polygon_vertices( polygon );
                    v++, local_vertex_count++ ) {
                    index_t vertex = this->surface_mesh_.polygon_vertex(
                        ElementLocalVertex( polygon, v ) );
                    next_local_vertex_around_vertex[local_vertex_count] =
                        vertex2polygon_local_vertex[vertex];
                    vertex2polygon_local_vertex[vertex] = local_vertex_count;
                }
            }

            local_vertex_count = 0;
            for( index_t polygon : polygons_to_connect ) {
                for( index_t v = 0;
                    v < this->surface_mesh_.nb_polygon_vertices( polygon );
                    v++, local_vertex_count++ ) {
                    if( !this->surface_mesh_.is_edge_on_border(
                        PolygonLocalEdge( polygon, v ) ) ) {
                        continue;
                    }
                    index_t vertex = this->surface_mesh_.polygon_vertex(
                        ElementLocalVertex( polygon, v ) );
                    index_t next_vertex = this->surface_mesh_.polygon_vertex(
                        this->surface_mesh_.next_polygon_vertex(
                            ElementLocalVertex( polygon, v ) ) );
                    for( index_t local_vertex =
                        vertex2polygon_local_vertex[next_vertex];
                        local_vertex != NO_ID; local_vertex =
                            next_local_vertex_around_vertex[local_vertex] ) {
                        if( local_vertex == local_vertex_count ) {
                            continue;
                        }
                        index_t adj_polygon = polygon_vertices[local_vertex].polygon_;
                        index_t adj_local_vertex =
                            polygon_vertices[local_vertex].local_vertex_;
                        index_t adj_next_vertex = this->surface_mesh_.polygon_vertex(
                            this->surface_mesh_.next_polygon_vertex(
                                ElementLocalVertex( adj_polygon,
                                    adj_local_vertex ) ) );
                        if( adj_next_vertex == vertex ) {
                            this->set_polygon_adjacent( polygon, v, adj_polygon );
                            this->set_polygon_adjacent( adj_polygon,
                                adj_local_vertex, polygon );
                            break;
                        }
                    }
                }
            }
        }

        void permute_polygons( const std::vector< index_t >& permutation )
        {
            do_permute_polygons( permutation );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Deletes a set of polygons.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that are
         * no longer incident to any entity are deleted.
         */
        void delete_polygons(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices )
        {
            do_delete_polygons( to_delete );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_polygon_linked_objects();
        }

        /*!@}
         * \name SurfaceMesh algorithms
         * @{
         */
        /**
         * \brief Removes the connected components that have an area
         *  smaller than a given threshold.
         * \param[in] min_area the connected components with an
         *  area smaller than this threshold are removed
         * \param[in] min_polygons the connected components with
         *  less than \param min_polygons polygons are removed
         */
        virtual void remove_small_connected_components(
            double min_area,
            index_t min_polygons ) = 0;
        virtual void triangulate(
            const SurfaceMeshBase< DIMENSION >& surface_in ) = 0;
        /*!@}
         */
        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        void remove_isolated_vertices()
        {
            std::vector< bool > to_delete( surface_mesh_.nb_vertices(), true );
            for( index_t p : range( surface_mesh_.nb_polygons() ) ) {
                for( index_t v : range( surface_mesh_.nb_polygon_vertices( p ) ) ) {
                    index_t vertex_id = surface_mesh_.polygon_vertex(
                        ElementLocalVertex( p, v ) );
                    to_delete[vertex_id] = false;
                }
            }
            this->delete_vertices( to_delete );
        }
    protected:
        SurfaceMeshBuilder( SurfaceMeshBase< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), surface_mesh_( mesh )
        {
        }

    private:
        /*!
         * @brief Deletes the NNSearch on polygons
         */
        void delete_polygon_nn_search()
        {
            surface_mesh_.nn_search_.reset();
        }

        /*!
         * @brief Deletes the AABB on polygons
         */
        void delete_polygon_aabb()
        {
            surface_mesh_.polygon_aabb_.reset();
        }

        void clear_vertex_linked_objects() override
        {
            this->delete_vertex_nn_search();
            clear_polygon_linked_objects();
        }

        void clear_polygon_linked_objects()
        {
            delete_polygon_aabb();
            delete_polygon_nn_search();
        }
        /*!
         * brief create polygons
         * @param[in] polygons is the vector of vertex index for each polygon
         * @param[in] polygon_ptr is the vector addressing the first polygon vertex for each polygon.
         */
        virtual void do_create_polygons(
            const std::vector< index_t >& polygons,
            const std::vector< index_t >& polygon_ptr ) = 0;
        /*!
         * \brief Creates a polygon
         * \param[in] vertices a const reference to a vector that
         *  contains the vertices
         * \return the index of the created polygon
         */
        virtual index_t do_create_polygon(
            const std::vector< index_t >& vertices ) = 0;
        /*!
         * \brief Creates a contiguous chunk of triangles
         * \param[in] nb_triangles number of triangles to create
         * \return the index of the first triangle
         */
        virtual index_t do_create_triangles( index_t nb_triangles ) = 0;
        /*!
         * \brief Creates a contiguous chunk of quads
         * \param[in] nb_quads number of quads to create
         * \return the index of the first quad
         */
        virtual index_t do_create_quads( index_t nb_quads ) = 0;
        /*!
         * @brief Sets a vertex of a polygon by local vertex index.
         * @param[in] polygon_id index of the polygon, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the polygon.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of the
         * polygon \param polygon_id. Index between 0 and @function nb() - 1.
         */
        virtual void do_set_polygon_vertex(
            index_t polygon_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0;
        /*!
         * @brief Sets an adjacent polygon by both its polygon \param polygon_id
         * and its local edge index \param edge_id.
         * @param[in] polygon_id the polygon index
         * @param[in] edge_id the local index of an edge in polygon \p polygon_id
         * @param[in] specifies the polygon adjacent to \param polygon_id along edge
         * \param edge_id or GEO::NO_FACET if the parameter \param edge_id is
         * on the border.
         */
        virtual void do_set_polygon_adjacent(
            index_t polygon_id,
            index_t edge_id,
            index_t specifies ) = 0;
        /*!
         * @brief Removes all the polygons and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_polygons( bool keep_attributes, bool keep_memory ) = 0;

        virtual void do_permute_polygons(
            const std::vector< index_t >& permutation ) = 0;
        /*!
         * @brief Deletes a set of polygons.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         */
        virtual void do_delete_polygons( const std::vector< bool >& to_delete ) = 0;

    protected:
        SurfaceMeshBase< DIMENSION >& surface_mesh_;
    };

    CLASS_DIMENSION_ALIASES( SurfaceMeshBuilder );

    template< index_t DIMENSION >
    using SurfaceMeshBuilderFactory = Factory< MeshType, SurfaceMeshBuilder< DIMENSION >, SurfaceMesh< DIMENSION >& >;

    CLASS_DIMENSION_ALIASES( SurfaceMeshBuilderFactory );

    template< index_t DIMENSION >
    class VolumeMeshBuilder: public MeshBaseBuilder< DIMENSION > {
    ringmesh_disable_copy( VolumeMeshBuilder );
        static_assert( DIMENSION == 3, "DIMENSION template should be 3" );
    public:
        virtual ~VolumeMeshBuilder() = default;

        static std::unique_ptr< VolumeMeshBuilder< DIMENSION > > create_builder(
            VolumeMesh< DIMENSION >& mesh );

        /*!
         * @brief Creates a contiguous chunk of cells of the same type.
         * @param[in] nb_cells number of cells to create
         * @param[in] type type of the cells to create, one of TETRAEDRON, HEXAEDRON,
         * CellType::PRISM, CellType::PYRAMID, CellType::UNCLASSIFIED.
         * @return the first created cell.
         */
        index_t create_cells( index_t nb_cells, CellType type )
        {
            index_t index = do_create_cells( nb_cells, type );
            clear_cell_linked_objects();
            return index;
        }
        /*
         * \brief Copies a tets mesh into this Mesh.
         * \details Cells adjacence are not computed.
         *   cell and corner attributes are zeroed.
         * \param[in] tets cells to vertex links
         * (using vector::swap).
         */
        void assign_cell_tet_mesh( const std::vector< index_t >& tets )
        {
            do_assign_cell_tet_mesh( tets );
            clear_cell_linked_objects();
        }
        /*!
         * @brief Sets a vertex of a cell by local vertex index.
         * @param[in] cell_id index of the cell, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the cell.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the global index of the vertex \param
         * local_vertex_id in the cell \param cell_id. Index between 0 and
         * @function nb() - 1.
         */
        void set_cell_vertex(
            index_t cell_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            do_set_cell_vertex( cell_id, local_vertex_id, vertex_id );
            clear_cell_linked_objects();
        }
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0.. @function nb() - 1
         * \param[in] vertex_index specifies the vertex that corner
         * \param corner_index is incident to
         */
        void set_cell_corner_vertex_index(
            index_t corner_index,
            index_t vertex_index )
        {
            do_set_cell_corner_vertex_index( corner_index, vertex_index );
            clear_cell_linked_objects();
        }
        /*!
         * \brief Sets the cell adjacent
         * \param[in] cell_index index of the cell
         * \param[in] facet_index local index of the cell facet
         * \param[in] cell_adjacent adjacent value to set
         */
        void set_cell_adjacent(
            index_t cell_index,
            index_t facet_index,
            index_t cell_adjacent )
        {
            do_set_cell_adjacent( cell_index, facet_index, cell_adjacent );
        }

        /*!
         * @brief Retrieve the adjacencies
         */
        virtual void connect_cells() = 0;

        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_cells( bool keep_attributes, bool keep_memory )
        {
            do_clear_cells( keep_attributes, keep_memory );
            clear_cell_linked_objects();
        }
        /*!
         * @brief Applies a permutation to the entities and their attributes.
         * On exit, permutation is modified (used for internal bookkeeping).
         * Applying a permutation permutation is equivalent to:
         * <code>
         *  for( i = 0 ; i < permutation.size() ; i++) {
         *      data2[i] = data[permutation[i]]
         *       }
         *  data = data2 ;
         *  </code>
         */
        void permute_cells( const std::vector< index_t >& permutation )
        {
            do_permute_cells( permutation );
            clear_cell_linked_objects();
        }
        /*!
         * @brief Deletes a set of cells.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that are
         * no longer incident to any entity are deleted.
         */
        void delete_cells(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices )
        {
            do_delete_cells( to_delete );
            if( remove_isolated_vertices ) {
                this->remove_isolated_vertices();
            }
            clear_cell_linked_objects();
        }

        void remove_isolated_vertices()
        {
            std::vector< bool > to_delete( volume_mesh_.nb_vertices(), true );
            for( index_t c : range( volume_mesh_.nb_cells() ) ) {
                for( index_t v : range( volume_mesh_.nb_cell_vertices( c ) ) ) {
                    index_t vertex_id = volume_mesh_.cell_vertex(
                        ElementLocalVertex( c, v ) );
                    to_delete[vertex_id] = false;
                }
            }
            this->delete_vertices( to_delete );
        }
    protected:
        VolumeMeshBuilder( VolumeMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), volume_mesh_( mesh )
        {
        }

    private:
        /*!
         * @brief Deletes the NNSearch on cells
         */
        void delete_cell_nn_search()
        {
            volume_mesh_.cell_nn_search_.reset();
            volume_mesh_.cell_facet_nn_search_.reset();
        }

        /*!
         * @brief Deletes the AABB on cells
         */
        void delete_cell_aabb()
        {
            volume_mesh_.cell_aabb_.reset();
        }

        void clear_vertex_linked_objects() override
        {
            this->delete_vertex_nn_search();
            clear_cell_linked_objects();
        }

        void clear_cell_linked_objects()
        {
            delete_cell_aabb();
            delete_cell_nn_search();
        }

        /*!
         * @brief Creates a contiguous chunk of cells of the same type.
         * @param[in] nb_cells number of cells to create
         * @param[in] type type of the cells to create, one of TETRAEDRON, HEXAEDRON,
         * CellType::PRISM, CellType::PYRAMID, CellType::UNCLASSIFIED.
         * @return the first created cell.
         */
        virtual index_t do_create_cells( index_t nb_cells, CellType type ) = 0;
        /*
         * \brief Copies a tets mesh into this Mesh.
         * \details Cells adjacence are not computed.
         *   cell and corner attributes are zeroed.
         * \param[in] tets cells to vertex links
         * (using vector::swap).
         */
        virtual void do_assign_cell_tet_mesh(
            const std::vector< index_t >& tets ) = 0;
        /*!
         * @brief Sets a vertex of a cell by local vertex index.
         * @param[in] cell_id index of the cell, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the cell.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the global index of the vertex \param
         * local_vertex_id in the cell \param cell_id. Index between 0 and
         * @function nb() - 1.
         */
        virtual void do_set_cell_vertex(
            index_t cell_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0;
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0.. @function nb() - 1
         * \param[in] vertex_index specifies the vertex that corner
         * \param corner_index is incident to
         */
        virtual void do_set_cell_corner_vertex_index(
            index_t corner_index,
            index_t vertex_index ) = 0;
        /*!
         * \brief Sets the cell adjacent
         * \param[in] cell_index index of the cell
         * \param[in] facet_index local index of the cell facet
         * \param[in] cell_adjacent adjacent value to set
         */
        virtual void do_set_cell_adjacent(
            index_t cell_index,
            index_t facet_index,
            index_t cell_adjacent ) = 0;
        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_cells( bool keep_attributes, bool keep_memory ) = 0;
        /*!
         * @brief Applies a permutation to the entities and their attributes.
         * On exit, permutation is modified (used for internal bookkeeping).
         * Applying a permutation permutation is equivalent to:
         * <code>
         *  for( i = 0 ; i < permutation.size() ; i++) {
         *      data2[i] = data[permutation[i]]
         *       }
         *  data = data2 ;
         *  </code>
         */
        virtual void do_permute_cells(
            const std::vector< index_t >& permutation ) = 0;
        /*!
         * @brief Deletes a set of cells.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         */
        virtual void do_delete_cells( const std::vector< bool >& to_delete ) = 0;

    protected:
        VolumeMesh< DIMENSION >& volume_mesh_;
    };

    using VolumeMeshBuilder3D = VolumeMeshBuilder< 3 >;

    template< index_t DIMENSION >
    using VolumeMeshBuilderFactory = Factory< MeshType, VolumeMeshBuilder< DIMENSION >, VolumeMesh< DIMENSION >& >;

    using VolumeMeshBuilderFactory3D = VolumeMeshBuilderFactory< 3 >;
}
