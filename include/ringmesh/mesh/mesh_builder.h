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

#include <memory>
#include <numeric>

#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/basic/factory.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMesh );

    struct CellLocalFacet;
    struct EdgeLocalVertex;
    struct ElementLocalVertex;
    struct PolygonLocalEdge;
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class MeshBaseBuilder
    {
        ringmesh_disable_copy_and_move( MeshBaseBuilder );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~MeshBaseBuilder() = default;
        /*!
         * \name general methods
         * @{
         */
        /*!
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @param[in] copy_attributes if true, all attributes are copied.
         * @return a modifiable reference to the point that corresponds to the
         * vertex.
         */
        void copy( const MeshBase< DIMENSION >& rhs, bool copy_attributes );

        void load_mesh( const std::string& filename )
        {
            clear( false, false );
            do_load_mesh( filename );
            update_attributes_size();
        }
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear( bool keep_attributes, bool keep_memory );

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
        void set_vertex( index_t v_id, const vecn< DIMENSION >& vertex );

        /*!
         * @brief Creates a new vertex.
         * @return the index of the created vertex
         */
        index_t create_vertex();

        /*!
         * @brief Creates a new vertex.
         * @param[in] coords a pointer to @function dimension() coordinate.
         * @return the index of the created vertex
         */
        index_t create_vertex( const vecn< DIMENSION >& vertex );

        /*!
         * @brief Creates a contiguous chunk of vertices.
         * @param[in] nb number of sub-entities to create.
         * @return the index of the first created vertex
         */
        index_t create_vertices( index_t nb );

        /*!
         * @brief set vertex coordinates from a std::vector of coordinates
         * @param[in] point_coordinates a set of x, y (, z) coordinates
         */
        void assign_vertices( const std::vector< double >& point_coordinates );

        /*!
         * @brief Deletes a set of vertices.
         * @param[in] to_delete     a vector of size @function nb(). If
         * to_delete[e] is true,
         * then entity e will be destroyed, else it will be kept.
         */
        void delete_vertices( const std::vector< bool >& to_delete );

        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_vertices( bool keep_attributes, bool keep_memory );

        void permute_vertices( const std::vector< index_t >& permutation );
        /*!@}
         */

        static std::unique_ptr< MeshBaseBuilder< DIMENSION > > create_builder(
            MeshBase< DIMENSION >& mesh );

    protected:
        explicit MeshBaseBuilder( MeshBase< DIMENSION >& mesh )
            : mesh_base_( mesh )
        {
        }

        void delete_vertex_nn_search();

        /*!
         * @brief Deletes the NNSearch on vertices
         */
        virtual void clear_vertex_linked_objects() = 0;

        void update_vertex_attributes_size();

    private:
        virtual void update_attributes_size() = 0;

        virtual void do_load_mesh( const std::string& filename ) = 0;
        /*!
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @param[in] copy_attributes if true, all attributes are copied.
         * @return a modifiable reference to the point that corresponds to the
         * vertex.
         */
        virtual void do_copy(
            const MeshBase< DIMENSION >& rhs, bool copy_attributes ) = 0;
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear( bool keep_attributes, bool keep_memory ) = 0;
        /*!
         * @brief Sets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @param[in] vertex the vertex coordinates
         * @return reference to the point that corresponds to the vertex.
         */
        virtual void do_set_vertex(
            index_t v_id, const vecn< DIMENSION >& vertex ) = 0;
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
         * @param[in] to_delete     a vector of size @function nb(). If
         * to_delete[e] is true,
         * then entity e will be destroyed, else it will be kept.
         */
        virtual void do_delete_vertices(
            const std::vector< bool >& to_delete ) = 0;
        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_vertices(
            bool keep_attributes, bool keep_memory ) = 0;
        virtual void do_permute_vertices(
            const std::vector< index_t >& permutation ) = 0;

    protected:
        MeshBase< DIMENSION >& mesh_base_;
    };

    ALIAS_2D_AND_3D( MeshBaseBuilder );

    template < index_t DIMENSION >
    class PointSetMeshBuilder : public MeshBaseBuilder< DIMENSION >
    {
    public:
        static std::unique_ptr< PointSetMeshBuilder< DIMENSION > >
            create_builder( PointSetMesh< DIMENSION >& mesh );

        void remove_isolated_vertices()
        {
            // All vertices are isolated in a Mesh0D
        }

    protected:
        explicit PointSetMeshBuilder( PointSetMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), pointset_mesh_( mesh )
        {
        }

    private:
        void update_attributes_size() final
        {
            this->update_vertex_attributes_size();
        }

        void clear_vertex_linked_objects() final
        {
            this->delete_vertex_nn_search();
        }

    protected:
        PointSetMesh< DIMENSION >& pointset_mesh_;
    };

    ALIAS_2D_AND_3D( PointSetMeshBuilder );

    template < index_t DIMENSION >
    using PointSetMeshBuilderFactory = Factory< MeshType,
        PointSetMeshBuilder< DIMENSION >,
        PointSetMesh< DIMENSION >& >;

    ALIAS_2D_AND_3D( PointSetMeshBuilderFactory );

    template < index_t DIMENSION >
    class LineMeshBuilder : public MeshBaseBuilder< DIMENSION >
    {
    public:
        static std::unique_ptr< LineMeshBuilder > create_builder(
            LineMesh< DIMENSION >& mesh );

        /*!
         * @brief Create a new edge.
         * @param[in] v1_id index of the starting vertex.
         * @param[in] v2_id index of the ending vertex.
         */
        void create_edge( index_t v1_id, index_t v2_id );

        /*!
         * \brief Creates a contiguous chunk of edges
         * \param[in] nb_edges number of edges to create
         * \return the index of the first edge
         */
        index_t create_edges( index_t nb_edges );

        /*!
         * @brief Sets a vertex of a edge by local vertex index.
         * @param[in] edge_local_vertex index of the edge and local index of the
         * vertex in the edge.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of
         * edge
         * \param edge_id. Index between 0 and @function nb() - 1.
         */
        void set_edge_vertex(
            const EdgeLocalVertex& edge_local_vertex, index_t vertex_id );

        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices
         * that are no longer incident to any entity are deleted.
         */
        void delete_edges( const std::vector< bool >& to_delete,
            bool remove_isolated_vertices );

        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_edges( bool keep_attributes, bool keep_memory );

        void permute_edges( const std::vector< index_t >& permutation );

        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        void remove_isolated_vertices();

    protected:
        explicit LineMeshBuilder( LineMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), line_mesh_( mesh )
        {
        }

        void update_edge_attributes_size()
        {
            line_mesh_.edge_attribute_manager().resize( line_mesh_.nb_edges() );
        }

    private:
        void update_attributes_size() final
        {
            this->update_vertex_attributes_size();
            update_edge_attributes_size();
        }
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
         * @param[in] edge_local_vertex index of the edge and local index of the
         * vertex in the edge.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of
         * edge
         * \param edge_id. Index between 0 and @function nb() - 1.
         */
        virtual void do_set_edge_vertex(
            const EdgeLocalVertex& edge_local_vertex, index_t vertex_id ) = 0;
        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         */
        virtual void do_delete_edges(
            const std::vector< bool >& to_delete ) = 0;
        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_edges(
            bool keep_attributes, bool keep_memory ) = 0;
        virtual void do_permute_edges(
            const std::vector< index_t >& permutation ) = 0;

    protected:
        LineMesh< DIMENSION >& line_mesh_;
    };

    ALIAS_2D_AND_3D( LineMeshBuilder );

    template < index_t DIMENSION >
    using LineMeshBuilderFactory = Factory< MeshType,
        LineMeshBuilder< DIMENSION >,
        LineMesh< DIMENSION >& >;

    ALIAS_2D_AND_3D( LineMeshBuilderFactory );

    template < index_t DIMENSION >
    class SurfaceMeshBuilder : public MeshBaseBuilder< DIMENSION >
    {
    public:
        static std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >
            create_builder( SurfaceMesh< DIMENSION >& mesh );

        /*!@}
         * \name Polygon related methods
         * @{
         */
        /*!
         * brief create polygons
         * @param[in] polygons is the vector of vertex index for each polygon
         * @param[in] polygon_ptr is the vector addressing the first polygon
         * vertex for each polygon.
         */
        void create_polygons( const std::vector< index_t >& polygons,
            const std::vector< index_t >& polygon_ptr )
        {
            for( auto p : range( polygon_ptr.size() - 1 ) )
            {
                index_t first{ polygon_ptr[p] };
                index_t last{ polygon_ptr[p + 1] };
                index_t nb_to_copy{ last - first };
                std::vector< index_t > polygon_vertices( nb_to_copy );
                for( auto i : range( nb_to_copy ) )
                {
                    polygon_vertices[i] = polygons[first + i];
                }
                do_create_polygon( polygon_vertices );
            }
            clear_polygon_linked_objects();
            update_polygon_attributes_size();
        }
        /*!
         * \brief Creates a polygon
         * \param[in] vertices a const reference to a vector that
         *  contains the vertices
         * \return the index of the created polygon
         */
        index_t create_polygon( const std::vector< index_t >& vertices )
        {
            auto index = do_create_polygon( vertices );
            clear_polygon_linked_objects();
            update_polygon_attributes_size();
            return index;
        }
        /*!
         * \brief Creates a contiguous chunk of triangles
         * \param[in] nb_triangles number of triangles to create
         * \return the index of the first triangle
         */
        index_t create_triangles( index_t nb_triangles )
        {
            auto index = do_create_triangles( nb_triangles );
            clear_polygon_linked_objects();
            update_polygon_attributes_size();
            return index;
        }
        /*!
         * \brief Creates a contiguous chunk of quads
         * \param[in] nb_quads number of quads to create
         * \return the index of the first quad
         */
        index_t create_quads( index_t nb_quads )
        {
            auto index = do_create_quads( nb_quads );
            clear_polygon_linked_objects();
            update_polygon_attributes_size();
            return index;
        }
        /*!
         * @brief Sets a vertex of a polygon by local vertex index.
         * @param[in] polygon_local_edge the polygon index and local index of an
         * edge.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of
         * the
         * polygon \param polygon_id. Index between 0 and @function nb() - 1.
         */
        void set_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex, index_t vertex_id )
        {
            do_set_polygon_vertex( polygon_local_vertex, vertex_id );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Sets an adjacent polygon by both its polygon \param polygon_id
         * and its local edge index \param edge_id.
         * @param[in] polygon_local_edge the polygon index and local index of an
         * edge.
         * @param[in] specifies the polygon adjacent to \param polygon_id along
         * edge
         * \param edge_id or GEO::NO_FACET if the parameter \param edge_id is
         * on the border.
         */
        void set_polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge, index_t specifies )
        {
            do_set_polygon_adjacent( polygon_local_edge, specifies );
        }
        /*!
         * @brief Removes all the polygons and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_polygons( bool keep_attributes, bool keep_memory )
        {
            do_clear_polygons( keep_attributes, keep_memory );
            surface_mesh_.polygon_attributes_manager_.clear( keep_attributes );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Retrieve the adjacencies of polygons
         */
        void connect_polygons()
        {
            std::vector< index_t > polygons_to_connect(
                surface_mesh_.nb_polygons() );
            std::iota(
                polygons_to_connect.begin(), polygons_to_connect.end(), 0 );
            connect_polygons( polygons_to_connect );
        }
        void connect_polygons(
            const std::vector< index_t >& polygons_to_connect )
        {
            // Initialization of the number of local vertices
            index_t nb_local_vertices{ 0 };
            for( auto polygon : polygons_to_connect )
            {
                nb_local_vertices +=
                    this->surface_mesh_.nb_polygon_vertices( polygon );
            }

            // Initialization of the polygon vertices
            std::vector< ElementLocalVertex > polygon_vertices;
            polygon_vertices.reserve( nb_local_vertices );
            for( auto polygon : polygons_to_connect )
            {
                for( auto v : range(
                         this->surface_mesh_.nb_polygon_vertices( polygon ) ) )
                {
                    polygon_vertices.emplace_back( polygon, v );
                }
            }

            std::vector< index_t > next_local_vertex_around_vertex(
                nb_local_vertices, NO_ID );
            std::vector< index_t > vertex2polygon_local_vertex(
                this->surface_mesh_.nb_vertices(), NO_ID );
            index_t local_vertex_count{ 0 };
            for( auto polygon : polygons_to_connect )
            {
                for( index_t v{ 0 };
                     v < this->surface_mesh_.nb_polygon_vertices( polygon );
                     v++, local_vertex_count++ )
                {
                    auto vertex =
                        this->surface_mesh_.polygon_vertex( { polygon, v } );
                    next_local_vertex_around_vertex[local_vertex_count] =
                        vertex2polygon_local_vertex[vertex];
                    vertex2polygon_local_vertex[vertex] = local_vertex_count;
                }
            }

            local_vertex_count = 0;
            for( auto polygon : polygons_to_connect )
            {
                for( index_t v = 0;
                     v < this->surface_mesh_.nb_polygon_vertices( polygon );
                     v++, local_vertex_count++ )
                {
                    if( !this->surface_mesh_.is_edge_on_border(
                            { polygon, v } ) )
                    {
                        continue;
                    }
                    auto vertex =
                        this->surface_mesh_.polygon_vertex( { polygon, v } );
                    auto next_vertex = this->surface_mesh_.polygon_vertex(
                        this->surface_mesh_.next_polygon_vertex(
                            { polygon, v } ) );
                    for( auto local_vertex =
                             vertex2polygon_local_vertex[next_vertex];
                         local_vertex != NO_ID;
                         local_vertex =
                             next_local_vertex_around_vertex[local_vertex] )
                    {
                        if( local_vertex == local_vertex_count )
                        {
                            continue;
                        }
                        auto adj_polygon =
                            polygon_vertices[local_vertex].element_id;
                        auto adj_local_vertex =
                            polygon_vertices[local_vertex].local_vertex_id;
                        auto adj_next_vertex =
                            this->surface_mesh_.polygon_vertex(
                                this->surface_mesh_.next_polygon_vertex(
                                    { adj_polygon, adj_local_vertex } ) );
                        if( adj_next_vertex == vertex )
                        {
                            this->set_polygon_adjacent(
                                { polygon, v }, adj_polygon );
                            this->set_polygon_adjacent(
                                { adj_polygon, adj_local_vertex }, polygon );
                            break;
                        }
                    }
                }
            }
        }

        void permute_polygons( const std::vector< index_t >& permutation )
        {
            do_permute_polygons( permutation );
            surface_mesh_.polygon_attributes_manager_.apply_permutation(
                permutation );
            clear_polygon_linked_objects();
        }
        /*!
         * @brief Deletes a set of polygons.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that
         * are
         * no longer incident to any entity are deleted.
         */
        void delete_polygons( const std::vector< bool >& to_delete,
            bool remove_isolated_vertices )
        {
            do_delete_polygons( to_delete );
            if( remove_isolated_vertices )
            {
                this->remove_isolated_vertices();
            }
            clear_polygon_linked_objects();
            update_polygon_attributes_size();
        }

        /*!@}
         */
        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        void remove_isolated_vertices();

    protected:
        explicit SurfaceMeshBuilder( SurfaceMeshBase< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), surface_mesh_( mesh )
        {
        }

        void update_polygon_attributes_size()
        {
            surface_mesh_.polygon_attribute_manager().resize(
                surface_mesh_.nb_polygons() );
        }

    private:
        void update_attributes_size() final
        {
            this->update_vertex_attributes_size();
            update_polygon_attributes_size();
        }

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
         * @param[in] polygon_local_vertex the polygon index and the local index
         * of a vertex in the polygon.
         * @param[in] vertex_id specifies the vertex between 0 and the number
         * of vertex in polygon.
         */
        virtual void do_set_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex,
            index_t vertex_id ) = 0;
        /*!
         * @brief Sets an adjacent polygon by both its polygon \param polygon_id
         * and its local edge index \param edge_id.
         * @param[in] polygon_local_edge the polygon index and the local index
         * of an edge.
         * @param[in] specifies the polygon adjacent to \param polygon_id along
         * edge
         * \param edge_id or GEO::NO_FACET if the parameter \param edge_id is
         * on the border.
         */
        virtual void do_set_polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge, index_t specifies ) = 0;
        /*!
         * @brief Removes all the polygons and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_polygons(
            bool keep_attributes, bool keep_memory ) = 0;

        virtual void do_permute_polygons(
            const std::vector< index_t >& permutation ) = 0;
        /*!
         * @brief Deletes a set of polygons.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         */
        virtual void do_delete_polygons(
            const std::vector< bool >& to_delete ) = 0;

    protected:
        SurfaceMeshBase< DIMENSION >& surface_mesh_;
    };

    ALIAS_2D_AND_3D( SurfaceMeshBuilder );

    template < index_t DIMENSION >
    using SurfaceMeshBuilderFactory = Factory< MeshType,
        SurfaceMeshBuilder< DIMENSION >,
        SurfaceMesh< DIMENSION >& >;

    ALIAS_2D_AND_3D( SurfaceMeshBuilderFactory );

    template < index_t DIMENSION >
    class VolumeMeshBuilder : public MeshBaseBuilder< DIMENSION >
    {
        static_assert( DIMENSION == 3, "DIMENSION template should be 3" );

    public:
        static std::unique_ptr< VolumeMeshBuilder< DIMENSION > > create_builder(
            VolumeMesh< DIMENSION >& mesh );

        /*!
         * @brief Creates a contiguous chunk of cells of the same type.
         * @param[in] nb_cells number of cells to create
         * @param[in] type type of the cells to create, one of TETRAEDRON,
         * HEXAEDRON,
         * CellType::PRISM, CellType::PYRAMID, CellType::UNCLASSIFIED.
         * @return the first created cell.
         */
        index_t create_cells( index_t nb_cells, CellType type )
        {
            index_t index = do_create_cells( nb_cells, type );
            clear_cell_linked_objects();
            update_cell_attributes_size();
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
            update_cell_attributes_size();
        }
        /*!
         * @brief Sets a vertex of a cell by local vertex index.
         * @param[in] cell_local_vertex index of the cell, and local index of
         * the vertex in the cell.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the global index of the vertex \param
         * local_vertex_id in the cell \param cell_id. Index between 0 and
         * @function nb() - 1.
         */
        void set_cell_vertex(
            const ElementLocalVertex& cell_local_vertex, index_t vertex_id )
        {
            do_set_cell_vertex( cell_local_vertex, vertex_id );
            clear_cell_linked_objects();
        }
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0.. @function nb() - 1
         * \param[in] vertex_index specifies the vertex that corner
         * \param corner_index is incident to
         */
        void set_cell_corner_vertex_index(
            index_t corner_index, index_t vertex_index )
        {
            do_set_cell_corner_vertex_index( corner_index, vertex_index );
            clear_cell_linked_objects();
        }
        /*!
         * \brief Sets the cell adjacent
         * \param[in] cell_local_facet index of the cell, and local index of the
         * cell facet
         * \param[in] cell_adjacent adjacent value to set
         */
        void set_cell_adjacent(
            const CellLocalFacet& cell_local_facet, index_t cell_adjacent )
        {
            do_set_cell_adjacent( cell_local_facet, cell_adjacent );
        }

        /*!
         * @brief Retrieve the adjacencies
         */
        virtual void connect_cells() = 0;

        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_cells( bool keep_attributes, bool keep_memory )
        {
            do_clear_cells( keep_attributes, keep_memory );
            volume_mesh_.cell_attributes_manager_.clear( keep_attributes );
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
            volume_mesh_.cell_attributes_manager_.apply_permutation(
                permutation );
            clear_cell_linked_objects();
        }
        /*!
         * @brief Deletes a set of cells.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that
         * are
         * no longer incident to any entity are deleted.
         */
        void delete_cells( const std::vector< bool >& to_delete,
            bool remove_isolated_vertices )
        {
            do_delete_cells( to_delete );
            if( remove_isolated_vertices )
            {
                this->remove_isolated_vertices();
            }
            clear_cell_linked_objects();
            update_cell_attributes_size();
        }

        void remove_isolated_vertices();

    protected:
        explicit VolumeMeshBuilder( VolumeMesh< DIMENSION >& mesh )
            : MeshBaseBuilder< DIMENSION >( mesh ), volume_mesh_( mesh )
        {
        }

        void update_cell_attributes_size();

    private:
        void update_attributes_size() final
        {
            this->update_vertex_attributes_size();
            update_cell_attributes_size();
        }
        /*!
         * @brief Deletes the NNSearch on cells
         */
        void delete_cell_nn_search();

        /*!
         * @brief Deletes the AABB on cells
         */
        void delete_cell_aabb();

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
         * @param[in] type type of the cells to create, one of TETRAEDRON,
         * HEXAEDRON,
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
         * @param[in] cell_local_vertex index of the cell,and local index of the
         * vertex in the cell.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the global index of the vertex \param
         * local_vertex_id in the cell \param cell_id. Index between 0 and
         * @function nb() - 1.
         */
        virtual void do_set_cell_vertex(
            const ElementLocalVertex& cell_local_vertex,
            index_t vertex_id ) = 0;
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0.. @function nb() - 1
         * \param[in] vertex_index specifies the vertex that corner
         * \param corner_index is incident to
         */
        virtual void do_set_cell_corner_vertex_index(
            index_t corner_index, index_t vertex_index ) = 0;
        /*!
         * \brief Sets the cell adjacent
         * \param[in] cell_local_facet index of the cell, and local index of the
         * cell facet
         * \param[in] cell_adjacent adjacent value to set
         */
        virtual void do_set_cell_adjacent(
            const CellLocalFacet& cell_local_facet, index_t cell_adjacent ) = 0;
        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are
         * destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void do_clear_cells(
            bool keep_attributes, bool keep_memory ) = 0;
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
         * If to_delete[e] is true, then entity e will be destroyed, else it
         * will be kept.
         */
        virtual void do_delete_cells(
            const std::vector< bool >& to_delete ) = 0;

    protected:
        VolumeMesh< DIMENSION >& volume_mesh_;
    };

    using VolumeMeshBuilder3D = VolumeMeshBuilder< 3 >;

    template < index_t DIMENSION >
    using VolumeMeshBuilderFactory = Factory< MeshType,
        VolumeMeshBuilder< DIMENSION >,
        VolumeMesh< DIMENSION >& >;

    using VolumeMeshBuilderFactory3D = VolumeMeshBuilderFactory< 3 >;
} // namespace RINGMesh
