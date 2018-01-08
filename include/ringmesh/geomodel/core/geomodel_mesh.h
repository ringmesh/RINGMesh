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

#include <ringmesh/geomodel/core/common.h>

#include <ringmesh/basic/pimpl.h>

#include <ringmesh/geomodel/core/entity_type.h>

/*!
 * @file ringmesh/geomodel_mesh.h
 * @brief Classes to manage globally the indexing of mesh entities of a GeoModel
 * @author Arnaud Botella and Jeanne Pellerin
 */

namespace GEO
{
    class AttributesManager;
} // namespace GEO

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEdges );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshWells );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshPolygons );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshPolygonsBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshCells );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshSet );
    FORWARD_DECLARATION_DIMENSION_CLASS( NNSearch );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineAABBTree );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceAABBTree );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeAABBTree );

    ALIAS_3D( GeoModel );
    ALIAS_3D( GeoModelMesh );
    ALIAS_3D( PointSetMesh );
    ALIAS_3D( SurfaceMesh );

    struct CellLocalFacet;
    struct ElementLocalVertex;
    struct PolygonLocalEdge;
} // namespace RINGMesh

namespace RINGMesh
{
    /*! @todo Move this global variables in a function */

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshCommon
    {
        ringmesh_disable_copy_and_move( GeoModelMeshCommon );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        bool is_initialized() const
        {
            return is_initialized_;
        }

    protected:
        GeoModelMeshCommon(
            GeoModelMesh< DIMENSION >& gmm, GeoModel< DIMENSION >& geomodel );

        virtual ~GeoModelMeshCommon() = default;

        void set_mesh( MeshBase< DIMENSION >* mesh )
        {
            mesh_base_ = mesh;
        }

        void save_mesh( const std::string& filename ) const
        {
            mesh_base_->save_mesh( filename );
        }

        void set_is_initialized( bool value ) const
        {
            is_initialized_ = value;
        }

    protected:
        /// Attached GeoModelMesh
        GeoModelMesh< DIMENSION >& gmm_;
        /// Attached GeoModel
        GeoModel< DIMENSION >& geomodel_;
        /// Attached MeshBase
        MeshBase< DIMENSION >* mesh_base_;

    private:
        /// initilization flag
        mutable bool is_initialized_{ false };
    };

    struct GMEVertex
    {
        GMEVertex( gmme_id t, index_t vertex_id_in )
            : gmme( std::move( t ) ), v_index( vertex_id_in )
        {
        }
        bool operator==( const GMEVertex& rhs ) const
        {
            return gmme == rhs.gmme && v_index == rhs.v_index;
        }
        bool is_defined() const
        {
            return gmme.is_defined() && v_index != NO_ID;
        }
        /// GeoModelEntity index in the GeoModel that owns it
        gmme_id gmme;
        /// Index of the vertex in the GeoModelEntity
        index_t v_index{ NO_ID };
    };

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshVerticesBase
        : public GeoModelMeshCommon< DIMENSION >
    {
        ringmesh_disable_copy_and_move( GeoModelMeshVerticesBase );

    public:
        friend class GeoModelMeshWells< DIMENSION >;
        friend class GeoModelMeshEdges< DIMENSION >;
        friend class GeoModelMeshPolygonsBase< DIMENSION >;
        friend class GeoModelMeshCells< DIMENSION >;

        ~GeoModelMeshVerticesBase();

        GEO::AttributesManager& attribute_manager() const;

        /*!
         * Tests if the mesh vertices need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const;

        /*!
         * @brief Number of non colocated vertices stored.
         */
        index_t nb() const;

        /*!
         * @brief Coordinates of a vertex of the GeoModel
         * @pre v < nb()
         */
        const vecn< DIMENSION >& vertex( index_t v ) const;

        /*!
         * @brief Returns the index of the given vertex in the geomodel
         * @param[in] p input point coordinates
         * @return index of the vertex in the geomodel if found
         * (distance < epsilon), otherwise NO_ID
         */
        index_t index( const vecn< DIMENSION >& p ) const;

        /*!
         * @brief Get the GeoModelMesh index of a GeoModelMeshEntity vertex from
         * its
         * index in that GeoModelMeshEntity
         * @param[in] mesh_entity GeoModelMeshEntity that vertex belongs to
         * @param[in] entity_vertex_index index of the query vertex in the
         * GeoModelMeshEntity
         * @return if found the vertex index in the geomodel, else NO_ID.
         */
        index_t geomodel_vertex_id(
            const gmme_id& mesh_entity, index_t entity_vertex_index = 0 ) const;

        /*!
         * @brief Get the GeoModelMesh index of a GeoModelMeshEntity vertex from
         * its
         * index in that GeoModelMeshEntity
         * @param[in] mesh_entity GeoModelMeshEntity that vertex belongs to
         * @param[in] element_local_vertex local vertex in the element
         * @return if found the vertex index in the geomodel, else NO_ID.
         */
        index_t geomodel_vertex_id( const gmme_id& mesh_entity,
            const ElementLocalVertex& element_local_vertex ) const;

        /*!
         * @brief Get the GeoModelMeshEntity vertices from its index in the
         * GeoModelMesh
         * @param[in] mesh_entity Unique id to a GeoModelMeshEntity
         * @param[in] geomodel_vertex_id Index of the query vertex in the
         * geomodel
         * @return Corresponding GeoModelMeshEntity vertices
         */
        std::vector< index_t > mesh_entity_vertex_id(
            const gmme_id& mesh_entity, index_t geomodel_vertex_id ) const;

        /*!
         * @brief Get the vertices in GeoModelEntity corresponding to the given
         * unique vertex
         * @param[in] vertex Vertex index in the geomodel
         * @return Corresponding GeoModelMeshEntity vertices
         */
        const std::vector< GMEVertex >& gme_vertices( index_t v ) const;

        /*!
         * @brief Get the vertex indices in the specified MeshEntity type
         * corresponding to the given unique vertex
         */
        std::vector< GMEVertex > gme_type_vertices(
            const MeshEntityType& entity_type, index_t vertex ) const;

        /*!
         * @brief Set the point coordinates of a vertex
         * @param[in] vertex Index of the vertex
         * @param[in] point New coordinates
         */
        void set_point( index_t v, const vecn< DIMENSION >& point );

        void update_vertex_mapping( const gmme_id& entity_id,
            index_t entity_vertex_index,
            index_t geomodel_vertex_index );

        /*!
         * @brief Clear the vertices - clear the gme_vertices_ -
         *        clear global vertex information in the all BMME
         * @warning Not stable - crashes if attributes are still bound
         */
        void clear() const;

        void unbind_geomodel_vertex_map( const gmme_id& mesh_entity_id );

        void bind_geomodel_vertex_map( const gmme_id& mesh_entity_id );

        const NNSearch< DIMENSION >& nn_search() const;

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated() const;

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the geomodel_vertx_id in the GeoModelMeshEntity
         * of the GeoModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with which information should be merged.
         *  It is recycled to give the mapping between old and new vertex
         * indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) const;

    private:
        /*!
         * @brief Initialize the vertices from the vertices
         *        of the GeoModel Corners, Lines, Surfaces and Regions
         * @details Fills the mesh_->vertices, computes the vertex mapping and
         *         delete colocated vertices
         */
        void initialize() const;

    protected:
        GeoModelMeshVerticesBase( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< PointSetMesh< DIMENSION > >& mesh );

        /*!
         *@brief return the number of all vertices
         * it is computed summing all entities.nb().
         *@note colocated vertices are counted twice or more.
         */
        virtual index_t nb_total_vertices() const;
        virtual index_t fill_vertices() const;
        void fill_vertices_for_entity_type(
            const GeoModel< DIMENSION >& geomodel,
            const MeshEntityType& entity_type,
            index_t& count ) const;

    protected:
        /// Attached Mesh
        std::unique_ptr< PointSetMesh< DIMENSION > >& mesh_;

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };

    ALIAS_2D_AND_3D( GeoModelMeshVerticesBase );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshVertices final
        : public GeoModelMeshVerticesBase< DIMENSION >
    {
    public:
        GeoModelMeshVertices( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< PointSetMesh< DIMENSION > >& mesh );
    };

    template <>
    class geomodel_core_api GeoModelMeshVertices< 3 > final
        : public GeoModelMeshVerticesBase< 3 >
    {
    public:
        GeoModelMeshVertices( GeoModelMesh3D& gmm,
            GeoModel3D& gm,
            std::unique_ptr< PointSetMesh3D >& mesh );

        void clear() const;
        index_t nb_total_vertices() const override;
        index_t fill_vertices() const override;
    };

    ALIAS_2D_AND_3D( GeoModelMeshVertices );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshPolygonsBase
        : public GeoModelMeshCommon< DIMENSION >
    {
        ringmesh_disable_copy_and_move( GeoModelMeshPolygonsBase );
        static const std::string surface_att_name;
        static const std::string polygon_surface_att_name;

    public:
        friend class GeoModelMeshBase< DIMENSION >;
        friend class GeoModelMesh< DIMENSION >;

        virtual ~GeoModelMeshPolygonsBase();

        GEO::AttributesManager& attribute_manager() const;

        /*!
         * Test if the mesh polygons are initialized
         */
        void test_and_initialize() const;

        /*!
         * @brief Number of polygons stored.
         */
        index_t nb() const;

        /*!
         * Get the number of vertices in the polygon
         * @param[in] polygon the polygon index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t polygon ) const;
        /*!
         * Get the vertex index of a vertex in a polygon
         * in the GeoModelMesh
         * @param[in] p the polygon index
         * @param[in] vertex the local vertex index [0, nb_vertices_in_polygon[
         * @return the vertex index
         */
        index_t vertex( const ElementLocalVertex& polygon_local_vertex ) const;
        /*!
         * Get the adjacent polygon index in the GeoModelMesh
         * @param[in] polygon_local_edge in the polygon
         * @return the adjacent polygon index
         */
        index_t adjacent( const PolygonLocalEdge& polygon_local_edge ) const;
        /*!
         * Get the surface index in the GeoModel according the polygon
         * index in the GeoModelMesh
         * @param[in] polygon the polygon index
         * @return the surface index
         */
        index_t surface( index_t polygon ) const;
        /*!
         * Get the polygon index in the GeoModelMesh restricted to
         * the surface owing the polygon
         * @param[in] p the polygon index
         * @return the polygon index varying from 0 to nb_polygons
         * in the surface owing \p p
         */
        index_t index_in_surface( index_t polygon ) const;
        /*!
         * Get the polygon index in the GeoModelMesh restricted to
         * the surface owing the polygon and its type
         * @param[in] polygon the polygon index
         * @return the type of the polygon \p p
         * and the polygon index varying from 0 to nb_polygons
         * of the corresponding type of \p p in the owing surface.
         */
        std::tuple< PolygonType, index_t > type( index_t polygon ) const;

        /*!
         * Get the number of polygons of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of polygons
         */
        index_t nb_polygons( PolygonType type = PolygonType::UNDEFINED ) const;
        /*!
         * Get the number of polygons of the corresponding type
         * in the given surface of the GeoModel
         * @param[in] surface the surface index
         * @param[in] type the corresponding type
         * @return the number of polygons
         */
        index_t nb_polygons(
            index_t surface, PolygonType type = PolygonType::UNDEFINED ) const;
        /*!
         * Get the polygon index in the GeoModelMesh
         * @param[in] surface the surface index owing the polygon
         * @param[in] polygon the polygon index varying from 0 to the number of
         * polygons
         * of type \p type in the surface \p surface.
         * @warning \p polygon is NOT a polygon id
         * of the surface \p surface.
         * It is pth polygon of type \p type in the internal storage of the
         * GeoModelMeshPolygons (see
         * GeoModelMeshPolygons::surface_polygon_ptr_).
         * @note to find the polygon id of the GeoModelMeshPolygons from a
         * surface
         * and a polygon id of this surface, you need to perform a search using
         * NNSearch and the barycenter of the polygon for instance.
         * @param[in] type it can specify the polygon type used. For example, if
         * type = QUAD
         * then \p polygon represents the fth quad in the surface \p surface and
         * \p polygon can vary from 0
         * to nb_quads( s ).
         * If \p type is PolygonType::ALL, all the polygon types are
         * taken into account.
         * @return the polygon index
         */
        index_t polygon( index_t surface,
            index_t polygon,
            PolygonType type = PolygonType::UNDEFINED ) const;

        /*!
         * Get the number of triangles in the GeoModelMesh
         * @return the number of triangles
         */
        index_t nb_triangle() const;
        /*!
         * Get the number of triangles in the given surface
         * @param[in] surface the surface index
         * @return the number of triangles
         */
        index_t nb_triangle( index_t surface ) const;
        /*!
         * Get the polygon index in the GeoModelMesh corresponding
         * to the asked triangle in the surface
         * @param[in] surface the surface index
         * @param[in] triangle the triangleth triangle index varying from 0 to
         * nb_triangles( surface )
         * @return the polygon index
         */
        index_t triangle( index_t surface, index_t triangle ) const;

        /*!
         * Get the number of quads in the GeoModelMesh
         * @return the number of quads
         */
        index_t nb_quad() const;
        /*!
         * Get the number of quads in the given surface
         * @param[in] surface the surface index
         * @return the number of quads
         */
        index_t nb_quad( index_t surface ) const;
        /*!
         * Get the polygon index in the GeoModelMesh corresponding
         * to the asked quad in the surface
         * @param[in] surface the surface index
         * @param[in] quad the quadth quad index varying from 0 to nb_quads(
         * surface )
         * @return the polygon index
         */
        index_t quad( index_t surface, index_t quad ) const;

        /*!
         * Get the number of unclassified polygons in the GeoModelMesh
         * @return the number of unclassified polygons
         */
        index_t nb_unclassified_polygon() const;
        /*!
         * Get the number of polygons in the given surface
         * @param[in] surface the surface index
         * @return the number of polygons
         */
        index_t nb_unclassified_polygon( index_t surface ) const;
        /*!
         * Get the polygon index in the GeoModelMesh corresponding
         * to the asked polygon in the surface
         * @param[in] surface the surface index
         * @param[in] polygon the polygonth polygon index varying from 0 to
         * nb_polygons( s )
         * @return the polygon index
         */
        index_t unclassified_polygon( index_t surface, index_t polygon ) const;

        /*!
         * Clear the polygons of the GeoModelMesh
         */
        void clear();

        /*!
         * Get the center of the given polygon
         * @param[in] polygon the polygon index
         */
        vecn< DIMENSION > center( index_t polygon ) const;
        /*!
         * Get the area of the polygon
         * @param[in] p the polygon index
         */
        double area( index_t polygon ) const;

        const NNSearch< DIMENSION >& nn_search() const;

        /*!
         * @brief return the AABB tree for the polygons of the mesh
         */
        const SurfaceAABBTree< DIMENSION >& aabb() const;

    private:
        /*!
         * Initialize the polygons of the GeoModelMesh
         * and sort them per surface and polygon type
         * Example for a mesh with two surfaces and only triangles and quads
         * [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
         * |          surface 0           |             surface 1           |
         */
        void initialize();
        /*!
         * Resize edge data: surface_id_ and polygon_id_
         */
        void resize_polygon_data( index_t nb_polygons );
        /*!
         * Clear edge data: surface_id_ and polygon_id_
         */
        void clear_polygon_data();

        /*!
         * @brief Removes polygon adjacencies along lines
         */
        void disconnect_along_lines();

        /*!
         * @brief Sorts the polygons by surface and type
         *  Permute polygons to sort them per surface and per type
         * Example for a mesh with two surfaces and only triangles and quads
         * [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
         * |          surface 0           |             surface 1           |
         */
        void sort_polygons();

    protected:
        GeoModelMeshPolygonsBase( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh );

    protected:
        /// Attached Mesh
        std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh_;

        /// Vector storing the surface index per polygon
        std::vector< index_t > surface_id_;
        /// Vector storing the polygon index in surface per polygon
        std::vector< index_t > polygon_id_;

        /*!
         * Vector storing the index of the starting polygon index
         * for a given surface and a given polygon type.
         * For example:
         *    the 2nd quad index of the surface index S will be found here:
         *    surface_polygon_ptr_[ALL*S + QUAD] + 2
         */
        std::vector< index_t > surface_polygon_ptr_;

        /// Number of triangles in the GeoModelMesh
        index_t nb_triangles_{ 0 };
        /// Number of quads in the GeoModelMesh
        index_t nb_quads_{ 0 };
        /// Number of unclassified polygons in the GeoModelMesh
        index_t nb_unclassified_polygons_{ 0 };
    };

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshPolygons final
        : public GeoModelMeshPolygonsBase< DIMENSION >
    {
    public:
        GeoModelMeshPolygons( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh );
    };

    template <>
    class geomodel_core_api GeoModelMeshPolygons< 3 > final
        : public GeoModelMeshPolygonsBase< 3 >
    {
    public:
        GeoModelMeshPolygons( GeoModelMesh3D& gmm,
            GeoModel3D& gm,
            std::unique_ptr< SurfaceMesh3D >& mesh );

        /*!
         * Get the normal of the polygon
         * @param[in] p the polygon index
         */
        vec3 normal( index_t p ) const;
    };

    ALIAS_2D_AND_3D( GeoModelMeshPolygons );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshEdges final
        : public GeoModelMeshCommon< DIMENSION >
    {
        ringmesh_disable_copy_and_move( GeoModelMeshEdges );
        static const std::string line_att_name;
        static const std::string edge_line_att_name;

    public:
        friend class GeoModelMeshBase< DIMENSION >;
        friend class GeoModelMesh< DIMENSION >;

        virtual ~GeoModelMeshEdges();

        GEO::AttributesManager& attribute_manager() const;

        /*!
         * Test if the mesh edges are initialized
         */
        void test_and_initialize() const;

        /*!
         * @brief Number of edges stored.
         */
        index_t nb() const;

        /*!
         * Get the vertex index of a vertex in a edge
         * in the GeoModelMesh
         * @param[in] edge the edge index
         * @param[in] vertex the local vertex index
         * @return the vertex index
         */
        index_t vertex( const ElementLocalVertex& edge_local_vertex ) const;

        /*!
         * Get the line index in the GeoModel according the edge
         * index in the GeoModelMesh
         * @param[in] edge the edge index
         * @return the line index
         */
        index_t line( index_t edge ) const;
        /*!
         * Get the edge index in the GeoModelMesh restricted to
         * the line owing the edge
         * @param[in] edge the edge index
         * @return the edge index varying from 0 to nb_edge()
         * in the line owing \p e
         */
        index_t index_in_line( index_t edge ) const;

        /*!
         * Get the number of edges in the given line
         * @param[in] edge the edge index
         * @return the number of edges
         */
        index_t nb_edges( index_t line ) const;
        /*!
         * Get the edge index in the GeoModelMesh corresponding
         * to the asked edge in the line
         * @param[in] l the line index
         * @param[in] edge the eth edge index varying from 0 to nb_edges( l )
         * @return the edge index
         */
        index_t edge( index_t line, index_t edge ) const;
        /*!
         * Clear the edges of the GeoModelMesh
         */
        void clear();

        /*!
         * Get the center of the given edge
         * @param[in] p the edge index
         */
        vecn< DIMENSION > center( index_t edge ) const;
        /*!
         * Get the length of the edge
         * @param[in] p the edge index
         */
        double length( index_t edge ) const;

        const NNSearch< DIMENSION >& nn_search() const;

        /*!
         * @brief return the AABB tree for the edges of the mesh
         */
        const LineAABBTree< DIMENSION >& aabb() const;

    private:
        /*!
         * Initialize the edges of the GeoModelMesh
         * and sort them per line
         * Example for a mesh with two lines
         * [EDGE, EDGE, ... ,EDGE , EDGE, EDGE, ... , EDGE]
         * |        line 0        |        line 1         |
         */
        void initialize();
        /*!
         * Resize edge data: line_id_ and edge_id_
         */
        void resize_edge_data();
        /*!
         * Clear edge data: line_id_ and edge_id_
         */
        void clear_edge_data();

    protected:
        GeoModelMeshEdges( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< LineMesh< DIMENSION > >& mesh );

    protected:
        /// Attached Mesh
        std::unique_ptr< LineMesh< DIMENSION > >& mesh_;

        /// Vector storing the line index per edge
        std::vector< index_t > line_id_;
        /// Vector storing the edge index in line per edge
        std::vector< index_t > edge_id_;

        /*!
         * Vector storing the index of the starting edge index
         * for a given line and a given edge type.
         * For example:
         *    the 2nd edge index of the line index L will be found here:
         *    line_edge_ptr_[L] + 2
         */
        std::vector< index_t > line_edge_ptr_;

        /// Number of edges in the GeoModelMesh
        index_t nb_edges_{ 0 };
    };

    ALIAS_2D_AND_3D( GeoModelMeshEdges );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshWells final
        : public GeoModelMeshCommon< DIMENSION >
    {
    public:
        explicit GeoModelMeshWells( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< LineMesh< DIMENSION > >& mesh );

        GEO::AttributesManager& attribute_manager() const;

        /*!
         * Tests if the mesh edges needs to be initialized and initialize it
         */
        void test_and_initialize() const;

        /*!
         * Gets the number of wells
         * @return the corresponding number
         */
        index_t nb_wells() const;
        /*!
         * Gets the number of edges in the MacroMesh
         * @return the corresponding number
         */
        index_t nb_edges() const;
        /*!
         * Gets the number of edges of a Well
         * @param[in] w the well index
         * @return the corresponding number
         */
        index_t nb_edges( index_t well ) const;
        /*!
         * Gets the vertex index of the GeoModelMesh
         * @param[in] w the well index
         * @param[in] edge the edge index in the well (from 0 to nb_edges in the
         * well)
         * @param[in] vertex the vertex index of the edge (0 or 1 )
         * @return the global vertex index
         */
        index_t vertex( index_t well, index_t edge, index_t vertex ) const;
        /*!
         * Clear the mesh edges
         */
        void clear();

        /*!
         * Initialize the mesh edges
         */
        void initialize();

        /*!
         * @brief return the AABB tree for the edges of the mesh
         */
        const LineAABBTree< DIMENSION >& aabb() const;

    private:
        /// Attached Mesh
        std::unique_ptr< LineMesh< DIMENSION > >& mesh_;

        /*!
         * Vector storing the index of the starting edge index
         * for a given well
         */
        std::vector< index_t > well_ptr_;
    };

    ALIAS_2D_AND_3D( GeoModelMeshWells );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshCells final
        : public GeoModelMeshCommon< DIMENSION >
    {
        static const std::string region_att_name;
        static const std::string cell_region_att_name;

    public:
        friend class GeoModelMeshBase< DIMENSION >;
        friend class GeoModelMesh< DIMENSION >;

        /*!
         * Several modes for vertex duplication algorithm:
         *  - NONE = no duplication
         *  - FAULT = duplication along faults
         *  - HORIZON = duplication along horizons
         *  - ALL = duplication along faults and horizons
         */
        enum DuplicateMode
        {
            NONE,
            FAULT,
            HORIZON,
            ALL
        };

        GeoModelMeshCells( GeoModelMesh< DIMENSION >& gmm,
            GeoModel< DIMENSION >& gm,
            std::unique_ptr< VolumeMesh< DIMENSION > >& mesh );

        GEO::AttributesManager& attribute_manager() const;

        /*!
         * Test if the mesh cells are duplicated
         */
        bool is_duplication_initialized() const;

        /*!
         * Test if the mesh cells need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const;

        /*!
         * @brief Number of cells stored.
         */
        index_t nb() const;
        /*!
         * Gets the number of duplicated points by the DuplicateMode algorithm
         * @return the corresponding number of duplications
         */
        index_t nb_duplicated_vertices() const;
        /*!
         * Gets the total number of vertices (mesh.vertices.nb() +
         * nb_duplicated_vertices())
         * @return the corresponding number of vertices
         */
        index_t nb_total_vertices() const;
        /*!
         * Check if the corner in a cell is duplicated,
         * if so give the duplicated vertex index
         * @param[in] c the cell index in the GeoModelMesh
         * @param[in] vertex the local vertex index in the cell \p c (0 to
         * nb_vertices( c ))
         * @return the duplicated vertex index (0 to nb_duplicated_vertices())
         * if the corner is duplicated, else NO_ID.
         */
        index_t duplicated_corner_index(
            const ElementLocalVertex& cell_local_vertex ) const;
        /*!
         * Get the vertex index in the GeoModelMesh corresponding
         * to the given duplicated vertex index
         * @param[in] duplicate_vertex_index the duplicated vertex index
         * @return the vertex index
         */
        index_t duplicated_vertex( index_t duplicate_vertex_index ) const;

        /*!
         * Get the number of vertices in the cell
         * @param[in] c the cell index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t cell ) const;
        /*!
         * Get the vertex index of a vertex in a cell
         * in the GeoModelMesh
         * @param[in] c the cell index
         * @param[in] vertex the local vertex index [0, nb_vertices_in_cell[
         * @return the vertex index
         */
        index_t vertex( const ElementLocalVertex& cell_local_vertex ) const;
        /*!
         * Get the number of edges in the cell
         * @param[in] cell the cell index
         */
        index_t nb_edges( index_t cell ) const;
        /*!
         * Get the number of facets in the cell
         * @param[in] cell the cell index
         */
        index_t nb_facets( index_t cell ) const;
        /*!
         * Get the number of facets in the cell
         */
        index_t nb_facet_vertices(
            const CellLocalFacet& cell_local_facet ) const;
        /*!
         * \brief Gets a cell vertex by local facet index and local
         *  vertex index in the edge
         * \param[in] cell_local_facet the local facet in a cell
         * \param[in] local_vertex the local index in the cell facet
         * \return vertex \p lv of facet \p lf in cell \p c
         */
        index_t facet_vertex( const CellLocalFacet& cell_local_facet,
            index_t local_vertex ) const;
        /*!
         * \brief Gets a cell vertex by local edge index and local
         *  vertex index in the edge
         * \param[in] cell the cell, in 0..nb()-1
         * \param[in] local_edge the local edge index, in 0..nb_edges(c)-1
         * \param[in] local_vertex the local index in the edge, one of 0,1
         * \return vertex \p local_vertex of edge \p local_edge in cell \p c
         */
        index_t edge_vertex(
            index_t cell, index_t local_edge, index_t local_vertex ) const;
        /*!
         * Get the adjacent cell index in the GeoModelMesh
         * @param[in] cell the cell index
         * @param[in] facet the facet index
         * @return the adjacent cell index
         */
        index_t adjacent( index_t cell, index_t facet ) const;
        /*!
         * Get the region index in the GeoModel according the cell
         * index in the GeoModelMesh
         * @param[in] cell the cell index
         * @return the region index
         */
        index_t region( index_t cell ) const;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell
         * @param[in] cell the cell index
         * @return the cell index varying from 0 to nb_cells
         * in the region owing \p cell
         */
        index_t index_in_region( index_t cell ) const;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell and its type
         * @param[in] cell the cell index
         * @param[out] index the cell index varying from 0 to nb_cells
         * of the corresponding type of \p c in the owing region
         * @return the type of the cell \p c
         */
        CellType type( index_t cell ) const;

        /*!
         * Get the number of cells of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells( CellType type = CellType::UNDEFINED ) const;
        /*!
         * Get the number of cells of the corresponding type
         * in the given region of the GeoModel
         * @param[in] region the region index
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells(
            index_t region, CellType type = CellType::UNDEFINED ) const;
        /*!
         * Get the cell index in the GeoModelMesh
         * @param[in] region the region index owing the cell
         * @param[in] cell the cell index varying from 0 to number of cells
         * of type \p type in the region \p region.
         * @warning \p cell is NOT a cell id of the region \p region,
         * It is cellth cell of type \p type in the internal storage of the
         * GeoModelMeshCells (see GeoModelMeshCells::region_cell_ptr_).
         * @note to find the cell id of the GeoModelMeshCells from a region
         * and a cell id of this region, you need to perform a search using
         * NNSearch and the barycenter of the cell for instance.
         * @param[in] type it can specify the cell type used. For example,
         * if type = HEXAEDRON then \p cell represents the cellth hex in the
         * region \p region and \p cell can vary from 0 to nb_hex( r ).
         * If \p type is CellType::UNDEFINED, all the cell types are
         * taken into account.
         * @return the cell index
         */
        index_t cell( index_t region,
            index_t cell,
            CellType type = CellType::UNDEFINED ) const;

        /*!
         * Get the number of tets in the GeoModelMesh
         * @return the number of tets
         */
        index_t nb_tet() const;
        /*!
         * Get the number of tets in the given region
         * @param[in] r the region index
         * @return the number of tets
         */
        index_t nb_tet( index_t region ) const;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked tet in the region
         * @param[in] region the region index
         * @param[in] tet the tetth tet index varying from 0 to nb_tet( region )
         * @return the cell index
         */
        index_t tet( index_t region, index_t tet ) const;

        /*!
         * Get the number of hexs in the GeoModelMesh
         * @return the number of hexs
         */
        index_t nb_hex() const;
        /*!
         * Get the number of hexs in the given region
         * @param[in] region the region index
         * @return the number of hexs
         */
        index_t nb_hex( index_t region ) const;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked hex in the region
         * @param[in] region the region index
         * @param[in] hex the hexth hex index varying from 0 to nb_hex( region )
         * @return the cell index
         */
        index_t hex( index_t region, index_t hex ) const;

        /*!
         * Get the number of prisms in the GeoModelMesh
         * @return the number of prisms
         */
        index_t nb_prism() const;
        /*!
         * Get the number of prisms in the given region
         * @param[in] region the region index
         * @return the number of prisms
         */
        index_t nb_prism( index_t region ) const;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked prism in the region
         * @param[in] region the region index
         * @param[in] prism the prismth prism index varying from 0 to nb_prism(
         * region )
         * @return the cell index
         */
        index_t prism( index_t region, index_t prism ) const;

        /*!
         * Get the number of pyramids in the GeoModelMesh
         * @return the number of pyramids
         */
        index_t nb_pyramid() const;
        /*!
         * Get the number of pyramids in the given region
         * @param[in] region the region index
         * @return the number of pyramids
         */
        index_t nb_pyramid( index_t region ) const;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked pyramid in the region
         * @param[in] r the region index
         * @param[in] p the pth pyramid index varying from 0 to nb_pyramid( r )
         * @return the cell index
         */
        index_t pyramid( index_t region, index_t pyramid ) const;

        /*!
         * Get the number of connectors in the GeoModelMesh
         * @return the number of connectors
         */
        index_t nb_connector() const;
        /*!
         * Get the number of connectors in the given region
         * @param[in] region the region index
         * @return the number of connectors
         */
        index_t nb_connector( index_t region ) const;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked connector in the region
         * @param[in] region the region index
         * @param[in] connector the connectorth connector index varying from 0
         * to nb_connector( r )
         * @return the cell index
         */
        index_t connector( index_t region, index_t connector ) const;

        /*!
         * Clear the mesh cells
         */
        void clear();
        /*!
         * Remove the duplication of the mesh cell facets
         */
        void clear_duplication();

        /*!
         * Determine if a cell facet is on a surface. If so, fill the \p action
         * with the surface index and the surface side encountered
         * @param[in] cell the cell index
         * @param[in] facet the facet index
         * @param[out] colocated_facet_index the facet index colocalised with
         * the cell facet
         * @param[out] side the side of the facet \p facet.
         * true = side of the facet normal, false = the other side
         * @return true is the cell facet is on a surface
         */
        bool is_cell_facet_on_surface( index_t cell,
            index_t facet_index,
            index_t& colocated_facet_index,
            bool& side ) const;

        /*!
         * Get the center of the given cell
         * @param[in] cell the cell index
         */
        vec3 barycenter( index_t cell ) const;
        /*!
         * Get the volume of the cell
         * @param[in] cell the cell index
         */
        double volume( index_t cell ) const;

        const NNSearch< DIMENSION >& cell_nn_search() const;

        const NNSearch< DIMENSION >& cell_facet_nn_search() const;

        /*!
         * @brief return the AABB tree for the cells of the mesh
         */
        const VolumeAABBTree< DIMENSION >& aabb() const;

    private:
        /// enum to characterize the action to do concerning a surface
        /// Action concerns the vertices of a Surface and not the Surface
        enum ActionOnSurface
        {
            /// do nothing
            SKIP = -2,
            /// need to be duplicated (don't know which side yet)
            TO_PROCESS = -1,
            /// need to duplicate the side opposite to the facet normal
            NEG_SIDE = 0,
            /// need to duplicate the side following the facet normal
            POS_SIDE = 1
        };
        /// Action to do according a surface index
        using action_on_surface = std::pair< index_t, ActionOnSurface >;

        /*!
         * @brief Initialize the  cells from the cells
         *        of the GeoModel Region cells
         */
        void initialize();

        /*!
         * Resize region_id and cell_id
         */
        void resize_cell_data();
        /*!
         * Clear cell_id; region_id and polygon id_
         */
        void clear_cell_data();

        /*!
         * Test if the mesh cell are duplicated according
         * the duplication mode, if not duplicate them.
         */
        void test_and_initialize_duplication() const;
        /*!
         * Duplicate the mesh cell along some surfaces defined
         * by the duplication mode
         */
        void initialize_duplication();
        /*!
         * Test if we need to duplicate mesh cell along the given
         * surface according the duplicate mode
         * @param[in] surface_id the surface index in the GeoModel
         */
        bool is_surface_to_duplicate( index_t surface_id ) const;

        /*!
         * Determine the actions to do according the action_on_surfaces
         * encountered during the propagation around a vertex (initialize())
         * @param[in] surfaces the action_on_surfaces encountered
         * @param[in,out] info the global information on what to do for each
         * surface.
         * This information is updated in this function according the
         * encountered action_on_surfaces
         * @return true if the corners should be duplicated
         */
        bool are_corners_to_duplicate(
            const std::vector< action_on_surface >& surfaces,
            std::vector< ActionOnSurface >& info );
        /*!
         * Test if the mesh cell facet vector is filled with
         * the colocalised facet. If not fill it.
         */
        void test_and_initialize_cell_facet() const;
        /*!
         * Initialize the mesh cell facet vector of colocalised facet.
         */
        void initialize_cell_facet();

        void sort_cells();

    private:
        /// Attached Mesh
        std::unique_ptr< VolumeMesh< DIMENSION > >& mesh_;

        /// Vector storing the region index per cell
        std::vector< index_t > region_id_;
        /// Vector storing the cell index in region per cell
        std::vector< index_t > cell_id_;

        /*!
         * Vector storing the index of the starting cell index
         * for a given region and a given cell type.
         * For example:
         *    the 2nd hex index of the region index R will be found here:
         *    surface_polygon_ptr_[CellType::UNDEFINED*R + CellType::HEXAEDRON]
         * + 2
         */
        std::vector< index_t > region_cell_ptr_;

        /// Number of tet in the GeoModelMesh
        index_t nb_tets_{ 0 };
        /// Number of hex in the GeoModelMesh
        index_t nb_hexs_{ 0 };
        /// Number of prism in the GeoModelMesh
        index_t nb_prisms_{ 0 };
        /// Number of pyramid in the GeoModelMesh
        index_t nb_pyramids_{ 0 };
        /// Number of connector in the GeoModelMesh
        index_t nb_connectors_{ 0 };

        /// Current duplicate mode applied on the mesh
        DuplicateMode mode_{ NONE };
        /*!
         * @brief Vector of duplicated vertices
         * @details Each value is a duplicated vertex, the index corresponds to
         * vertex index in mesh.vertices.
         */
        std::vector< index_t > duplicated_vertex_indices_;

        /*!
         * @brief Vector storing the colocalised polygon index per cell facet
         * @detail If a cell facet is on a surface, the vector is equal to
         * the index of the corresponding polygon.
         */
        std::vector< index_t > polygon_id_;
    };

    ALIAS_2D_AND_3D( GeoModelMeshCells );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshBase
    {
        ringmesh_disable_copy_and_move( GeoModelMeshBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelMeshBase();

        const GeoModel< DIMENSION >& geomodel() const
        {
            return geomodel_;
        }

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated_vertices();

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the geomodel_vertex_id in the GeoModelMeshEntity
         * of the GeoModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with which information should be merged.
         *  It is recycled to give the mapping between old and new vertex
         * indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete );

        /*!
         * @brief Remove all invalid GMEVertex and delete the vertices
         * that are not anymore in any GeoModelEntity
         */
        void erase_invalid_vertices();

        void change_point_set_mesh_data_structure( const MeshType& type );
        void change_line_mesh_data_structure( const MeshType& type );
        void change_surface_mesh_data_structure( const MeshType& type );

    protected:
        /*! Attached GeoModel */
        const GeoModel< DIMENSION >& geomodel_;

        /// Mesh storing all the elements
        std::unique_ptr< MeshSet< DIMENSION > > mesh_set_;

    protected:
        GeoModelMeshBase(
            GeoModelMesh< DIMENSION >& gmm, GeoModel< DIMENSION >& geomodel );

    public:
        GeoModelMeshVertices< DIMENSION > vertices;
        GeoModelMeshEdges< DIMENSION > edges;
        GeoModelMeshWells< DIMENSION > wells;
        GeoModelMeshPolygons< DIMENSION > polygons;
    };

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMesh final
        : public GeoModelMeshBase< DIMENSION >
    {
    public:
        explicit GeoModelMesh( GeoModel< DIMENSION >& geomodel );
    };

    template <>
    class geomodel_core_api GeoModelMesh< 3 > final
        : public GeoModelMeshBase< 3 >
    {
        ringmesh_disable_copy_and_move( GeoModelMesh );

    public:
        explicit GeoModelMesh( GeoModel3D& geomodel );
        virtual ~GeoModelMesh();

        /*! @}
         * \name Transfers of attributes
         * @{
         */
        void transfer_attributes_from_gmm_to_gm_regions() const;
        void transfer_attributes_from_gm_regions_to_gmm() const;
        void transfer_cell_attributes_from_gmm_to_gm_regions() const;
        void transfer_cell_attributes_from_gm_regions_to_gmm() const;
        void transfer_vertex_attributes_from_gmm_to_gm_regions() const;
        void transfer_vertex_attributes_from_gm_regions_to_gmm() const;

        /*! @}
         * \name Vertex duplication
         * @{
         */
        /*!
         * Access the DuplicateMode
         * @return the current DuplicateMode
         */
        GeoModelMeshCells3D::DuplicateMode duplicate_mode() const
        {
            return mode_;
        }
        /*!
         * Set a new DuplicateMode
         * @param[in] mode the new DuplicateMode for the GeoModelMesh
         */
        void set_duplicate_mode(
            const GeoModelMeshCells3D::DuplicateMode& mode ) const
        {
            if( mode_ == mode )
                return;
            mode_ = mode;
            const_cast< GeoModelMesh3D* >( this )->cells.clear_duplication();
        }

        /*! @}
         * \name Data structure change
         * @{
         */
        void change_volume_mesh_data_structure( const MeshType& type );

    private:
        /// Optional duplication mode to compute the duplication of cells on
        /// surfaces
        mutable GeoModelMeshCells3D::DuplicateMode mode_{
            GeoModelMeshCells3D::NONE
        };

    public:
        GeoModelMeshCells3D cells;
    };

    ALIAS_2D_AND_3D( GeoModelMesh );

} // namespace RINGMesh
