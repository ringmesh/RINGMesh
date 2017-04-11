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

#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/geomodel/entity_type_manager.h>

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>

/*!
 * @file ringmesh/geomodel_mesh.h
 * @brief Classes to manage globally the indexing of mesh entities of a GeoModel
 * @author Arnaud Botella and Jeanne Pellerin
 */

/*! @todo URGENT : Rename all parameters of all the function in this file. 
 *                 A lot of comments could then be removed.  
 */

namespace RINGMesh {
    class GeoModelMesh ;
    class GeoModelMeshVertices ;
    class GeoModel ;
    class GeoModelEntity ;
    class GeoModelMeshEntity ;
}

namespace RINGMesh {

    /*! @todo Move this global variables in a function */
    static const std::string surface_att_name = "region" ;
    static const std::string region_att_name = "region" ;
    static const std::string cell_region_att_name = "cell_region" ;
    static const std::string facet_surface_att_name = "facet_surface" ;

    class RINGMESH_API GeoModelMeshBase {
    protected:
        GeoModelMeshBase( GeoModelMesh& gmm, GeoModel& gm ) ;

        void set_mesh( MeshBase* mesh )
        {
            mesh_base_ = mesh ;
        }

        void save_mesh( const std::string& filename ) const
        {
            mesh_base_->save_mesh( filename ) ;
        }
    protected:
        /// Attached GeoModelMesh
        GeoModelMesh& gmm_ ;
        /// Attached GeoModel
        GeoModel& gm_ ;
        /// Attached MeshBase
        MeshBase* mesh_base_ ;
    } ;

    class RINGMESH_API GeoModelMeshVertices: public GeoModelMeshBase {
    ringmesh_disable_copy( GeoModelMeshVertices ) ;
    public:
        friend class GeoModelMeshEdges ;
        friend class GeoModelMeshFacets ;
        friend class GeoModelMeshCells ;

        GeoModelMeshVertices( GeoModelMesh& gmm, GeoModel& gm ) ;

        ~GeoModelMeshVertices() ;

        GEO::AttributesManager& attribute_manager() const
        {
            return mesh_->vertex_attribute_manager() ;
        }

        /*!
         * Tests if the mesh vertices are initialized
         */
        bool is_initialized() const ;
        /*!
         * Tests if the mesh vertices need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;

        /*!
         * @brief Number of vertices stored.
         */
        index_t nb() const ;

        /*!
         * @brief Coordinates of a vertex of the GeoModel
         * @pre v < nb()
         */
        const vec3& vertex( index_t v ) const ;

        /*!
         * @brief Returns the index of the given vertex in the geomodel
         * @param[in] p input point coordinates
         * @return index of the vertex in the geomodel if found
         * (distance < epsilon), otherwise NO_ID
         */
        index_t index( const vec3& p ) const ;

        /*!
         * @brief Get the GeoModelMesh index of a GeoModelMeshEntity vertex from its
         * index in that GeoModelMeshEntity
         * @param[in] mesh_entity GeoModelMeshEntity that vertex belongs to
         * @param[in] entity_vertex_index index of the query vertex in the GeoModelMeshEntity
         * @return if found the vertex index in the geomodel, else NO_ID.
         */
        index_t geomodel_vertex_id(
            const gmme_id& mesh_entity,
            index_t entity_vertex_index = 0 ) const ;

        /*!
         * @brief Get the GeoModelMesh index of a GeoModelMeshEntity vertex from its
         * index in that GeoModelMeshEntity
         * @param[in] mesh_entity GeoModelMeshEntity that vertex belongs to
         * @param[in] entity_mesh_element_index index of the mesh element that vertex belongs to
         * @param[in] vertex_local_index local index of the query vertex in the mesh element
         * @return if found the vertex index in the geomodel, else NO_ID.
         */
        index_t geomodel_vertex_id(
            const gmme_id& mesh_entity,
            index_t entity_mesh_element_index,
            index_t vertex_local_index ) const ;

        /*!
         * @brief Get the GeoModelMeshEntity vertices from its index in the GeoModelMesh
         * @param[in] mesh_entity Unique id to a GeoModelMeshEntity
         * @param[in] geomodel_vertex_id Index of the query vertex in the geomodel
         * @param[out] mesh_entity_vertex_ids Corresponding GeoModelMeshEntity vertices
         */
        void mesh_entity_vertex_id(
            const gmme_id& mesh_entity,
            index_t geomodel_vertex_id,
            std::vector< index_t >& mesh_entity_vertex_ids ) const ;

        /*!
         * @brief Get the vertices in GeoModelEntity corresponding to the given unique vertex
         * @param[in] v Vertex index in the geomodel
         * @param[out] gme_vertices Result
         */
        void gme_vertices(
            index_t v,
            std::vector< GMEVertex >& gme_vertices ) const ;

        /*!
         * @brief Get the vertex indices in the specified MeshEntity type
         * corresponding to the given unique vertex
         */
        void gme_type_vertices(
            const MeshEntityType& entity_type,
            index_t v,
            std::vector< GMEVertex >& gme_vertices ) const ;

        /*!
         * @brief To use when building the geomodel by first adding its vertices
         * @return the first index of created vertices
         * @warning The client is responsible for setting the mapping between the points
         * of the GeoModelEntity and the unique vertex
         */
        index_t add_vertex( const vec3& point ) ;
        index_t add_vertices( const std::vector< vec3 >& points ) ;

        /*!
         * @brief Set the point coordinates of all the vertices that
         *        share this unique vertex, including the unique vertex itself.
         * @param[in] v Index of the vertex
         * @param[in] point New coordinates
         */
        void update_point( index_t v, const vec3& point ) ;

        void update_vertex_mapping(
            const gmme_id& entity_id,
            index_t entity_vertex_index,
            index_t geomodel_vertex_index ) ;

        /*!
         * @brief Clear the vertices - clear the gme_vertices_ -
         *        clear global vertex information in the all BMME
         * @warning Not stable - crashes if attributes are still bound
         */
        void clear() ;

        void unbind_geomodel_vertex_map( const gmme_id& mesh_entity_id ) ;

        void bind_geomodel_vertex_map( const gmme_id& mesh_entity_id ) ;

        const NNSearch& nn_search() const
        {
            test_and_initialize() ;
            return mesh_->vertices_nn_search() ;
        }

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated() ;

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the geomodel_vertx_id in the GeoModelMeshEntity
         * of the GeoModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with which information should be merged.
         *  It is recycled to give the mapping between old and new vertex indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) ;

    private:
        void fill_vertices(
            const GeoModel& M,
            const MeshEntityType& entity_type,
            index_t& count ) ;

        /*!
         * @brief Initialize the vertices from the vertices
         *        of the GeoModel Corners, Lines, Surfaces and Regions
         * @details Fills the mesh_.vertices, computes the vertex mapping and
         *         delete colocated vertices
         */
        void initialize() ;

    private:
        /*!
         * Class which manages the mapping informations between vertices
         * of GeoModelMeshEntites (entity_index) and GeoModelMeshVertices (global index)
         */
        class RINGMESH_API GeoModelVertexMapper {
        ringmesh_disable_copy( GeoModelVertexMapper ) ;
        public:
            GeoModelVertexMapper(
                GeoModelMeshVertices& geomodel_vertices,
                const GeoModel& geomodel ) ;

            /*!
             * \name Query
             * @{
             */

            /*!
             * @brief Returns the index of a GeoModelMeshEntity vertex in the geomodel
             * global indexing
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             * @param[in] mesh_entity_vertex_index Index of query vertex in the
             * GeoModelMeshEntity indexing
             * @returns Model index of the GeoModelMeshEntity vertex
             */
            index_t geomodel_vertex_index(
                const gmme_id& mesh_entity_id,
                index_t mesh_entity_vertex_index ) const ;

            /*!
             * @brief Returns all the corresponding vertices in GeoModelMeshEntities
             * to a given geomodel vertex
             * @param[in] v Model vertex index
             * @returns All the corresponding vertices in their local indexing
             */
            const std::vector< GMEVertex >& mesh_entity_vertex_indices(
                index_t v ) const ;

            /*!
             * @brief Returns all the corresponding vertices in GeoModelMeshEntities
             * of a specific type to a given geomodel vertex
             * @param[in] v Model vertex index
             * @param[in] mesh_entity_type Type of GeoModelMeshEntity
             * @param[out] result corresponding vertices in GeoModelMeshEntities
             * of a specific type
             */
            void mesh_entity_vertex_indices(
                index_t v,
                const MeshEntityType& mesh_entity_type,
                std::vector< GMEVertex >& result ) const ;

            /*!
             * @brief Returns all the corresponding vertices to a geomodel vertex
             * in a specific GeoModelMeshEntities
             * @param[in] v Model vertex index
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             * @param[out] result corresponding vertices in the GeoModelMeshEntity
             * @returns All the corresponding vertices in their local indexing
             */
            void mesh_entity_vertex_indices(
                index_t v,
                const gmme_id& mesh_entity_id,
                std::vector< index_t >& result ) const ;

            const GEO::Attribute< index_t >& vertex_map(
                const gmme_id& mesh_entity_id ) const ;

            GEO::Attribute< index_t >& vertex_map( const gmme_id& mesh_entity_id ) ;

            /*! @}
             * \name Updating
             * @{
             */

            /*!
             * @brief Sets the geomodel vertex mapping value of a given vertex
             * in a GeoModelMeshEntity
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             * @param[in] mesh_entity_vertex_index Index of query vertex in the
             * GeoModelMeshEntity indexing
             * @param[in] geomodel_entity_vertex_index Model vertex index to map with
             */
            void set_vertex_map_value(
                const gmme_id& mesh_entity_id,
                index_t mesh_entity_vertex_index,
                index_t geomodel_entity_vertex_index ) ;

            void add_to_gme_vertices(
                const GMEVertex& gme_vertex,
                index_t geomodel_vertex_index ) ;

            /*!
             * @brief Updates all the vertex maps with regards to the global indexing
             * changes
             * @param[in] old2new Map between actual geomodel indexing and wanted
             * geomodel indexing. Its size is equal to the number of geomodel vertices.
             */
            void update_mesh_entity_maps_and_gmes(
                const std::vector< index_t >& old2new ) ;

            /*! @}
             * \name Initialization
             * @{
             */

            /*!
             * @brief Resizes the GME_Vertex vectors
             * @param[in] nb Size of the vector
             */
            void resize_geomodel_vertex_gmes( const index_t nb )
            {
                gme_vertices_.resize( nb ) ;
            }

            /*!
             * @brief Clears and resizes the GME_Vertex vectors
             * @param[in] nb Size of the vector
             */
            void clear_and_resize_geomodel_vertex_gmes( const index_t nb )
            {
                gme_vertices_.clear() ;
                resize_geomodel_vertex_gmes( nb ) ;
            }

            void bind_all_mesh_entity_vertex_maps() ;

            /*! @}
             * \name Clearing
             * @{
             */

            /*!
             * @brief Clears all the information about vertex mapping (attribute maps
             * and vectors of GME_Vertices
             */
            void clear() ;

            /*!
             * @brief Clears the GME_Vertices about one geomodel vertex
             */
            void clear_geomodel_vertex_gmes( index_t v )
            {
                ringmesh_assert( v < gme_vertices_.size() ) ;
                gme_vertices_[v].clear() ;
            }

            void unbind_vertex_map( const gmme_id& mesh_entity_id ) ;

            GEO::Attribute< index_t >& bind_vertex_map(
                const gmme_id& mesh_entity_id ) ;

            /*!
             * @}
             */

        private:
            /*!
             * @brief Initializes the given GeoModelMeshEntity vertex map
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             */
            void initialize_mesh_entity_vertex_map( const gmme_id& mesh_entity_id ) ;

            /*!
             * @brief Tests if the given GeoModelMeshEntity vertex map is initialized.
             * If not, initializes it.
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             * @return True is the map was initialized, false if not.
             */
            bool test_and_initialize_mesh_entity_vertex_map(
                const gmme_id& mesh_entity_id ) ;

            /*!
             * @brief Tests if the given GeoModelMeshEntity vertex map exists.
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             * @return True is the map exists, false if not.
             */
            bool is_mesh_entity_vertex_map_initialized(
                const gmme_id& mesh_entity_id ) const ;

            /*!
             * @brief Unbinds all the GeoModelMeshEntity vertex maps
             */
            void clear_all_mesh_entity_vertex_map() ;

            void resize_all_mesh_entity_vertex_maps() ;

            /*!
             * @brief Returns the vertex attribute of a GeoModelMeshEntity
             * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
             */
            GEO::AttributesManager& mesh_entity_vertex_attribute_manager(
                const gmme_id& mesh_entity_id ) const ;

        private:
            GeoModelMeshVertices& geomodel_vertices_ ;
            const GeoModel& geomodel_ ;

            /// Vertex maps
            AttributeVector< index_t > corner_vertex_maps_ ;
            AttributeVector< index_t > line_vertex_maps_ ;
            AttributeVector< index_t > surface_vertex_maps_ ;
            AttributeVector< index_t > region_vertex_maps_ ;
            std::map< MeshEntityType, AttributeVector< index_t >* > vertex_maps_ ;

            /// GeoModelEntity Vertices for each geomodel vertex
            std::vector< std::vector< GMEVertex > > gme_vertices_ ;
        } ;

    private:
        /// Attached Mesh
        std::unique_ptr< Mesh0D > mesh_ ;
        /// Mapper from/to GeoModelMeshEntity vertices
        GeoModelVertexMapper vertex_mapper_ ;
    } ;

    class RINGMESH_API GeoModelMeshFacets: public GeoModelMeshBase {
    ringmesh_disable_copy( GeoModelMeshFacets ) ;
    public:
        friend class GeoModelMesh ;

        enum FacetType {
            TRIANGLE, QUAD, POLYGON, ALL, NO_FACET
        } ;

        GeoModelMeshFacets( GeoModelMesh& gmm, GeoModel& gm ) ;
        ~GeoModelMeshFacets() ;

        GEO::AttributesManager& attribute_manager() const
        {
            return mesh_->facet_attribute_manager() ;
        }

        /*!
         * Test if the mesh facets are initialized
         */
        bool is_initialized() const ;
        void test_and_initialize() const ;

        /*!
         * @brief Number of facets stored.
         */
        index_t nb() const ;

        /*!
         * Get the number of vertices in the facet
         * @param[in] f the facet index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t f ) const ;
        /*!
         * Get the vertex index of a vertex in a facet
         * in the GeoModelMesh
         * @param[in] f the facet index
         * @param[in] v the local vertex index [0, nb_vertices_in_facet[
         * @return the vertex index
         */
        index_t vertex( index_t f, index_t v ) const ;
        /*!
         * Get the adjacent facet index in the GeoModelMesh
         * @param[in] f the facet index
         * @param[in] e the edge index
         * @return the adjacent facet index
         */
        index_t adjacent( index_t f, index_t e ) const ;
        /*!
         * Get the surface index in the GeoModel according the facet
         * index in the GeoModelMesh
         * @param[in] f the facet index
         * @return the surface index
         */
        index_t surface( index_t f ) const ;
        /*!
         * Get the facet index in the GeoModelMesh restricted to
         * the surface owing the facet
         * @param[in] f the facet index
         * @return the facet index varying from 0 to nb_facets
         * in the surface owing \p f
         */
        index_t index_in_surface( index_t f ) const ;
        /*!
         * Get the facet index in the GeoModelMesh restricted to
         * the surface owing the facet and its type
         * @param[in] f the facet index
         * @param[out] index the facet index varying from 0 to nb_facets
         * of the corresponding type of \p f in the owing surface
         * @return the type of the facet \p f
         */
        FacetType type( index_t f, index_t& index ) const ;

        /*!
         * Get the number of facets of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of facets
         */
        index_t nb_facets( FacetType type = ALL ) const ;
        /*!
         * Get the number of facets of the corresponding type
         * in the given surface of the GeoModel
         * @param[in] s the surface index
         * @param[in] type the corresponding type
         * @return the number of facets
         */
        index_t nb_facets( index_t s, FacetType type = ALL ) const ;
        /*!
         * Get the facet index in the GeoModelMesh
         * @param[in] s the surface index owing the facet
         * @param[in] f the facet index varying from 0 to the number of facets
         * of type \p type in the surface \p s.
         * @warning \p f is NOT a facet id
         * of the surface \p s.
         * It is fth facet of type \p type in the internal storage of the
         * GeoModelMeshFacets (see GeoModelMeshFacets::surface_facet_ptr_).
         * @note to find the facet id of the GeoModelMeshFacets from a surface
         * and a facet id of this surface, you need to perform a search using
         * NNSearch and the barycenter of the facet for instance.
         * @param[in] type it can specify the facet type used. For example, if type = QUAD
         * then \p f represents the fth quad in the surface \p s and \p f can vary from 0
         * to nb_quads( s ).
         * If \p type is FacetType::ALL, all the facet types are
         * taken into account.
         * @return the facet index
         */
        index_t facet( index_t s, index_t f, FacetType type = ALL ) const ;

        /*!
         * Get the number of triangles in the GeoModelMesh
         * @return the number of triangles
         */
        index_t nb_triangle() const ;
        /*!
         * Get the number of triangles in the given surface
         * @param[in] s the surface index
         * @return the number of triangles
         */
        index_t nb_triangle( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked triangle in the surface
         * @param[in] s the surface index
         * @param[in] t the tth triangle index varying from 0 to nb_triangles( s )
         * @return the facet index
         */
        index_t triangle( index_t s, index_t t ) const ;

        /*!
         * Get the number of quads in the GeoModelMesh
         * @return the number of quads
         */
        index_t nb_quad() const ;
        /*!
         * Get the number of quads in the given surface
         * @param[in] s the surface index
         * @return the number of quads
         */
        index_t nb_quad( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked quad in the surface
         * @param[in] s the surface index
         * @param[in] q the qth quad index varying from 0 to nb_quads( s )
         * @return the facet index
         */
        index_t quad( index_t s, index_t q ) const ;

        /*!
         * Get the number of polygons in the GeoModelMesh
         * @return the number of polygons
         */
        index_t nb_polygon() const ;
        /*!
         * Get the number of polygons in the given surface
         * @param[in] s the surface index
         * @return the number of polygons
         */
        index_t nb_polygon( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked polygon in the surface
         * @param[in] s the surface index
         * @param[in] p the pth polygon index varying from 0 to nb_polygons( s )
         * @return the facet index
         */
        index_t polygon( index_t s, index_t p ) const ;

        /*!
         * Clear the facets of the GeoModelMesh
         */
        void clear() ;

        /*!
         * Get the center of the given facet
         * @param[in] f the facet index
         */
        vec3 center( index_t f ) const ;
        /*!
         * Get the area of the facet
         * @param[in] f the facet index
         */
        double area( index_t f ) const ;
        /*!
         * Get the normal of the facet
         * @param[in] f the facet index
         */
        vec3 normal( index_t f ) const ;

        const NNSearch& nn_search() const
        {
            test_and_initialize() ;
            return mesh_->facets_nn_search() ;
        }

        /*!
         * @brief return the AABB tree for the facets of the mesh
         */
        const AABBTree2D& aabb() const ;

    private:
        /*!
         * Initialize the facets of the GeoModelMesh
         * and sort them per surface and facet type
         * Example for a mesh with two surfaces and only triangles and quads
         * [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
         * |          surface 0           |             surface 1           |
         */
        void initialize() ;

        /*!
         * Bind attribute to the facets attribute manager
         */
        void bind_attribute() ;
        /*!
         * Unbind attribute to the facets attribute manager
         */
        void unbind_attribute() ;
        /*!
         * @brief Removes facet adjacencies along lines
         */
        void disconnect_along_lines() ;

    private:
        /// Attached Mesh
        std::unique_ptr< Mesh2D > mesh_ ;

        /// Attribute storing the surface index per facet
        GEO::Attribute< index_t > surface_id_ ;
        /// Attribute storing the facet index in surface per facet
        GEO::Attribute< index_t > facet_id_ ;

        /*!
         * Vector storing the index of the starting facet index
         * for a given surface and a given facet type.
         * For example:
         *    the 2nd quad index of the surface index S will be found here:
         *    surface_facet_ptr_[ALL*S + QUAD] + 2
         */
        std::vector< index_t > surface_facet_ptr_ ;

        /// Number of triangles in the GeoModelMesh
        index_t nb_triangle_ ;
        /// Number of quads in the GeoModelMesh
        index_t nb_quad_ ;
        /// Number of polygons in the GeoModelMesh
        index_t nb_polygon_ ;
    } ;

    class RINGMESH_API GeoModelMeshEdges: public GeoModelMeshBase {
    ringmesh_disable_copy( GeoModelMeshEdges ) ;
    public:
        GeoModelMeshEdges( GeoModelMesh& gmm, GeoModel& gm ) ;
        ~GeoModelMeshEdges() ;

        GEO::AttributesManager& attribute_manager() const
        {
            return mesh_->edge_attribute_manager() ;
        }

        /*!
         * Test if the mesh edges are initialized
         */
        bool is_initialized() const ;
        /*!
         * Tests if the mesh edges needs to be initialized and initialize it
         */
        void test_and_initialize() const ;

        /*!
         * Gets the number of wells
         * @return the corresponding number
         */
        index_t nb_wells() const ;
        /*!
         * Gets the number of edges in the MacroMesh
         * @return the corresponding number
         */
        index_t nb_edges() const ;
        /*!
         * Gets the number of edges of a Well
         * @param[in] w the well index
         * @return the corresponding number
         */
        index_t nb_edges( index_t w ) const ;
        /*!
         * Gets the vertex index of the GeoModelMesh
         * @param[in] w the well index
         * @param[in] e the edge index in the well (from 0 to nb_edges in the well)
         * @param[in] v the vertex index of the edge (0 or 1 )
         * @return the global vertex index
         */
        index_t vertex( index_t w, index_t e, index_t v ) const ;
        /*!
         * Clear the mesh edges
         */
        void clear() ;

        /*!
         * Initialize the mesh edges
         */
        void initialize() ;

        /*!
         * @brief return the AABB tree for the edges of the mesh
         */
        const AABBTree1D& aabb() const ;

    private:
        /// Attached Mesh
        std::unique_ptr< Mesh1D > mesh_ ;

        /*!
         * Vector storing the index of the starting edge index
         * for a given well
         */
        std::vector< index_t > well_ptr_ ;
    } ;

    class RINGMESH_API GeoModelMeshCells: public GeoModelMeshBase {
    ringmesh_disable_copy( GeoModelMeshCells ) ;
    public:
        friend class GeoModelMesh ;

        /*!
         * Several modes for vertex duplication algorithm:
         *  - NONE = no duplication
         *  - FAULT = duplication along faults
         *  - HORIZON = duplication along horizons
         *  - ALL = duplication along faults and horizons
         */
        enum DuplicateMode {
            NONE, FAULT, HORIZON, ALL
        } ;

        GeoModelMeshCells( GeoModelMesh& gmm, GeoModel& gm ) ;

        GEO::AttributesManager& attribute_manager() const
        {
            return mesh_->cell_attribute_manager() ;
        }
        /*!
         * Test if the mesh cells are initialized
         */
        bool is_initialized() const ;
        /*!
         * Test if the mesh cells are duplicated
         */
        bool is_duplication_initialized() const ;

        /*!
         * Test if the mesh cells need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;

        /*!
         * @brief Number of cells stored.
         */
        index_t nb() const ;
        /*!
         * Gets the number of duplicated points by the DuplicateMode algorithm
         * @return the corresponding number of duplications
         */
        index_t nb_duplicated_vertices() const ;
        /*!
         * Gets the total number of vertices (mesh.vertices.nb() + nb_duplicated_vertices())
         * @return the corresponding number of vertices
         */
        index_t nb_total_vertices() const ;
        /*!
         * Check if the corner in a cell is duplicated,
         * if so give the duplicated vertex index
         * @param[in] c the cell index in the GeoModelMesh
         * @param[in] v the local vertex index in the cell \p c (0 to nb_vertices( c ))
         * @param[out] duplicate_vertex_index the duplicated vertex index (0 to nb_duplicated_vertices())
         * @return true if the corner is duplicated
         */
        bool is_corner_duplicated(
            index_t c,
            index_t v,
            index_t& duplicate_vertex_index ) const ;
        /*!
         * Get the vertex index in the GeoModelMesh corresponding
         * to the given duplicated vertex index
         * @param[in] duplicate_vertex_index the duplicated vertex index
         * @return the vertex index
         */
        index_t duplicated_vertex( index_t duplicate_vertex_index ) const ;

        /*!
         * Get the number of vertices in the cell
         * @param[in] c the cell index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t c ) const ;
        /*!
         * Get the vertex index of a vertex in a cell
         * in the GeoModelMesh
         * @param[in] c the cell index
         * @param[in] v the local vertex index [0, nb_vertices_in_cell[
         * @return the vertex index
         */
        index_t vertex( index_t c, index_t v ) const ;
        /*!
         * Get the number of edges in the cell
         * @param[in] c the cell index
         */
        index_t nb_edges( index_t c ) const ;
        /*!
         * Get the number of facets in the cell
         * @param[in] c the cell index
         */
        index_t nb_facets( index_t c ) const ;
        /*!
         * Get the number of facets in the cell
         * @param[in] c the cell index
         * @param[in] lf the cell facet index
         */
        index_t nb_facet_vertices( index_t c, index_t lf ) const ;
        /*!
         * \brief Gets a cell vertex by local facet index and local
         *  vertex index in the edge
         * \param[in] c the cell, in 0..nb()-1
         * \param[in] lf the local facet index, in 0..nb_facets(c)-1
         * \param[in] lv the local index in the cell facet
         * \return vertex \p lv of facet \p lf in cell \p c
         */
        index_t facet_vertex( index_t c, index_t lf, index_t lv ) const ;
        /*!
         * \brief Gets a cell vertex by local edge index and local
         *  vertex index in the edge
         * \param[in] c the cell, in 0..nb()-1
         * \param[in] le the local edge index, in 0..nb_edges(c)-1
         * \param[in] lv the local index in the edge, one of 0,1
         * \return vertex \p lv of edge \p le in cell \p c
         */
        index_t edge_vertex( index_t c, index_t le, index_t lv ) const ;
        /*!
         * Get the adjacent cell index in the GeoModelMesh
         * @param[in] c the cell index
         * @param[in] f the edge index
         * @return the adjacent cell index
         */
        index_t adjacent( index_t c, index_t f ) const ;
        /*!
         * Get the region index in the GeoModel according the cell
         * index in the GeoModelMesh
         * @param[in] c the cell index
         * @return the region index
         */
        index_t region( index_t c ) const ;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell
         * @param[in] c the cell index
         * @return the cell index varying from 0 to nb_cells
         * in the region owing \p c
         */
        index_t index_in_region( index_t c ) const ;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell and its type
         * @param[in] c the cell index
         * @param[out] index the cell index varying from 0 to nb_cells
         * of the corresponding type of \p c in the owing region
         * @return the type of the cell \p c
         */
        GEO::MeshCellType type( index_t c ) const ;

        /*!
         * Get the number of cells of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells( GEO::MeshCellType type = GEO::MESH_NB_CELL_TYPES ) const ;
        /*!
         * Get the number of cells of the corresponding type
         * in the given region of the GeoModel
         * @param[in] r the region index
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells(
            index_t r,
            GEO::MeshCellType type = GEO::MESH_NB_CELL_TYPES ) const ;
        /*!
         * Get the cell index in the GeoModelMesh
         * @param[in] r the region index owing the cell
         * @param[in] c the cell index varying from 0 to number of cells
         * of type \p type in the region \p r.
         * @warning \p c is NOT a cell id of the region \p r,
         * It is cth cell of type \p type in the internal storage of the
         * GeoModelMeshCells (see GeoModelMeshCells::region_cell_ptr_).
         * @note to find the cell id of the GeoModelMeshCells from a region
         * and a cell id of this region, you need to perform a search using
         * NNSearch and the barycenter of the cell for instance.
         * @param[in] type it can specify the cell type used. For example,
         * if type = GEO::MESH_HEX then \p c represents the cth hex in the
         * region \p r and \p c can vary from 0 to nb_hex( r ).
         * If \p type is GEO::MESH_NB_CELL_TYPES, all the cell types are
         * taken into account.
         * @return the cell index
         */
        index_t cell( index_t r, index_t c, GEO::MeshCellType type =
            GEO::MESH_NB_CELL_TYPES ) const ;

        /*!
         * Get the number of tets in the GeoModelMesh
         * @return the number of tets
         */
        index_t nb_tet() const ;
        /*!
         * Get the number of tets in the given region
         * @param[in] r the region index
         * @return the number of tets
         */
        index_t nb_tet( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked tet in the region
         * @param[in] r the region index
         * @param[in] t the tth tet index varying from 0 to nb_tet( r )
         * @return the cell index
         */
        index_t tet( index_t r, index_t t ) const ;

        /*!
         * Get the number of hexs in the GeoModelMesh
         * @return the number of hexs
         */
        index_t nb_hex() const ;
        /*!
         * Get the number of hexs in the given region
         * @param[in] r the region index
         * @return the number of hexs
         */
        index_t nb_hex( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked hex in the region
         * @param[in] r the region index
         * @param[in] h the hth hex index varying from 0 to nb_hex( r )
         * @return the cell index
         */
        index_t hex( index_t r, index_t h ) const ;

        /*!
         * Get the number of prisms in the GeoModelMesh
         * @return the number of prisms
         */
        index_t nb_prism() const ;
        /*!
         * Get the number of prisms in the given region
         * @param[in] r the region index
         * @return the number of prisms
         */
        index_t nb_prism( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked prism in the region
         * @param[in] r the region index
         * @param[in] p the pth prism index varying from 0 to nb_prism( r )
         * @return the cell index
         */
        index_t prism( index_t r, index_t p ) const ;

        /*!
         * Get the number of pyramids in the GeoModelMesh
         * @return the number of pyramids
         */
        index_t nb_pyramid() const ;
        /*!
         * Get the number of pyramids in the given region
         * @param[in] r the region index
         * @return the number of pyramids
         */
        index_t nb_pyramid( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked pyramid in the region
         * @param[in] r the region index
         * @param[in] p the pth pyramid index varying from 0 to nb_pyramid( r )
         * @return the cell index
         */
        index_t pyramid( index_t r, index_t p ) const ;

        /*!
         * Get the number of connectors in the GeoModelMesh
         * @return the number of connectors
         */
        index_t nb_connector() const ;
        /*!
         * Get the number of connectors in the given region
         * @param[in] r the region index
         * @return the number of connectors
         */
        index_t nb_connector( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked connector in the region
         * @param[in] r the region index
         * @param[in] c the cth connector index varying from 0 to nb_connector( r )
         * @return the cell index
         */
        index_t connector( index_t r, index_t c ) const ;

        /*!
         * Clear the mesh cells
         */
        void clear() ;
        /*!
         * Remove the duplication of the mesh cell facets
         */
        void clear_duplication() ;

        /*!
         * Determine if a cell facet is on a surface. If so, fill the \p action
         * with the surface index and the surface side encountered
         * @param[in] c the cell index
         * @param[in] f the facet index
         * @param[out] facet the facet index colocalised with the cell facet
         * @param[out] side the side of the facet \p facet.
         * true = side of the facet normal, false = the other side
         * @return true is the cell facet is on a surface
         */
        bool is_cell_facet_on_surface(
            index_t c,
            index_t f,
            index_t& facet,
            bool& side ) const ;

        /*!
         * Get the center of the given cell
         * @param[in] c the cell index
         */
        vec3 barycenter( index_t c ) const ;
        /*!
         * Get the volume of the cell
         * @param[in] c the cell index
         */
        double volume( index_t c ) const ;

        const NNSearch& cell_nn_search() const
        {
            test_and_initialize() ;
            return mesh_->cells_nn_search() ;
        }
        const NNSearch& cell_facet_nn_search() const
        {
            test_and_initialize() ;
            return mesh_->cell_facets_nn_search() ;
        }

        /*!
         * @brief return the AABB tree for the cells of the mesh
         */
        const AABBTree3D& aabb() const ;

    private:
        /// enum to characterize the action to do concerning a surface
        /// Action concerns the vertices of a Surface and not the Surface
        enum ActionOnSurface {
            /// do nothing
            SKIP = -2,
            /// need to be duplicated (don't know which side yet)
            TO_PROCESS = -1,
            /// need to duplicate the side opposite to the facet normal
            NEG_SIDE = 0,
            /// need to duplicate the side following the facet normal
            POS_SIDE = 1
        } ;
        /// Action to do according a surface index
        using action_on_surface = std::pair< index_t, ActionOnSurface > ;

        /*!
         * @brief Initialize the  cells from the cells
         *        of the GeoModel Region cells
         * @details Fills the mesh_.cells
         */
        void initialize() ;

        /*!
         * Bind attribute to the cells attribute manager
         */
        void bind_attribute() ;
        /*!
         * Unbind attribute to the cells attribute manager
         */
        void unbind_attribute() ;

        /*!
         * Test if the mesh cell are duplicated according
         * the duplication mode, if not duplicate them.
         */
        void test_and_initialize_duplication() const ;
        /*!
         * Duplicate the mesh cell along some surfaces defined
         * by the duplication mode
         */
        void initialize_duplication() ;
        /*!
         * Test if we need to duplicate mesh cell along the given
         * surface according the duplicate mode
         * @param[in] s the surface index in the GeoModel
         */
        bool is_surface_to_duplicate( index_t s ) const ;

        /*!
         * Determine the actions to do according the action_on_surfaces
         * encountered during the propagation around a vertex (initialize())
         * @param[in] surfaces the action_on_surfaces encountered
         * @param[in,out] info the global information on what to do for each surface.
         * This information is updated in this function according the encountered action_on_surfaces
         * @return true if the corners should be duplicated
         */
        bool are_corners_to_duplicate(
            const std::vector< action_on_surface >& surfaces,
            std::vector< ActionOnSurface >& info ) ;
        /*!
         * Test if the mesh cell facet attribute is filled with
         * the colocalised facet. If not fill it.
         */
        void test_and_initialize_cell_facet() const ;
        /*!
         * Initialize the mesh cell facet attribute of colocalised facet.
         */
        void initialize_cell_facet() ;

    private:
        /// Attached Mesh
        std::unique_ptr< Mesh3D > mesh_ ;

        /// Attribute storing the region index per cell
        GEO::Attribute< index_t > region_id_ ;
        /// Attribute storing the cell index in region per cell
        GEO::Attribute< index_t > cell_id_ ;

        /*!
         * Vector storing the index of the starting cell index
         * for a given region and a given cell type.
         * For example:
         *    the 2nd hex index of the region index R will be found here:
         *    surface_facet_ptr_[GEO::MESH_NB_CELL_TYPES*R + HEX] + 2
         */
        std::vector< index_t > region_cell_ptr_ ;

        /// Number of tet in the GeoModelMesh
        index_t nb_tet_ ;
        /// Number of hex in the GeoModelMesh
        index_t nb_hex_ ;
        /// Number of prism in the GeoModelMesh
        index_t nb_prism_ ;
        /// Number of pyramid in the GeoModelMesh
        index_t nb_pyramid_ ;
        /// Number of connector in the GeoModelMesh
        index_t nb_connector_ ;

        /// Current duplicate mode applied on the mesh
        DuplicateMode mode_ ;
        /*!
         * @brief Vector of duplicated vertices
         * @details Each value is a duplicated vertex, the index corresponds to
         * vertex index in mesh.vertices.
         */
        std::vector< index_t > duplicated_vertex_indices_ ;

        /*!
         * @brief Attribute storing the colocalised facet index per cell facet
         * @detail If a cell facet is on a surface, the attribute is equal to
         * the index of the corresponding facet.
         */
        GEO::Attribute< index_t > facet_id_ ;
    } ;

    class RINGMESH_API GeoModelMesh {
    ringmesh_disable_copy( GeoModelMesh ) ;
    public:
        GeoModelMesh( GeoModel& geomodel ) ;
        ~GeoModelMesh() ;

        const GeoModel& geomodel() const
        {
            return geomodel_ ;
        }

        /*!
         * @brief Transfer attributes from the GeoModelMesh to the
         * GeoModel
         */
        void transfert_attributes() const ;
        void transfert_attributes_from_gmm_to_gm() const ;

        /*!
         * @brief Transfer attributes from the GeoModelMeshCell to the
         * GeoModel
         */
        void transfert_cell_attributes() const ;
        void transfert_cell_attributes_from_gmm_to_gm() const ;
        /*!
         * @brief Transfer attributes from the GeoModelMeshVertices to the
         * GeoModel
         */
        void transfert_vertex_attributes() const ;
        void transfert_vertex_attributes_from_gmm_to_gm() const ;

        /*!
         * Access the DuplicateMode
         * @return the current DuplicateMode
         */
        GeoModelMeshCells::DuplicateMode duplicate_mode() const
        {
            return mode_ ;
        }
        /*!
         * Set a new DuplicateMode
         * @param[in] mode the new DuplicateMode for the GeoModelMesh
         */
        void set_duplicate_mode( const GeoModelMeshCells::DuplicateMode& mode ) const
        {
            if( mode_ == mode ) return ;
            mode_ = mode ;
            const_cast< GeoModelMesh* >( this )->cells.clear_duplication() ;
        }

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated_vertices() ;

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the geomodel_vertex_id in the GeoModelMeshEntity
         * of the GeoModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with which information should be merged.
         *  It is recycled to give the mapping between old and new vertex indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) ;

        /*!
         * @brief Remove all invalid GMEVertex and delete the vertices
         * that are not anymore in any GeoModelEntity
         */
        void erase_invalid_vertices() ;

    private:
        /*! Attached GeoModel */
        const GeoModel& geomodel_ ;

        /// Optional duplication mode to compute the duplication of cells on surfaces
        mutable GeoModelMeshCells::DuplicateMode mode_ ;

    public:
        GeoModelMeshVertices vertices ;
        GeoModelMeshEdges edges ;
        GeoModelMeshFacets facets ;
        GeoModelMeshCells cells ;
    } ;

}
