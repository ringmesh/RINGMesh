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

/*!
 * @file Declaration of GeoModelMeshEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/entity_type.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>

namespace RINGMesh {
    class GeoModel;
    class GeoModelGeologicalEntity;
}

namespace RINGMesh {
    /*!
     * @brief Abstract base class for GeoModelMeshEntity.
     * @details The GeoModelMeshEntity geometrical representation
     * is stored as a RINGMesh::Mesh. We defines generic access to
     * the RINGMesh::Mesh geometry.
     */
    class RINGMESH_API GeoModelMeshEntity: public GeoModelEntity {
    ringmesh_disable_copy( GeoModelMeshEntity );
        friend class GeoModelMeshEntityAccess;
        friend class GeoModelMeshEntityConstAccess;

    public:
        virtual ~GeoModelMeshEntity();

        virtual MeshEntityType type_name() const = 0;

        gmme_id gmme() const
        {
            return gmme_id( type_name(), id_ );
        }

        MeshEntityType mesh_entity_type() const
        {
            return gmme().type();
        }
        /*!
         * @brief Global validity of the entity
         */
        virtual bool is_valid() const final
        {
            return is_mesh_valid() && are_geomodel_vertex_indices_valid();
        }
        /*!
         * Check that required information for the TYPE is defined
         *    and that reverse information is stored by the corresponding
         *    entities
         */
        virtual bool is_connectivity_valid() const;

        /*!
         *  Check that for each parent type, this entity has one and only
         *  one parent of that type.
         */
        bool is_parent_connectivity_valid() const;

        /*!
         * \name Boundary relationship functions
         * @{
         */
        index_t nb_boundaries() const
        {
            return static_cast< index_t >( boundaries_.size() );
        }
        const gmme_id& boundary_gmme( index_t x ) const;
        const GeoModelMeshEntity& boundary( index_t x ) const;

        index_t nb_incident_entities() const
        {
            return static_cast< index_t >( incident_entities_.size() );
        }
        const gmme_id& incident_entity_gmme( index_t x ) const;
        const GeoModelMeshEntity& incident_entity( index_t x ) const;

        /*!
         * @brief Check if one entity is twice in the boundary
         */
        bool has_inside_border() const;
        /*!
         * @brief Check if this entity is an inside border of rhs
         * @details That can be Surface stopping in a Region, or Line stopping in a Surface.
         * @param[in] rhs The entity to test
         */
        bool is_inside_border( const GeoModelMeshEntity& rhs ) const;

        /*! @}
         * \name Parents/children relationship functions
         * @{
         */
        bool has_parent() const
        {
            return nb_parents() != 0;
        }
        /*!
         * @brief Check if the entity has a parent of the given type
         */
        bool has_parent( const GeologicalEntityType& parent_type ) const
        {
            return could_be_undefined_parent_gmge( parent_type ).is_defined();
        }
        index_t nb_parents() const
        {
            return static_cast< index_t >( parents_.size() );
        }
        const gmge_id& parent_gmge( index_t id ) const;
        /*!
         * @brief Returns the gmge_id of the parent of the given type.
         * @pre The code assumes that this entity has a parent of the given type
         * \p parent_type_name, else it will crash in debug and return an invalid
         * entity in release.
         * The method has_parent(GeologicalEntityType) is a safe way to check
         * if the type \p parent_type_name is a valid parent type for the
         * GeoModelMeshEntity.
         * @param[in] parent_type_name the asking parent type
         *
         */
        gmge_id parent_gmge( const GeologicalEntityType& parent_type ) const;
        const GeoModelGeologicalEntity& parent( index_t id ) const;
        const GeoModelGeologicalEntity& parent(
            const GeologicalEntityType& parent_type ) const;

        /*!
         * @}
         */

        void save( const std::string& filename ) const
        {
            mesh_->save_mesh( filename );
        }
        /*!
         * @brief Return the NNSearch for the Entity vertices.
         */

        const NNSearch& vertex_nn_search() const
        {
            return mesh_->vertex_nn_search();
        }

        /*!
         * \name Local access to the GeoModelMeshEntity geometry
         * @{
         */
        index_t nb_vertices() const
        {
            return mesh_->nb_vertices();
        }
        /*!
         * @brief Coordinates of the \p vertex_index.
         */
        const vec3& vertex( index_t vertex_index ) const
        {
            return mesh_->vertex( vertex_index );
        }

        /*!
         * @brief Get the number constitutive elements of the mesh
         * @details Constitutive elements are those of the dimension of the object.
         * segments for lines, polygons for surfaces and cells for volumetric meshes.
         */
        virtual index_t nb_mesh_elements() const = 0;
        /*!
         * @brief Number of vertices of a constitutive element of the mesh
         */
        virtual index_t nb_mesh_element_vertices(
            index_t mesh_element_index ) const = 0;
        /*!
         * @brief Convert the index in a mesh element to an index in the Entity.
         * @param[in] mesh_element_index Index of a constitutive element of the mesh
         * @param[in] vertex_local_index Local index of a vertex in the constitutive
         * mesh element
         * @return the global index of the vertex in the GeoModelMeshEntity
         */
        virtual index_t mesh_element_vertex_index(
            index_t mesh_element_index,
            index_t vertex_local_index ) const = 0;
        /*!
         * @brief Coordinates of a vertex of a mesh element.
         * @param[in] mesh_element_index Index of a constitutive element of the mesh
         * @param[in] vertex_local_index Local index of a vertex in the constitutive
         * mesh element
         * @return the vertex coordinates in the GeoModelMeshEntity
         */
        const vec3& mesh_element_vertex(
            index_t mesh_element_index,
            index_t vertex_local_index ) const
        {
            return vertex(
                mesh_element_vertex_index( mesh_element_index, vertex_local_index ) );
        }

        /*! @}
         * \name Geometrical request on Entity
         * @{
         */
        virtual double mesh_element_size( index_t mesh_element_index ) const = 0;
        virtual vec3 mesh_element_barycenter( index_t mesh_element_index ) const = 0;
        virtual vec3 entity_barycenter() const
        {
            vec3 result( 0., 0., 0. );
            for( index_t v = 0; v < nb_vertices(); v++ ) {
                result += vertex( v );
            }
            ringmesh_assert( nb_vertices() > 0 );
            return result / static_cast< double >( nb_vertices() );
        }
        virtual double size() const
        {
            double size = 0.0;
            for( index_t i = 0; i < nb_mesh_elements(); ++i ) {
                size += mesh_element_size( i );
            }
            return size;
        }

        /*! @}
         */
        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_->vertex_attribute_manager();
        }

    protected:
        GeoModelMeshEntity(
            const GeoModel& geomodel,
            index_t id,
            const std::string& name = "No_name" )
            : GeoModelEntity( geomodel, id, name ), mesh_( nullptr )
        {
        }
        virtual void copy_mesh_entity( const GeoModelMeshEntity& from )
        {
            copy_name( from );
            id_ = from.id_;
            boundaries_ = from.boundaries_;
            incident_entities_ = from.incident_entities_;
            parents_ = from.parents_;
        }
        virtual bool is_index_valid() const final;
        virtual bool is_mesh_valid() const
        {
            return mesh_ != nullptr;
        }

        void set_mesh( std::shared_ptr< MeshBase > mesh )
        {
            ringmesh_assert( mesh != nullptr );
            mesh_ = std::move( mesh );
        }

        /*!
         * All entities in the boundary must have this in their
         *  incident_entity vector
         */
        bool is_boundary_connectivity_valid() const;
        /*!
         * All entities must be at least in the boundary of another entity
         * and all entities in the incident_entity must have this entity in their
         * boundary vector
         */
        bool is_incident_entity_connectivity_valid() const;
        /*!
         * @brief Check that geomodel vertex indices are consistent
         * with what is stored at the GeoModel level.
         * @todo Review: this dependancy is bad, and puts us in a lot of trouble [JP]
         */
        bool are_geomodel_vertex_indices_valid() const;

        void unbind_vertex_mapping_attribute() const;
        void bind_vertex_mapping_attribute() const;

        virtual void change_mesh_data_structure( const MeshType type ) = 0;
    private:
        gmge_id defined_parent_gmge(
            const GeologicalEntityType& parent_type ) const;

        gmge_id could_be_undefined_parent_gmge(
            const GeologicalEntityType& parent_type ) const;
    protected:

        /// Boundary relations of this entity
        std::vector< index_t > boundaries_;

        /// Incident-entity relations of this entity
        std::vector< index_t > incident_entities_;

        /// Parents relations of this entity
        std::vector< index_t > parents_;
    private:
        /// The RINGMesh::Mesh giving the geometry of this entity
        std::shared_ptr< MeshBase > mesh_;
    };

    /*!
     * @brief A GeoModelEntity of type CORNER
     * @details It is a unique point.
     */
    class RINGMESH_API Corner: public GeoModelMeshEntity {
    public:
        friend class GeoModelMeshEntityAccess;
        friend class GeoModelMeshEntityConstAccess;

        virtual ~Corner()
        {
            unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static()
        {
            return MeshEntityType( "Corner" );
        }

        virtual MeshEntityType type_name() const override
        {
            return type_name_static();
        }

        /*!
         * @brief Checks if this entity define the geomodel external boundary
         * @details Test if the entity is in the Surfaces defining the universe
         */
        virtual bool is_on_voi() const final;

        /*!
         * @return 0, no mesh_element are defined for corners.
         */
        virtual index_t nb_mesh_elements() const final
        {
            return 0;
        }
        /*!
         * @return the number of vertices of the Corner
         */
        virtual index_t nb_mesh_element_vertices( index_t mesh_element = 0 ) const override
        {
            ringmesh_unused( mesh_element );
            index_t nb_vertices = mesh0d_->nb_vertices();
            ringmesh_assert( nb_vertices < 2 );
            return nb_vertices;
        }
        const Line& incident_entity( index_t x ) const;

        /*! @}
         * \name Geometrical request on Corner
         * @{
         */
        virtual double mesh_element_size( index_t mesh_element = 0 ) const override
        {
            ringmesh_unused( mesh_element );
            return 0.0;
        }
        virtual double size() const
        {
            return 0.0;
        }
        virtual vec3 mesh_element_barycenter( index_t mesh_element = 0 ) const override
        {
            ringmesh_unused( mesh_element );
            return vertex( 0 );
        }

        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const PointSetMesh& low_level_mesh_storage() const
        {
            return *mesh0d_;
        }
    protected:
        /*! @brief Creates a Corner.
         *  A point is added to its Mesh.
         */
        Corner( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( PointSetMesh::create_mesh( type ) );
        }

        /*!
         * @brief Get the index of the unique vertex constituting of the Corner.
         * @return 0.
         */
        virtual index_t mesh_element_vertex_index(
            index_t mesh_element = 0,
            index_t vertex_index = 0 ) const override
        {
            ringmesh_unused( mesh_element );
            ringmesh_unused( vertex_index );
            return 0;
        }

        /*!
         * @brief Check that the Corner mesh is a unique point
         */
        virtual bool is_mesh_valid() const final;

    private:

        void update_mesh_storage_type( std::unique_ptr< PointSetMesh > mesh )
        {
            mesh0d_ = std::move( mesh );
            GeoModelMeshEntity::set_mesh( mesh0d_ );
        }
        virtual void change_mesh_data_structure( const MeshType type ) override;

    private:
        std::shared_ptr< PointSetMesh > mesh0d_;
    };

    /*!
     * @brief A GeoModelEntity of type LINE
     *
     * @details This must be one connected component (one part) of
     * a 1-manifold (Line with no T intersections).
     */
    class RINGMESH_API Line: public GeoModelMeshEntity {
    public:
        friend class GeoModelMeshEntityAccess;

        virtual ~Line()
        {
            unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static()
        {
            return MeshEntityType( "Line" );
        }

        virtual MeshEntityType type_name() const override
        {
            return type_name_static();
        }

        virtual bool is_on_voi() const final;

        const Corner& boundary( index_t x ) const;

        const Surface& incident_entity( index_t x ) const;
        virtual bool is_connectivity_valid() const final;

        const LineAABBTree& edge_aabb() const
        {
            return mesh1d_->edge_aabb();
        }

        /*!
         * @brief Return the NNSearch for the edges of the line
         * @details The barycenter of the edges is used.
         */
        const NNSearch& edge_nn_search() const
        {
            return mesh1d_->edge_nn_search();
        }

        /*!
         * Get the number of edges of the Line
         */
        virtual index_t nb_mesh_elements() const final
        {
            return mesh1d_->nb_edges();
        }

        /*!
         * @return The number of vertices per edge: 2.
         */
        virtual index_t nb_mesh_element_vertices( index_t mesh_element = 0 ) const final
        {
            ringmesh_unused( mesh_element );
            return 2;
        }

        /*!
         * @brief Get the index of a vertex in the Line from the
         * @param vertex_index in a given edge @param edge_index.
         */
        virtual index_t mesh_element_vertex_index(
            index_t edge_index,
            index_t vertex_index ) const final
        {
            ringmesh_assert( edge_index < nb_mesh_elements() );
            ringmesh_assert( vertex_index < 2 );
            return mesh1d_->edge_vertex( edge_index, vertex_index );
        }

        /*!
         * @brief A Line is closed if its two extremities are identitical.
         */
        bool is_closed() const
        {
            ringmesh_assert( nb_boundaries() == 2 );
            return ( boundary_gmme( 0 ).is_defined() )
                && ( boundary_gmme( 0 ) == boundary_gmme( 1 ) );
        }

        /*!
         * @brief Gets the length of an edge
         */
        virtual double mesh_element_size( index_t edge_index ) const final
        {
            ringmesh_assert( edge_index < nb_mesh_elements() );
            return mesh1d_->edge_length( edge_index );
        }

        /*!
         * @brief Gets the barycenter of an edge
         */
        virtual vec3 mesh_element_barycenter( index_t edge_index ) const final
        {
            ringmesh_assert( edge_index < nb_mesh_elements() );
            return mesh1d_->edge_barycenter( edge_index );
        }

        bool is_first_corner_first_vertex() const;

        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const LineMesh& low_level_mesh_storage() const
        {
            return *mesh1d_;
        }
    protected:
        Line( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( LineMesh::create_mesh( type ) );
        }

        /*!
         * @brief Check that the mesh of the Line is valid
         * @details Check that
         *  - the GEO::Mesh has more than 1 vertex - more than 1 edge - no polygon - no cell.
         *  - global indices of vertices in the geomodel are in a valid range
         *  - each vertex is in 2 edges except extremities that are in 1 edge
         *
         * Does not check:
         *  - Self-intersection - I suppose there are no segment - segment intersection (JP)
         *  - Duplicated edge - most probably ruled out with the duplicated vertex test (JP)
         *  - Duplicated vertex (verified at GeoModel level)
         */
        virtual bool is_mesh_valid() const final;

    private:
        void update_mesh_storage_type( std::unique_ptr< LineMesh > mesh )
        {
            mesh1d_ = std::move( mesh );
            GeoModelMeshEntity::set_mesh( mesh1d_ );
        }
        virtual void change_mesh_data_structure( const MeshType type ) override;

    private:
        std::shared_ptr< LineMesh > mesh1d_;
    };

    /*!
     * @brief A GeoModelEntity of type SURFACE
     *
     * @details One connected component (part) of a 2-manifold surface
     * (all edges of the polygons are in at most 2 polygons)
     */
    class RINGMESH_API Surface: public GeoModelMeshEntity {
    public:
        friend class GeoModelMeshEntityAccess;

        virtual ~Surface()
        {
            unbind_vertex_mapping_attribute();
        }

        virtual MeshEntityType type_name() const override
        {
            return type_name_static();
        }

        static MeshEntityType type_name_static()
        {
            return MeshEntityType( "Surface" );
        }

        virtual bool is_on_voi() const final;
        const Line& boundary( index_t x ) const;

        const Region& incident_entity( index_t x ) const;

        bool is_simplicial() const
        {
            return mesh2d_->polygons_are_simplicies();
        }

        const SurfaceAABBTree& polygon_aabb() const
        {
            return mesh2d_->polygon_aabb();
        }

        /*!
         * @brief Return the NNSearch for the polygons of the surface
         * @details The barycenter of the polygons is used.
         */
        const NNSearch& polygon_nn_search() const
        {
            return mesh2d_->polygon_nn_search();
        }

        GEO::AttributesManager& polygon_attribute_manager() const
        {
            return mesh2d_->polygon_attribute_manager();
        }

        /*!
         * \name Accessors to Surface polygons, edges and vertices
         * @{
         */
        /*!
         * Number of polygons of the Surface.
         */
        virtual index_t nb_mesh_elements() const final
        {
            return mesh2d_->nb_polygons();
        }

        /*!
         * Number of vertices of a polygon
         */
        virtual index_t nb_mesh_element_vertices( index_t polygon_index ) const final
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            return mesh2d_->nb_polygon_vertices( polygon_index );
        }

        /*!
         * @brief Index of the vertex in the Surface
         * from its index in a polygon of the mesh.
         */
        virtual index_t mesh_element_vertex_index(
            index_t polygon_index,
            index_t vertex_index ) const final
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            ringmesh_assert(
                vertex_index < nb_mesh_element_vertices( polygon_index ) );
            return mesh2d_->polygon_vertex( polygon_index, vertex_index );
        }

        /*!
         * @brief Gets the next vertex index in a polygon.
         */
        index_t next_polygon_vertex_index(
            index_t polygon_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            ringmesh_assert(
                vertex_index < nb_mesh_element_vertices( polygon_index ) );
            return mesh2d_->next_polygon_vertex( polygon_index, vertex_index );
        }

        /*!
         * @brief Gets the previous vertex index in a polygon.
         */
        index_t prev_polygon_vertex_index(
            index_t polygon_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            ringmesh_assert(
                vertex_index < nb_mesh_element_vertices( polygon_index ) );
            return mesh2d_->prev_polygon_vertex( polygon_index, vertex_index );
        }

        /*!
         * @brief Gets the polygon adjacent along an edge of a polygon.
         * @param polygon_index in the polygon
         * @param edge_index in the edge
         * @note The edge index is assumed to be the index of the vertex
         * at which it is starting.
         */
        index_t polygon_adjacent_index(
            index_t polygon_index,
            index_t edge_index ) const
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            ringmesh_assert(
                edge_index < nb_mesh_element_vertices( polygon_index ) );
            return mesh2d_->polygon_adjacent( polygon_index, edge_index );
        }

        /*!
         * @brief Get the previous edge on the border
         * @details The returned border edge is the previous in the way of polygon edges
         * orientation.
         * @param[in] p Input polygon index
         * @param[in] e Edge index in the polygon
         * @param[out] prev_p Previous polygon index
         * @param[out] prev_e Previous edge index in the polygon
         *
         * @pre the surface must be correctly oriented and
         * the given polygon edge must be on border
         * @warning the edge index is in fact the index of the vertex where the edge starts.
         */
        void prev_on_border(
            index_t p,
            index_t e,
            index_t& prev_p,
            index_t& prev_e ) const
        {
            return mesh2d_->prev_on_border( p, e, prev_p, prev_e );
        }

        /*!
         * @brief Get the next edge on the border
         * @details The returned border edge is the next in the way of polygon edges
         * orientation.
         * @param[in] p Input polygon index
         * @param[in] e Edge index in the polygon
         * @param[out] next_p Next polygon index
         * @param[out] next_e Next edge index in the polygon
         *
         * @pre the given polygon edge must be on border
         * @warning the edge index is in fact the index of the vertex where the edge starts.
         */
        void next_on_border(
            index_t p,
            index_t e,
            index_t& next_p,
            index_t& next_e ) const
        {
            return mesh2d_->next_on_border( p, e, next_p, next_e );
        }

        /*!
         * @brief Get the vertex index in a polygon @param polygon_index from its
         * index in the Surface @param surface_vertex_index
         * @return NO_ID or index of the vertex in the polygon
         */
        index_t vertex_index_in_polygon(
            index_t polygon_index,
            index_t surface_vertex_index ) const
        {
            return mesh2d_->vertex_index_in_polygon( polygon_index,
                surface_vertex_index );
        }
        /*!
         * @brief Get the first polygon of the surface that has an edge linking the two vertices (ids in the surface)
         *
         * @param[in] in0 Index of the first vertex in the surface
         * @param[in] in1 Index of the second vertex in the surface
         * @return NO_ID or the index of the polygon
         */
        index_t polygon_from_surface_vertex_ids( index_t in0, index_t in1 ) const
        {
            return mesh2d_->polygon_from_vertex_ids( in0, in1 );
        }

        /*!
         * @brief Determines the polygons around a vertex
         * @param[in] surf_vertex_id Index of the vertex in the surface
         * @param[in] border_only If true only facets on the border are considered
         * @param[in] f0 (Optional) Index of one facet containing the vertex @p surf_vertex_id
         * @return Indices of the facets containing @param surf_vertex_id
         * @note If a facet containing the vertex is given, facets around this
         * vertex is search by propagation. Else, a first facet is found by brute
         * force algorithm, and then the other by propagation
         */
        std::vector< index_t > polygons_around_vertex(
            index_t surf_vertex_id,
            bool border_only,
            index_t first_polygon = NO_ID ) const /// TODO [BC] keep the default NO_ID. Not put in Region::cells_around_vertex...
        {
            return mesh2d_->polygons_around_vertex( surf_vertex_id, border_only,
                first_polygon );
        }

        /*! @}
         * \name Geometrical request on polygons
         * @{
         */
        /*!
         * @return Normal to the polygon
         */
        vec3 polygon_normal( index_t polygon_index ) const
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            return mesh2d_->polygon_normal( polygon_index );
        }

        /*!
         * @brief Computes the normal of the surface at the vertex location
         * it computes the average value of polygon normal neighbors
         * @param[in] vertex_id the vertex index
         * @param[in] p0 index of a polygon that contain the vertex \param vertex_id
         * @return the normal of the surface at the given vertex
         */
        vec3 normal_at_vertex( index_t vertex_id, index_t p0 = NO_ID ) const
        {
            ringmesh_assert( vertex_id < nb_vertices() );
            return mesh2d_->normal_at_vertex( vertex_id, p0 );
        }
        /*!
         * @return Polygon barycenter.
         */
        virtual vec3 mesh_element_barycenter( index_t polygon_index ) const final
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            return mesh2d_->polygon_barycenter( polygon_index );
        }

        /*!
         * @return Area of a polygon.
         */
        virtual double mesh_element_size( index_t polygon_index ) const final
        {
            ringmesh_assert( polygon_index < nb_mesh_elements() );
            return mesh2d_->polygon_area( polygon_index );
        }

        /*!
         * @brief Compute closest vertex in a polygon of a Surface to a point
         * @param[in] polygon_index Polygon index
         * @param[in] query_point Coordinates of the point to which distance is measured
         * @return Index of the vertex of @param polygon_index closest to @param query_point
         */
        index_t closest_vertex_in_polygon(
            index_t polygon_index,
            const vec3& query_point ) const
        {
            return mesh2d_->closest_vertex_in_polygon( polygon_index, query_point );
        }

        /*!
         * Is the edge starting with the given vertex of the polygon on a border of the Surface?
         */
        bool is_on_border( index_t polygon_index, index_t vertex_index ) const
        {
            return mesh2d_->is_edge_on_border( polygon_index, vertex_index );
        }

        /*!
         * Is one of the edges of the polygon on the border of the surface?
         */
        bool is_on_border( index_t polygon_index ) const
        {
            return mesh2d_->is_polygon_on_border( polygon_index );
        }
        /*! @}
         */

        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const SurfaceMesh& low_level_mesh_storage() const
        {
            return *mesh2d_;
        }
    protected:
        Surface( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( SurfaceMesh::create_mesh( type ) );
        }

        /*!
         * @brief Check that the mesh of the Surface is valid
         * @details Check that
         *  - the GEO::Mesh has more than 2 vertices, at least 1 polygon, no cells.
         *  - global indices of vertices in the geomodel are in a valid range
         *  - no degenerate polygon
         *  - one connected component
         *
         *  Some tests are not performed here but globally on the GeoModel
         *  - intersection of polygons
         *  - non-manifold edges
         *  - duplicated vertices are on a boundary Line ending in the Surface
         *
         *
         *  Some tests are not performed
         *  - non-manifold points
         *  - surface orientability
         *  - planarity of polygonal polygons
         *
         * @todo Check that there is no duplicated polygon
         */
        virtual bool is_mesh_valid() const final;

    private:
        void update_mesh_storage_type( std::unique_ptr< SurfaceMesh > mesh )
        {
            mesh2d_ = std::move( mesh );
            GeoModelMeshEntity::set_mesh( mesh2d_ );
        }
        virtual void change_mesh_data_structure( const MeshType type ) override;
    private:
        std::shared_ptr< SurfaceMesh > mesh2d_;
    };

    /*!
     * @brief A GeoModelEntity of type REGION
     *
     * @details A Region a volumetric connected component of the geomodel
     * defined by a set of surfaces.
     * The Region can be only defined by its boundary Surfaces.
     * Its volumetric mesh is optional.
     */
    class RINGMESH_API Region: public GeoModelMeshEntity {
    public:
        friend class GeoModelMeshEntityAccess;

        virtual ~Region()
        {
            unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static()
        {
            return MeshEntityType( "Region" );
        }

        virtual MeshEntityType type_name() const override
        {
            return type_name_static();
        }

        virtual bool is_on_voi() const final;
        const Surface& boundary( index_t x ) const;

        virtual bool is_connectivity_valid() const final;

        bool is_meshed() const
        {
            return mesh3d_->nb_cells() > 0;
        }

        bool is_simplicial() const
        {
            return mesh3d_->cells_are_simplicies();
        }

        const VolumeAABBTree& cell_aabb() const
        {
            return mesh3d_->cell_aabb();
        }

        const NNSearch& cell_facets_nn_search() const
        {
            return mesh3d_->cell_facet_nn_search();
        }

        /*!
         * @brief Return the NNSearch for the cells of the region
         * @details The barycenter of the cells is used.
         */
        const NNSearch& cell_nn_search() const
        {
            return mesh3d_->cell_nn_search();
        }

        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh3d_->cell_attribute_manager();
        }

        /*!
         * \name Accessors to Region cells, facets, edges and vertices
         * @{
         */

        /*!
         * Get the number of cells of the Region.
         */
        virtual index_t nb_mesh_elements() const final
        {
            return mesh3d_->nb_cells();
        }

        /*!
         * Get the number of vertex in the cell \param cell_index of the Region.
         */
        virtual index_t nb_mesh_element_vertices( index_t cell_index ) const final
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->nb_cell_vertices( cell_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        /*!
         * @brief Index of a vertex in the Region from its index in a cell
         */
        virtual index_t mesh_element_vertex_index(
            index_t cell_index,
            index_t vertex_index ) const final
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert(
                    vertex_index < nb_mesh_element_vertices( cell_index ) );
                return mesh3d_->cell_vertex( cell_index, vertex_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        /*!
         * Get the type of a given cell.
         * @pre The Region must be meshed
         */
        CellType cell_type( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->cell_type( cell_index );
            }
            ringmesh_assert_not_reached;
            return CellType::UNDEFINED;
        }

        index_t nb_cell_edges( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->nb_cell_edges( cell_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        index_t nb_cell_facets( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->nb_cell_facets( cell_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        index_t nb_cell_facet_vertices(
            index_t cell_index,
            index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                return mesh3d_->nb_cell_facet_vertices( cell_index, facet_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        index_t cell_edge_vertex_index(
            index_t cell_index,
            index_t edge_index,
            index_t vertex_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( edge_index < nb_cell_edges( cell_index ) );
                ringmesh_assert(
                    vertex_index < nb_mesh_element_vertices( cell_index ) );
                return mesh3d_->cell_edge_vertex( cell_index, edge_index,
                    vertex_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        index_t cell_facet_vertex_index(
            index_t cell_index,
            index_t facet_index,
            index_t vertex_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                ringmesh_assert(
                    vertex_index < nb_mesh_element_vertices( cell_index ) );
                return mesh3d_->cell_facet_vertex( cell_index, facet_index,
                    vertex_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        index_t cell_adjacent_index( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                return mesh3d_->cell_adjacent( cell_index, facet_index );
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        /*! @}
         * \name Geometrical request on Region Entity
         * @{
         */
        bool is_cell_facet_on_border( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                return mesh3d_->cell_adjacent( cell_index, facet_index )
                    == GEO::NO_CELL;
            }
            ringmesh_assert_not_reached;
            return false;
        }

        vec3 cell_facet_barycenter( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                return mesh3d_->cell_facet_barycenter( cell_index, facet_index );
            }
            ringmesh_assert_not_reached;
            return vec3();
        }

        vec3 cell_facet_normal( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
                return mesh3d_->cell_facet_normal( cell_index, facet_index );
            }
            ringmesh_assert_not_reached;
            return vec3();
        }

        /*!
         * @brief Volume of a cell
         */
        virtual double mesh_element_size( index_t cell_index ) const final
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->cell_volume( cell_index );
            }
            ringmesh_assert_not_reached;
            return 0;
        }
        /*!
         * @brief Compute the volume of the Region
         */
        virtual double size() const final
        {
            double result = 0.;
            for( index_t i = 0; i < nb_boundaries(); i++ ) {
                const Surface& surface = dynamic_cast< const Surface& >( boundary(
                    i ) );

                for( index_t t = 0; t < surface.nb_mesh_elements(); t++ ) {
                    const vec3& p0 = surface.mesh_element_vertex( t, 0 );
                    for( index_t v = 1;
                        v + 1 < surface.nb_mesh_element_vertices( t ); ++v ) {
                        double cur_volume = ( dot( p0,
                            cross( surface.mesh_element_vertex( t, v ),
                                surface.mesh_element_vertex( t, v + 1 ) ) ) )
                            / static_cast< double >( 6 );
                        side( i ) ? result -= cur_volume : result += cur_volume;
                    }
                }
            }
            return fabs( result );
        }
        /*!
         * @brief Get the center of the cell \param cell_index
         */
        virtual vec3 mesh_element_barycenter( index_t cell_index ) const final
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() );
                return mesh3d_->cell_barycenter( cell_index );
            }
            ringmesh_assert_not_reached;
            return vec3();
        }

        std::vector< index_t > cells_around_vertex(
            index_t vertex_id,
            index_t cell_hint ) const
        {
            return mesh3d_->cells_around_vertex( vertex_id, cell_hint );
        }

        void compute_region_volumes_per_cell_type(
            double& tet_volume,
            double& pyramid_volume,
            double& prism_volume,
            double& hex_volume,
            double& poly_volume ) const;

        bool side( index_t i ) const
        {
            return sides_[i];
        }
        /*! @}
         */

        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const VolumeMesh& low_level_mesh_storage() const
        {
            return *mesh3d_;
        }
    protected:
        Region( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( VolumeMesh::create_mesh( type ) );
        }

        virtual bool is_mesh_valid() const final;

    private:
        void update_mesh_storage_type( std::unique_ptr< VolumeMesh > mesh )
        {
            mesh3d_ = std::move( mesh );
            GeoModelMeshEntity::set_mesh( mesh3d_ );
        }
        virtual void change_mesh_data_structure( const MeshType type ) override;

        virtual void copy_mesh_entity( const GeoModelMeshEntity& from ) final
        {
            const Region& region_from = dynamic_cast< const Region& >( from );
            GeoModelMeshEntity::copy_mesh_entity( from );
            sides_ = region_from.sides_;
        }
    protected:
        /*! Additional information to store oriented boundary Surfaces
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_;
    private:
        std::shared_ptr< VolumeMesh > mesh3d_;
    };

    class GeoModelMeshEntityConstAccess {
    ringmesh_disable_copy( GeoModelMeshEntityConstAccess );
        friend class GeoModelBuilderGeometry;
        friend class GeoModelBuilderTopology;

    private:
        GeoModelMeshEntityConstAccess( const GeoModelMeshEntity& gme )
            : gmme_( gme )
        {
        }

        const std::shared_ptr< MeshBase >& mesh() const
        {
            return gmme_.mesh_;
        }

        const std::vector< index_t >& incident_entity_relation_ids() const
        {
            return gmme_.incident_entities_;
        }

        const std::vector< index_t >& boundary_relation_ids() const
        {
            return gmme_.boundaries_;
        }
    private:
        const GeoModelMeshEntity& gmme_;
    };

    class GeoModelMeshEntityAccess {
    ringmesh_disable_copy( GeoModelMeshEntityAccess );
        friend class GeoModelBuilderTopology;
        friend class GeoModelBuilderGeometry;
        friend class GeoModelBuilderGeology;
        friend class GeoModelBuilderInfo;
        friend class GeoModelBuilderRemoval;

    private:
        GeoModelMeshEntityAccess( GeoModelMeshEntity& gme )
            : gmme_( gme )
        {
        }

        std::string& modifiable_name()
        {
            return gmme_.name_;
        }

        index_t& modifiable_index()
        {
            return gmme_.id_;
        }

        std::vector< index_t >& modifiable_boundaries()
        {
            return gmme_.boundaries_;
        }

        std::vector< index_t >& modifiable_incident_entities()
        {
            return gmme_.incident_entities_;
        }

        std::vector< bool >& modifiable_sides()
        {
            ringmesh_assert( gmme_.type_name() == Region::type_name_static() );
            return dynamic_cast< Region& >( gmme_ ).sides_;
        }

        std::vector< index_t >& modifiable_parents()
        {
            return gmme_.parents_;
        }

        std::shared_ptr< MeshBase >& modifiable_mesh()
        {
            return gmme_.mesh_;
        }

        void copy( const GeoModelMeshEntity& from )
        {
            gmme_.copy_mesh_entity( from );
        }
        void change_mesh_data_structure( const MeshType type );

        template< typename ENTITY >
        static std::unique_ptr< ENTITY > create_entity(
            const GeoModel& geomodel,
            index_t id,
            const MeshType type )
        {
            return std::unique_ptr< ENTITY >( new ENTITY( geomodel, id, type ) );
        }

    private:
        GeoModelMeshEntity& gmme_;
    };
}
