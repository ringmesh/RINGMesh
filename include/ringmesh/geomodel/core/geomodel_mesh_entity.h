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

/*!
 * @file Declaration of GeoModelMeshEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella
 */

#pragma once

#include <ringmesh/geomodel/core/common.h>

#include <ringmesh/geomodel/core/geomodel_entity.h>

namespace GEO
{
    class AttributesManager;
} // namespace GEO

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntityAccess );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntityConstAccess );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometryBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometry );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoveBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemove );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderInfo );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( Corner );
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    FORWARD_DECLARATION_DIMENSION_CLASS( Line );
    FORWARD_DECLARATION_DIMENSION_CLASS( Region );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( NNSearch );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineAABBTree );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceAABBTree );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeAABBTree );

    ALIAS_2D_AND_3D( GeoModel );

    class GeologicalEntityType;
    class MeshEntityType;
    struct gmge_id;
    struct gmme_id;
    struct ElementLocalVertex;
    struct PolygonLocalEdge;
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Abstract base class for GeoModelMeshEntity.
     * @details The GeoModelMeshEntity geometrical representation
     * is stored as a RINGMesh::Mesh. We defines generic access to
     * the RINGMesh::Mesh geometry.
     */
    template < index_t DIMENSION >
    class geomodel_core_api GeoModelMeshEntity
        : public GeoModelEntity< DIMENSION >
    {
        ringmesh_disable_copy_and_move( GeoModelMeshEntity );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelMeshEntityAccess< DIMENSION >;
        friend class GeoModelMeshEntityConstAccess< DIMENSION >;

    public:
        virtual ~GeoModelMeshEntity();

        virtual MeshEntityType type_name() const = 0;

        gmme_id gmme() const;
        MeshEntityType mesh_entity_type() const;
        /*!
         * @brief Global validity of the entity
         */
        bool is_valid() const final
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
        const GeoModelMeshEntity< DIMENSION >& boundary( index_t x ) const;

        index_t nb_incident_entities() const
        {
            return static_cast< index_t >( incident_entities_.size() );
        }
        const gmme_id& incident_entity_gmme( index_t x ) const;
        const GeoModelMeshEntity< DIMENSION >& incident_entity(
            index_t x ) const;

        /*!
         * @brief Check if one entity is twice in the boundary
         */
        bool has_inside_border() const;
        /*!
         * @brief Check if this entity is an inside border of rhs
         * @details That can be Surface stopping in a Region, or Line stopping
         * in a Surface.
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
        bool has_parent( const GeologicalEntityType& parent_type ) const;

        index_t nb_parents() const
        {
            return static_cast< index_t >( parents_.size() );
        }

        const gmge_id& parent_gmge( index_t id ) const;

        /*!
         * @brief Returns the gmge_id of the parent of the given type.
         * @pre The code assumes that this entity has a parent of the given type
         * \p parent_type_name, else it will crash in debug and return an
         * invalid
         * entity in release.
         * The method has_parent(GeologicalEntityType) is a safe way to check
         * if the type \p parent_type_name is a valid parent type for the
         * GeoModelMeshEntity.
         * @param[in] parent_type_name the asking parent type
         *
         */
        const gmge_id parent_gmge(
            const GeologicalEntityType& parent_type ) const;
        const GeoModelGeologicalEntity< DIMENSION >& parent(
            index_t parent_index ) const;
        const GeoModelGeologicalEntity< DIMENSION >& parent(
            const GeologicalEntityType& parent_type ) const;

        /*!
         * @}
         */

        void save( const std::string& filename ) const;
        /*!
         * @brief Return the NNSearch for the Entity vertices.
         */

        const NNSearch< DIMENSION >& vertex_nn_search() const;

        /*!
         * \name Local access to the GeoModelMeshEntity geometry
         * @{
         */

        index_t nb_vertices() const;

        /*!
         * @brief Coordinates of the \p vertex_index.
         */
        const vecn< DIMENSION >& vertex( index_t vertex_index ) const;

        /*!
         * @brief Get the number constitutive elements of the mesh
         * @details Constitutive elements are those of the dimension of the
         * object.
         * segments for lines, polygons for surfaces and cells for volumetric
         * meshes.
         */
        virtual index_t nb_mesh_elements() const = 0;
        /*!
         * @brief Number of vertices of a constitutive element of the mesh
         */
        virtual index_t nb_mesh_element_vertices(
            index_t mesh_element_index ) const = 0;
        /*!
         * @brief Convert the index in a mesh element to an index in the Entity.
         * @param[in] mesh_element_index Index of a constitutive element of the
         * mesh
         * @param[in] vertex_local_index Local index of a vertex in the
         * constitutive
         * mesh element
         * @return the global index of the vertex in the GeoModelMeshEntity
         */
        virtual index_t mesh_element_vertex_index(
            const ElementLocalVertex& element_local_vertex ) const = 0;
        /*!
         * @brief Coordinates of a vertex of a mesh element.
         * @param[in] mesh_element_index Index of a constitutive element of the
         * mesh
         * @param[in] vertex_local_index Local index of a vertex in the
         * constitutive
         * mesh element
         * @return the vertex coordinates in the GeoModelMeshEntity
         */
        const vecn< DIMENSION >& mesh_element_vertex(
            const ElementLocalVertex& element_local_vertex ) const
        {
            return vertex( mesh_element_vertex_index( element_local_vertex ) );
        }

        /*! @}
         * \name Geometrical request on Entity
         * @{
         */

        virtual double mesh_element_size(
            index_t mesh_element_index ) const = 0;
        virtual vecn< DIMENSION > mesh_element_barycenter(
            index_t mesh_element_index ) const = 0;
        vecn< DIMENSION > entity_barycenter() const
        {
            vecn< DIMENSION > result;
            for( auto v : range( nb_vertices() ) )
            {
                result += vertex( v );
            }
            ringmesh_assert( nb_vertices() > 0 );
            return result / static_cast< double >( nb_vertices() );
        }

        virtual double size() const
        {
            double size = 0.0;
            for( auto i : range( nb_mesh_elements() ) )
            {
                size += mesh_element_size( i );
            }
            return size;
        }

        /*! @}
         */
        GEO::AttributesManager& vertex_attribute_manager() const;

    protected:
        GeoModelMeshEntity( const GeoModel< DIMENSION >& geomodel, index_t id )
            : GeoModelEntity< DIMENSION >( geomodel, id )
        {
        }

        virtual void copy_mesh_entity(
            const GeoModelMeshEntity< DIMENSION >& from )
        {
            this->copy_name( from );
            this->id_ = from.index();
            boundaries_ = from.boundaries_;
            incident_entities_ = from.incident_entities_;
            parents_ = from.parents_;
        }

        virtual bool is_index_valid() const final;
        virtual bool is_mesh_valid() const
        {
            return mesh_ != nullptr;
        }

        void set_mesh( std::shared_ptr< MeshBase< DIMENSION > > mesh )
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
         * and all entities in the incident_entity must have this entity in
         * their
         * boundary vector
         */
        bool is_incident_entity_connectivity_valid() const;
        /*!
         * @brief Check that geomodel vertex indices are consistent
         * with what is stored at the GeoModel level.
         * @todo Review: this dependancy is bad, and puts us in a lot of trouble
         * [JP]
         */
        bool are_geomodel_vertex_indices_valid() const;

        void unbind_vertex_mapping_attribute() const;
        void bind_vertex_mapping_attribute() const;

        virtual void change_mesh_data_structure( const MeshType& type ) = 0;

    private:
        gmge_id defined_parent_gmge(
            const GeologicalEntityType& parent_type ) const;

        gmge_id could_be_undefined_parent_gmge(
            const GeologicalEntityType& parent_type ) const;

    protected:
        /// Boundary relations of this entity
        std::vector< index_t > boundaries_{};

        /// Incident-entity relations of this entity
        std::vector< index_t > incident_entities_{};

        /// Parents relations of this entity
        std::vector< index_t > parents_{};

    private:
        /// The RINGMesh::Mesh giving the geometry of this entity
        std::shared_ptr< MeshBase< DIMENSION > > mesh_{};
    };
    ALIAS_2D_AND_3D( GeoModelMeshEntity );

    /*!
     * @brief A GeoModelEntity of type CORNER
     * @details It is a unique point.
     */
    template < index_t DIMENSION >
    class geomodel_core_api Corner final
        : public GeoModelMeshEntity< DIMENSION >
    {
        ringmesh_disable_copy_and_move( Corner );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        friend class GeoModelMeshEntityAccess< DIMENSION >;
        friend class GeoModelMeshEntityConstAccess< DIMENSION >;

        virtual ~Corner()
        {
            this->unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static();

        MeshEntityType type_name() const final;

        /*!
         * @brief Checks if this entity define the geomodel external boundary
         * @details Test if the entity is in the Surfaces defining the universe
         */
        bool is_on_voi() const final;

        /*!
         * @brief For corners, the number of mesh elements is always equal to 1.
         */
        index_t nb_mesh_elements() const final;

        /*!
         * @return the number of vertices of the Corner
         */
        index_t nb_mesh_element_vertices(
            index_t mesh_element = 0 ) const final;

        const Line< DIMENSION >& incident_entity( index_t x ) const;

        /*! @}
         * \name Geometrical request on Corner
         * @{
         */
        double mesh_element_size( index_t mesh_element = 0 ) const final
        {
            ringmesh_unused( mesh_element );
            return 0.0;
        }

        double size() const final
        {
            return 0.0;
        }

        vecn< DIMENSION > mesh_element_barycenter(
            index_t mesh_element = 0 ) const final
        {
            ringmesh_unused( mesh_element );
            ringmesh_assert( this->nb_vertices() > 0 );
            return this->vertex( 0 );
        }

        /*!
         * @brief Get the low level mesh data structure
         */
        const PointSetMesh< DIMENSION >& mesh() const
        {
            return *point_set_mesh_;
        }

    protected:
        /*! @brief Creates a Corner.
         *  A point is added to its Mesh.
         */
        Corner( const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const MeshType& type );

        /*!
         * @brief Get the index of the unique vertex constituting of the Corner.
         * @return 0.
         */
        index_t mesh_element_vertex_index(
            const ElementLocalVertex& element_local_vertex ) const final
        {
            ringmesh_unused( element_local_vertex );
            return 0;
        }

        /*!
         * @brief Check that the Corner mesh is a unique point
         */
        bool is_mesh_valid() const final;

    private:
        void update_mesh_storage_type(
            std::unique_ptr< PointSetMesh< DIMENSION > > mesh );

        void change_mesh_data_structure( const MeshType& type ) final;

    private:
        std::shared_ptr< PointSetMesh< DIMENSION > > point_set_mesh_{};
    };
    ALIAS_2D_AND_3D( Corner );

    MeshEntityType geomodel_core_api corner_type_name_static();

    /*!
     * @brief A GeoModelEntity of type LINE
     *
     * @details This must be one connected component (one part) of
     * a 1-manifold (Line with no T intersections).
     */
    template < index_t DIMENSION >
    class geomodel_core_api Line final : public GeoModelMeshEntity< DIMENSION >
    {
        ringmesh_disable_copy_and_move( Line );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        friend class GeoModelMeshEntityAccess< DIMENSION >;

        virtual ~Line()
        {
            this->unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static();

        MeshEntityType type_name() const final;

        bool is_on_voi() const final;

        const Corner< DIMENSION >& boundary( index_t x ) const;

        const Surface< DIMENSION >& incident_entity( index_t x ) const
        {
            return static_cast< const Surface< DIMENSION >& >(
                GeoModelMeshEntity< DIMENSION >::incident_entity( x ) );
        }

        bool is_connectivity_valid() const final;

        const LineAABBTree< DIMENSION >& edge_aabb() const;

        /*!
         * @brief Return the NNSearch for the edges of the line
         * @details The barycenter of the edges is used.
         */
        const NNSearch< DIMENSION >& edge_nn_search() const;

        /*!
         * Get the number of edges of the Line
         */
        index_t nb_mesh_elements() const;

        /*!
         * @return The number of vertices per edge: 2.
         */
        index_t nb_mesh_element_vertices( index_t mesh_element = 0 ) const final
        {
            ringmesh_unused( mesh_element );
            return 2;
        }

        /*!
         * @brief Get the index of a vertex in the Line from the
         * @param vertex_index in a given edge @param edge_index.
         */
        index_t mesh_element_vertex_index(
            const ElementLocalVertex& element_local_vertex ) const final;

        /*!
         * @brief A Line is closed if its two extremities are identitical.
         */
        bool is_closed() const
        {
            ringmesh_assert( this->nb_boundaries() == 2 );
            return ( this->boundary_gmme( 0 ).is_defined() )
                   && ( this->boundary_gmme( 0 ) == this->boundary_gmme( 1 ) );
        }

        /*!
         * @brief Gets the length of an edge
         */
        double mesh_element_size( index_t edge_index ) const final;
        /*!
         * @brief Gets the barycenter of an edge
         */
        vecn< DIMENSION > mesh_element_barycenter(
            index_t edge_index ) const final;

        bool is_first_corner_first_vertex() const;

        /*!
         * @brief Get the low level mesh data structure
         */
        const LineMesh< DIMENSION >& mesh() const
        {
            return *line_mesh_;
        }

    protected:
        Line( const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const MeshType& type );

        /*!
         * @brief Check that the mesh of the Line is valid
         * @details Check that
         *  - the GEO::Mesh has more than 1 vertex - more than 1 edge - no
         * polygon - no cell.
         *  - global indices of vertices in the geomodel are in a valid range
         *  - each vertex is in 2 edges except extremities that are in 1 edge
         *
         * Does not check:
         *  - Self-intersection - I suppose there are no segment - segment
         * intersection (JP)
         *  - Duplicated edge - most probably ruled out with the duplicated
         * vertex test (JP)
         *  - Duplicated vertex (verified at GeoModel level)
         */
        bool is_mesh_valid() const final;

    private:
        void update_mesh_storage_type(
            std::unique_ptr< LineMesh< DIMENSION > > mesh );

        void change_mesh_data_structure( const MeshType& type ) final;

    private:
        std::shared_ptr< LineMesh< DIMENSION > > line_mesh_{};
    };
    ALIAS_2D_AND_3D( Line );

    MeshEntityType geomodel_core_api line_type_name_static();

    /*!
     * @brief A GeoModelEntity of type SURFACE
     *
     * @details One connected component (part) of a 2-manifold surface
     * (all edges of the polygons are in at most 2 polygons)
     */
    template < index_t DIMENSION >
    class geomodel_core_api SurfaceBase : public GeoModelMeshEntity< DIMENSION >
    {
        ringmesh_disable_copy_and_move( SurfaceBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        friend class GeoModelMeshEntityAccess< DIMENSION >;

        virtual ~SurfaceBase()
        {
            this->unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static();

        MeshEntityType type_name() const final;

        const Line< DIMENSION >& boundary( index_t x ) const;

        const Region< DIMENSION >& incident_entity( index_t x ) const;

        bool is_simplicial() const;

        const SurfaceAABBTree< DIMENSION >& polygon_aabb() const;

        /*!
         * @brief Return the NNSearch for the polygons of the surface
         * @details The barycenter of the polygons is used.
         */
        const NNSearch< DIMENSION >& polygon_nn_search() const;

        GEO::AttributesManager& polygon_attribute_manager() const;

        /*!
         * \name Accessors to Surface polygons, edges and vertices
         * @{
         */
        /*!
         * Number of polygons of the Surface.
         */
        index_t nb_mesh_elements() const final;

        /*!
         * Number of vertices of a polygon
         */
        index_t nb_mesh_element_vertices( index_t polygon_index ) const final;

        /*!
         * @brief Index of the vertex in the Surface
         * from its index in a polygon of the mesh.
         */
        index_t mesh_element_vertex_index(
            const ElementLocalVertex& element_local_vertex ) const final;

        /*!
         * @brief Gets the polygon adjacent along an edge of a polygon.
         * @param polygon_index in the polygon
         * @param edge_index in the edge
         * @note The edge index is assumed to be the index of the vertex
         * at which it is starting.
         */
        index_t polygon_adjacent_index(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*! @}
         * @{
         */
        /*!
         * @return Polygon barycenter.
         */
        vecn< DIMENSION > mesh_element_barycenter(
            index_t polygon_index ) const final;

        /*!
         * @return Area of a polygon.
         */
        double mesh_element_size( index_t polygon_index ) const final;

        /*!
         * @brief Get the low level mesh data structure
         */
        const SurfaceMesh< DIMENSION >& mesh() const
        {
            return *surface_mesh_;
        }

    protected:
        SurfaceBase( const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const MeshType& type );

        /*!
         * @brief Check that the mesh of the Surface is valid
         * @details Check that
         *  - the GEO::Mesh has more than 2 vertices, at least 1 polygon, no
         * cells.
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
        bool is_mesh_valid_base() const;

    private:
        void update_mesh_storage_type(
            std::unique_ptr< SurfaceMesh< DIMENSION > > mesh );

        void change_mesh_data_structure( const MeshType& type ) final;

    private:
        std::shared_ptr< SurfaceMesh< DIMENSION > > surface_mesh_{};
    };

    MeshEntityType geomodel_core_api surface_type_name_static();

    template < index_t DIMENSION >
    class geomodel_core_api Surface final : public SurfaceBase< DIMENSION >
    {
    };

    template <>
    class geomodel_core_api Surface< 2 > final : public SurfaceBase< 2 >
    {
        friend class GeoModelMeshEntityAccess< 2 >;

    private:
        Surface( const GeoModel2D& geomodel, index_t id, const MeshType type )
            : SurfaceBase< 2 >( geomodel, id, type )
        {
        }

        bool is_mesh_valid() const final;

    public:
        bool is_on_voi() const final;
        bool side( index_t i ) const
        {
            return sides_[i];
        }

        /*!
         * @brief Tells whether the Surface2D is meshed or not
         */
        bool is_meshed() const;

    private:
        /*! Additional information to store oriented boundary Lines
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_{};
    };

    template <>
    class geomodel_core_api Surface< 3 > final : public SurfaceBase< 3 >
    {
        friend class GeoModelMeshEntityAccess< 3 >;

    private:
        Surface( const GeoModel3D& geomodel, index_t id, const MeshType type )
            : SurfaceBase< 3 >( geomodel, id, type )
        {
        }

        bool is_mesh_valid() const final;

    public:
        bool is_on_voi() const final;
        const Region< 3 >& incident_entity( index_t x ) const;
    };
    ALIAS_2D_AND_3D( Surface );

    /*!
     * @brief A GeoModelEntity of type REGION
     *
     * @details A Region a volumetric connected component of the geomodel
     * defined by a set of surfaces.
     * The Region can be only defined by its boundary Surfaces.
     * Its volumetric mesh is optional.
     */
    template < index_t DIMENSION >
    class geomodel_core_api Region final
        : public GeoModelMeshEntity< DIMENSION >
    {
        ringmesh_disable_copy_and_move( Region );
        ringmesh_template_assert_3d( DIMENSION );

    public:
        friend class GeoModelMeshEntityAccess< DIMENSION >;

        virtual ~Region()
        {
            this->unbind_vertex_mapping_attribute();
        }

        static MeshEntityType type_name_static();

        MeshEntityType type_name() const final;

        bool is_on_voi() const final;

        const Surface< DIMENSION >& boundary( index_t x ) const;

        bool is_connectivity_valid() const final;

        bool is_meshed() const;

        bool is_simplicial() const;

        const VolumeAABBTree< DIMENSION >& cell_aabb() const;

        /*!
         * @brief Return the NNSearch for the cells of the region
         * @details The barycenter of the cells is used.
         */
        const NNSearch< DIMENSION >& cell_nn_search() const;

        GEO::AttributesManager& cell_attribute_manager() const;

        /*!
         * \name Accessors to Region cells, facets, edges and vertices
         * @{
         */

        /*!
         * Get the number of cells of the Region.
         */
        index_t nb_mesh_elements() const final;

        /*!
         * Get the number of vertex in the cell \param cell_index of the Region.
         */
        index_t nb_mesh_element_vertices( index_t cell_index ) const;

        /*!
         * @brief Index of a vertex in the Region from its index in a cell
         */
        index_t mesh_element_vertex_index(
            const ElementLocalVertex& element_local_vertex ) const final;

        /*!
         * Get the type of a given cell.
         * @pre The Region must be meshed
         */
        CellType cell_type( index_t cell_index ) const;

        index_t nb_cell_edges( index_t cell_index ) const;

        index_t nb_cell_facets( index_t cell_index ) const;

        index_t nb_cell_facet_vertices(
            index_t cell_index, index_t facet_index ) const;

        index_t cell_edge_vertex_index( index_t cell_index,
            index_t edge_index,
            index_t vertex_index ) const;

        index_t cell_facet_vertex_index( index_t cell_index,
            index_t facet_index,
            index_t vertex_index ) const;

        index_t cell_adjacent_index(
            index_t cell_index, index_t facet_index ) const;

        ElementLocalVertex find_cell_from_colocated_vertex_if_any(
            const vecn< DIMENSION >& vertex_vec ) const;

        /*! @}
         * \name Geometrical request on Region Entity
         * @{
         */

        /*!
         * @brief Volume of a cell
         */
        double mesh_element_size( index_t cell_index ) const final;

        /*!
         * @brief Compute the volume of the Region
         */
        double size() const final;

        /*!
         * @brief Get the center of the cell \param cell_index
         */
        vecn< DIMENSION > mesh_element_barycenter(
            index_t cell_index ) const final;

        std::vector< index_t > cells_around_vertex(
            index_t vertex_id, index_t cell_hint ) const;

        bool side( index_t i ) const
        {
            return sides_[i];
        }
        /*! @}
         */

        /*!
         * @brief Get the low level mesh data structure
         */
        const VolumeMesh< DIMENSION >& mesh() const
        {
            return *volume_mesh_;
        }

    private:
        Region( const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const MeshType& type );

        bool is_mesh_valid() const final;

        void update_mesh_storage_type(
            std::unique_ptr< VolumeMesh< DIMENSION > > mesh );

        void change_mesh_data_structure( const MeshType& type ) final;

        void copy_mesh_entity(
            const GeoModelMeshEntity< DIMENSION >& from ) final
        {
            const auto& region_from =
                dynamic_cast< const Region< DIMENSION >& >( from );
            GeoModelMeshEntity< DIMENSION >::copy_mesh_entity( from );
            sides_ = region_from.sides_;
        }

    protected:
        /*! Additional information to store oriented boundary Surfaces
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_{};

    private:
        std::shared_ptr< VolumeMesh< DIMENSION > > volume_mesh_{};
    };

    MeshEntityType geomodel_core_api region_type_name_static();

    ALIAS_3D( Region );

} // namespace RINGMesh
