/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#ifndef __RINGMESH_GEO_MODEL_MESH_ENTITY__
#define __RINGMESH_GEO_MODEL_MESH_ENTITY__

#include <ringmesh/basic/common.h>

#include <string>
#include <vector>

#include <ringmesh/geomodel/geo_model_entity.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelGeologicalEntity ;

    /*!
     * @brief Abstract base class for GeoModelMeshEntity.
     * @details The GeoModelMeshEntity geometrical representation
     * is stored as a RINGMesh::Mesh. We defines generic access to
     * the RINGMesh::Mesh geometry.
     */
    class RINGMESH_API GeoModelMeshEntity: public GeoModelEntity {
    ringmesh_disable_copy( GeoModelMeshEntity ) ;
    public:
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
        friend class GeoModelRepair ;

        virtual ~GeoModelMeshEntity() ;

        typedef std::string EntityType ;

        static const EntityType default_entity_type_name()
        {
            return "No_entity_type" ;
        }

        virtual const EntityType type_name() const
        {
            return default_entity_type_name() ;
        }
        /*!
         * @brief Global validity of the entity
         */
        virtual bool is_valid() const
        {
            return is_identification_valid() && is_connectivity_valid()
                && is_mesh_valid() && are_geomodel_vertex_indices_valid() ;
        }
        virtual bool is_connectivity_valid() const ;

        /*!
         * \name Boundary relationship functions
         * @{
         */
        index_t nb_boundaries() const
        {
            return static_cast< index_t >( boundaries_.size() ) ;
        }
        const gme_t& boundary_gme( index_t x ) const
        {
            ringmesh_assert( x < nb_boundaries() ) ;
            return boundaries_[x] ;
        }
        const GeoModelMeshEntity& boundary( index_t x ) const ;

        index_t nb_in_boundary() const
        {
            return static_cast< index_t >( in_boundary_.size() ) ;
        }
        const gme_t& in_boundary_gme( index_t x ) const
        {
            ringmesh_assert( x < nb_in_boundary() ) ;
            return in_boundary_[x] ;
        }
        const GeoModelMeshEntity& in_boundary( index_t x ) const ;

        /*!
         * @brief Check if one entity is twice in the boundary
         */
        bool has_inside_border() const ;
        /*!
         * @brief Check if this entity is an inside border of rhs
         * @details That can be Surface stopping in a Region, or Line stopping in a Surface.
         * @param[in] rhs The entity to test
         */
        bool is_inside_border( const GeoModelMeshEntity& rhs ) const ;

        /*! @}
         * \name Parents/children relationship functions
         * @{
         */
        bool has_parent() const
        {
            return nb_parents() != 0 ;
        }
        /*!
         * @brief Check if the entity has a parent of the given type
         */
        bool has_parent( const EntityType& parent_type_name ) const
        {
            return parent_gme( parent_type_name ).is_defined() ;
        }
        index_t nb_parents() const
        {
            return static_cast< index_t >( parents_.size() ) ;
        }
        const gme_t& parent_gme( index_t id ) const
        {
            ringmesh_assert( id < nb_parents() ) ;
            return parents_[id] ;
        }
        /*!
         * @brief Returns the gme_t of the parent of the given type.
         * @note If this entity has no parent of the given type,
         * it will return an undefined gme_t (with no type and no id).
         * You should check on the returned gme_t.
         * @param[in] parent_type_name the asking parent type
         */
        const gme_t& parent_gme( const EntityType& parent_type_name ) const ;
        const GeoModelGeologicalEntity& parent( index_t id ) const ;
        const GeoModelGeologicalEntity& parent(
            const EntityType& parent_type_name ) const ;

        /*!
         * @}
         */

        /*! @todo To remove when GFX Mesh is encapsulated */
        const GEO::Mesh& gfx_mesh() const
        {
            return mesh_->gfx_mesh() ;
        }

        void save( const std::string& filename ) const
        {
            mesh_->save_mesh( filename ) ;
        }
        /*!
         * @brief Return the colocater for the Entity vertices.
         */

        const ColocaterANN& vertex_colocater_ann() const
        {
            return mesh_->vertices_colocater_ann() ;
        }

        /*!
         * \name Local access to the GeoModelMeshEntity geometry
         * @{
         */
        index_t nb_vertices() const
        {
            return mesh_->nb_vertices() ;
        }
        /*!
         * @brief Coordinates of the \p vertex_index.
         */
        const vec3& vertex( index_t vertex_index ) const
        {
            return mesh_->vertex( vertex_index ) ;
        }

        /*!
         * @brief Get the number constitutive elements of the mesh
         * @details Constitutive elements are those of the dimension of the object.
         * segments for lines, facets for surfaces and cells for volumetric meshes.
         */
        virtual index_t nb_mesh_elements() const = 0 ;
        /*!
         * @brief Number of vertices of a constitutive element of the mesh
         */
        virtual index_t nb_mesh_element_vertices(
            index_t mesh_element_index ) const = 0 ;
        /*!
         * @brief Convert the index in a mesh element to an index in the Entity.
         * @param[in] mesh_element_index Index of a constitutive element of the mesh
         * @param[in] vertex_local_index Local index of a vertex in the constitutive
         * mesh element
         * @return the global index of the vertex in the GeoModelMeshEntity
         */
        virtual index_t mesh_element_vertex_index(
            index_t mesh_element_index,
            index_t vertex_local_index ) const = 0 ;
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
                mesh_element_vertex_index( mesh_element_index, vertex_local_index ) ) ;
        }

        /*! @}
         * \name Geometrical request on Entity
         * @{
         */
        virtual double mesh_element_size( index_t mesh_element_index ) const = 0 ;
        virtual vec3 mesh_element_barycenter(
            index_t mesh_element_index ) const = 0 ;
        virtual vec3 entity_barycenter() const
        {
            vec3 result( 0., 0., 0. ) ;
            for( index_t v = 0; v < nb_vertices(); v++ ) {
                result += vertex( v ) ;
            }
            ringmesh_assert( nb_vertices() > 0 ) ;
            return result / static_cast< double >( nb_vertices() ) ;
        }
        virtual double size() const
        {
            double size = 0.0 ;
            for( index_t i = 0; i < nb_mesh_elements(); ++i ) {
                size += mesh_element_size( i ) ;
            }
            return size ;
        }

        /*! @}
         */
        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_->vertex_attribute_manager() ;
        }

    protected:
        GeoModelMeshEntity(
            const GeoModel& geomodel,
            index_t id,
            const std::string& name = "No_name",
            GEOL_FEATURE geological_feature = NO_GEOL )
            : GeoModelEntity( geomodel, id, name, geological_feature ), mesh_( NULL )
        {
        }
        virtual void copy( const GeoModelEntity& from )
        {
            GME::copy( from ) ;
            const GeoModelMeshEntity& mesh_entity_from =
                dynamic_cast< const GeoModelMeshEntity& >( from ) ;
            boundaries_ = mesh_entity_from.boundaries_ ;
            in_boundary_ = mesh_entity_from.in_boundary_ ;
            parents_ = mesh_entity_from.parents_ ;
        }
        virtual bool is_index_valid() const ;
        virtual bool is_mesh_valid() const
        {
            return mesh_ != NULL ;
        }

        void set_mesh( MeshBase* mesh )
        {
            ringmesh_assert( mesh != NULL ) ;
            mesh_ = mesh ;
        }

        bool is_boundary_connectivity_valid() const ;
        bool is_in_boundary_connectivity_valid() const ;
        bool is_parent_connectivity_valid() const ;
        /*!
         * @brief Check that geomodel vertex indices are consistent
         * with what is stored at the GeoModel level.
         * @todo Review: this dependancy is bad, and puts us in a lot of trouble [JP]
         */
        bool are_geomodel_vertex_indices_valid() const ;

    protected:

        /// Entities on the boundary of this entity
        std::vector< gme_t > boundaries_ ;

        /// Entities in which boundary this entity is
        std::vector< gme_t > in_boundary_ ;

        /// The optional GeoModelGeologicalEntities 
        /// (groups of GeoModelMeshEntity this entity belongs to)
        std::vector< gme_t > parents_ ;
    private:
        /// The RINGMesh::Mesh giving the geometry of this entity
        MeshBase* mesh_ ;
    } ;

    /*!
     * @brief A GeoModelEntity of type CORNER
     * @details It is a unique point.
     */
    class RINGMESH_API Corner: public GeoModelMeshEntity {
    public:
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;

        virtual ~Corner()
        {
        }

        static const EntityType type_name_static()
        {
            return "Corner" ;
        }

        virtual const EntityType type_name() const
        {
            return type_name_static() ;
        }

        virtual bool is_on_voi() const ;

        /*!
         * @return 0, no mesh_element are defined for corners.
         */
        virtual index_t nb_mesh_elements() const
        {
            return 0 ;
        }
        /*!
         * @return 1 the number of vertices of the Corner
         */
        virtual index_t nb_mesh_element_vertices( index_t mesh_element = 0 ) const
        {
            ringmesh_unused( mesh_element ) ;
            return 1 ;
        }

        /*! @}
         * \name Geometrical request on Corner
         * @{
         */
        virtual double mesh_element_size( index_t mesh_element = 0 ) const
        {
            ringmesh_unused( mesh_element ) ;
            return 0.0 ;
        }
        virtual double size() const
        {
            return 0.0 ;
        }
        virtual vec3 mesh_element_barycenter( index_t mesh_element = 0 ) const
        {
            ringmesh_unused( mesh_element ) ;
            return vertex( 0 ) ;
        }

        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        Mesh0D& low_level_mesh_storage()
        {
            return *mesh0d_ ;
        }
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const Mesh0D& low_level_mesh_storage() const
        {
            return *mesh0d_ ;
        }
    protected:
        /*! @brief Creates a Corner.
         *  A point is added to its Mesh.
         */
        Corner( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( Mesh0D::create_mesh( type ) ) ;
            Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder( *mesh0d_ ) ;
            builder->create_vertex() ;

            id_.type = type_name_static() ;
        }

        /*!
         * @brief Get the index of the unique vertex constituting of the Corner.
         * @return 0.
         */
        virtual index_t mesh_element_vertex_index(
            index_t mesh_element = 0,
            index_t vertex_index = 0 ) const
        {
            ringmesh_unused( mesh_element ) ;
            ringmesh_unused( vertex_index ) ;
            return 0 ;
        }

        virtual bool is_mesh_valid() const ;

        void update_mesh_storage_type( Mesh0D* mesh )
        {
            mesh0d_ = mesh ;
            GeoModelMeshEntity::set_mesh( mesh0d_ ) ;
        }


    private:
        Mesh0D* mesh0d_ ;
    } ;

    /*!
     * @brief A GeoModelEntity of type LINE
     *
     * @details This must be one connected component (one part) of
     * a 1-manifold (Line with no T intersections).
     */
    class RINGMESH_API Line: public GeoModelMeshEntity {
    public:
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
        friend class GeoModelRepair ;

        virtual ~Line()
        {
        }

        static const EntityType type_name_static()
        {
            return "Line" ;
        }

        virtual const EntityType type_name() const
        {
            return type_name_static() ;
        }

        virtual bool is_on_voi() const ;

        virtual bool is_connectivity_valid() const ;

        /*!
         * @brief Return the colocater for the edges of the line
         * @details The barycenter of the edges is used.
         */
        const ColocaterANN& edge_colocater_ann() const
        {
            return mesh1d_->edges_colocater_ann() ;
        }

        /*!
         * Get the number of edges of the Line
         */
        virtual index_t nb_mesh_elements() const
        {
            return mesh1d_->nb_edges() ;
        }

        /*!
         * @return The number of vertices per edge: 2.
         */
        virtual index_t nb_mesh_element_vertices( index_t mesh_element = 0 ) const
        {
            ringmesh_unused( mesh_element ) ;
            return 2 ;
        }

        /*!
         * @brief Get the index of a vertex in the Line from the
         * @param vertex_index in a given edge @param edge_index.
         */
        virtual index_t mesh_element_vertex_index(
            index_t edge_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( edge_index < nb_mesh_elements() ) ;
            ringmesh_assert( vertex_index < 2 ) ;
            return mesh1d_->edge_vertex( edge_index, vertex_index ) ;
        }

        /*!
         * @brief A Line is closed if its two extremities are identitical.
         */
        bool is_closed() const
        {
            ringmesh_assert( nb_boundaries() == 2 ) ;
            return ( boundaries_[0].is_defined() )
                && ( boundaries_[0] == boundaries_[1] ) ;
        }

        /*!
         * @brief Gets the length of an edge
         */
        virtual double mesh_element_size( index_t edge_index ) const
        {
            ringmesh_assert( edge_index < nb_mesh_elements() ) ;
            return mesh1d_->edge_length( edge_index ) ;
        }

        /*!
         * @brief Gets the barycenter of an edge
         */
        virtual vec3 mesh_element_barycenter( index_t edge_index ) const
        {
            ringmesh_assert( edge_index < nb_mesh_elements() ) ;
            return 0.5
                * ( mesh_element_vertex( edge_index, 0 )
                    + mesh_element_vertex( edge_index, 1 ) ) ;
        }

        bool is_first_corner_first_vertex() const ;
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        Mesh1D& low_level_mesh_storage()
        {
            return *mesh1d_ ;
        }
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const Mesh1D& low_level_mesh_storage() const
        {
            return *mesh1d_ ;
        }
    protected:
        Line( const GeoModel& geomodel, index_t id, const MeshType type ) ;

        virtual bool is_mesh_valid() const ;

        void update_mesh_storage_type( Mesh1D* mesh )
        {
            mesh1d_ = mesh ;
            GeoModelMeshEntity::set_mesh( mesh1d_ ) ;
        }

    private:
        Mesh1D* mesh1d_ ;
    } ;

    /*!
     * @brief A GeoModelEntity of type SURFACE
     *
     * @details One connected component (part) of a 2-manifold surface
     * (all edges of the facets are in at most 2 facets)
     */
    class RINGMESH_API Surface: public GeoModelMeshEntity {
    public:
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
        friend class GeoModelRepair ;

        virtual ~Surface()
        {
        }

        virtual const EntityType type_name() const
        {
            return type_name_static() ;
        }

        static const EntityType type_name_static()
        {
            return "Surface" ;
        }

        virtual bool is_on_voi() const ;

        bool is_simplicial() const
        {
            return mesh2d_->facets_are_simplicies() ;
        }

        const AABBTree2D& facets_aabb() const
        {
            return mesh2d_->facets_aabb() ;
        }

        /*!
         * @brief Return the colocater for the facets of the surface
         * @details The barycenter of the facets is used.
         */
        const ColocaterANN& facet_colocater_ann() const
        {
            return mesh2d_->facets_colocater_ann() ;
        }

        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh2d_->facet_attribute_manager() ;
        }

        /*!
         * \name Accessors to Surface facets, edges and vertices
         * @{
         */
        /*!
         * Number of facets of the Surface.
         */
        virtual index_t nb_mesh_elements() const
        {
            return mesh2d_->nb_facets() ;
        }

        /*!
         * Number of vertices of a facet
         */
        virtual index_t nb_mesh_element_vertices( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh2d_->nb_facet_vertices( facet_index ) ;
        }

        /*!
         * @brief Index of the vertex in the Surface
         * from its index in a facet of the mesh.
         */
        virtual index_t mesh_element_vertex_index(
            index_t facet_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( vertex_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh2d_->facet_vertex( facet_index, vertex_index ) ;
        }

        /*!
         * @brief Gets the next vertex index in a facet.
         */
        index_t next_facet_vertex_index(
            index_t facet_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( vertex_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh2d_->next_facet_vertex( facet_index, vertex_index ) ;
        }

        /*!
         * @brief Gets the previous vertex index in a facet.
         */
        index_t prev_facet_vertex_index(
            index_t facet_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( vertex_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh2d_->prev_facet_vertex( facet_index, vertex_index ) ;
        }

        /*!
         * @brief Gets the facet adjacent along an edge of a facet.
         * @param edge_index in the facet
         * @note The edge index is assumed to be the index of the vertex
         * at which it is starting.
         */
        index_t facet_adjacent_index( index_t facet_index, index_t edge_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( edge_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh2d_->facet_adjacent( facet_index, edge_index ) ;
        }

        /*!
         * @brief Get the previous edge on the border
         * @details The returned border edge is the previous in the way of facet edges
         * orientation.
         * @param[in] f Input facet index
         * @param[in] e Edge index in the facet
         * @param[out] prev_f Previous facet index
         * @param[out] prev_e Previous edge index in the facet
         *
         * @pre the surface must be correctly oriented and
         * the given facet edge must be on border
         */
        void prev_on_border(
            index_t f,
            index_t e,
            index_t& prev_f,
            index_t& prev_e ) const ;

        /*!
         * @brief Get the next edge on the border
         * @details The returned border edge is the next in the way of facet edges
         * orientation.
         * @param[in] f Input facet index
         * @param[in] e Edge index in the facet
         * @param[out] next_f Next facet index
         * @param[out] next_e Next edge index in the facet
         *
         * @pre the given facet edge must be on border
         */
        void next_on_border(
            index_t f,
            index_t e,
            index_t& next_f,
            index_t& next_e ) const ;

        /*!
         * @brief Get the vertex index in a facet @param facet_index from its
         * index in the Surface @param surface_vertex_index
         * @return NO_ID or index of the vertex in the facet
         */
        index_t vertex_index_in_facet(
            index_t facet_index,
            index_t surface_vertex_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            for( index_t v = 0; v < nb_mesh_element_vertices( facet_index ); v++ ) {
                if( mesh_element_vertex_index( facet_index, v )
                    == surface_vertex_index ) {
                    return v ;
                }
            }
            return NO_ID ;
        }
        index_t facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const ;

        /*!
         * @brief Determines the facets around a vertex
         * @param[in] surf_vertex_id Index of the vertex in the surface
         * @param[in] result Indices of the facets containing @param P
         * @param[in] border_only If true only facets on the border are considered
         * @param[in] f0 (Optional) Index of one facet containing the vertex @param P
         * @return The number of facets found
         * @note If a facet containing the vertex is given, facets around this
         * vertex is search by propagation. Else, a first facet is found by brute
         * force algorithm, and then the other by propagation
         * @todo Try to use a AABB tree to remove @param first_facet. [PA]
         */
        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet = NO_ID ) const ;

        /*! @}
         * \name Geometrical request on facets
         * @{
         */
        /*!
         * @return Normal to the facet
         */
        vec3 facet_normal( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh2d_->facet_normal( facet_index ) ;
        }

        /*!
         * @return Facet barycenter.
         */
        virtual vec3 mesh_element_barycenter( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh2d_->facet_barycenter( facet_index ) ;
        }

        /*!
         * @return Area of a facet.
         */
        virtual double mesh_element_size( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh2d_->facet_area( facet_index ) ;
        }

        index_t closest_vertex_in_facet(
            index_t facet_index,
            const vec3& to_point ) const ;

        /*!
         * Is the edge starting with the given vertex of the facet on a border of the Surface?
         */
        bool is_on_border( index_t facet_index, index_t vertex_index ) const
        {
            return facet_adjacent_index( facet_index, vertex_index ) == NO_ID ;
        }

        /*!
         * Is one of the edges of the facet on the border of the surface?
         */
        bool is_on_border( index_t facet_index ) const
        {
            for( index_t v = 0; v < mesh2d_->nb_facet_vertices( facet_index );
                v++ ) {
                if( is_on_border( facet_index, v ) ) {
                    return true ;
                }
            }
            return false ;
        }
        /*! @}
         */
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        Mesh2D& low_level_mesh_storage()
        {
            return *mesh2d_ ;
        }
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const Mesh2D& low_level_mesh_storage() const
        {
            return *mesh2d_ ;
        }
    protected:
        Surface( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( Mesh2D::create_mesh( type ) ) ;
            id_.type = type_name_static() ;
        }

        virtual bool is_mesh_valid() const ;

        void update_mesh_storage_type( Mesh2D* mesh )
        {
            mesh2d_ = mesh ;
            GeoModelMeshEntity::set_mesh( mesh2d_ ) ;
        }
    private:
        Mesh2D* mesh2d_ ;
    } ;

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
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
        friend class GeoModelRepair ;

        virtual ~Region()
        {
        }

        static const EntityType type_name_static()
        {
            return "Region" ;
        }

        virtual const EntityType type_name() const
        {
            return type_name_static() ;
        }

        virtual bool is_on_voi() const ;

        virtual bool is_connectivity_valid() const ;

        bool is_meshed() const
        {
            return mesh3d_->nb_cells() > 0 ;
        }

        bool is_simplicial() const
        {
            return mesh3d_->cells_are_simplicies() ;
        }

        const AABBTree3D& cells_aabb() const
        {
            return mesh3d_->cells_aabb() ;
        }

        /*!
         * @brief Return the colocater for the cells of the region
         * @details The barycenter of the cells is used.
         */
        const ColocaterANN& cell_colocater_ann() const
        {
            return mesh3d_->cells_colocater_ann() ;
        }

        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh3d_->cell_attribute_manager() ;
        }

        /*!
         * \name Accessors to Region cells, facets, edges and vertices
         * @{
         */

        /*!
         * Get the number of cells of the Region.
         */
        virtual index_t nb_mesh_elements() const
        {
            return mesh3d_->nb_cells() ;
        }

        /*!
         * Get the number of vertex in the cell \param cell_index of the Region.
         */
        virtual index_t nb_mesh_element_vertices( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->nb_cell_vertices( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        /*!
         * @brief Index of a vertex in the Region from its index in a cell
         */
        virtual index_t mesh_element_vertex_index(
            index_t cell_index,
            index_t vertex_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( vertex_index < nb_mesh_element_vertices( cell_index ) ) ;
                return mesh3d_->cell_vertex( cell_index, vertex_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        /*!
         * Get the type of a given cell.
         * @pre The Region must be meshed
         */
        GEO::MeshCellType cell_type( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->cell_type( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return GEO::MESH_NB_CELL_TYPES ;
        }

        index_t nb_cell_edges( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->nb_cell_edges( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        index_t nb_cell_facets( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->nb_cell_facets( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        index_t nb_cell_facet_vertices(
            index_t cell_index,
            index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh3d_->nb_cell_facet_vertices( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        index_t cell_edge_vertex_index(
            index_t cell_index,
            index_t edge_index,
            index_t vertex_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( edge_index < nb_cell_edges( cell_index ) ) ;
                ringmesh_assert( vertex_index < nb_mesh_element_vertices( cell_index ) ) ;
                return mesh3d_->cell_edge_vertex( cell_index, edge_index,
                    vertex_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        index_t cell_facet_vertex_index(
            index_t cell_index,
            index_t facet_index,
            index_t vertex_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                ringmesh_assert( vertex_index < nb_mesh_element_vertices( cell_index ) ) ;
                return mesh3d_->cell_facet_vertex( cell_index, facet_index,
                    vertex_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        index_t cell_adjacent_index( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh3d_->cell_adjacent( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        /*! @}
         * \name Geometrical request on Region Entity
         * @{
         */
        bool is_cell_facet_on_border( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh3d_->cell_adjacent( cell_index, facet_index )
                    == GEO::NO_CELL ;
            }
            ringmesh_assert_not_reached ;
            return false ;
        }

        vec3 cell_facet_barycenter( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh3d_->cell_facet_barycenter( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return vec3() ;
        }

        vec3 cell_facet_normal( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh3d_->cell_facet_normal( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return vec3() ;
        }

        /*!
         * @brief Volume of a cell
         */
        virtual double mesh_element_size( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->cell_volume( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return 0 ;
        }
        /*!
         * @brief Compute the volume of the Region
         */
        virtual double size() const
        {
            double result = 0. ;
            for( index_t i = 0; i < nb_boundaries(); i++ ) {
                const Surface& surface = dynamic_cast< const Surface& >( boundary(
                    i ) ) ;

                for( index_t t = 0; t < surface.nb_mesh_elements(); t++ ) {
                    const vec3& p0 = surface.mesh_element_vertex( t, 0 ) ;
                    for( index_t v = 1;
                        v + 1 < surface.nb_mesh_element_vertices( t ); ++v ) {
                        double cur_volume = ( dot( p0,
                            cross( surface.mesh_element_vertex( t, v ),
                                surface.mesh_element_vertex( t, v + 1 ) ) ) )
                            / static_cast< double >( 6 ) ;
                        side( i ) ? result -= cur_volume : result += cur_volume ;
                    }
                }
            }
            return fabs( result ) ;
        }
        /*!
         * @brief Get the center of the cell \param cell_index
         */
        virtual vec3 mesh_element_barycenter( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh3d_->cell_barycenter( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return vec3() ;
        }

        index_t cells_around_vertex(
            index_t vertex_id,
            std::vector< index_t >& result,
            index_t cell_hint ) const ;

        void compute_region_volumes_per_cell_type(
            double& tet_volume,
            double& pyramid_volume,
            double& prism_volume,
            double& hex_volume,
            double& poly_volume ) const ;

        bool side( index_t i ) const
        {
            return sides_[i] ;
        }
        /*! @}
         */
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        Mesh3D& low_level_mesh_storage()
        {
            return *mesh3d_ ;
        }
        /*!
         * @brief Get the low level mesh data structure
         * @warn This function is for ADVANCED user only. If you use it,
         * you are responsible for low level mesh consistency.
         */
        const Mesh3D& low_level_mesh_storage() const
        {
            return *mesh3d_ ;
        }
    protected:
        Region( const GeoModel& geomodel, index_t id, const MeshType type )
            : GeoModelMeshEntity( geomodel, id )
        {
            update_mesh_storage_type( Mesh3D::create_mesh( type ) ) ;
            id_.type = type_name_static() ;
        }

        void copy( const GeoModelEntity& from )
        {
            const Region& region_from = dynamic_cast< const Region& >( from ) ;
            GeoModelMeshEntity::copy( from ) ;
            sides_ = region_from.sides_ ;
        }

        virtual bool is_mesh_valid() const ;

        void update_mesh_storage_type( Mesh3D* mesh )
        {
            mesh3d_ = mesh ;
            GeoModelMeshEntity::set_mesh( mesh3d_ ) ;
        }

    protected:
        /*! Additional information to store oriented boundary Surfaces
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_ ;
    private:
        Mesh3D* mesh3d_ ;
    } ;
}

#endif
