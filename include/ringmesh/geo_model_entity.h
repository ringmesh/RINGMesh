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
 * @file Declaration of GeoModelEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#ifndef __RINGMESH_GEO_MODEL_ENTITY__
#define __RINGMESH_GEO_MODEL_ENTITY__

#include <ringmesh/common.h>

#include <string>
#include <vector>

#include <ringmesh/mesh.h>

namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {

    /*!
     * @brief Generic class describing one entity of a GeoModel
     * 
     * @todo We need to think once again about this design
     * to be able to ease the addition of optional elements to the GeoModel [JP]
     */
    class RINGMESH_API GeoModelEntity {
    ringmesh_disable_copy( GeoModelEntity ) ;
        friend class GeoModelEditor ;
    public:
        /*!
         * @brief Geological feature types for GeoModelEntity
         * @todo Read all possible geological features used in RESQML.
         */
        enum GEOL_FEATURE {
            /// All geological features 
            ALL_GEOL,
            /// Default value - No geological feature defined
            NO_GEOL,
            /// Stratigraphic surface - an horizon
            STRATI,
            /*!
             * Unconformity
             * @todo: Distinguish between as erosive, baselap or intrusive
             * unconformities --GC
             */
            UNCONFORMITY,
            /// A normal fault
            NORMAL_FAULT,
            /// A reverse fault 
            REVERSE_FAULT,
            /// An unspecified fault 
            FAULT,
            /// Volume Of Interest
            VOI
        } ;

        /*!
         * @brief Type of the GeoModelEntity
         * @details When no type is defined NO_TYPE should be used
         * There are two categories of entities
         *   - low-level entities (CORNER, LINE, SURFACE, REGION) which have a 
         *     geometry and connectivity relationships
         *   - high-level entities (CONTACT, INTERFACE, LAYER) 
         *     that are constituted of low-level entities
         * TYPE is used extensively to manage entities, iterate on them, etc.
         *  
         * @warning DO NOT MODIFY THIS ENUM SINCE ALL THE CODE RELIES ON IT. 
         * 
         * @todo Add fault blocks.
         */
        enum TYPE {
            /// Points at LINE extremities
            CORNER = 0,
            /// One connected component of the intersection of at least 2 SURFACE 
            LINE,
            /// One 2-manifold connected component
            SURFACE,
            /// One volumetric region defined by its boundary SURFACE and
            /// is optionally meshed
            REGION,
            /// A group of LINE, intersection of at least 2 INTERFACE
            CONTACT,
            /// A group of SURFACE, delimit maximum 2 LAYER 
            INTERFACE,
            /// A group of REGION 
            LAYER,
            /// Default TYPE
            NO_TYPE,
            /// Any type - to access generically entities
            ALL_TYPES
        } ;

        /*! 
         * @brief Unique identification of a GeoModelEntity in a GeoModel
         * @details Stores the TYPE of the entity and its index in the GeoModel.
         *          Default values are NO_TYPE and NO_ID
         * @todo Should we change this name? it looks like index_t but does not enforce
         *       the programming guidelines [JP]
         */
        struct gme_t {
            gme_t()
                : type( NO_TYPE ), index( NO_ID )
            {
            }
            gme_t( TYPE t, index_t id )
                : type( t ), index( id )
            {
            }
            bool operator!=( const gme_t& rhs ) const
            {
                return type != rhs.type || index != rhs.index ;
            }
            bool operator==( const gme_t& rhs ) const
            {
                return type == rhs.type && index == rhs.index ;
            }
            /*!
             * @brief Sort GME identifiers
             * @details Compare first types, then compare indices, 
             *          beginning with NO_ID indices. 
             * @note In a sorted vector v of gme_t one can find the first surface with
             *       std::lower_bound( v.begin(), v.end(), gme_t( SURFACE, NO_ID ) ) ;
             */
            bool operator<( const gme_t& rhs ) const
            {
                if( type != rhs.type ) {
                    return type < rhs.type ;
                } else {
                    if( index == NO_ID ) return true ;
                    if( rhs.index == NO_ID ) return false ;
                    return index < rhs.index ;
                }
            }
            friend std::ostream& operator<<( std::ostream& os, const gme_t& in )
            {
                os << GeoModelEntity::type_name( in.type ) << " " << in.index ;
                return os ;
            }

            bool is_defined() const
            {
                return type != NO_TYPE && type != ALL_TYPES && index != NO_ID ;
            }
            /*!
             * TYPE of the GeoModelEntity
             */
            TYPE type ;
            /*!
             * Index of the entity in the GeoModel
             */
            index_t index ;
        } ;

        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type(
            const std::vector< GEOL_FEATURE >& types ) ;
        static std::string geol_name( GEOL_FEATURE ) ;
        static bool is_fault( GEOL_FEATURE T )
        {
            return T == FAULT || T == REVERSE_FAULT || T == NORMAL_FAULT ;
        }
        static bool is_stratigraphic_limit( GEOL_FEATURE T )
        {
            return T == STRATI || T == UNCONFORMITY ;
        }

        /*!
         * \name Key functions to access relationships between TYPEs 
         * @{
         */
        static std::string type_name( TYPE t ) ;

        static TYPE parent_type( TYPE t ) ;
        static TYPE child_type( TYPE t ) ;
        static TYPE boundary_type( TYPE t ) ;
        static TYPE in_boundary_type( TYPE t ) ;
        static index_t dimension( TYPE t ) ;
        static bool has_mesh( TYPE t ) ;

        static bool parent_allowed( TYPE t )
        {
            return parent_type( t ) != NO_TYPE ;
        }
        static bool child_allowed( TYPE t )
        {
            return child_type( t ) != NO_TYPE ;
        }
        static bool boundary_allowed( TYPE t )
        {
            return boundary_type( t ) != NO_TYPE ;
        }
        static bool in_boundary_allowed( TYPE t )
        {
            return in_boundary_type( t ) != NO_TYPE ;
        }

        /*!@}
         */

        virtual ~GeoModelEntity()
        {
        }

        /*!@}
         * \name Validity checks
         * @{
         */

        /*!
         * @brief Global validity check of the GME
         */
        virtual bool is_valid() const
        {
            return is_connectivity_valid() ;
        }

        /*!
         * @brief Basic checks on the minimum required information 
         * @details Required connectivity information depends on the TYPE.   
         */
        bool is_connectivity_valid() const ;

        /*!
         * \name Accessors to basic information
         * @{
         */

        const GeoModel& model() const
        {
            return model_ ;
        }
        bool has_name() const
        {
            return !name().empty() ;
        }
        const std::string& name() const
        {
            return name_ ;
        }
        const gme_t& gme_id() const
        {
            return id_ ;
        }
        index_t index() const
        {
            return gme_id().index ;
        }
        TYPE type() const
        {
            return gme_id().type ;
        }
        bool has_geological_feature() const
        {
            return geol_feature_ != NO_GEOL ;
        }
        GEOL_FEATURE geological_feature() const
        {
            return geol_feature_ ;
        }
        bool is_on_voi() const ;

        /*!@}
         * \name Connectivity - boundary and in_boundary
         * @todo Change in_boundary to incident_entities? in_boundary is obscure? [JP]
         * @{
         */
        index_t nb_boundaries() const
        {
            return static_cast< index_t >( boundaries_.size() ) ;
        }
        const gme_t& boundary_gme( index_t x ) const
        {
            return boundaries_[x] ;
        }
        const GeoModelEntity& boundary( index_t x ) const ;

        index_t nb_in_boundary() const
        {
            return static_cast< index_t >( in_boundary_.size() ) ;
        }
        const gme_t& in_boundary_gme( index_t x ) const
        {
            return in_boundary_[x] ;
        }
        const GeoModelEntity& in_boundary( index_t x ) const ;

        bool is_inside_border( const GeoModelEntity& e ) const ;
        bool has_inside_border() const ;

        /*!@}
         * \name Parent - children relationships
         * @{
         */
        bool has_parent() const
        {
            return parent_id().is_defined() ;
        }
        const gme_t& parent_id() const
        {
            return parent_ ;
        }
        const GeoModelEntity& parent() const ;

        index_t nb_children() const
        {
            return static_cast< index_t >( children_.size() ) ;
        }
        const gme_t& child_id( index_t x ) const
        {
            return children_[x] ;
        }
        const GeoModelEntity& child( index_t x ) const ;

    protected:
        /*!
         * @brief Constructs a GeoModelEntity
         * Client code should only create GeoModelEntities through
         * GeoModelEditor derived classes.
         *
         * @param[in] model Constant reference to the parent model of this entity.
         * @param[in] entity_type Type of the entity to create
         * @param[in] id Index of the entity in the corresponding vector in the model
         * @param[in] name Name of the entity, empty by default.
         * @param[in] geological_feature Feature of the entity, none by default.
         */
        GeoModelEntity(
            const GeoModel& model,
            TYPE entity_type,
            index_t id,
            const std::string& name = "",
            GEOL_FEATURE geological_feature = NO_GEOL )
            :
                model_( model ),
                id_( entity_type, id ),
                name_( name ),
                geol_feature_( geological_feature )
        {
        }

    protected:
        /// Reference to the GeoModel owning this entity
        const GeoModel& model_ ;

        /// Unique identifier of the GeoModelEntity in the model
        gme_t id_ ;

        /// Name of the entity - default is an empty string
        std::string name_ ;

        /// Geological feature of this object - default is NO_GEOL
        GEOL_FEATURE geol_feature_ ;

        /// Entities on the boundary of this entity - see boundary_type( TYPE )
        std::vector< gme_t > boundaries_ ;

        /// Entities in which boundary this entity is - see in_boundary_type( TYPE )
        std::vector< gme_t > in_boundary_ ;

        /// Parent identification - see parent_type( TYPE )
        gme_t parent_ ;

        /// Entities constituting this one - see child_type( TYPE )
        std::vector< gme_t > children_ ;
    } ;

    typedef GeoModelEntity GME ;

    /*!
     * @brief Vertex in a GeoModelEntity
     */
    struct GMEVertex {
        GMEVertex( GME::gme_t t, index_t vertex_id_in )
            : gme_id( t ), v_id( vertex_id_in )
        {
        }
        GMEVertex()
            : gme_id(), v_id( NO_ID )
        {
        }
        bool operator<( const GMEVertex& rhs ) const
        {
            if( gme_id != rhs.gme_id ) {
                return gme_id < rhs.gme_id ;
            } else {
                return v_id < rhs.v_id ;
            }
        }
        bool operator==( const GMEVertex& rhs ) const
        {
            return gme_id == rhs.gme_id && v_id == rhs.v_id ;
        }
        bool is_defined() const
        {
            return gme_id.is_defined() && v_id != NO_ID ;
        }
        /// GeoModelEntity index in the GeoModel that owns it
        GME::gme_t gme_id ;
        /// Index of the vertex in the GeoModelEntity
        index_t v_id ;
    } ;

    /*!
     * @brief Abstract base class for GeoModelMeshEntity.
     * @details The GeoModelMeshEntity geometrical representation
     * is stored as a RINGMesh::Mesh. We defines generic access to
     * the RINGMesh::Mesh geometry. We also provide functions to link
     * the GeoModelMeshEntity with GeoModelMesh.
     */
    class RINGMESH_API GeoModelMeshEntity: public GeoModelEntity {
    ringmesh_disable_copy( GeoModelMeshEntity ) ;
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
        friend class GeoModelRepair ;

    public:
        GeoModelMeshEntity(
            const GeoModel& model,
            TYPE entity_type,
            index_t id,
            const std::string& name = "",
            GEOL_FEATURE geological_feature = NO_GEOL )
            :
                GeoModelEntity( model, entity_type, id, name, geological_feature ),
                mesh_( model, 3, false )
        {
            model_vertex_id_.bind( mesh_.vertex_attribute_manager(),
                model_vertex_id_att_name() ) ;
        }

        virtual ~GeoModelMeshEntity() ;

        /*! @todo To remove when GFX Mesh is encapsulated */
        const GEO::Mesh& gfx_mesh() const
        {
            return mesh_.gfx_mesh() ;
        }

        /*!
         * @brief Global validity of the entity
         */
        virtual bool is_valid() const
        {
            return is_connectivity_valid() && is_mesh_valid() ;
            /* @todo Test and add the model vertex validity test            
             * (no time right now JP)
             * are_model_vertex_indices_valid() ;
             */
        }
        void save(
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags ) const
        {
            mesh_.save_mesh( filename, ioflags ) ;
        }
        /*!
         * @brief Return the colocater for the Entity vertices.
         */
        const ColocaterANN& vertex_colocater_ann() const
        {
            return mesh_.colotater_ann( ColocaterANN::VERTICES ) ;
        }

        /*!
         * \name Local access to the GeoModelMeshEntity geometry
         * @{
         */
        /*!
         * @brief Number of vertices 
         */
        index_t nb_vertices() const
        {
            return mesh_.nb_vertices() ;
        }
        /*!
         * @brief Coordinates of the @param vertex_index.
         */
        const vec3& vertex( index_t vertex_index ) const
        {
            return mesh_.vertex( vertex_index ) ;
        }

        /*!
         * @brief Access to the vertex from the corner index.
         * @todo Review: this depends on the mesh implementation. To remove [JP].
         */
        virtual index_t vertex_index( index_t corner_index ) const = 0 ;

        /*!
         * @brief Get the number constitutive elements of the mesh 
         * @details Constitutive elements are those of the dimension of the object.
         * segments for lines, facets for surfaces and cells for volumetric meshes.
         */
        virtual index_t nb_mesh_elements() const = 0 ;
        /*!
         * @brief Number of vertices of a consitutive element of the mesh
         */
        virtual index_t nb_mesh_element_vertices( index_t mesh_element_index ) const = 0 ;
        /*!
         * @brief Convert the index in a mesh element to an index in the Entity.
         * @todo Review: I find this name confusing, maybe we should change it [JP]
         */
        virtual index_t mesh_element_vertex_index(
            index_t mesh_element_index,
            index_t vertex_index ) const = 0 ;
        /*!
         * @brief Coordinates of a vertex of a mesh element.
         */
        const vec3& mesh_element_vertex(
            index_t mesh_element_index,
            index_t vertex_index ) const
        {
            return vertex( mesh_element_vertex_index( mesh_element_index, vertex_index ) ) ;
        }
        
        /*! @}
         * \name Geometrical request on Entity
         * @{
         */
        virtual double mesh_element_size( index_t mesh_element_index ) const = 0 ;
        virtual double size() const
        {
            double size = 0.0 ;
            for( index_t i = 0; i < nb_mesh_elements(); ++i ) {
                size += mesh_element_size( i ) ;
            }
            return size ;
        }
        virtual vec3 mesh_element_center( index_t mesh_element_index ) const = 0 ;
        virtual vec3 center() const
        {
            vec3 result( 0., 0., 0. ) ;
            for( index_t v = 0; v < nb_vertices(); v++ ) {
                result += vertex( v ) ;
            }
            return result / static_cast< double >( nb_vertices() ) ;
        }

        /*! @}
         */
        /*!
         * \name Linking to GeoModelMesh indexing
         * @{
         */

        /*!
         * @brief Name of the attribute storing the global index at the GeoModel level.
         */
        static const std::string model_vertex_id_att_name() ;

        /*!
         * @brief Get the global index at the GeoModel from the index in the entity.
         */
        index_t model_vertex_id( index_t gmme_vertex_index = 0 ) const
        {
            ringmesh_assert( gmme_vertex_index < nb_vertices() ) ;
            return model_vertex_id_[gmme_vertex_index] ;
        }
        /*!
         * @brief Get the global index at the GeoModel from a 
         * index in a mesh element.
         */
        index_t model_vertex_id(
            index_t mesh_element_index,
            index_t vertex_index ) const
        {
            return model_vertex_id(
                mesh_element_vertex_index( mesh_element_index, vertex_index ) ) ;
        }
        /*!
         * @brief Get the first vertex index of the Entity 
         *  that corresponds to the vertex @param model_vertex_id.
         * @details Returns NO_ID if no matching point is found.
         *
         * @param model_vertex_id Index of the vertex in the GeoModel
         * @todo Change that name.
         */
        index_t gmme_vertex_index_from_model( index_t model_vertex_id ) const ;
        /*!
         * @ Get all the vertices of the Entity which corresponds
         * to the same point in the GeoModel.
         */
        std::vector< index_t > gme_vertex_indices( index_t model_vertex_id ) const ;

        /*!
         * @}
         * \name Attribute management
         * @{
         */
        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_.vertex_attribute_manager() ;
        }
        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh_.facet_attribute_manager() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_.cell_attribute_manager() ;
        }
        void bind_attributes() ;
        void unbind_attributes() ;
        /*! @}
         */
    protected:
        /*!
         * @brief Check if the mesh stored is valid.
         */
        virtual bool is_mesh_valid() const = 0 ;
        /*!
         * @brief Check that model vertex indinces are consistent 
         * with what is stored at the GeoModel level.
         * @todo Review: this dependancy is bad, and puts us in a lot of trouble [JP]
         */
        bool are_model_vertex_indices_valid() const ;

    protected:
        /// RINGMesh::Mesh implementing the geometry of this entity 
        Mesh mesh_ ;
        /*! Attribute on the Mesh vertices storing the index of
         *  the vertex in the GeoModel owning this entity
         */
        GEO::Attribute< index_t > model_vertex_id_ ;
    } ;

    /*!
     * @brief A GeoModelEntity of type CORNER
     * @details It is a unique point.
     */
    class RINGMESH_API Corner: public GeoModelMeshEntity {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        /*! @brief Creates a Corner. 
         *  A point is added to its Mesh.
         */
        Corner( const GeoModel& model, index_t id )
            : GeoModelMeshEntity( model, CORNER, id )
        {
            MeshBuilder builder( mesh_ ) ;
            builder.create_vertex() ;
        }

        ~Corner()
        {
        }

        virtual index_t vertex_index( index_t corner_index = 0 ) const {
            ringmesh_unused( corner_index ) ;
            return 0 ;
        }
        /*!
         * @brief Get the index of the unique vertex constituting of the Corner.
         * @return 0.
         */
        virtual index_t mesh_element_vertex_index(
            index_t /*mesh_element*/, index_t /*vertex_index*/ ) const 
        {
            return 0 ;
        }
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
        virtual index_t nb_mesh_element_vertices( index_t /*mesh_element*/ ) const
        {
            return 1 ;
        }

        /*! @}
         * \name Geometrical request on Corner
         * @{
         */
        virtual double mesh_element_size( index_t /*mesh_element*/ ) const
        {
            return 0.0 ;
        }
        virtual double size() const
        {
            return 0.0 ;
        }
        virtual vec3 mesh_element_center( index_t /*mesh_element*/ ) const
        {
            return vertex( 0 ) ;
        }

    protected:
        virtual bool is_mesh_valid() const ;
    } ;

    /*!
     * @brief A GeoModelEntity of type LINE
     *
     * @details This must be one connected component (one part) of 
     * a 1-manifold (Line with no T intersections).
     */
    class RINGMESH_API Line: public GeoModelMeshEntity {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        Line( const GeoModel& model, index_t id ) ;

        ~Line()
        {
        }

        virtual index_t vertex_index (index_t corner_index) const {
            return mesh_.edge_vertex(corner_index/2, corner_index%2) ;
        }
        /*!
         * Get the number of edges of the Line
         */
        virtual index_t nb_mesh_elements() const
        {
            return mesh_.nb_edges() ;
        }
        /*!
         * @Return The number of vertices per edge: 2.
         */
        virtual index_t nb_mesh_element_vertices( index_t /*mesh_element*/ ) const
        {
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
            return mesh_.edge_vertex( edge_index, vertex_index ) ;
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

        /*! @}
         * \name Geometrical request on Line
         * @{
         */
        /*!
         * @brief Gets the length of an edge
         */
        virtual double mesh_element_size( index_t edge_index ) const
        {
            ringmesh_assert( edge_index < nb_mesh_elements() ) ;
            return mesh_.edge_length( edge_index ) ;
        }
        /*!
         * @brief Gets the barycenter of an edge
         */
        virtual vec3 mesh_element_center( index_t edge_index ) const
        {
            ringmesh_assert( edge_index < nb_mesh_elements() ) ;
            return 0.5
                * ( mesh_element_vertex( edge_index, 0 )
                    + mesh_element_vertex( edge_index, 1 ) ) ;
        }

    private:
        virtual bool is_mesh_valid() const ;
    } ;

    /*!
     * @brief A GeoModelEntity of type SURFACE
     *
     * @details One connected component (part) of a 2-manifold surface 
     * (all edges of the facets are in at most 2 facets)
     */
    class RINGMESH_API Surface: public GeoModelMeshEntity {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    
    public:
        static const index_t NO_ADJACENT = index_t( -1 ) ;

        Surface( const GeoModel& model, index_t id )
            : GeoModelMeshEntity( model, SURFACE, id )
        {
        }

        ~Surface()
        {
        }

        bool is_simplicial() const
        {
            return mesh_.facets_are_simplicies() ;
        }
        const GEO::MeshFacetsAABB& facets_aabb() const
        {
            return mesh_.facets_aabb() ;
        }
        /*!
         * @brief Return the a colocater for the facets of the surface
         * @details The barycenter of the facets is used.
         */
        const ColocaterANN& facet_colocater_ann() const
        {
            return mesh_.colotater_ann( ColocaterANN::FACETS ) ;
        }

        /*!
         * \name Accessors to Surface facets, edges and vertices
         * @{
         */
        virtual index_t vertex_index( index_t corner_index ) const {
            return mesh_.facet_corner_vertex( corner_index ) ;
        }
        /*!
         * Number of facets of the Surface.
         */
        virtual index_t nb_mesh_elements() const
        {
            return mesh_.nb_facets() ;
        }
        /*!
         * Number of vertices of a facet 
         */
        virtual index_t nb_mesh_element_vertices( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh_.nb_facet_vertices( facet_index ) ;
        }
        /*!
         * @brief Index of the vertex in the Surface 
         * from its index in a facet of the mesh.
         * 
         * @todo Review Shouldn't we remove the virtual - These different functions are
         * not meant to be overriden. [JP]
         */
        virtual index_t mesh_element_vertex_index(
            index_t facet_index,
            index_t vertex_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( vertex_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh_.facet_vertex( facet_index, vertex_index ) ;
        }

        /*!
         * @brief Index of the first corner of a facet.
         * @todo Review: To remove? this intrinsically linked to the GEO::Mesh implementation [JP]       
         */
        index_t facet_begin( index_t facet_index ) const
        {
            return mesh_.facet_begin( facet_index ) ;
        }
        /*!
         * @brief Index of the last corner of a facet.
         * @todo Review: To remove? this intrinsically linked to the GEO::Mesh implementation [JP]
         */
        index_t facet_end( index_t facet_id ) const
        {
            return mesh_.facet_end( facet_id ) ;
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
            return mesh_.next_facet_vertex( facet_index, vertex_index ) ;
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
            return mesh_.prev_facet_vertex( facet_index, vertex_index ) ;
        }
        /*!
         * @brief Gets the facet adjacent along an edge of a facet.
         * @param edge_index in the facet
         */
        index_t facet_adjacent_index( index_t facet_index, index_t edge_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            ringmesh_assert( edge_index < nb_mesh_element_vertices( facet_index ) ) ;
            return mesh_.facet_adjacent( facet_index, edge_index ) ;
        }
        void next_on_border(
            index_t f,
            index_t from,
            index_t v,
            index_t& next_f,
            index_t& v_in_next,
            index_t& to ) const ;
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
            for( index_t v = 0; v < nb_mesh_element_vertices( facet_index ); v++ ) {
                if( mesh_element_vertex_index( facet_index, v )
                    == surface_vertex_index ) {
                    return v ;
                }
            }
            return NO_ID ;
        }
        index_t facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const ;

        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only ) const ;
        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet ) const ;

        /*! @}
         * \name Linking to GeoModelMesh indexing
         * @{
         */
        void edge_from_model_vertex_ids(
            index_t i0,
            index_t i1,
            index_t& f,
            index_t& e ) const ;
        void oriented_edge_from_model_vertex_ids(
            index_t i0,
            index_t i1,
            index_t& facet,
            index_t& edge ) const ;
        /*!
         * @brief Convert model vertex index to an index in a facet
         * @param[in] f Index of the facet
         * @param[in] model_v_id Index of the vertex in the GeoModel
         * @return NO_ID or index of the vertex in the facet
         */
        index_t facet_id_from_model( index_t f, index_t model_vertex_index ) const
        {
            for( index_t v = 0; v < nb_mesh_element_vertices( f ); v++ ) {
                if( model_vertex_id( f, v ) == model_vertex_index ) {
                    return v ;
                }
            }
            return NO_ID ;
        }
        index_t facet_from_model_vertex_ids( index_t i0, index_t i1 ) const ;

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
            return mesh_.facet_normal( facet_index ) ;
        }
        /*!
         * @return Facet barycenter.
         */
        virtual vec3 mesh_element_center( index_t facet_index ) const
        {            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh_.facet_barycenter( facet_index ) ;
        }
        bool facet_is_triangle( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh_.is_triangle( facet_index ) ;
        }

        /*!
         * @return Area of a facet.
         */
        virtual double mesh_element_size( index_t facet_index ) const
        {
            ringmesh_assert( facet_index < nb_mesh_elements() ) ;
            return mesh_.facet_area( facet_index ) ;
        }
        index_t closest_vertex_in_facet(
            index_t facet_index,
            const vec3& to_point ) const ;
        /*! 
         * Is the edge starting the given vertex of the facet on a border of the Surface?
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
            for( index_t v = 0; v < mesh_.nb_facet_vertices( facet_index ); v++ ) {
                if( is_on_border( facet_index, v ) ) {
                    return true ;
                }
            }
            return false ;
        }
        /*! @}
         */
    private:
        virtual bool is_mesh_valid() const ;
    } ;

    /*!
     * @brief A GeoModelEntity of type REGION
     *
     * @details A Region a volumetric connected component of the model defined
     * by a set of surfaces.
     * The Region can be only defined by its boundary Surfaces. 
     * Its volumetric mesh is optional.
     */
    class RINGMESH_API Region: public GeoModelMeshEntity {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    
    public:
        Region( const GeoModel& model, index_t id )
            : GeoModelMeshEntity( model, REGION, id )
        {
        }
        Region(
            const GeoModel& model,
            index_t id,
            const std::string& name,
            GEOL_FEATURE geological_feature )
            : GeoModelMeshEntity( model, REGION, id, name, geological_feature )
        {
        }
        ~Region()
        {
        }

        bool is_meshed() const
        {
            return mesh_.nb_cells() > 0 ;
        }
        bool is_simplicial() const
        {
            return mesh_.cells_are_simplicies() ;
        }
        const GEO::MeshCellsAABB& cells_aabb() const
        {
            return mesh_.cells_aabb() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_.cell_attribute_manager() ;
        }
        bool side( index_t i ) const
        {
            return sides_[i] ;
        }

        virtual index_t vertex_index (index_t corner_index) const {
            return mesh_.cell_corner_vertex(corner_index) ;
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
            return mesh_.nb_cells() ;
        }
        index_t nb_cell_edges( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh_.nb_cell_edges( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }
        index_t nb_cell_facets( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh_.nb_cell_facets( cell_index ) ;
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
                ringmesh_assert( facet_index < nb_cell_facets(cell_index) ) ;
                return mesh_.nb_cell_facet_vertices( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }
        /*!
         * Get the number of vertex in the cell \param cell_index of the Region.
         */
        virtual index_t nb_mesh_element_vertices( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh_.nb_cell_vertices( cell_index ) ;
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
                return mesh_.cell_vertex( cell_index, vertex_index ) ;
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
                return mesh_.cell_edge_vertex( cell_index, edge_index, vertex_index ) ;
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
                return mesh_.cell_facet_vertex( cell_index, facet_index, vertex_index ) ;
            }
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }
        index_t cell_adjacent_index( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh_.cell_adjacent( cell_index, facet_index ) ;
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
                return mesh_.cell_adjacent( cell_index, facet_index ) == GEO::NO_CELL ;
            }
            ringmesh_assert_not_reached ;
            return false ;
        }
        GEO::MeshCellType cell_type( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh_.cell_type( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return GEO::MESH_NB_CELL_TYPES ;
        }
        vec3 cell_facet_barycenter( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh_.cell_facet_barycenter( cell_index, facet_index ) ;
            }
            ringmesh_assert_not_reached ;
            return vec3() ;
        }
        vec3 cell_facet_normal( index_t cell_index, index_t facet_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                ringmesh_assert( facet_index < nb_cell_facets( cell_index ) ) ;
                return mesh_.cell_facet_normal( cell_index, facet_index ) ;
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
                return mesh_.cell_volume( cell_index ) ;
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
                    for( index_t v = 1; v + 1 < surface.nb_mesh_element_vertices( t );
                        ++v ) {
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
        virtual vec3 mesh_element_center( index_t cell_index ) const
        {
            if( is_meshed() ) {
                ringmesh_assert( cell_index < nb_mesh_elements() ) ;
                return mesh_.cell_barycenter( cell_index ) ;
            }
            ringmesh_assert_not_reached ;
            return vec3() ;
        }
        void compute_region_volumes_per_cell_type(
            double& tet_volume,
            double& pyramid_volume,
            double& prism_volume,
            double& hex_volume,
            double& poly_volume ) const ;

        /*! @}
         */

        /*!
         * \todo Is connectivity valid should be virtual and we should 
         * reimplement here to check 
         *     if( nb_boundaries() != sides_.size() ) {
         *       GEO::Logger::err( "GeoModelEntity" )
         *       << gme_id() << " boundary sides are invalid "
         *       << std::endl ;
         *       return false ;
         *     }

         */

    private:
        virtual bool is_mesh_valid() const ;

    private:
        /*! Additional information to store oriented boundary Surfaces
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_ ;
    } ;

} // namespace

#endif
