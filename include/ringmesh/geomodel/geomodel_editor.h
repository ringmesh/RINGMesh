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

#ifndef __RINGMESH_GEOMODEL_EDITOR__
#define __RINGMESH_GEOMODEL_EDITOR__

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

/*!
 * @file Declaration of GeoModelEditor class.
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
     * @brief Basic edition of a GeoModel
     *        Topological level. NO GEOMETRY HERE.
     *        GeoModelBuilder deals with the dirty geometry
     *
     * @todo Separate the high level functions from the low-level ones 
     * that typically provides access to attributes of the GeoModel or Entities.
     * This class is a monster.
     */
    class RINGMESH_API GeoModelEditor {
    public:
        typedef std::string EntityType ;

        GeoModelEditor( GeoModel& M ) ;
        virtual ~GeoModelEditor() ;

        void copy_macro_topology( const GeoModel& from ) ;

        /*!
         *@brief Set the name of the geomodel
         */
        void set_geomodel_name( const std::string& name ) ;

        /*! @}
         * \name Creation - Deletion - Access to GeoModelEntities.
         * @{
         */
        template< typename T >
        gme_t create_mesh_entity( const MeshType type = "" )
        {
            const EntityType entity_type = T::type_name_static() ;
            index_t nb_entities( geomodel().nb_mesh_entities( entity_type ) ) ;
            index_t new_id( nb_entities ) ;
            T* new_entity = new T( geomodel(), new_id, type ) ;
            modifiable_mesh_entities( entity_type ).push_back( new_entity ) ;
            return new_entity->gme_id() ;

        }

        /*!
         * @brief Transfer general mesh information from one mesh data structure to another one
         * @param[in] id the GeoModelMeshEntity id to operate on
         * @param[in] type the new mesh data structure type
         */
        void change_mesh_data_structure( const gme_t& id, const MeshType type ) ;

        gme_t create_geological_entity( const EntityType& type ) ;

        // ---- Duplicate of protected functions. Dangerous open-bar access.
        GeoModelMeshEntity& mesh_entity( const gme_t& id )
        {
            return modifiable_mesh_entity( id ) ;
        }
        GeoModelMeshEntity& mesh_entity(
            const EntityType& entity_type,
            index_t entity_index )
        {
            return mesh_entity( gme_t( entity_type, entity_index ) ) ;
        }
        GeoModelGeologicalEntity& geological_entity( const gme_t& id )
        {
            return modifiable_geological_entity( id ) ;
        }
        GeoModelGeologicalEntity& geological_entity(
            const EntityType& entity_type,
            index_t entity_index )
        {
            return geological_entity( gme_t( entity_type, entity_index ) ) ;
        }

        /*! @}
         * \name Filling GeoModelEntity attributes.
         * @{
         */

        /*!
         * @brief Complete missing information in GeoModelEntities
         * boundaries - in_boundary - parent - children
         * @details For all 7 types of entities, check what information is available
         * for the first one and fill the entities of the same type accordingly
         * THIS MEANS that the all the entities of the same type have been initialized with
         * the same information
         */
        void complete_entity_connectivity() ;

        /*!
         * @brief Fill the boundaries of all entities of the given type
         * @details If the boundary entities do not have any in_boundary
         * information, nothing is done.
         */
        void fill_mesh_entities_boundaries( const EntityType& type ) ;

        /*!
         * @brief Fill the in_boundary vector of all entities of the given type
         * @details If the in_boundary entities do not have any boundary
         * information, nothing is done, and geomodel construction will eventually fail.
         */
        void fill_mesh_entities_in_boundaries( const EntityType& type ) ;

        /*!
         * @brief Fill the parent of all entities of the given type
         * @details If the parents do not have any child nothing is done.
         */
        void fill_mesh_entities_parent( const EntityType& type ) ;

        /*!
         * @brief Fill the children of all entities of the given type
         * @details If the children entities do not have any parent information
         * nothing is done.
         */
        void fill_geological_entities_children( const EntityType& type ) ;

        void complete_mesh_entities_geol_feature_from_first_parent(
            const EntityType& type ) ;
        void complete_geological_entities_geol_feature_from_first_child(
            const EntityType& type ) ;

        void set_entity_name( const gme_t& t, const std::string& name )
        {
            if( geomodel().is_mesh_entity_type( t.type ) ) {
                mesh_entity( t ).name_ = name ;
            } else {
                geological_entity( t ).name_ = name ;
            }
        }
        void set_entity_geol_feature( const gme_t& t, GME::GEOL_FEATURE geol )
        {
            if( geomodel().is_mesh_entity_type( t.type ) ) {
                mesh_entity( t ).geol_feature_ = geol ;
            } else {
                geological_entity( t ).geol_feature_ = geol ;
            }
        }

        void add_mesh_entity_boundary(
            const gme_t& t,
            index_t boundary_id,
            bool side = false )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const EntityType& b_type = geomodel().entity_type_manager().boundary_type(
                t.type ) ;
            gme_t boundary( b_type, boundary_id ) ;
            mesh.boundaries_.push_back( boundary ) ;

            if( t.type == Region::type_name_static() ) {
                dynamic_cast< Region& >( mesh ).sides_.push_back( side ) ;
            }
        }

        void set_mesh_entity_boundary(
            const gme_t& t,
            index_t id,
            index_t boundary_id,
            bool side = false )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const EntityType& b_type = geomodel().entity_type_manager().boundary_type(
                t.type ) ;
            gme_t boundary( b_type, boundary_id ) ;
            mesh.boundaries_[id] = boundary ;

            if( t.type == Region::type_name_static() ) {
                dynamic_cast< Region& >( mesh ).sides_[id] = side ;
            }
        }

        void add_universe_boundary( index_t boundary_id, bool side )
        {
            gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            geomodel().universe_.boundary_surfaces_.push_back( boundary ) ;
            geomodel().universe_.boundary_surface_sides_.push_back( side ) ;
        }

        void set_universe_boundary( index_t id, index_t boundary_id, bool side )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            ringmesh_assert( id < geomodel().universe_.nb_boundaries() ) ;
            gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            geomodel().universe_.boundary_surfaces_[id] = boundary ;
            geomodel().universe_.boundary_surface_sides_[id] = side ;
        }

        void add_mesh_entity_in_boundary( const gme_t& t, index_t in_boundary_id )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const EntityType& in_b_type =
                geomodel().entity_type_manager().in_boundary_type( t.type ) ;
            gme_t in_boundary( in_b_type, in_boundary_id ) ;
            mesh.in_boundary_.push_back( in_boundary ) ;
        }

        void set_mesh_entity_in_boundary(
            const gme_t& t,
            index_t id,
            index_t in_boundary_id )
        {
            /// No check on the validity of the index of the entity in_boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert( id < mesh.nb_in_boundary() ) ;
            const EntityType& in_b_type =
                geomodel().entity_type_manager().in_boundary_type( t.type ) ;
            gme_t in_boundary( in_b_type, in_boundary_id ) ;
            mesh.in_boundary_[id] = in_boundary ;
        }

        void add_mesh_entity_parent( const gme_t& t, const gme_t& parent_index )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            mesh.parents_.push_back( parent_index ) ;
        }

        void set_mesh_entity_parent(
            const gme_t& t,
            index_t id,
            const gme_t& parent_index )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& entity = mesh_entity( t ) ;
            ringmesh_assert( id < entity.nb_parents() ) ;
            entity.parents_[id] = parent_index ;
        }

        void add_geological_entity_child( const gme_t& t, index_t child_id )
        {
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            const EntityType& child_type = geomodel().entity_type_manager().child_type(
                t.type ) ;
            gme_t child( child_type, child_id ) ;
            entity.children_.push_back( child ) ;
        }

        void set_geological_entity_child(
            const gme_t& t,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            const EntityType& child_type = geomodel().entity_type_manager().child_type(
                t.type ) ;
            gme_t child( child_type, child_id ) ;
            entity.children_[id] = child ;
        }

        void remove_mesh_entities( const std::set< gme_t >& entities ) ;

        void remove_geological_entities( const std::set< gme_t >& entities ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_entities( std::set< gme_t >& entities ) const ;

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gme_t >& entities_to_remove ) ;

        const GeoModel& geomodel() const
        {
            return geomodel_ ;
        }

    protected:
        /*!
         * @ brief The geomodel under construction
         */
        GeoModel& geomodel()
        {
            return geomodel_ ;
        }

        bool create_mesh_entities(
            const EntityType& type,
            index_t nb_additional_entities ) ;

        template< typename ENTITY >
        bool create_mesh_entities(
            index_t nb_additionnal_entities,
            const MeshType type = "" )
        {
            const EntityType entity_type = ENTITY::type_name_static() ;
            std::vector< GeoModelMeshEntity* >& store = modifiable_mesh_entities(
                entity_type ) ;
            index_t old_size = static_cast< index_t >( store.size() ) ;
            index_t new_size = old_size + nb_additionnal_entities ;
            store.resize( new_size, nil ) ;
            for( index_t i = old_size; i < new_size; i++ ) {
                store[i] = new ENTITY( geomodel(), i, type ) ;
            }
            return true ;
        }
        bool create_geological_entities( const EntityType& type, index_t nb ) ;

        index_t find_or_create_geological_entity_type( const EntityType& type )
        {
            index_t type_index = geomodel_.geological_entity_type_index( type ) ;
            if( type_index == NO_ID ) {
                type_index = create_geological_entity_type( type ) ;
            }
            return type_index ;
        }
        index_t create_geological_entity_type( const EntityType& type ) ;

        GeoModelGeologicalEntity* create_geological_entity(
            const EntityType& type,
            index_t index_in_geomodel ) ;
        void delete_mesh_entity( const EntityType& type, index_t index )
        {
            std::vector< GeoModelMeshEntity* >& store = modifiable_mesh_entities(
                type ) ;
            delete store[index] ;
            store[index] = nil ;
        }
        void delete_geological_entity( const EntityType& type, index_t index )
        {
            std::vector< GeoModelGeologicalEntity* >& store =
                modifiable_geological_entities( type ) ;
            delete store[index] ;
            store[index] = nil ;
        }
        EntityTypeManager& entity_type_manager()
        {
            return geomodel().entity_type_manager_ ;
        }
        const EntityTypeManager& entity_type_manager() const
        {
            return geomodel().entity_type_manager_ ;
        }

        std::vector< GeoModelMeshEntity* >& modifiable_mesh_entities(
            const EntityType& type )
        {
            // Beurk again 
            return const_cast< std::vector< GeoModelMeshEntity* >& >( geomodel().mesh_entities(
                type ) ) ;
        }
        std::vector< GeoModelGeologicalEntity* >& modifiable_geological_entities(
            const EntityType& type )
        {
            return const_cast< std::vector< GeoModelGeologicalEntity* >& >( geomodel().geological_entities(
                type ) ) ;
        }
        GeoModelEntity& modifiable_entity( const gme_t& id )
        {
            if( geomodel().is_mesh_entity_type( id.type ) ) {
                return modifiable_mesh_entity( id ) ;
            } else if( geomodel().is_geological_entity_type( id.type ) ) {
                return modifiable_geological_entity( id ) ;
            } else {
                ringmesh_assert_not_reached;
                gme_t first_surface( Surface::type_name_static(), 0 ) ;
                return modifiable_mesh_entity( first_surface ) ;
            }
        }
        GeoModelMeshEntity& modifiable_mesh_entity( const gme_t& id )
        {
            return *modifiable_mesh_entities( id.type )[id.index] ;
        }
        GeoModelGeologicalEntity& modifiable_geological_entity( const gme_t& id )
        {
            return *modifiable_geological_entities( id.type )[id.index] ;
        }
        Universe& universe()
        {
            return geomodel_.universe_ ;
        }

        // The more I think about it the more this design
        // for editing the GeoModel and its Entities appears bad.
        // AM.
        void set_entity_index( GeoModelEntity& E, index_t new_index_in_geomodel )
        {
            E.id_.index = new_index_in_geomodel ;
        }
        void set_boundary_sign( Region& R, index_t boundary_index, bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() ) ;
            R.sides_[boundary_index] = new_side ;
        }

        std::vector< gme_t >& modifiable_children( GeoModelGeologicalEntity& E )
        {
            return E.children_ ;
        }
        std::vector< gme_t >& modifiable_boundaries( GeoModelMeshEntity& E )
        {
            return E.boundaries_ ;
        }
        std::vector< gme_t >& modifiable_boundaries( Universe& U )
        {
            return U.boundary_surfaces_ ;
        }
        std::vector< gme_t >& modifiable_in_boundaries( GeoModelMeshEntity& E )
        {
            return E.in_boundary_ ;
        }
        std::vector< gme_t >& modifiable_parents( GeoModelMeshEntity& E )
        {
            return E.parents_ ;
        }
        std::vector< bool >& modifiable_sides( Region& R )
        {
            return R.sides_ ;
        }
        std::vector< bool >& modifiable_sides( Universe& U )
        {
            return U.boundary_surface_sides_ ;
        }

    private:
        template< typename T >
        void complete_mesh_entity_connectivity() ;

        template< typename T >
        void copy_mesh_entity_topology( const GeoModel& from ) ;
        void copy_geological_entity_topology(
            const GeoModel& from,
            const EntityType& type ) ;

    private:
        /*! The geomodel edited
        */
        GeoModel& geomodel_ ;
        /*! Parameter to forbid element creation. Crucial to control
        *  building of the geomodel and detect errors in find_or_create functions
        */
        bool create_entity_allowed_ ;
    } ;

}

#endif
