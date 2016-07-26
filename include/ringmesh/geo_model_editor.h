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

#ifndef __RINGMESH_GEO_MODEL_EDITOR__
#define __RINGMESH_GEO_MODEL_EDITOR__

#include <ringmesh/common.h>

#include <set>

#include <ringmesh/geo_model_entity.h> 
#include <ringmesh/geo_model_geological_entity.h> 
#include <ringmesh/geo_model.h>

/*!
 * @file Declaration of GeoModelEditor class.
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
     * @brief Basic edition of a GeoModel
     *        Topological level. NO GEOMETRY HERE.
     *        GeoModelBuilder deals with the dirty geometry
     * @todo Force editor to use EntityRelationships
     */
    class RINGMESH_API GeoModelEditor {
    public:
        GeoModelEditor( GeoModel& M ) ;
        virtual ~GeoModelEditor() ;

        void copy_macro_topology( const GeoModel& from ) ;    

        /*!
         *@brief Set the name of the model
         */
        void set_model_name( const std::string& name ) ;

        /*! 
         * @brief When model_topolgy is fixed, any attempt to modify by it 
         *        will throw an assertion/exception
         */
        void forbid_entity_creation()
        {
            create_entity_allowed_ = false ;
        }

        /*! @}
         * \name Creation - Deletion - Access to GeoModelEntities.
         * @{
         */
        GME::gme_t create_entity( const std::string& type ) ;

        template < typename T >
        void add_entity_to_model( T* entity )
        {
            const std::string entity_type = T::type_name_static() ;
            modifiable_mesh_entities( entity_type ).push_back( entity ) ;
        }

        template < typename T >
        GME::gme_t create_mesh_entity()
        {
            const std::string entity_type = T::type_name_static() ;
            index_t nb_entities( model().nb_entities( entity_type ) ) ;
            index_t new_id( nb_entities ) ;
            T* new_entity = new T( model(), new_id ) ;            
            add_entity_to_model( new_entity ) ;
            return new_entity->gme_id() ;
        }

        GME::gme_t create_geological_entity( const std::string& type ) ;

        /*! @}
         * @brief Access to modifiable entities of the edited GeoModel
         * @{
         */
        GeoModelEntity& entity( const GME::gme_t id )
        {
            if( model().is_mesh_entity_type( id.type ) ) {
                return mesh_entity( id ) ;
            } else {
                return geological_entity( id );
            }
        }
        //GeoModelEntity& entity( const std::string& entity_type, index_t entity_index ) ;
        GeoModelMeshEntity& mesh_entity( const GME::gme_t& id )
        {
            return modifiable_mesh_entity( id ) ;
        }
        GeoModelMeshEntity& mesh_entity(
            const std::string& entity_type,
            index_t entity_index )
        {
            return mesh_entity( GME::gme_t( entity_type, entity_index ) ) ;
        }
        GeoModelGeologicalEntity& geological_entity( const GME::gme_t& id )
        {
            return modifiable_geological_entity( id ) ;
        }
        GeoModelGeologicalEntity& geological_entity(
            const std::string& entity_type,
            index_t entity_index )
        {
            return geological_entity( GME::gme_t( entity_type, entity_index ) ) ;
        }

        /*! @}
         * \name Filling GeoModelEntity attributes.
         * @{
         */

        /*!
         * @brief Complete missing information in GeoModelEntities
         * boundaries - in_boundary - parent - children
         */
        void complete_entity_connectivity() ;

        /*!
         * @brief Fill the boundaries of all entities of the given type
         * @details If the boundary entities do not have any in_boundary
         * information, nothing is done.
         */
        void fill_mesh_entities_boundaries( const std::string& type ) ;

        /*!
         * @brief Fill the in_boundary vector of all entities of the given type
         * @details If the in_boundary entities do not have any boundary
         * information, nothing is done, and model construction will eventually fail.
         */
        void fill_mesh_entities_in_boundaries( const std::string& type ) ;

        /*!
         * @brief Fill the parent of all entities of the given type
         * @details If the parents do not have any child nothing is done.
         */
        void fill_mesh_entities_parent( const std::string& type ) ;

        /*!
         * @brief Fill the children of all entities of the given type
         * @details If the children entities do not have any parent information
         * nothing is done.
         */
        void fill_geological_entities_children( const std::string& type ) ;

        void complete_mesh_entities_geol_feature_from_first_parent(
            const std::string& type ) ;
        void complete_geological_entities_geol_feature_from_first_child(
            const std::string& type ) ;

        void set_entity_name( const GME::gme_t& t, const std::string& name )
        {
            entity( t ).name_ = name ;
        }
        void set_entity_geol_feature(
            const GME::gme_t& t,
            GME::GEOL_FEATURE geol )
        {
            entity( t ).geol_feature_ = geol ;
        }
        
        void add_mesh_entity_boundary(
            const GME::gme_t& t,
            index_t boundary_id,
            bool side = false )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const std::string& b_type = model().entity_relationships().boundary_type(
                t.type ) ;
            GME::gme_t boundary( b_type, boundary_id ) ;
            mesh.boundaries_.push_back( boundary ) ;

            if( t.type == Region::type_name_static() ) {
                dynamic_cast< Region& >( mesh ).sides_.push_back( side ) ;
            }
        }

        void set_mesh_entity_boundary(
            const GME::gme_t& t,
            index_t id,
            index_t boundary_id,
            bool side = false )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const std::string& b_type = model().entity_relationships().boundary_type(
                t.type ) ;
            GME::gme_t boundary( b_type, boundary_id ) ;
            mesh.boundaries_[id] = boundary ;

            if( t.type == Region::type_name_static() ) {
                dynamic_cast< Region& >( mesh ).sides_[id] = side ;
            }
        }

        void add_universe_boundary(
            index_t boundary_id,
            bool side )
        {
            GME::gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            model().universe_.boundary_surfaces_.push_back( boundary ) ;
            model().universe_.boundary_surface_sides_.push_back( side ) ;
        }

        void set_universe_boundary(
            index_t id,
            index_t boundary_id,
            bool side )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            ringmesh_assert( id < model().universe_.nb_boundaries() ) ;
            GME::gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            model().universe_.boundary_surfaces_[id] = boundary ;
            model().universe_.boundary_surface_sides_[id] = side ;
        }

        void add_mesh_entity_in_boundary(
            const GME::gme_t& t,
            index_t in_boundary_id )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            const std::string& in_b_type =
                model().entity_relationships().in_boundary_type( t.type ) ;
            GME::gme_t in_boundary( in_b_type, in_boundary_id ) ;
            mesh.in_boundary_.push_back( in_boundary ) ;
        }

        void set_mesh_entity_in_boundary(
            const GME::gme_t& t,
            index_t id,
            index_t in_boundary_id )
        {
            /// No check on the validity of the index of the entity in_boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert( id < mesh.nb_in_boundary() ) ;
            const std::string& in_b_type =
                model().entity_relationships().in_boundary_type( t.type ) ;
            GME::gme_t in_boundary( in_b_type, in_boundary_id ) ;
            mesh.in_boundary_[id] = in_boundary ;
        }

        void add_mesh_entity_parent(
            const GME::gme_t& t,
            const GME::gme_t& parent_index )
        {
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            mesh.parents_.push_back( parent_index ) ;
        }

        void set_mesh_entity_parent(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& parent_index )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert( id < mesh.nb_parents() ) ;
            mesh.parents_[id] = parent_index ;
        }

        void add_geological_entity_child(
            const GME::gme_t& t,
            index_t child_id )
        {
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            const std::string& child_type =
                model().entity_relationships().child_type( t.type ) ;
            GME::gme_t child( child_type, child_id ) ;
            entity.children_.push_back( child ) ;
        }

        void set_geological_entity_child(
            const GME::gme_t& t,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            const std::string& child_type =
                model().entity_relationships().child_type( t.type ) ;
            GME::gme_t child( child_type, child_id ) ;
            entity.children_[id] = child ;
        }

        void remove_entities( const std::set< GME::gme_t >& entities ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_entities( std::set< GME::gme_t >& entities ) const ;

        /*!
         * @brief Not implemented
         * @note The code was copied and wrong. To reimplement if needed. [JP]
         */
        void remove_entities_and_dependencies(
            const std::set< GME::gme_t >& entities_to_remove ) ;

        const GeoModel& model() const
        {
            return model_ ;
        }

    protected:
        /*!
        * @ brief The model under construction
        */
        GeoModel& model()
        {
            return model_ ;
        }
        
        /*!  
         * @brief Throws an assertion if the input GeoModelEntity id is invalid
         *        in the model under construction
         */
        void assert_gme_id_validity( GME::gme_t id )
        {
            return model().assert_gme_id_validity( id ) ;
        }

        EntityRelationships& entity_relationships()
        {
            return model().entity_relationships_ ;
        }
              
        std::vector< GeoModelEntity* >& modifiable_entities(
            const std::string& type )
        {
            return const_cast< std::vector< GeoModelEntity* >& >
                ( model().entities( type ) ) ;
        }
        std::vector< GeoModelMeshEntity* >& modifiable_mesh_entities(
            const std::string& type )
        {
            // Beurk again 
            return const_cast< std::vector< GeoModelMeshEntity* >& >
                ( model().mesh_entities( type ) ) ;
        }      
        std::vector< GeoModelGeologicalEntity* >& modifiable_geological_entities(
            const std::string& type )
        {
            return const_cast< std::vector< GeoModelGeologicalEntity* >& >(
                model().geological_entities( type ) ) ;
        }                           
        GeoModelEntity& modifiable_entity( const GME::gme_t& id )
        {
            if( model().is_mesh_entity_type( id.type ) ) {
                return modifiable_mesh_entity( id ) ;
            } else if( model().is_geological_entity_type( id.type ) ) {
                return modifiable_geological_entity(id) ;
            } else {
                ringmesh_assert_not_reached ;
                GME::gme_t first_surface( Surface::type_name_static(), 0 ) ;
                return modifiable_mesh_entity( first_surface ) ;
            }
        }
        GeoModelMeshEntity& modifiable_mesh_entity(
            const GME::gme_t& id )
        {
            model_.assert_gme_id_validity( id ) ;
            return *modifiable_mesh_entities( id.type )[id.index] ;            
        }
        GeoModelGeologicalEntity& modifiable_geological_entity( const GME::gme_t& id )
        {
            model_.assert_gme_id_validity( id ) ;
            return *modifiable_geological_entities( id.type )[id.index] ;            
        }
        Universe& universe()
        {
            return model_.universe_ ;
        }

        index_t create_mesh_entities( const std::string& type, index_t nb ) ;
        index_t create_geological_entities( const std::string& type, index_t nb ) ;

        index_t create_geological_entity_type( const std::string& type ) ;
        index_t delete_geological_entity_type( const std::string& type ) ;        

        void assert_entity_creation_allowed()
        {
            ringmesh_assert( create_entity_allowed_ ) ;
        }

        // The more I think about it the more this design 
        // for editing the GeoModel and its Entities appears bad.
        void set_entity_index( GeoModelEntity& E, index_t new_index_in_geomodel )
        {        
            E.id_.index = new_index_in_geomodel ;
        }

        void set_boundary_sign( Region& R, index_t boundary_index, bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() ) ;
            R.sides_[boundary_index] = new_side ;
        }

        std::vector< GME::gme_t >& modifiable_children( GeoModelGeologicalEntity& E )
        {
            return E.children_ ;
        }
        std::vector< GME::gme_t >& modifiable_boundaries( GeoModelMeshEntity& E )
        {
            return E.boundaries_ ;
        }
        std::vector< GME::gme_t >& modifiable_boundaries( Universe& U )
        {
            return U.boundary_surfaces_ ;
        }
        std::vector< GME::gme_t >& modifiable_in_boundaries( GeoModelMeshEntity& E )
        {
            return E.in_boundary_ ;
        }
        std::vector< GME::gme_t >& modifiable_parents( GeoModelMeshEntity& E )
        {
            return E.parents_;
        }
        std::vector< bool >& modifiable_sides( Region& R )
        {
            return R.sides_ ;
        }
        std::vector< bool >& modifiable_sides( Universe& U )
        {
            return U.boundary_surface_sides_ ;
        }
        void delete_entity( const std::string& type, index_t index )
        {
            std::vector< GeoModelEntity* >& store = modifiable_entities( type ) ;
            delete store[index] ;
            store[index] = nil ;
        }

    private:
        template< typename E >
        void complete_mesh_entity_connectivity() ;

        template< typename E >
        void copy_mesh_entity_topology( const GeoModel& from ) ;

        GeoModelGeologicalEntity* new_geological_entity(
            const std::string& type,
            index_t id ) ;

    private:
        /*! The model edited 
         */
        GeoModel& model_ ;
        /*! Parameter to forbid element creation. Crucial to control
         *  building of the model and detect errors in find_or_create functions
         */
        bool create_entity_allowed_ ;
    } ;

}


#endif
