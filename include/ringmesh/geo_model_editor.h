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
        GeoModelEditor( GeoModel& M )
            : model_(M), create_entity_allowed_(true)
        {
        }

        ~GeoModelEditor()
        {
        }

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

        GME::gme_t create_mesh_entity( const std::string& type ) ;
        GME::gme_t create_geological_entity( const std::string& type ) ;

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         * @pre The id must refer to a valid entity.
         * @note Stupid override of a function of GeoModel
         */
        GeoModelMeshEntity& mesh_entity( const GME::gme_t& id )
        {
            return model_.modifiable_mesh_entity( id ) ;
        }

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         */
        GeoModelMeshEntity& mesh_entity(
            const std::string& entity_type,
            index_t entity_index )
        {
            return mesh_entity( GME::gme_t( entity_type, entity_index ) ) ;
        }
        /*!
         * @brief Reference to a modifiable entity of the model
         * @pre The id must refer to a valid entity of the model
         * @note Stupid override of a function of GeoModel
         */
        GeoModelGeologicalEntity& geological_entity( const GME::gme_t& id )
        {
            return model_.modifiable_geological_entity( id ) ;
        }
        /*!
         * @brief Reference to a modifiable meshed entity of the model
         */
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

        void set_mesh_entity_name( const GME::gme_t& t, const std::string& name )
        {
            mesh_entity( t ).name_ = name ;
        }
        void set_geological_entity_name(
            const GME::gme_t& t,
            const std::string& name )
        {
            geological_entity( t ).name_ = name ;
        }

        void set_mesh_entity_geol_feature(
            const GME::gme_t& t,
            GME::GEOL_FEATURE geol )
        {
            mesh_entity( t ).geol_feature_ = geol ;
        }
        void set_geological_entity_geol_feature(
            const GME::gme_t& t,
            GME::GEOL_FEATURE geol )
        {
            geological_entity( t ).geol_feature_ = geol ;
        }

        void add_mesh_entity_boundary(
            const GME::gme_t& t,
            const GME::gme_t& boundary,
            bool side = false )
        {
            ringmesh_assert( boundary.is_defined() ) ;
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert( mesh.boundary_type() == boundary.type ) ;
            mesh.boundaries_.push_back( boundary ) ;

            if( t.type == Region::type_name_ ) {
                dynamic_cast< Region& >( mesh ).sides_.push_back( side ) ;
            }
        }

        void set_mesh_entity_boundary(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& boundary,
            bool side = false )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert( mesh.boundary_type() == boundary.type ) ;
            ringmesh_assert( id < mesh.nb_boundaries() ) ;
            mesh.boundaries_[id] = boundary ;

            if( t.type == Region::type_name_ ) {
                dynamic_cast< Region& >( mesh ).sides_[id] = side ;
            }
        }

        void add_universe_boundary(
            const GME::gme_t& boundary,
            bool side )
        {
            ringmesh_assert( boundary.is_defined() ) ;
            model().universe_.boundary_surfaces_.push_back( boundary ) ;
            model().universe_.boundary_surface_sides_.push_back( side ) ;
        }

        void set_universe_boundary(
            index_t id,
            const GME::gme_t& boundary,
            bool side )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            ringmesh_assert( id < model().universe_.nb_boundaries() ) ;
            model().universe_.boundary_surfaces_[id] = boundary ;
            model().universe_.boundary_surface_sides_[id] = side ;
        }

        void add_mesh_entity_in_boundary(
            const GME::gme_t& t,
            const GME::gme_t& in_boundary )
        {
            ringmesh_assert( in_boundary.is_defined() ) ;
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert(
                mesh.in_boundary_type() == in_boundary.type ) ;
            mesh.in_boundary_.push_back( in_boundary ) ;
        }

        void set_mesh_entity_in_boundary(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& in_boundary )
        {
            /// No check on the validity of the index of the entity in_boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh = mesh_entity( t ) ;
            ringmesh_assert(
                mesh.in_boundary_type() == in_boundary.type ) ;
            ringmesh_assert( id < mesh.nb_in_boundary() ) ;
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
            const GME::gme_t& child_index )
        {
            ringmesh_assert( child_index.is_defined() ) ;
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            ringmesh_assert( entity.child_type() == child_index.type ) ;
            entity.children_.push_back( child_index ) ;
        }

        void set_geological_entity_child(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& child_index )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity& entity = geological_entity( t ) ;
            ringmesh_assert( entity.child_type() == child_index.type ) ;
            entity.children_[id] = child_index ;
        }

        void remove_entities( const std::set< GME::gme_t >& entities ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_entities( std::set< GME::gme_t >& entities ) const ;

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

        EntityRelationships& entity_relationships()
        {
            return model().entity_relationships_ ;
        }

        void delete_elements( std::vector< std::vector< index_t > >& to_erase ) ;

        index_t create_mesh_entities( const std::string& type, index_t nb ) ;
        index_t create_geological_entities( const std::string& type, index_t nb ) ;

        index_t create_geological_entity_type( const std::string& type ) ;
        void erase_invalid_entity_references( const GME::gme_t& E_id ) ;

        void assert_entity_creation_allowed()
        {
            ringmesh_assert( create_entity_allowed_ ) ;
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
        GeoModel& model_ ;
        bool create_entity_allowed_ ;
    } ;
}


#endif
