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

        GME::gme_t create_entity( GME::TYPE e_type ) ;

        void create_geomodel_entities( GME::TYPE entity_type, index_t nb_entities )
        {
            for( index_t i = 0; i < nb_entities; ++i ) {
                create_entity( entity_type ) ;
            }
        }

        /*! @}
         * \name Creation - Deletion - Access to GeoModelEntities.
         * @{
         */

        /*!
         * @brief Reference to a modifiable entity of the model
         * @pre The id must refer to a valid entity of the model
         * @note Stupid override of a function of GeoModel
         */
        GeoModelEntity& entity( const GME::gme_t& id ) ;

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         * @pre The id must refer to a valid entity.
         * @note Stupid override of a function of GeoModel
         */
        GeoModelMeshEntity& mesh_entity( const GME::gme_t& id )
        {
            ringmesh_assert( GME::has_mesh( id.type ) ) ;
            return dynamic_cast< GeoModelMeshEntity& >( entity( id ) ) ;
        }

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         */
        GeoModelMeshEntity& mesh_entity(
            GeoModelEntity::TYPE entity_type,
            index_t entity_index )
        {
            return mesh_entity( GME::gme_t( entity_type, entity_index ) ) ;
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
        void fill_entities_boundaries( GME::TYPE type ) ;

        /*!
         * @brief Fill the in_boundary vector of all entities of the given type
         * @details If the in_boundary entities do not have any boundary
         * information, nothing is done, and model construction will eventually fail.
         */
        void fill_entities_in_boundaries( GME::TYPE type ) ;

        /*!
         * @brief Fill the parent of all entities of the given type
         * @details If the parents do not have any child nothing is done.
         */
        void fill_entities_parent( GME::TYPE type ) ;

        /*!
         * @brief Fill the children of all entities of the given type
         * @details If the children entities do not have any parent information
         * nothing is done.
         */
        void fill_entities_children( GME::TYPE type ) ;

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

        void add_entity_boundary(
            const GME::gme_t& t,
            const GME::gme_t& boundary,
            bool side = false )
        {
            ringmesh_assert( boundary.is_defined() ) ;
            ringmesh_assert( GME::boundary_type( t.type ) == boundary.type ) ;
            entity( t ).boundaries_.push_back( boundary ) ;

            if( t.type == GME::REGION ) {
                dynamic_cast< Region& >( entity( t ) ).sides_.push_back( side ) ;
            }
        }

        void set_entity_boundary(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& boundary,
            bool side = false )
        {
            /// No check on the validity of the index of the entity boundary
            /// NO_ID is used to flag entities to delete
            ringmesh_assert( GME::boundary_type( t.type ) == boundary.type ) ;
            ringmesh_assert( id < entity( t ).nb_boundaries() ) ;
            entity( t ).boundaries_[id] = boundary ;

            if( t.type == GME::REGION ) {
                dynamic_cast< Region& >( entity( t ) ).sides_[id] = side ;
            }
        }

        void add_entity_in_boundary(
            const GME::gme_t& t,
            const GME::gme_t& in_boundary )
        {
            ringmesh_assert( in_boundary.is_defined() ) ;
            ringmesh_assert(
                GME::in_boundary_type( t.type ) == in_boundary.type ) ;
            entity( t ).in_boundary_.push_back( in_boundary ) ;
        }

        void set_entity_in_boundary(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& in_boundary )
        {
            /// No check on the validity of the index of the entity in_boundary
            /// NO_ID is used to flag entities to delete
            ringmesh_assert(
                GME::in_boundary_type( t.type ) == in_boundary.type ) ;
            ringmesh_assert( id < entity( t ).nb_in_boundary() ) ;
            entity( t ).in_boundary_[id] = in_boundary ;
        }

        void set_entity_parent(
            const GME::gme_t& t,
            const GME::gme_t& parent_index )
        {
            ringmesh_assert( GME::parent_type( t.type ) == parent_index.type ) ;
            entity( t ).parent_ = parent_index ;
        }

        void add_entity_child(
            const GME::gme_t& t,
            const GME::gme_t& child_index )
        {
            ringmesh_assert( child_index.is_defined() ) ;
            ringmesh_assert( GME::child_type( t.type ) == child_index.type ) ;
            entity( t ).children_.push_back( child_index ) ;
        }

        void set_entity_child(
            const GME::gme_t& t,
            index_t id,
            const GME::gme_t& child_index )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            ringmesh_assert( GME::child_type( t.type ) == child_index.type ) ;
            entity( t ).children_[id] = child_index ;
        }

        void set_region_side( Region& region, index_t surf_boundary, bool side )
        {
            region.sides_[surf_boundary] = side ;
        }

        /*!
         * @brief Set an entity of the model.
         * @details It is on purpose that entity validity is not checked.
         *          This way nil pointers can be set for a further entity removal.
         * @param id Id card of the entity to modify. The ownership of the previous entity
         *           is given up by the GeoModel.
         * @param E Entity to set. The ownership is transferred to the GeoModel.
         */
        void set_entity( const GME::gme_t& id, GeoModelEntity* E ) const ;

        void remove_entities( const std::set< GME::gme_t >& entities ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_entities( std::set< GME::gme_t >& entities ) const ;

        void remove_entities_and_dependencies(
            const std::set< GME::gme_t >& entities_to_remove ) ;

    protected:
        /*!
        * @ brief The model under construction
        */
        const GeoModel& model() const
        {
            return model_ ;
        }
        GeoModel& model()
        {
            return model_ ;
        }

        void delete_elements( std::vector< std::vector< index_t > >& to_erase ) ;

        /*!
         * @brief Adds \p nb new entities of type \p type
         * @details Adds \p nb new entities at the end of the
         * vector of entities of type \p type.
         * @param[in] type the type of entity to add.
         * @param[in] nb the number of new entities.
         * @return the old number of entities of the type \p type.
         * It corresponds to the index of the first new entity.
         */
        index_t create_entities( GME::TYPE type, index_t nb ) ;

        void erase_invalid_entity_references( const GME::gme_t& E_id ) ;

        void assert_entity_creation_allowed()
        {
            ringmesh_assert( create_entity_allowed_ ) ;
        }

    private:
        void copy_entity_topology(
            const GME::gme_t& lhs_id,
            const GeoModel& from,
            const GME::gme_t& rhs_id ) ;

        /*!
         * @brief Creates an empty entity of the right type in the GeoModel
         */
        GME* new_entity( GME::TYPE type, index_t id ) ;
        GME* new_entity( GME::TYPE type ) ;

    private:
        GeoModel& model_ ;
        bool create_entity_allowed_ ;
    } ;
}

#endif
