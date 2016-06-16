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

#ifndef __RINGMESH_GEO_MODEL__
#define __RINGMESH_GEO_MODEL__

#include <ringmesh/common.h>

#include <ringmesh/geo_model_entity.h>
#include <ringmesh/geo_model_mesh.h>

#include <vector>

/*!
 * @file ringmesh/geo_model.h
 * @brief Class representing a geological structural model: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {
    class WellGroup ;
}

namespace RINGMesh {

    /*!
     * @brief The class to describe a geological model represented 
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    class RINGMESH_API GeoModel {
        ringmesh_disable_copy( GeoModel ) ;
        friend class GeoModelBuilder ;
        friend class GeoModelEditor ;
        friend class GeoModelRepair ;

    public:
        static const index_t NO_ID = index_t( -1 ) ;

        /*!
         * @brief Constructs an empty GeoModel
         */
        GeoModel() ;

        /*!
         * @brief Deletes all the GeoModelEntities of the GeoModel
         */
        virtual ~GeoModel() ;

        /*
         * @todo Implement a real copy_constructor and operator= [JP]
         */
        void copy( const GeoModel& from ) ;

        /*!
         * @brief Name of the model
         */
        const std::string& name() const
        {
            return name_ ;
        }

        /*!
         * \name Generic GeoModelEntity accessors
         * @{
         */

        /*!
         * @brief Returns the number of entities of the given type
         * @details Default value is 0
         * @param[in] type the entity type
         */
        index_t nb_entities( GME::TYPE type ) const
        {
            switch( type ) {
                case GME::CORNER:
                    return nb_corners();
                case GME::LINE:
                    return nb_lines();
                case GME::SURFACE:
                    return nb_surfaces();
                case GME::REGION:
                    return nb_regions();
                case GME::CONTACT:
                    return nb_contacts(); 
                case GME::INTERFACE:
                    return nb_interfaces();
                case GME::LAYER:
                    return nb_layers();
                case GME::ALL_TYPES:
                    ringmesh_assert( !nb_entities_per_type_.empty() ) ;
                    return nb_entities_per_type_.back() ;
                default:
                    ringmesh_assert_not_reached ;
                    return 0 ;
            }
        }

        /*!
         * @brief Returns a const reference the identified GeoModelEntity
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         * @pre Entity identification is valid.
         */
        const GeoModelEntity& entity( GME::gme_t id ) const
        {
            return *entity_ptr( id ) ;
        }
        
        /*!
         * Convenient overload of entity( GME::gme_t id )
         */
        const GeoModelEntity& entity(
            GME::TYPE entity_type,
            index_t entity_index ) const
        {
            return entity( GME::gme_t( entity_type, entity_index ) ) ;
        }
   
        /*!
         * @brief Generic access to a meshed entity
         * @pre Type of the entity is CORNER, LINE, SURFACE, or REGION
         */
        const GeoModelMeshEntity& mesh_entity( GME::gme_t id ) const
        {
            ringmesh_assert( GME::has_mesh( id.type ) ) ;
            return dynamic_cast< const GeoModelMeshEntity& >( entity( id ) ) ;
        }

        /*!
         * Convenient overload of mesh_entity( GME::gme_t id )
         */
        const GeoModelMeshEntity& mesh_entity(
            GME::TYPE entity_type,
            index_t entity_index ) const
        {
            return mesh_entity( GME::gme_t( entity_type, entity_index ) ) ;
        }

        /*! @}
         * \name Specialized accessors.
         * @{
         */
        index_t nb_corners() const
        {
            return static_cast< index_t >( corners_.size() ) ;
        }
        index_t nb_lines() const
        {
            return static_cast< index_t >( lines_.size() ) ;
        }
        index_t nb_surfaces() const
        {
            return static_cast< index_t >( surfaces_.size() ) ;
        }
        index_t nb_regions() const
        {
            return static_cast< index_t >( regions_.size() ) ;
        }
        index_t nb_contacts() const
        {
            return static_cast< index_t >( contacts_.size() ) ;
        }
        index_t nb_interfaces() const
        {
            return static_cast< index_t >( interfaces_.size() ) ;
        }
        index_t nb_layers() const
        {
            return static_cast< index_t >( layers_.size() ) ;
        }

        const Corner& corner( index_t index ) const
        {
            return *corners_.at( index ) ;
        }
        const Line& line( index_t index ) const
        {
            return *lines_.at( index ) ;
        }
        const Surface& surface( index_t index ) const
        {                                  
            return *surfaces_.at( index ) ;
        }
        const Region& region( index_t index ) const
        {
            return *regions_.at( index ) ;
        }
        const GeoModelEntity& contact( index_t index ) const
        {
            return entity( GME::CONTACT, index ) ;
        }
        const GeoModelEntity& one_interface( index_t index ) const
        {
            return entity( GME::INTERFACE, index ) ;
        }
        const GeoModelEntity& layer( index_t index ) const
        {
            return entity( GME::LAYER, index ) ;
        }
        const Region& universe() const
        {
            return universe_ ;
        }

        /*!
         * @}
         */
        /* @todo Wells have nothing at all to do here [JP] */
        void set_wells( const WellGroup* wells ) ;
        const WellGroup* wells() const
        {
            return wells_ ;
        }

    private:
        /*!
         * @brief Convert a global GME index into a typed index
         * @details Relies on the nb_entities_per_type_ vector that
         *          must be up to date at all times
         *          See the GeoModelBuilder::end_model() function
         * @param[in] global A GME id of TYPE - ALL_TYPES
         * @return A GME id of an entity of the model, or a invalid one if nothing found
         */
        inline GME::gme_t global_to_typed_id( const GME::gme_t& global ) const
        {
            ringmesh_assert( global.type == GME::ALL_TYPES ) ;

            index_t t = NO_ID ;
            for( index_t i = 1; i < nb_entities_per_type_.size(); i++ ) {
                if( global.index >= nb_entities_per_type_[i - 1]
                    && global.index < nb_entities_per_type_[i] ) {
                    t = i - 1 ;
                    break ;
                }
            }
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            if( T < GME::NO_TYPE ) {
                index_t i = global.index - nb_entities_per_type_[t] ;
                return GME::gme_t( T, i ) ;
            } else {
                ringmesh_assert_not_reached ;
                return GME::gme_t() ;
            }
        }

        /*!
         * @brief Generic accessor to the storage of entities of the given type
         * @pre The type must be valid NO_TYPE or ALL_TYPES will throw an assertion
         */
        inline std::vector< GME* >& modifiable_entities( GME::TYPE type )
        {
            return const_cast< std::vector< GME* >& >( entities( type ) ) ;
        }

        /*!
         * @brief Generic accessor to the storage of entities of the given type
         * @pre The type must be valid. NO_TYPE or ALL_TYPES will throw an assertion
         */
        const std::vector< GME* >& entities( GME::TYPE type ) const
        {
            // The following casts are really nasty, I know.
            // But we need this generic access to vectors of GME* [JP]
            switch( type ) {
                case GME::CORNER:
                    return *(std::vector<GME*> *)&corners_ ;
                case GME::LINE:
                    return *(std::vector<GME*> *)&lines_ ;
                case GME::SURFACE:
                    return *(std::vector<GME*> *)&surfaces_ ;
                case GME::REGION:
                    return *(std::vector<GME*> *)&regions_ ;
                case GME::CONTACT:
                    return contacts_ ;
                case GME::INTERFACE:
                    return interfaces_ ;
                case GME::LAYER:
                    return layers_ ;
                default:
                    ringmesh_assert_not_reached ;
                    return *(std::vector<GME*> *)&surfaces_ ;
            }
        }

        /*!
         * @brief Modifiable pointer to an entity of the model
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         *
         * @pre Entity identification is valid.
         */
        const GeoModelEntity* entity_ptr( const GME::gme_t& id ) const
        {
            if( id.type == GME::REGION && id.index == NO_ID ) {
                return &universe_ ;
            } else {
                if( id.type < GME::NO_TYPE ) {
                    ringmesh_assert( id.index < nb_entities( id.type ) ) ;
                    return entities( id.type )[id.index] ;
                } else if( id.type == GME::ALL_TYPES ) {
                    return entity_ptr( global_to_typed_id( id ) ) ;
                } else {
                    ringmesh_assert_not_reached ;
                    return &universe_ ;
                }
            }
        }

        /*!
         * @brief Reference to a modifiable entity of the model
         * @pre The id must refer to a valid entity of the model
         */
        GeoModelEntity& modifiable_entity( const GME::gme_t& id )
        {
            return const_cast< GeoModelEntity& >( *entity_ptr( id ) ) ;
        }

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         * @pre Assert in debug model that the given id refers to a meshed entity.
         *      The id must refer to a valid entity.
         */
        inline GeoModelMeshEntity& modifiable_mesh_entity(
            const GME::gme_t& id ) //const
        {
            ringmesh_assert( GME::has_mesh( id.type ) ) ;
            return dynamic_cast< GeoModelMeshEntity& >( modifiable_entity( id ) ) ;
        }

        /*!
         * @brief Clears and fills the model nb_entities_per_type_ vector
         * @details See global entity access with GeoModel::entity( GME::TYPE, index_t )
         */
        void init_global_model_entity_access()
        {
            nb_entities_per_type_.clear() ;

            index_t count = 0 ;
            nb_entities_per_type_.push_back( count ) ;
            for( index_t type = GME::CORNER; type < GME::NO_TYPE; type++ ) {
                count += nb_entities( (GME::TYPE) type ) ;
                nb_entities_per_type_.push_back( count ) ;
            }
        }

    public:
        GeoModelMesh mesh ;

    private:
        // Name of the model
        std::string name_ ;

        /*!
         * \name Mandatory entities of the model
         * @{
         */
        std::vector< Corner* > corners_ ;
        std::vector< Line* > lines_ ;
        std::vector< Surface* > surfaces_ ;   
        std::vector< Region* > regions_ ;

        /*!
         * The Region defining the model extension
         * @todo Put it as the last item in regions_ and do not
         *       forget to create it in the Builder
         */
        Region universe_ ;

        /*! @}
         * \name Optional geological entities
         * @{
         */

        /*!
         * @brief Entities of type CONTACT
         */
        std::vector< GeoModelEntity* > contacts_ ;
        /*!
         * @brief Entities of type INTERFACE
         */
        std::vector< GeoModelEntity* > interfaces_ ;
        /*!
         * @brief Entities of type LAYER
         */
        std::vector< GeoModelEntity* > layers_ ;

        /*!
         * @}
         */

        /*
         * @brief Global access to GeoModelEntities.
         * It MUST be updated if one entity is added.
         * @warning It must be up to date at all times
         */
        std::vector< index_t > nb_entities_per_type_ ;

        /*! Optional WellGroup associated with the model
         * @todo Give a more general name - this could be anything [JP]
         * @todo Does it really have to be an attribute of GeoModel ? [JP]
         */
        const WellGroup* wells_ ;
    } ;

}

#endif
