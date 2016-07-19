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

#include <geogram/basic/factory.h>

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
         * @brief Returns the number of mesh entities of the given type
         * @details Default value is 0
         * @param[in] type the mesh entity type
         */
        index_t nb_mesh_entities( const std::string& type ) const
        {
            if( type == Corner::type_name_ ) {
                return nb_corners();
            } else if( type == Line::type_name_ ) {
                return nb_lines();
            } else if( type == Surface::type_name_ ) {
                return nb_surfaces();
            } else if( type == Region::type_name_ ) {
                return nb_regions();
            } else {
                ringmesh_assert_not_reached ;
                return 0 ;
            }
        }


        /*!
         * @brief Returns the number of geological entities of the given type
         * @details Default value is 0
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entities( const std::string& type ) const
        {
            index_t index = geological_entity_type( type ) ;
            if( index == NO_ID ) return 0 ;
            return geological_entities_[index].size() ;
        }

        /*!
         * @brief Returns the index of the geological entity type storage
         * @details Default value is NO_ID
         * @param[in] type the geological entity type
         */
        index_t geological_entity_type( const std::string& type ) const
        {
            return find( geological_entity_types_, type ) ;

        }

        /*!
         * @brief Returns a const reference the identified GeoModelEntity
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         * @pre Entity identification is valid.
         */
        const GeoModelEntity& geological_entity( GME::gme_t id ) const
        {
            return *geological_entity_ptr( id ) ;
        }
        
        /*!
         * Convenient overload of entity( GME::gme_t id )
         */
        const GeoModelEntity& geological_entity(
            const std::string& entity_type,
            index_t entity_index ) const
        {
            return geological_entity( GME::gme_t( entity_type, entity_index ) ) ;
        }
   
        /*!
         * @brief Generic access to a meshed entity
         * @pre Type of the entity is CORNER, LINE, SURFACE, or REGION
         */
        const GeoModelMeshEntity& mesh_entity( GME::gme_t id ) const
        {
            return dynamic_cast< const GeoModelMeshEntity& >( mesh_entity( id ) ) ;
        }

        /*!
         * Convenient overload of mesh_entity( GME::gme_t id )
         */
        const GeoModelMeshEntity& mesh_entity(
            const std::string& entity_type,
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
        const Universe& universe() const
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
         * @brief Generic accessor to the storage of mesh entities of the given type
         */
        inline std::vector< GeoModelMeshEntity* >& modifiable_mesh_entities(
            const std::string& type )
        {
            return const_cast< std::vector< GeoModelMeshEntity* >& >( mesh_entities(
                type ) ) ;
        }
        /*!
         * @brief Generic accessor to the storage of geological entities of the given type
         */
        inline std::vector< GeoModelGeologicalEntity* >& modifiable_geological_entities(
            const std::string& type )
        {
            return const_cast< std::vector< GeoModelGeologicalEntity* >& >( mesh_entities(
                type ) ) ;
        }

        /*!
         * @brief Generic accessor to the storage of mesh entities of the given type
         */
        const std::vector< GeoModelMeshEntity* >& mesh_entities(
            const std::string& type ) const
        {
            if( type == Corner::type_name_ ) {
                return *(std::vector< GME* > *) &corners_ ;
            } else if( type == Line::type_name_ ) {
                return *(std::vector< GME* > *) &lines_ ;
            } else if( type == Surface::type_name_ ) {
                return *(std::vector< GME* > *) &surfaces_ ;
            } else if( type == Region::type_name_ ) {
                return *(std::vector< GME* > *) &regions_ ;
            } else {
                ringmesh_assert_not_reached ;
                return *(std::vector< GME* > *) &surfaces_ ;
            }
        }


        /*!
         * @brief Generic accessor to the storage of geologcial entities of the given type
         */
        const std::vector< GeoModelGeologicalEntity* >& geological_entities(
            const std::string& type ) const
        {
            index_t index = geological_entity_type( type ) ;
            ringmesh_assert( index != NO_ID ) ;
            return geological_entities_[index] ;
        }

        /*!
         * @brief Modifiable pointer to a mesh entity of the model
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         *
         * @pre Entity identification is valid.
         */
        const GeoModelMeshEntity* mesh_entity_ptr(
            const GME::gme_t& id ) const
        {
            return mesh_entities( id.type )[id.index] ;
        }

        /*!
         * @brief Modifiable pointer to an entity of the model
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         *
         * @pre Entity identification is valid.
         */
        const GeoModelGeologicalEntity* geological_entity_ptr(
            const GME::gme_t& id ) const
        {
            return geological_entities( id.type )[id.index] ;
        }

        /*!
         * @brief Reference to a modifiable meshed entity of the model
         * @pre Assert in debug model that the given id refers to a meshed entity.
         *      The id must refer to a valid entity.
         */
        inline GeoModelMeshEntity& modifiable_mesh_entity(
            const GME::gme_t& id )
        {
            return dynamic_cast< GeoModelMeshEntity& >( mesh_entity_ptr( id ) ) ;
        }

        /*!
         * @brief Reference to a modifiable entity of the model
         * @pre The id must refer to a valid entity of the model
         */
        GeoModelGeologicalEntity& modifiable_geological_entity( const GME::gme_t& id )
        {
            return const_cast< GeoModelGeologicalEntity& >( *geological_entity_ptr( id ) ) ;
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
         */
        Universe universe_ ;

        /*! @}
         * \name Optional geological entities
         * @{
         */

        /*!
         * @brief Geological entities
         */
        std::vector< std::vector< GeoModelGeologicalEntity* > > geological_entities_ ;
        /*!
         * @brief Geological entity types
         */
        std::vector< std::string > geological_entity_types_ ;

        /*!
         * @}
         */

        /*! Optional WellGroup associated with the model
         * @todo Give a more general name - this could be anything [JP]
         * @todo Does it really have to be an attribute of GeoModel ? [JP]
         */
        const WellGroup* wells_ ;
    } ;

    GEO::Factory0< GeoModelGeologicalEntity > GeoModelGeologicalEntityFactory ;
#define ringmesh_register_GeoModelGeologicalEntity_creator( type, name ) \
    geo_register_creator( GeoModelGeologicalEntityFactory, type, name )

}

#endif
