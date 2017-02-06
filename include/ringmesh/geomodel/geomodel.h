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

#ifndef __RINGMESH_GEOMODEL__
#define __RINGMESH_GEOMODEL__

#include <ringmesh/basic/common.h>

#include <geogram/basic/factory.h>

#include <ringmesh/basic/algorithm.h>

#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh.h>

#include <vector>

/*!
 * @file ringmesh/geomodel.h
 * @brief Class representing a geological structural geomodel: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {
    class WellGroup ;
    class GeoModelGeologicalEntity ;
    class GeoModelMeshEntity ;
    class Corner ;
    class Surface ;
    class Line ;
    class Region ;
}

namespace RINGMesh {
    /*!  
     * @brief Manages the type relationship between GeoModelEntities
     * Each GeoModel owns one instance of it.
     */
    class RINGMESH_API EntityTypeManager {
    ringmesh_disable_copy( EntityTypeManager ) ;
    public:
        friend class GeoModelBuilderTopology ;
        EntityTypeManager()
        {
        }

        typedef std::string EntityType ;
        typedef std::string MeshEntityType ;
        typedef std::string GeologicalEntityType ;
        typedef std::map< MeshEntityType, std::set< GeologicalEntityType > > MeshEntityToParents ;
        typedef std::map< GeologicalEntityType, MeshEntityType > GeologicalEntityToChild ;

        static bool is_defined_type( const EntityType& type )
        {
            return type != default_entity_type() ;
        }
        static const EntityType default_entity_type() ;
        bool is_valid_type( const EntityType& type ) const
        {
            return is_defined_type( type )
                && ( is_mesh_entity_type( type ) || is_geological_entity_type( type ) ) ;
        }

        // ---- Static members : access to relationships between MeshEntities
        // ---- Maybe they could be outside the class.
        static bool is_mesh_entity_type( const EntityType& type ) ;
        static bool is_corner( const EntityType& type ) ;
        static bool is_line( const EntityType& type ) ;
        static bool is_surface( const EntityType& type ) ;
        static bool is_region( const EntityType& type ) ;
        static const MeshEntityType& boundary_type( const MeshEntityType& type ) ;
        static const MeshEntityType& in_boundary_type( const MeshEntityType& type ) ;
        static const std::vector< MeshEntityType >& mesh_entity_types() ;
        static index_t nb_mesh_entity_types() ;

        // --- GeologicalEntity and MeshEntity relationships
        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const
        {
            MeshEntityToParents::const_iterator itr = child_to_parents_.find(
                child_type ) ;
            std::vector< GeologicalEntityType > result ;
            if( itr != child_to_parents_.end() ) {
                result.insert( result.begin(), itr->second.begin(),
                    itr->second.end() ) ;
            }
            return result ;
        }
        index_t nb_parent_types( const MeshEntityType& child_type ) const
        {
            return static_cast< index_t >( parent_types( child_type ).size() ) ;
        }
        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const
        {
            GeologicalEntityToChild::const_iterator itr = parent_to_child_.find(
                parent_type ) ;
            if( itr == parent_to_child_.end() ) {
                return default_entity_type() ;
            } else {
                return itr->second ;
            }
        }

        // --- Geological entity types 
        bool is_geological_entity_type( const EntityType& type ) const ;
        index_t nb_geological_entity_types() const
        {
            return static_cast< index_t >( geological_entity_types_.size() ) ;
        }
        const std::vector< GeologicalEntityType >& geological_entity_types() const
        {
            return geological_entity_types_ ;
        }
        const EntityType& geological_entity_type( index_t index ) const
        {
            return geological_entity_types_.at( index ) ;
        }
        index_t geological_entity_type_index( const EntityType& type ) const
        {
            return find( geological_entity_types_, type ) ;
        }

    private:
        void register_geological_entity_type(
            const GeologicalEntityType& geological_type_name )
        {
            if( find( geological_entity_types_, geological_type_name ) == NO_ID ) {
                geological_entity_types_.push_back( ( geological_type_name ) ) ;
            }
        }

        void register_relationship(
            const EntityType& parent_type_name,
            const EntityType& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name ;
            child_to_parents_[child_type_name].insert( parent_type_name ) ;
        }

    private:
        std::vector< GeologicalEntityType > geological_entity_types_ ;
        GeologicalEntityToChild parent_to_child_ ;
        MeshEntityToParents child_to_parents_ ;
    } ;

    /*!
     * @brief The class to describe a geological geomodel represented 
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    class RINGMESH_API GeoModel {
    ringmesh_disable_copy( GeoModel ) ;
        friend class GeoModelAccess ;

    public:
        typedef std::string EntityType ;

        /*!
         * @brief Constructs an empty GeoModel
         */
        GeoModel() ;
        /*!
         * @brief Deletes all the GeoModelEntities of the GeoModel
         */
        virtual ~GeoModel() ;

        /*!
         * @brief Gets the name of the GeoModel
         */
        const std::string& name() const
        {
            return geomodel_name_ ;
        }

        /*!
         * @brief Gets the EntityTypeManager associated to the GeoModel
         */
        const EntityTypeManager& entity_type_manager() const
        {
            return entity_type_manager_ ;
        }
        /*!
         * @brief Tests if the EntityType corresponds to a GeoModelMeshEntity
         * @param[in] type the EntityType to test
         * @return true is it is a GeoModelMeshEntity type
         */
        bool is_mesh_entity_type( const EntityType& type ) const
        {
            return EntityTypeManager::is_mesh_entity_type( type ) ;
        }
        /*!
         * @brief Tests if the EntityType corresponds to a GeoModelGeologicalEntity
         * @param[in] type the EntityType to test
         * @return true is it is a GeoModelGeologicalEntity type
         */
        bool is_geological_entity_type( const EntityType& type ) const
        {
            return entity_type_manager_.is_geological_entity_type( type ) ;
        }

        /*!
         * @brief Returns the number of mesh entities of the given type
         * @details Default value is 0
         * @param[in] type the mesh entity type
         */
        index_t nb_mesh_entities( const EntityType& type ) const ;

        /*!
         * @brief Returns the number of geological entities of the given type
         * @details Default value is 0
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entities( const EntityType& type ) const
        {
            if( !is_geological_entity_type( type ) ) {
                return 0 ;
            } else {
                return static_cast< index_t >( geological_entities( type ).size() ) ;
            }
        }
        /*!
         * @brief Returns the index of the geological entity type storage
         * @details Default value is NO_ID
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entity_types() const
        {
            return entity_type_manager_.nb_geological_entity_types() ;
        }

        const EntityType& geological_entity_type( index_t index ) const
        {
            return entity_type_manager_.geological_entity_type( index ) ;
        }
        /*!
         * @brief Returns a const reference the identified GeoModelGeologicalEntity
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         * @pre Entity identification is valid.
         */
        const GeoModelGeologicalEntity& geological_entity( gme_t id ) const
        {
            return *geological_entities( id.type )[id.index] ;
        }
        /*!
         * Convenient overload of entity( gme_t id )
         */
        const GeoModelGeologicalEntity& geological_entity(
            const EntityType& entity_type,
            index_t entity_index ) const
        {
            return geological_entity( gme_t( entity_type, entity_index ) ) ;
        }
        /*!
         * @brief Generic access to a meshed entity
         * @pre Type of the entity is CORNER, LINE, SURFACE, or REGION
         */
        const GeoModelMeshEntity& mesh_entity( gme_t id ) const ;
        /*!
         * Convenient overload of mesh_entity( gme_t id )
         */
        const GeoModelMeshEntity& mesh_entity(
            const EntityType& entity_type,
            index_t entity_index ) const
        {
            return mesh_entity( gme_t( entity_type, entity_index ) ) ;
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

        double epsilon() const ;

        double epsilon2() const
        {
            return epsilon() * epsilon() ;
        }

        double epsilon3() const
        {
            return epsilon2() * epsilon() ;
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

    public:
        GeoModelMesh mesh ;

    private:
        /*!
         * Access to the position of the entity of that type in storage.
         */
        index_t geological_entity_type_index( const EntityType& type ) const
        {
            return entity_type_manager_.geological_entity_type_index( type ) ;
        }
        /*!
         * @brief Generic accessor to the storage of mesh entities of the given type
         */
        const std::vector< GeoModelMeshEntity* >& mesh_entities(
            const EntityType& type ) const ;

        /*!
         * @brief Generic accessor to the storage of geological entities of the given type
         */
        const std::vector< GeoModelGeologicalEntity* >& geological_entities(
            const EntityType& type ) const
        {
            index_t entity_index = geological_entity_type_index( type ) ;
            return geological_entities( entity_index ) ;
        }

        const std::vector< GeoModelGeologicalEntity* >& geological_entities(
            index_t geological_entity_type_index ) const
        {
            ringmesh_assert( geological_entity_type_index != NO_ID ) ;
            return geological_entities_[geological_entity_type_index] ;
        }

    private:
        std::string geomodel_name_ ;
        mutable double epsilon_ ;

        EntityTypeManager entity_type_manager_ ;

        /*!
         * \name Mandatory entities of the geomodel
         * @{
         */
        std::vector< Corner* > corners_ ;
        std::vector< Line* > lines_ ;
        std::vector< Surface* > surfaces_ ;
        std::vector< Region* > regions_ ;

        /*!
         * The Universe defines the extension of the GeoModel
         */
        Universe universe_ ;

        /*!
         * @brief Geological entities. They are optional.
         * The EntityTypes are managed by the EntityTypeManager of the class.
         */
        std::vector< std::vector< GeoModelGeologicalEntity* > > geological_entities_ ;

        /*!
         * @}
         */

        /*! Optional WellGroup associated with the geomodel
         * @todo Move it out. It has nothing to do here. [JP]
         */
        const WellGroup* wells_ ;
    } ;

    class GeoModelAccess {
    ringmesh_disable_copy( GeoModelAccess ) ;
        friend class GeoModelBuilder ;
        friend class GeoModelBuilderGM ;
        friend class GeoModelBuilderTopology ;
        friend class GeoModelBuilderGeometry ;
        friend class GeoModelBuilderGeology ;
        friend class GeoModelBuilderRemoval ;
        friend class GeoModelBuilderRepair ;
        friend class GeoModelBuilderCopy ;
        friend class GeoModelBuilderInfo ;
        friend class GeoModelBuilderFromSurfaces ;

    private:
        GeoModelAccess( GeoModel& geomodel )
            : geomodel_( geomodel )
        {
        }

        std::string& modifiable_name()
        {
            return geomodel_.geomodel_name_ ;
        }

        EntityTypeManager& modifiable_entity_type_manager()
        {
            return geomodel_.entity_type_manager_ ;
        }

        std::vector< GeoModelMeshEntity* >& modifiable_mesh_entities(
            const EntityType& type )
        {
            return const_cast< std::vector< GeoModelMeshEntity* >& >( geomodel_.mesh_entities(
                type ) ) ;
        }

        GeoModelMeshEntity& modifiable_mesh_entity( const gme_t& id )
        {
            return *modifiable_mesh_entities( id.type )[id.index] ;
        }

        std::vector< std::vector< GeoModelGeologicalEntity* > >& modifiable_geological_entities()
        {
            return geomodel_.geological_entities_ ;
        }

        std::vector< GeoModelGeologicalEntity* >& modifiable_geological_entities(
            const EntityType& type )
        {
            return const_cast< std::vector< GeoModelGeologicalEntity* >& >( geomodel_.geological_entities(
                type ) ) ;
        }

        GeoModelGeologicalEntity& modifiable_geological_entity( const gme_t& id )
        {
            return *modifiable_geological_entities( id.type )[id.index] ;
        }

        Universe& modifiable_universe()
        {
            return geomodel_.universe_ ;
        }

        double& modifiable_epsilon()
        {
            return geomodel_.epsilon_ ;
        }

    private:
        GeoModel& geomodel_ ;
    } ;

    typedef GEO::Factory1< GeoModelGeologicalEntity, GeoModel > GeoModelGeologicalEntityFactory ;

#define ringmesh_register_GeoModelGeologicalEntity_creator( type ) \
    geo_register_creator( GeoModelGeologicalEntityFactory, type, type::type_name_static() )

}

#endif
