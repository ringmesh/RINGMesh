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
#include <ringmesh/algorithm.h>

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
     * @brief Manages the type relationship between GeoModelEntities
     * One instance owned by each GeoModel.
     */
    class RINGMESH_API EntityRelationships {        
    public:
        friend GeoModelEditor ;

        typedef std::string EntityType ;
        typedef std::string MeshEntityType ;
        typedef std::string GeologicalEntityType ;
        typedef std::map< MeshEntityType, std::set< GeologicalEntityType > > MeshEntityToParents ;
        typedef std::map< GeologicalEntityType, MeshEntityType > GeologicalEntityToChild ;

        static bool is_valid_type( const EntityType& type )
        {
            // Change the name of the function 
            // True does not mean that this type is valid for this model
            return type != default_entity_type() ;
        }               
        static const EntityType default_entity_type() ;

        // ---- static memeber for relationships between MeshEntities         
        static bool is_mesh_entity_type( const EntityType& type ) ;
        static const MeshEntityType& boundary_type( const MeshEntityType& type ) ;
        static const MeshEntityType& in_boundary_type( const MeshEntityType& type ) ;
        static const std::vector< MeshEntityType >& mesh_entity_types() ;
        static index_t nb_mesh_entity_types() ;

        const std::set< GeologicalEntityType >& parent_types( const MeshEntityType& child_type ) const 
        {
            MeshEntityToParents::const_iterator
                itr = child_to_parents_.find( child_type );
            ringmesh_assert( itr != child_to_parents_.end() ) ;
            return itr->second ;
        }
        index_t nb_parent_types( const MeshEntityType& child_type ) const
        {
            MeshEntityToParents::const_iterator itr =
                child_to_parents_.find( child_type ) ;
            if( itr == child_to_parents_.end() ) {
                return 0 ;
            }
            return itr->second.size() ;
        }
        const MeshEntityType& child_type( const GeologicalEntityType& parent_type ) const
        {
           GeologicalEntityToChild::const_iterator
                itr = parent_to_child_.find( parent_type );
           ringmesh_assert( itr != parent_to_child_.end() ) ;
           return itr->second ;
        }
        
        // --- Geological entities depend on the GeoModel               
        bool is_geological_entity_type( const EntityType& type ) const ;
        index_t nb_geological_entity_types() const
        {
            return geological_entity_types_.size() ;
        }
        const EntityType& geological_entity_type( index_t index ) const
        {
            ringmesh_assert( index < nb_geological_entity_types() ) ;
            return geological_entity_types_[index] ;
        }
        index_t geological_entity_type_index( const EntityType& type ) const
        {
            return find( geological_entity_types_, type ) ;
        }

    private:
        void register_relationship( const EntityType& parent_type_name,
            const EntityType& child_type_name )
        {
            register_child_type( parent_type_name, child_type_name ) ;
            register_parent_type( parent_type_name, child_type_name ) ;
        }
        void register_child_type( const EntityType& parent_type_name, 
            const EntityType& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name ;
        }
        void register_parent_type( const EntityType& parent_type_name,
            const EntityType& child_type_name )
        {
            child_to_parents_[child_type_name].insert( parent_type_name ) ;
        }

    private:
        GeologicalEntityToChild parent_to_child_ ;
        MeshEntityToParents child_to_parents_ ; 

        std::vector< GeologicalEntityType > geological_entity_types_ ;
    };



    /*!
     * @brief The class to describe a geological model represented 
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    class RINGMESH_API GeoModel {    
    public:
        friend class GeoModelBuilder ;
        friend class GeoModelEditor ;
        friend class GeoModelRepair ;

        typedef std::string EntityType ;

        /*!
         * @brief Constructs an empty GeoModel
         */
        GeoModel() ;

        /*!
         * @brief Deletes all the GeoModelEntities of the GeoModel
         */
        virtual ~GeoModel() ;   

        const std::string& name() const
        {
            return name_ ;
        }

        const EntityRelationships& entity_relationships() const
        {
            return entity_relationships_ ;
        }
        bool is_mesh_entity_type( const EntityType& type ) const
        {
            return EntityRelationships::is_mesh_entity_type( type ) ;
        }
        bool is_geological_entity_type( const EntityType& type ) const
        {
            return entity_relationships_.is_geological_entity_type( type ) ;
        }

        index_t nb_entities( const EntityType& type ) const
        {
            if( is_mesh_entity_type( type ) ) {
                return nb_mesh_entities( type );
            }
            else if( is_geological_entity_type( type ) ) {
                return nb_geological_entities( type ) ;
            } else {
                return 0 ;
            }
        }
        /*!
         * @brief Returns the number of mesh entities of the given type
         * @details Default value is 0
         * @param[in] type the mesh entity type
         */
        index_t nb_mesh_entities( const EntityType& type ) const
        {
            if( type == Corner::type_name_static() ) {
                return nb_corners();
            } else if( type == Line::type_name_static() ) {
                return nb_lines();
            } else if( type == Surface::type_name_static() ) {
                return nb_surfaces();
            } else if( type == Region::type_name_static() ) {
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
        index_t nb_geological_entities( const EntityType& type ) const
        {
            index_t index = geological_entity_type_index( type ) ;
            if( index == NO_ID ) return 0 ;
            return static_cast< index_t >( geological_entities_[index].size() ) ;
        }

        /*!
         * @brief Returns the index of the geological entity type storage
         * @details Default value is NO_ID
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entity_types() const
        {
            return entity_relationships_.nb_geological_entity_types() ;
        }
        const EntityType& geological_entity_type( index_t index ) const
        {
            return entity_relationships_.geological_entity_type( index ) ;
        }

/*        static index_t nb_mesh_entity_types() {
            return 4 ;
        }
        // Hard encoded mapping between for meshed entities. Always the same.
/*        index_t mesh_entity_type( const EntityType& type ) const
        {
            if( type == Corner::type_name_static() ) {
                return 0 ;
            } else if( type == Line::type_name_static() ) {
                return 1 ;
            } else if( type == Surface::type_name_static() ) {
                return 2 ; 
            } else if( type == Region::type_name_static() ) {
                return 3 ;
            } else {
                return NO_ID ;
            }
        }
/*        const EntityType mesh_entity_type( index_t id )
        {
            ringmesh_assert( id < nb_mesh_entity_types() ) ;            
            return mesh_entity_types_[id] ;
        }  
 */
        
        void assert_gme_id_validity( GME::gme_t id )
        {
            bool is_valid_type = is_mesh_entity_type( id.type )
                || is_geological_entity_type( id.type ) ;            
            ringmesh_assert( is_valid_type ) ;

            bool is_valid_index = id.index < nb_entities( id.type ) ;
            ringmesh_assert( is_valid_index ) ;
        }

        /*!
         * @brief Returns a const reference the identified GeoModelEntity
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         * @pre Entity identification is valid.
         */
        const GeoModelGeologicalEntity& geological_entity( GME::gme_t id ) const
        {
            ringmesh_assert( id.index < nb_geological_entities( id.type ) ) ;
            return *geological_entities( id.type )[id.index] ;
        }
        /*!
         * Convenient overload of entity( GME::gme_t id )
         */
        const GeoModelGeologicalEntity& geological_entity(
            const EntityType& entity_type,
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
            const EntityType& type = id.type ;
            index_t index = id.index ;
            if( type == Corner::type_name_static() ) {
                return corner( index ) ;
            } else if( type == Line::type_name_static() ) {
                return line( index ) ;
            } else if( type == Surface::type_name_static() ) {
                return surface( index ) ;
            } else if( type == Region::type_name_static() ) {
                return region( index ) ;
            }
            ringmesh_assert_not_reached ;
            return surface( 0 );
        }
        /*!
         * Convenient overload of mesh_entity( GME::gme_t id )
         */
        const GeoModelMeshEntity& mesh_entity(
            const EntityType& entity_type,
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
        // Maybe implement them properly one day [JP]
        // For now, prevent any copy as a GeoModel manages a lot of memory
        ringmesh_disable_copy( GeoModel ) ;

        /*!
         * Access to the position of the entity of that type
         * in the private storage
         */
        index_t geological_entity_type_index( const EntityType& type ) const
        {
            return entity_relationships_.geological_entity_type_index( type ) ;
        }

        // I do know that this casts are really ugly. But we still have a big big design issue [JP]
        const std::vector< GeoModelEntity* >& entities( const EntityType& type )
        {
            if( is_mesh_entity_type( type ) ) {
                return *(std::vector< GeoModelEntity* > *) (&mesh_entities( type )) ;
            } else if( is_geological_entity_type( type ) ) {
                return *(std::vector< GeoModelEntity* > *) (&geological_entities( type )) ;
            } else {
                ringmesh_assert_not_reached ;
                return *(std::vector< GeoModelEntity* > *) &surfaces_ ;
            }
        }
        
        /*!
         * @brief Generic accessor to the storage of mesh entities of the given type
         */
        const std::vector< GeoModelMeshEntity* >& mesh_entities(
            const EntityType& type ) const
        {
            if( type == Corner::type_name_static() ) {
                return *(std::vector< GeoModelMeshEntity* > *) &corners_ ;
            } else if( type == Line::type_name_static() ) {
                return *(std::vector< GeoModelMeshEntity* > *) &lines_ ;
            } else if( type == Surface::type_name_static() ) {
                return *(std::vector< GeoModelMeshEntity* > *) &surfaces_ ;
            } else if( type == Region::type_name_static() ) {
                return *(std::vector< GeoModelMeshEntity* > *) &regions_ ;
            } else {
                ringmesh_assert_not_reached ;
                return *(std::vector< GeoModelMeshEntity* > *) &surfaces_ ;
            }
        }

        /*!
         * @brief Generic accessor to the storage of geologcial entities of the given type
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

    public:
        GeoModelMesh mesh ;

    private:
        // Name of the model
        std::string name_ ;

        // Responsible for EntityTypes management
        EntityRelationships entity_relationships_ ;

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

        /*!
         * @brief Geological entities. They are optional.
         * The EntityTypes are managed by the EntityRelationships of the class.
         */
        std::vector< std::vector< GeoModelGeologicalEntity* > > geological_entities_ ;
        
        /*!
         * @}
         */

        /*! Optional WellGroup associated with the model
         * @todo Move it out. It has nothing to do here [JP]
         */
        const WellGroup* wells_ ;
    } ;

    typedef GEO::Factory1< GeoModelGeologicalEntity, GeoModel > GeoModelGeologicalEntityFactory ;

#define ringmesh_register_GeoModelGeologicalEntity_creator( type ) \
    geo_register_creator( GeoModelGeologicalEntityFactory, type, type::type_name_static() )

}

#endif
