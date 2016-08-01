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
 * @file Implementation of THE GeoModel
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geo_model.h>

#include <ringmesh/geo_model_mesh_entity.h>
#include <ringmesh/geo_model_geological_entity.h>


namespace RINGMesh {

    typedef std::map< EntityType, EntityType > EntityTypeMap;

    struct EntityTypeBoundaryMap
    {
        EntityTypeBoundaryMap()
        {
            register_boundary< Corner, GeoModelEntity >();
            register_boundary< Line, Corner >();
            register_boundary< Surface, Line >();
            register_boundary< Region, Surface >();
        }
        template< typename TYPE, typename BOUNDARY >
        void register_boundary()
        {
            map.insert( std::pair<EntityType, EntityType>( 
                TYPE::type_name_static(), BOUNDARY::type_name_static() ) ) ;
        }
        EntityTypeMap map ;
    };

    struct EntityTypeInBoundaryMap
    {
        EntityTypeInBoundaryMap()
        {
            register_in_boundary< Corner, Line >();
            register_in_boundary< Line, Surface >();
            register_in_boundary< Surface, Region >();
            register_in_boundary< Region, GeoModelEntity >();
        }
        template< typename TYPE, typename IN_BOUNDARY >
        void register_in_boundary()
        {
            map.insert( std::pair<EntityType, EntityType>(
                TYPE::type_name_static(), IN_BOUNDARY::type_name_static() ) ) ;
        }
        EntityTypeMap map ;
    };
   
    // Maybe these should be declared as static in the functions 
    // that use them.
    const EntityTypeBoundaryMap boundary_relationships ;
    const EntityTypeInBoundaryMap in_boundary_relationships ; 
   
     // Not the smartest but hopefully compiles in C++98   
    static const EntityType hard_encoded_mesh_entity_types_array[4] = {
        Corner::type_name_static(),
        Line::type_name_static(),
        Surface::type_name_static(),
        Region::type_name_static()} ;
    
    static const std::vector< EntityType > hard_encoded_mesh_entity_types(
        &hard_encoded_mesh_entity_types_array[0], &hard_encoded_mesh_entity_types_array[4] ) ;

    index_t EntityTypeManager::nb_mesh_entity_types()
    {
        return hard_encoded_mesh_entity_types.size() ;
    }
    bool EntityTypeManager::is_corner( const EntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[0] ;
    }
    bool EntityTypeManager::is_line( const EntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[1] ;
    }
    bool EntityTypeManager::is_surface( const EntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[2] ;
    }
    bool EntityTypeManager::is_region( const EntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[3] ;
    }    
    bool EntityTypeManager::is_mesh_entity_type( const EntityType& type )
    {
        return find( hard_encoded_mesh_entity_types, type ) != NO_ID ;
    }
    const EntityType EntityTypeManager::default_entity_type()
    {
        return "No_entity_type" ;
    }
    const std::vector< EntityType >& EntityTypeManager::mesh_entity_types()
    {
        return hard_encoded_mesh_entity_types ;
    }
    
    /*! @todo What is the cost of using such maps ?
     */
    const EntityType& EntityTypeManager::boundary_type( const EntityType& mesh_entity_type ) 
    {
        EntityTypeMap::const_iterator
            itr = boundary_relationships.map.find( mesh_entity_type );
        ringmesh_assert( itr != boundary_relationships.map.end() ) ;
        return itr->second ;
    }
    const EntityType& EntityTypeManager::in_boundary_type( const EntityType& mesh_entity_type )
    {
       EntityTypeMap::const_iterator
            itr = in_boundary_relationships.map.find( mesh_entity_type );
        ringmesh_assert( itr != in_boundary_relationships.map.end() ) ;
        return itr->second ;
    }
    
    bool EntityTypeManager::is_geological_entity_type( const EntityType& type ) const
    {
        return find( geological_entity_types_, type ) != NO_ID ;
    }

    // ------------------------------------------------------------------------//


    GeoModel::GeoModel()
        :
            mesh( *this ),
            universe_( *this ),
            wells_( nil )
    {
    }

    GeoModel::~GeoModel()
    {        
        for( index_t i = 0; i < corners_.size(); ++i ) {
            delete corners_[i] ;
        }
        for( index_t i = 0; i < lines_.size(); ++i ) {
            delete lines_[i] ;
        }
        for( index_t i = 0; i < surfaces_.size(); ++i ) {
            delete surfaces_[i] ;
        }
        for( index_t i = 0; i < regions_.size(); ++i ) {
            delete regions_[i] ;
        }

        for( index_t i = 0 ; i < geological_entities_.size(); ++i ){
            for( index_t j = 0 ; j < geological_entities_[i].size(); ++j ) {
                delete geological_entities_[i][j] ;
            }
        }
    } 

   const GeoModelEntity& GeoModel::entity( const EntityType& entity_type,
            index_t entity_index ) const
   {
       if( is_mesh_entity_type( entity_type ) ) {
           return static_cast<const GeoModelEntity&>(
               mesh_entity( entity_type, entity_index )) ;
       } else {
           return static_cast<const GeoModelEntity&>(
               geological_entity( entity_type, entity_index ));
       }
   }

   index_t GeoModel::nb_mesh_entities( const EntityType& type ) const
   {
       if( EntityTypeManager::is_corner( type ) ) {
           return nb_corners();
       } else if( EntityTypeManager::is_line( type ) ) {
           return nb_lines();
       } else if( EntityTypeManager::is_surface( type ) ) {
           return nb_surfaces();
       } else if( EntityTypeManager::is_region( type ) ) {
           return nb_regions();
       } else {
           ringmesh_assert_not_reached ;
           return 0 ;
       }
   }

   const GeoModelMeshEntity& GeoModel::mesh_entity( gme_t id ) const
   {
       const EntityType& type = id.type ;
       index_t index = id.index ;
       if( EntityTypeManager::is_corner( type ) ) {
           return corner( index ) ;
       } else if( EntityTypeManager::is_line( type ) ) {
           return line( index ) ;
       } else if( EntityTypeManager::is_surface( type ) ) {
           return surface( index ) ;
       } else if( EntityTypeManager::is_region( type ) ) {
           return region( index ) ;
       }
       ringmesh_assert_not_reached ;
       return surface( 0 );
   }

   // I do know that this casts are really ugly. But we still have a big big design issue [JP]
   const std::vector< GeoModelEntity* >& GeoModel::entities( const EntityType& type )
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

   const std::vector< GeoModelMeshEntity* >& GeoModel::mesh_entities(
       const EntityType& type ) const
   {
       if( EntityTypeManager::is_corner( type ) ) {
           return *(std::vector< GeoModelMeshEntity* > *) &corners_ ;
       } else if( EntityTypeManager::is_line( type ) ) {
           return *(std::vector< GeoModelMeshEntity* > *) &lines_ ;
       } else if( EntityTypeManager::is_surface( type ) ) {
           return *(std::vector< GeoModelMeshEntity* > *) &surfaces_ ;
       } else if( EntityTypeManager::is_region( type ) ) {
           return *(std::vector< GeoModelMeshEntity* > *) &regions_ ;
       } else {
           ringmesh_assert_not_reached ;
           return *(std::vector< GeoModelMeshEntity* > *) &surfaces_ ;
       }
   }


    /*!
     * Associates a WellGroup to the GeoModel
     * @param[in] wells the WellGroup
     * @todo Review : What is this for ?
     * @todo Extend to other object types.
     */
    void GeoModel::set_wells( const WellGroup* wells )
    {
        wells_ = wells ;
    }


} // namespace
