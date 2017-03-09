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

#include <ringmesh/geomodel/entity_type_manager.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

namespace RINGMesh {

    // Not the smartest but hopefully compiles in C++98
    static const MeshEntityType hard_encoded_mesh_entity_types_array[4] = {
        Corner::type_name_static(), Line::type_name_static(),
        Surface::type_name_static(), Region::type_name_static() } ;
    static const std::vector< MeshEntityType > hard_encoded_mesh_entity_types(
        &hard_encoded_mesh_entity_types_array[0],
        &hard_encoded_mesh_entity_types_array[4] ) ;

    MeshEntityTypeManager::MeshEntityTypeManager( EntityTypeManager& type_manager )
        :
            type_manager_( type_manager ),
            boundary_relationships_(),
            in_boundary_relationships_()
    {

    }

    bool MeshEntityTypeManager::is_corner( const MeshEntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[0] ;
    }
    bool MeshEntityTypeManager::is_line( const MeshEntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[1] ;
    }
    bool MeshEntityTypeManager::is_surface( const MeshEntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[2] ;
    }
    bool MeshEntityTypeManager::is_region( const MeshEntityType& type )
    {
        return type == hard_encoded_mesh_entity_types[3] ;
    }

    bool MeshEntityTypeManager::is_valid_type(const MeshEntityType& type ) {
        return find( hard_encoded_mesh_entity_types, type ) != NO_ID ;
    }

    const MeshEntityType& MeshEntityTypeManager::boundary_type(
        const MeshEntityType& mesh_entity_type )
    {
        MeshEntityTypeMap::const_iterator itr = boundary_relationships_.map.find(
            mesh_entity_type ) ;
        ringmesh_assert( itr != boundary_relationships_.map.end() ) ;
        return itr->second ;
    }
    const MeshEntityType& MeshEntityTypeManager::in_boundary_type(
        const MeshEntityType& mesh_entity_type )
    {
        MeshEntityTypeMap::const_iterator itr = in_boundary_relationships_.map.find(
            mesh_entity_type ) ;
        ringmesh_assert( itr != in_boundary_relationships_.map.end() ) ;
        return itr->second ;
    }

    const std::vector< MeshEntityType >& MeshEntityTypeManager::mesh_entity_types()
    {
        return hard_encoded_mesh_entity_types ;
    }
    index_t MeshEntityTypeManager::nb_mesh_entity_types()
    {
        return static_cast< index_t >( hard_encoded_mesh_entity_types.size() ) ;

    }

    GeologicalTypeManager::GeologicalTypeManager( EntityTypeManager& type_manager )
        : type_manager_( type_manager )
    {

    }

    index_t GeologicalTypeManager::nb_geological_entity_types() const
    {
        return static_cast< index_t >( geological_entity_types_.size() ) ;
    }
    const std::vector< GeologicalEntityType >& GeologicalTypeManager::geological_entity_types() const
    {
        return geological_entity_types_ ;
    }
    const GeologicalEntityType& GeologicalTypeManager::geological_entity_type(
        index_t index ) const
    {
        return geological_entity_types_.at( index ) ;
    }
    index_t GeologicalTypeManager::geological_entity_type_index(
        const GeologicalEntityType& type ) const
    {
        return find( geological_entity_types_, type ) ;
    }

    RelationshipManager::RelationshipManager( EntityTypeManager& type_manager )
        : type_manager_( type_manager )
    {

    }

    std::vector< GeologicalEntityType > RelationshipManager::parent_types(
        const MeshEntityType& child_type ) const
    {
        MeshEntityToParents::const_iterator itr = child_to_parents_.find(
            child_type ) ;
        std::vector< GeologicalEntityType > result ;
        if( itr != child_to_parents_.end() ) {
            result.insert( result.begin(), itr->second.begin(), itr->second.end() ) ;
        }
        return result ;
    }
    index_t RelationshipManager::nb_parent_types(
        const MeshEntityType& child_type ) const
    {
        return static_cast< index_t >( parent_types( child_type ).size() ) ;
    }
    const MeshEntityType RelationshipManager::child_type(
        const GeologicalEntityType& parent_type ) const
    {
        GeologicalEntityToChild::const_iterator itr = parent_to_child_.find(
            parent_type ) ;
        if( itr == parent_to_child_.end() ) {
            return DefaultMeshEntityType::default_entity_type() ;
        } else {
            return itr->second ;
        }
    }

}
