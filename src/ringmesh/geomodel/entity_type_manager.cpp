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

    template< index_t DIMENSION >
    MeshEntityTypeBoundaryMap< DIMENSION >::MeshEntityTypeBoundaryMap()
    {
        initialize_base();
    }

    template< index_t DIMENSION >
    void MeshEntityTypeBoundaryMap< DIMENSION >::initialize_base()
    {
        register_boundary( Corner< DIMENSION >::type_name_static(),
            ForbiddenMeshEntityType::type_name_static() );
        register_boundary( Line< DIMENSION >::type_name_static(),
            Corner< DIMENSION >::type_name_static() );
        register_boundary( Surface< DIMENSION >::type_name_static(),
            Line< DIMENSION >::type_name_static() );
    }

    template< >
    MeshEntityTypeBoundaryMap< 3 >::MeshEntityTypeBoundaryMap()
    {
        initialize_base();
        register_boundary( Region< 3 >::type_name_static(),
            Surface< 3 >::type_name_static() );
    }

    template< index_t DIMENSION >
    void MeshEntityTypeIncidentEntityMap< DIMENSION >::initialize_base()
    {
        register_incident_entity( Corner< DIMENSION >::type_name_static(),
            Line< DIMENSION >::type_name_static() );
        register_incident_entity( Line< DIMENSION >::type_name_static(),
            Surface< DIMENSION >::type_name_static() );
    }

    template< >
    MeshEntityTypeIncidentEntityMap< 3 >::MeshEntityTypeIncidentEntityMap()
    {
        initialize_base();
        register_incident_entity( Surface< 3 >::type_name_static(),
            Region< 3 >::type_name_static() );
        register_incident_entity( Region< 3 >::type_name_static(),
            ForbiddenMeshEntityType::type_name_static() );
    }

    template< >
    MeshEntityTypeIncidentEntityMap< 2 >::MeshEntityTypeIncidentEntityMap()
    {
        initialize_base();
        register_incident_entity( Surface< 2 >::type_name_static(),
            ForbiddenMeshEntityType::type_name_static() );
    }

    template< index_t DIMENSION >
    MeshEntityTypes< DIMENSION >::MeshEntityTypes()
    {
        initialize_base();
    }

    template< index_t DIMENSION >
    void MeshEntityTypes< DIMENSION >::initialize_base()
    {
        mesh_entity_types_.push_back( Corner< DIMENSION >::type_name_static() );
        mesh_entity_types_.push_back( Line< DIMENSION >::type_name_static() );
        mesh_entity_types_.push_back( Surface< DIMENSION >::type_name_static() );
    }

    template< >
    MeshEntityTypes< 3 >::MeshEntityTypes()
    {
        initialize_base();
        mesh_entity_types_.push_back( Region< 3 >::type_name_static() );
    }

    index_t GeologicalTypeManager::nb_geological_entity_types() const
    {
        return static_cast< index_t >( geological_entity_types_.size() );
    }

    const std::vector< GeologicalEntityType >& GeologicalTypeManager::geological_entity_types() const
    {
        return geological_entity_types_;
    }

    const GeologicalEntityType& GeologicalTypeManager::geological_entity_type(
        index_t index ) const
    {
        return geological_entity_types_.at( index );
    }

    index_t GeologicalTypeManager::geological_entity_type_index(
        const GeologicalEntityType& type ) const
    {
        return find( geological_entity_types_, type );
    }

    bool GeologicalTypeManager::is_valid_type(
        const GeologicalEntityType& type ) const
    {
        return contains( geological_entity_types_, type );
    }

    std::vector< GeologicalEntityType > RelationshipManager::parent_types(
        const MeshEntityType& child_type ) const
    {
        MeshEntityToParents::const_iterator itr = child_to_parents_.find(
            child_type );
        std::vector< GeologicalEntityType > result;
        if( itr != child_to_parents_.end() ) {
            result.insert( result.begin(), itr->second.begin(), itr->second.end() );
        }
        return result;
    }

    index_t RelationshipManager::nb_parent_types(
        const MeshEntityType& child_type ) const
    {
        return static_cast< index_t >( parent_types( child_type ).size() );
    }

    const MeshEntityType RelationshipManager::child_type(
        const GeologicalEntityType& parent_type ) const
    {
        GeologicalEntityToChild::const_iterator itr = parent_to_child_.find(
            parent_type );
        if( itr == parent_to_child_.end() ) {
            return ForbiddenMeshEntityType::type_name_static();
        } else {
            return itr->second;
        }
    }

    template< > MeshEntityTypes< 2 > RINGMESH_API MeshEntityTypeManagerBase< 2 >::mesh_entity_types_ { };
    template< > MeshEntityTypeBoundaryMap< 2 > RINGMESH_API MeshEntityTypeManagerBase<
        2 >::boundary_relationships_ { };
    template< > MeshEntityTypeIncidentEntityMap< 2 > RINGMESH_API MeshEntityTypeManagerBase<
        2 >::incident_entity_relationships_ { };
    template< > MeshEntityTypes< 3 > RINGMESH_API MeshEntityTypeManagerBase< 3 >::mesh_entity_types_ { };
    template< > MeshEntityTypeBoundaryMap< 3 > RINGMESH_API MeshEntityTypeManagerBase<
        3 >::boundary_relationships_ { };
    template< > MeshEntityTypeIncidentEntityMap< 3 > RINGMESH_API MeshEntityTypeManagerBase<
        3 >::incident_entity_relationships_ { };

    template class RINGMESH_API MeshEntityTypes< 2 > ;
    template class RINGMESH_API MeshEntityTypeManagerBase< 2 > ;
    template class RINGMESH_API MeshEntityTypeManager< 2 > ;
    template class RINGMESH_API MeshEntityTypeIncidentEntityMap< 2 > ;
    template class RINGMESH_API MeshEntityTypeBoundaryMap< 2 > ;

    template class RINGMESH_API MeshEntityTypes< 3 > ;
    template class RINGMESH_API MeshEntityTypeManagerBase< 3 > ;
    template class RINGMESH_API MeshEntityTypeManager< 3 > ;
    template class RINGMESH_API MeshEntityTypeIncidentEntityMap< 3 > ;
    template class RINGMESH_API MeshEntityTypeBoundaryMap< 3 > ;
}
