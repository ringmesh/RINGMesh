/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geomodel/builder/geomodel_builder_access.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_base.h>

/*!
 * @file ringmesh/geomodel/builder/geomodel_builder_access.cpp
 * @brief Implementation of the classes to access the GeoModel and its Entities
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    template <>
    std::vector< bool >& GeoModelMeshEntityAccess< 2 >::modifiable_sides()
    {
        ringmesh_assert( gmme_.type_name() == Surface2D::type_name_static() );
        return dynamic_cast< Surface2D& >( gmme_ ).sides_;
    }

    template <>
    std::vector< bool >& GeoModelMeshEntityAccess< 3 >::modifiable_sides()
    {
        ringmesh_assert( gmme_.type_name() == Region3D::type_name_static() );
        return dynamic_cast< Region3D& >( gmme_ ).sides_;
    }

    template < index_t DIMENSION >
    void GeoModelMeshEntityAccess< DIMENSION >::change_mesh_data_structure(
        const MeshType& type )
    {
        if( gmme_.mesh_->type_name() != type )
        {
            gmme_.unbind_vertex_mapping_attribute();
            gmme_.change_mesh_data_structure( type );
            gmme_.bind_vertex_mapping_attribute();
        }
    }
    template < index_t DIMENSION >
    std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > >
        GeoModelGeologicalEntityAccess< DIMENSION >::create_geological_entity(
            const GeologicalEntityType& type,
            const GeoModel< DIMENSION >& geomodel,
            index_t index_in_geomodel )
    {
        auto GMGE = GeoModelGeologicalEntityFactory< DIMENSION >::create(
            type, geomodel );
        GMGE->id_ = index_in_geomodel;
        return GMGE;
    }

    template < index_t DIMENSION >
    std::string& GeoModelAccess< DIMENSION >::modifiable_name()
    {
        return geomodel_.geomodel_name_;
    }

    template < index_t DIMENSION >
    EntityTypeManager< DIMENSION >&
        GeoModelAccess< DIMENSION >::modifiable_entity_type_manager()
    {
        return geomodel_.entity_type_manager_;
    }

    template < index_t DIMENSION >
    std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >&
        GeoModelAccess< DIMENSION >::modifiable_mesh_entities(
            const MeshEntityType& type )
    {
        return const_cast< std::vector<
            std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >& >(
            geomodel_.mesh_entities( type ) );
    }

    template < index_t DIMENSION >
    GeoModelMeshEntity< DIMENSION >&
        GeoModelAccess< DIMENSION >::modifiable_mesh_entity( const gmme_id& id )
    {
        return *modifiable_mesh_entities( id.type() )[id.index()];
    }

    template < index_t DIMENSION >
    std::vector< std::vector<
        std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > > >&
        GeoModelAccess< DIMENSION >::modifiable_geological_entities()
    {
        return geomodel_.geological_entities_;
    }

    template < index_t DIMENSION >
    std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >&
        GeoModelAccess< DIMENSION >::modifiable_geological_entities(
            const GeologicalEntityType& type )
    {
        return const_cast< std::vector<
            std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& >(
            geomodel_.geological_entities( type ) );
    }

    template < index_t DIMENSION >
    GeoModelGeologicalEntity< DIMENSION >&
        GeoModelAccess< DIMENSION >::modifiable_geological_entity(
            const gmge_id& id )
    {
        return *modifiable_geological_entities( id.type() )[id.index()];
    }

    template < index_t DIMENSION >
    double& GeoModelAccess< DIMENSION >::modifiable_epsilon()
    {
        return geomodel_.epsilon_;
    }

    template class geomodel_builder_api GeoModelMeshEntityAccess< 2 >;
    template class geomodel_builder_api GeoModelGeologicalEntityAccess< 2 >;
    template class geomodel_builder_api GeoModelAccess< 2 >;

    template class geomodel_builder_api GeoModelMeshEntityAccess< 3 >;
    template class geomodel_builder_api GeoModelGeologicalEntityAccess< 3 >;
    template class geomodel_builder_api GeoModelAccess< 3 >;

} // namespace RINGMesh
