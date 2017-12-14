/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#pragma once

#include <ringmesh/geomodel/builder/common.h>

/*!
 * @brief Classes to access the GeoModel and its Entities
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderCopy );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRepair );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGM );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometryBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometry );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoveBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemove );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderInfo );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBase );
    FORWARD_DECLARATION_DIMENSION_STRUCT( EntityTypeManager );

    class GeologicalEntityType;
    class MeshEntityType;
    ;
    struct gmge_id;
    struct gmme_id;
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelMeshEntityConstAccess
    {
        ringmesh_disable_copy_and_move( GeoModelMeshEntityConstAccess );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderGeometryBase< DIMENSION >;
        friend class GeoModelBuilderGeometry< DIMENSION >;
        friend class GeoModelBuilderTopologyBase< DIMENSION >;
        friend class GeoModelBuilderTopology< DIMENSION >;

    private:
        explicit GeoModelMeshEntityConstAccess(
            const GeoModelMeshEntity< DIMENSION >& gme )
            : gmme_( gme )
        {
        }

        ~GeoModelMeshEntityConstAccess() = default;

        const std::shared_ptr< MeshBase< DIMENSION > >& mesh() const
        {
            return gmme_.mesh_;
        }

        const std::vector< index_t >& incident_entity_relation_ids() const
        {
            return gmme_.incident_entities_;
        }

        const std::vector< index_t >& boundary_relation_ids() const
        {
            return gmme_.boundaries_;
        }

    private:
        const GeoModelMeshEntity< DIMENSION >& gmme_;
    };

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelMeshEntityAccess
    {
        ringmesh_disable_copy_and_move( GeoModelMeshEntityAccess );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderTopologyBase< DIMENSION >;
        friend class GeoModelBuilderTopology< DIMENSION >;
        friend class GeoModelBuilderGeometryBase< DIMENSION >;
        friend class GeoModelBuilderGeometry< DIMENSION >;
        friend class GeoModelBuilderGeology< DIMENSION >;
        friend class GeoModelBuilderInfo< DIMENSION >;
        friend class GeoModelBuilderRemoveBase< DIMENSION >;
        friend class GeoModelBuilderRemove< DIMENSION >;

    private:
        explicit GeoModelMeshEntityAccess(
            GeoModelMeshEntity< DIMENSION >& gme )
            : gmme_( gme )
        {
        }

        ~GeoModelMeshEntityAccess() = default;

        std::string& modifiable_name()
        {
            return gmme_.name_;
        }

        index_t& modifiable_index()
        {
            return gmme_.id_;
        }

        std::vector< index_t >& modifiable_boundaries()
        {
            return gmme_.boundaries_;
        }

        std::vector< index_t >& modifiable_incident_entities()
        {
            return gmme_.incident_entities_;
        }

        std::vector< bool >& modifiable_sides();

        std::vector< index_t >& modifiable_parents()
        {
            return gmme_.parents_;
        }

        std::shared_ptr< MeshBase< DIMENSION > >& modifiable_mesh()
        {
            return gmme_.mesh_;
        }

        void copy( const GeoModelMeshEntity< DIMENSION >& from )
        {
            gmme_.copy_mesh_entity( from );
        }

        void change_mesh_data_structure( const MeshType& type );

        template < template < index_t > class ENTITY >
        static std::unique_ptr< ENTITY< DIMENSION > > create_entity(
            const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const MeshType& type )
        {
            return std::unique_ptr< ENTITY< DIMENSION > >(
                new ENTITY< DIMENSION >( geomodel, id, type ) );
        }

    private:
        GeoModelMeshEntity< DIMENSION >& gmme_;
    };
    ALIAS_2D_AND_3D( GeoModelMeshEntityAccess );

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelGeologicalEntityAccess
    {
        ringmesh_disable_copy_and_move( GeoModelGeologicalEntityAccess );
        friend class GeoModelBuilderTopology< DIMENSION >;
        friend class GeoModelBuilderGeology< DIMENSION >;
        friend class GeoModelBuilderInfo< DIMENSION >;
        friend class GeoModelBuilderRemoveBase< DIMENSION >;

    private:
        explicit GeoModelGeologicalEntityAccess(
            GeoModelGeologicalEntity< DIMENSION >& gmge )
            : gmge_( gmge )
        {
        }
        ~GeoModelGeologicalEntityAccess() = default;

        std::string& modifiable_name()
        {
            return gmge_.name_;
        }

        index_t& modifiable_index()
        {
            return gmge_.id_;
        }

        typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE&
            modifiable_geol_feature()
        {
            return gmge_.geol_feature_;
        }

        std::vector< index_t >& modifiable_children()
        {
            return gmge_.children_;
        }

        static std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > >
            create_geological_entity( const GeologicalEntityType& type,
                const GeoModel< DIMENSION >& geomodel,
                index_t index_in_geomodel );

        void copy( const GeoModelGeologicalEntity< DIMENSION >& from )
        {
            gmge_.copy_geological_entity( from );
        }

    private:
        GeoModelGeologicalEntity< DIMENSION >& gmge_;
    };

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelAccess
    {
        ringmesh_disable_copy_and_move( GeoModelAccess );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;
        friend class GeoModelBuilderGM< DIMENSION >;
        friend class GeoModelBuilderTopologyBase< DIMENSION >;
        friend class GeoModelBuilderTopology< DIMENSION >;
        friend class GeoModelBuilderGeometryBase< DIMENSION >;
        friend class GeoModelBuilderGeometry< DIMENSION >;
        friend class GeoModelBuilderGeology< DIMENSION >;
        friend class GeoModelBuilderRemoveBase< DIMENSION >;
        friend class GeoModelBuilderRemove< DIMENSION >;
        friend class GeoModelBuilderRepair< DIMENSION >;
        friend class GeoModelBuilderCopy< DIMENSION >;
        friend class GeoModelBuilderInfo< DIMENSION >;

    private:
        explicit GeoModelAccess( GeoModel< DIMENSION >& geomodel )
            : geomodel_( geomodel )
        {
        }
        ~GeoModelAccess() = default;

        std::string& modifiable_name();

        EntityTypeManager< DIMENSION >& modifiable_entity_type_manager();

        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >&
            modifiable_mesh_entities( const MeshEntityType& type );

        GeoModelMeshEntity< DIMENSION >& modifiable_mesh_entity(
            const gmme_id& id );

        std::vector< std::vector<
            std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > > >&
            modifiable_geological_entities();

        std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >&
            modifiable_geological_entities( const GeologicalEntityType& type );

        GeoModelGeologicalEntity< DIMENSION >& modifiable_geological_entity(
            const gmge_id& id );

        double& modifiable_epsilon();

    private:
        GeoModel< DIMENSION >& geomodel_;
    };
} // namespace RINGMesh
