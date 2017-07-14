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

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @brief Builder tools to edit and build GeoModel topology
 * (connectivity, entity creation and deletion, ...).
 * @author Pierre Anquez
 */

namespace RINGMesh {
    class GeoModelBuilder;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderTopology {
    ringmesh_disable_copy( GeoModelBuilderTopology );
        friend class GeoModelBuilder;

    public:
        /*!
         * @brief Copy topological information from a geomodel
         * @details Copy all the geomodel entities and their relationship
         * ignoring their geometry
         * @param[in] from Model to copy the information from
         */
        void copy_topology( const GeoModel& from );

        /*!
         * @brief Add to the vector the entities which cannot exist if
         *        an entity in the set does not exist.
         * @return True if at least one entity was added.
         * @todo I do not think that it work for only else than region and layer
         * because if you remove something else, for instance a corner or a line
         * or a surface, the incident boundaries may still exist but you removed
         * one of its boundaries... to handle in the future [BC].
         */
        bool get_dependent_entities(
            std::set< gmme_id >& in_mesh_entities,
            std::set< gmge_id >& in_geological_entities ) const;

        template< typename ENTITY >
        gmme_id create_mesh_entity( const MeshType mesh_type = "" )
        {
            const MeshEntityType entity_type = ENTITY::type_name_static();
            index_t nb_entities( geomodel_.nb_mesh_entities( entity_type ) );
            index_t new_id( nb_entities );
            geomodel_access_.modifiable_mesh_entities( entity_type ).emplace_back(
                GeoModelMeshEntityAccess::create_entity< ENTITY >( geomodel_, new_id,
                    mesh_type ) );
            return geomodel_access_.modifiable_mesh_entities( entity_type ).back()->gmme();
        }

        bool create_mesh_entities(
            const MeshEntityType& type,
            index_t nb_additional_entities )
        {
            if( MeshEntityTypeManager::is_corner( type ) ) {
                return create_mesh_entities< Corner >( nb_additional_entities );
            } else if( MeshEntityTypeManager::is_line( type ) ) {
                return create_mesh_entities< Line >( nb_additional_entities );
            } else if( MeshEntityTypeManager::is_surface( type ) ) {
                return create_mesh_entities< Surface >( nb_additional_entities );
            } else if( MeshEntityTypeManager::is_region( type ) ) {
                return create_mesh_entities< Region >( nb_additional_entities );
            } else {
                ringmesh_assert_not_reached;
                return false;
            }
        }

        void remove_mesh_entity_boundary_relation(
            const gmme_id& incident_entity,
            const gmme_id& boundary );

        void add_mesh_entity_boundary_relation(
            const gmme_id& boundary,
            const gmme_id& incident_entity,
            bool side = false );

        void set_mesh_entity_boundary(
            const gmme_id& gme_id,
            index_t id,
            index_t boundary_id,
            bool side = false );

        void add_universe_boundary( index_t boundary_id, bool side );

        void set_universe_boundary( index_t id, index_t boundary_id, bool side );

        void set_mesh_entity_incident_entity(
            const gmme_id& gme_id,
            index_t id,
            index_t incident_entity_id );

        void delete_mesh_entity( const MeshEntityType& type, index_t index );

        /*!
         * @brief Finds or creates a corner at given coordinates.
         * @param[in] point Geometric location of the Corner
         * @return Index of the Corner
         */
        gmme_id find_or_create_corner( const vec3& point );
        gmme_id find_or_create_corner( index_t geomodel_point_id );

        /*!
         * @brief Finds or creates a line
         * @param[in] vertices Coordinates of the vertices of the line
         * @return Index of the Line
         */
        gmme_id find_or_create_line( const std::vector< vec3 >& vertices );

        /*!
         * @brief Finds or creates a line knowing its topological adjacencies
         */
        gmme_id find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gmme_id& first_corner,
            const gmme_id& second_corner );

        void compute_universe();

    private:
        GeoModelBuilderTopology( GeoModelBuilder& builder, GeoModel& geomodel );

        template< typename ENTITY >
        bool create_mesh_entities(
            index_t nb_additionnal_entities,
            const MeshType type = "" )
        {
            const MeshEntityType& entity_type = ENTITY::type_name_static();
            std::vector< std::unique_ptr< GeoModelMeshEntity > >& store =
                geomodel_access_.modifiable_mesh_entities( entity_type );
            index_t old_size = static_cast< index_t >( store.size() );
            index_t new_size = old_size + nb_additionnal_entities;
            store.reserve( new_size );
            for( index_t i = old_size; i < new_size; i++ ) {
                store.emplace_back(
                    GeoModelMeshEntityAccess::create_entity< ENTITY >( geomodel_, i,
                        type ) );
            }
            return true;
        }

        template< typename ENTITY >
        void copy_mesh_entity_topology( const GeoModel& from )
        {
            const MeshEntityType& type = ENTITY::type_name_static();
            create_mesh_entities< ENTITY >( from.nb_mesh_entities( type ) );

            RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < geomodel_.nb_mesh_entities( type ); ++e ) {
                gmme_id id( type, e );
                GeoModelMeshEntityAccess gmme_access(
                    geomodel_access_.modifiable_mesh_entity( id ) );
                gmme_access.copy( from.mesh_entity( id ) );
            }
        }

        index_t check_if_boundary_incident_entity_relation_already_exists(
                const gmme_id& incident_entity,
                const gmme_id& boundary ) ;

    private:
        GeoModelBuilder& builder_;
        GeoModel& geomodel_;
        GeoModelAccess geomodel_access_;
    };

}
