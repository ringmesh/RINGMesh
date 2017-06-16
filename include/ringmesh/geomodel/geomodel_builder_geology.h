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
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/geomodel/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/geomodel_builder_remove.h>
#include <ringmesh/geomodel/geomodel_builder_repair.h>
#include <ringmesh/geomodel/geomodel_builder_topology.h>

/*!
 * @brief Classes to build GeoModel geological entities
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModelBuilderBase;
    template< index_t DIMENSION > class GeoModelBuilder;
}

namespace RINGMesh {

    template< index_t DIMENSION >
    class GeoModelBuilderGeology final {
    ringmesh_disable_copy( GeoModelBuilderGeology );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;

    public:
        void copy_geology( const GeoModel< DIMENSION >& from );

        /*!
         * @brief Create and store a geological entity of the given type
         * @return The index of the created geological entity
         */
        gmge_id create_geological_entity( const GeologicalEntityType& type );

        bool create_geological_entities(
            const GeologicalEntityType& type,
            index_t nb );

        void set_geological_entity_geol_feature(
            const gmge_id& gmge_id,
            typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE geol_feature );

        void set_mesh_entity_parent(
            const gmme_id& child_gmme,
            index_t id,
            const gmge_id& parent_gmge )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity< DIMENSION >& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( child_gmme );
            ringmesh_assert( id < mesh_entity.nb_parents() );
            GeoModelMeshEntityAccess< DIMENSION > gmme_access( mesh_entity );
            index_t relationship_id = gmme_access.modifiable_parents()[id];
            RelationshipManager& manager =
                geomodel_access_.modifiable_entity_type_manager().relationship_manager;
            manager.set_parent_to_parent_child_relationship( relationship_id,
                parent_gmge );
        }

        void add_parent_children_relation(
            const gmge_id& parent,
            const gmme_id& children );

        void remove_parent_children_relation(
            const gmge_id& parent,
            const gmme_id& children );
        void set_geological_entity_child(
            const gmge_id& parent_gmge,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity< DIMENSION >& geol_entity =
                geomodel_access_.modifiable_geological_entity( parent_gmge );
            const MeshEntityType& child_type =
                geomodel_.entity_type_manager().relationship_manager.child_type(
                    parent_gmge.type() );
            gmme_id child( child_type, child_id );
            GeoModelGeologicalEntityAccess< DIMENSION > gmge_access( geol_entity );
            index_t relationship_id = gmge_access.modifiable_children()[id];
            RelationshipManager& manager =
                geomodel_access_.modifiable_entity_type_manager().relationship_manager;
            manager.set_child_to_parent_child_relationship( relationship_id, child );
        }

        void delete_geological_entity(
            const GeologicalEntityType& type,
            index_t index );

        /*!
         * @brief Build the Contacts
         * @details One contact is a group of lines shared by the same Interfaces
         */
        void build_contacts();

    protected:
        GeoModelBuilderGeology(
            GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

    private:
        index_t create_geological_entity_type( const GeologicalEntityType& type );

        index_t find_or_create_geological_entity_type(
            const GeologicalEntityType& type );

        void copy_geological_entity_topology(
            const GeoModel< DIMENSION >& from,
            const GeologicalEntityType& type );

        bool check_if_boundary_incident_entity_relation_already_exists(
            const gmge_id& parent,
            const gmme_id& children );

    private:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };
}
