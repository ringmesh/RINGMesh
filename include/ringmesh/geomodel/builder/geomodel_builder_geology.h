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

#pragma once

#include <ringmesh/geomodel/builder/common.h>
#include <ringmesh/geomodel/builder/geomodel_builder_access.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>

/*!
 * @brief Classes to build GeoModel geological entities
 * @author Jeanne Pellerin
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderGeology final
    {
        ringmesh_disable_copy_and_move( GeoModelBuilderGeology );
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
            const GeologicalEntityType& type, index_t nb_additional_entities );

        void set_geological_entity_geol_feature( const gmge_id& gmge_id,
            typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE
                geol_feature );

        void set_mesh_entity_parent( const gmme_id& child_gmme,
            index_t id,
            const gmge_id& parent_gmge,
            bool update_parent = false );

        void add_parent_children_relation(
            const gmge_id& parent, const gmme_id& children );

        void remove_parent_children_relation(
            const gmge_id& parent, const gmme_id& children );
        void set_geological_entity_child(
            const gmge_id& parent_gmge, index_t id, index_t child_id );

        void delete_geological_entity(
            const GeologicalEntityType& type, index_t index );

        /*!
         * @brief Build the Contacts
         * @details One contact is a group of lines shared by the same
         * Interfaces
         */
        void build_contacts();

    protected:
        GeoModelBuilderGeology( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );
        ~GeoModelBuilderGeology() = default;

    private:
        index_t create_geological_entity_type(
            const GeologicalEntityType& type );

        index_t find_or_create_geological_entity_type(
            const GeologicalEntityType& type );

        void copy_geological_entity_topology( const GeoModel< DIMENSION >& from,
            const GeologicalEntityType& type );

        bool check_if_children_to_parent_relation_already_exists(
            const gmge_id& parent, const gmme_id& children );

        bool check_if_parent_to_children_relation_already_exists(
            const gmge_id& parent, const gmme_id& children );

        void update_parent_entity_children( index_t relationship_id,
            const gmge_id& new_parent_gmge,
            const gmge_id& old_parent_gmge );

    private:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };
} // namespace RINGMesh
