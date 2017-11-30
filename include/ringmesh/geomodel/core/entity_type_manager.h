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

#include <ringmesh/geomodel/core/common.h>

#include <ringmesh/basic/pimpl.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );

    class GeologicalEntityType;
    class MeshEntityType;
    struct gmme_id;
    struct gmge_id;
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief this class contains only static methods to manage the type of the
     * GeoModelMeshEntity. It gives access to the number of meshed entities of
     * each
     * type and also their (in) boundary
     */
    template < index_t DIMENSION >
    class geomodel_core_api MeshEntityTypeManagerBase
    {
    public:
        ~MeshEntityTypeManagerBase();

        bool is_corner( const MeshEntityType& type ) const;

        bool is_line( const MeshEntityType& type ) const;

        bool is_surface( const MeshEntityType& type ) const;

        bool is_valid_type( const MeshEntityType& type ) const;

        const MeshEntityType& boundary_entity_type(
            const MeshEntityType& mesh_entity_type ) const;

        const MeshEntityType& incident_entity_type(
            const MeshEntityType& mesh_entity_type ) const;

        const std::vector< MeshEntityType >& mesh_entity_types() const;

        index_t nb_mesh_entity_types() const
        {
            return static_cast< index_t >( mesh_entity_types().size() );
        }

    protected:
        MeshEntityTypeManagerBase();

    protected:
        IMPLEMENTATION_MEMBER( impl_ );
    };

    template < index_t DIMENSION >
    class geomodel_core_api MeshEntityTypeManager
        : public MeshEntityTypeManagerBase< DIMENSION >
    {
    };

    template <>
    class geomodel_core_api MeshEntityTypeManager< 3 >
        : public MeshEntityTypeManagerBase< 3 >
    {
    public:
        bool is_region( const MeshEntityType& type ) const;
    };
    ALIAS_2D_AND_3D( MeshEntityTypeManager );

    /*!
     * @brief this class contains methods to manage the type of the
     * GeoModelGeologicalEntity. It gives access to the number of geological
     * entities of each
     * type and also give the opportunity to create and manage new one.
     */
    class geomodel_core_api GeologicalTypeManager
    {
        friend class GeoModelBuilderGeology< 2 >;
        friend class GeoModelBuilderGeology< 3 >;

    public:
        index_t nb_geological_entity_types() const;
        const std::vector< GeologicalEntityType >&
            geological_entity_types() const;
        const GeologicalEntityType& geological_entity_type(
            index_t index ) const;
        index_t geological_entity_type_index(
            const GeologicalEntityType& type ) const;
        bool is_valid_type( const GeologicalEntityType& type ) const;

    private:
        std::vector< GeologicalEntityType > geological_entity_types_;

    private:
        void register_geological_entity_type(
            const GeologicalEntityType& geological_type_name );
    };

    /*!
     * @brief this class contains methods to manage relations between Geological
     * and
     * Mesh entities. For instance:
     * A "Contact" can be the parent of one or more "Line"
     * An "Interface" can the parent of one or more "Surface"
     * A "Layer" can be the parent of one or more "Region"
     *
     */
    class geomodel_core_api RelationshipManager
    {
        friend class GeoModelBuilderGeology< 2 >;
        friend class GeoModelBuilderTopologyBase< 2 >;
        friend class GeoModelBuilderGeology< 3 >;
        friend class GeoModelBuilderTopologyBase< 3 >;

    public:
        RelationshipManager();
        ~RelationshipManager();

        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const;
        index_t nb_parent_types( const MeshEntityType& child_type ) const;
        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const;

        /*! @}
         * \name Access to the Boundary/Incident entity relations
         * @{
         */

        const gmme_id& boundary_gmme( index_t relation_id ) const;

        const gmme_id& incident_entity_gmme( index_t relation_id ) const;

        /*! @}
         * \name Access to the Parent / Child relations
         * @{
         */

        const gmge_id& parent_of_gmme( index_t relation_id ) const;

        const gmme_id& child_of_gmge( index_t relation_id ) const;

    private:
        void copy( const RelationshipManager& from );

        void register_geology_relationship(
            const GeologicalEntityType& parent_type_name,
            const MeshEntityType& child_type_name );

        /*! @}
         * \name Boundary Relationship manager
         * @{
         */

        index_t add_boundary_relationship(
            const gmme_id& incident_entity, const gmme_id& boundary );

        index_t find_boundary_relationship(
            const gmme_id& incident_entity, const gmme_id& boundary );

        void set_boundary_to_boundary_relationship(
            index_t relationship_id, const gmme_id& boundary );

        void set_incident_entity_to_boundary_relationship(
            index_t relationship_id, const gmme_id& incident_entity );

        /*! @}
         * \name Parent/Child Relationship manager
         * @{
         */

        index_t add_parent_child_relationship(
            const gmge_id& parent, const gmme_id& child );

        index_t find_parent_child_relationship(
            const gmge_id& parent, const gmme_id& child );

        void set_parent_to_parent_child_relationship(
            index_t relationship_id, const gmge_id& parent );

        void set_child_to_parent_child_relationship(
            index_t relationship_id, const gmme_id& child );

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };

    /*!
     * @brief Global entity manager which could be associated to a geomodel
     * to give access to different manager to handle the entity types
     */
    template < index_t DIMENSION >
    struct geomodel_core_api EntityTypeManager
    {
        MeshEntityTypeManager< DIMENSION > mesh_entity_manager;
        GeologicalTypeManager geological_entity_manager;
        RelationshipManager relationship_manager;
    };
} // namespace RINGMesh
