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
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );

    ALIAS_2D_AND_3D( GeoModelBuilder );
}

namespace RINGMesh {

    template< index_t DIMENSION >
    class GeoModelBuilderTopologyBase {
    ringmesh_disable_copy( GeoModelBuilderTopologyBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderBase< DIMENSION > ;
        friend class GeoModelBuilder< DIMENSION > ;

    public:
        virtual ~GeoModelBuilderTopologyBase() = default;
        /*!
         * @brief Copy topological information from a geomodel
         * @details Copy all the geomodel entities and their relationship
         * ignoring their geometry
         * @param[in] from Model to copy the information from
         */
        void copy_topology( const GeoModel< DIMENSION >& from );

        /*!
         * @brief Add to the vector the entities which cannot exist if
         *        an entity in the set does not exist.
         * @return True if at least one entity was added.
         */
        bool get_dependent_entities(
            std::set< gmme_id >& in_mesh_entities,
            std::set< gmge_id >& in_geological_entities ) const;

        virtual gmme_id create_mesh_entity(
            const MeshEntityType& type );

        virtual bool create_mesh_entities(
            const MeshEntityType& type,
            index_t nb_additional_entities );

        void remove_mesh_entity_boundary_relation(
            const gmme_id& incident_entity,
            const gmme_id& boundary );

        virtual void add_mesh_entity_boundary_relation(
            const gmme_id& boundary,
            const gmme_id& incident_entity,
            bool side = false );

        virtual void set_mesh_entity_boundary(
            const gmme_id& gme_id,
            index_t id,
            index_t boundary_id,
            bool side = false );

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
        gmme_id find_or_create_corner( const vecn< DIMENSION >& point );
        gmme_id find_or_create_corner( index_t geomodel_point_id );

        /*!
         * @brief Finds or creates a line
         * @param[in] vertices Coordinates of the vertices of the line
         * @return Index of the Line
         */
        gmme_id find_or_create_line(
            const std::vector< vecn< DIMENSION > >& vertices );

        /*!
         * @brief Finds or creates a line knowing its topological adjacencies
         */
        gmme_id find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gmme_id& first_corner,
            const gmme_id& second_corner );

    protected:
        GeoModelBuilderTopologyBase(
            GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

        template< template< index_t > class ENTITY >
        gmme_id create_mesh_entity( const MeshType mesh_type = "" );

        template< template< index_t > class ENTITY >
        bool create_mesh_entities(
            index_t nb_additionnal_entities,
            const MeshType type = "" );

        template< template< index_t > class ENTITY >
        void copy_mesh_entity_topology( const GeoModel< DIMENSION >& from )
        {
            const MeshEntityType& type = ENTITY< DIMENSION >::type_name_static();
            create_mesh_entities< ENTITY >( from.nb_mesh_entities( type ) );

            parallel_for( geomodel_.nb_mesh_entities( type ),
                [&type, &from, this] ( index_t i ) {
                    gmme_id id( type, i );
                    GeoModelMeshEntityAccess< DIMENSION > gmme_access(
                        geomodel_access_.modifiable_mesh_entity( id ) );
                    gmme_access.copy( from.mesh_entity( id ) );
            } );
        }

        index_t check_if_boundary_incident_entity_relation_already_exists(
            const gmme_id& incident_entity,
            const gmme_id& boundary );

        virtual void copy_all_mesh_entity_topology(
            const GeoModel< DIMENSION >& from );

    protected:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };

    ALIAS_2D_AND_3D( GeoModelBuilderTopologyBase );

    template< index_t DIMENSION >
    class GeoModelBuilderTopology: public GeoModelBuilderTopologyBase< DIMENSION > {
    };

    template< >
    class RINGMESH_API GeoModelBuilderTopology< 2 > : public GeoModelBuilderTopologyBase< 2 > {
        friend class GeoModelBuilderBase< 2 > ;
        friend class GeoModelBuilder< 2 > ;
    public:
        virtual ~GeoModelBuilderTopology() = default;

        void add_universe_boundary( index_t boundary_id, bool side );

        void set_universe_boundary( index_t id, index_t boundary_id, bool side );

        void compute_universe();

        void add_mesh_entity_boundary_relation(
            const gmme_id& boundary,
            const gmme_id& incident_entity,
            bool side = false ) override;

        void set_mesh_entity_boundary(
            const gmme_id& gmme,
            index_t id,
            index_t boundary_id,
            bool side = false ) override;

    private:
        GeoModelBuilderTopology(
            GeoModelBuilder2D& builder,
            GeoModel2D& geomodel )
            : GeoModelBuilderTopologyBase< 2 >( builder, geomodel )
        {
        }
    };

    template< >
    class RINGMESH_API GeoModelBuilderTopology< 3 > : public GeoModelBuilderTopologyBase< 3 > {
        friend class GeoModelBuilderBase< 3 > ;
        friend class GeoModelBuilder< 3 > ;
    public:
        virtual ~GeoModelBuilderTopology() = default;

        void add_universe_boundary( index_t boundary_id, bool side );

        void set_universe_boundary( index_t id, index_t boundary_id, bool side );

        void compute_universe();

        gmme_id create_mesh_entity(
            const MeshEntityType& type ) override;

        bool create_mesh_entities(
            const MeshEntityType& type,
            index_t nb_additional_entities ) override;

        void add_mesh_entity_boundary_relation(
            const gmme_id& boundary,
            const gmme_id& incident_entity,
            bool side = false ) override;

        void set_mesh_entity_boundary(
            const gmme_id& gmme,
            index_t id,
            index_t boundary_id,
            bool side = false ) override;

    private:
        GeoModelBuilderTopology(
            GeoModelBuilder3D& builder,
            GeoModel3D& geomodel )
            : GeoModelBuilderTopologyBase< 3 >( builder, geomodel )
        {
        }

        void copy_all_mesh_entity_topology( const GeoModel3D& from ) override;
    };

}
