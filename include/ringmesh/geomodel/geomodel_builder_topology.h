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
    class GeoModelBuilder ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderTopology {
    ringmesh_disable_copy( GeoModelBuilderTopology ) ;
        friend class GeoModelBuilder ;

    public:
        /*!
         * @brief Copy topological information from a geomodel
         * @details Copy all the geomodel entities and their relationship
         * ignoring their geometry
         * @param[in] from Model to copy the information from
         */
        void copy_topology( const GeoModel& from ) ;

        /*!
         * @brief Add to the vector the entities which cannot exist if
         *        an entity in the set does not exist.
         * @return True if at least one entity was added.
         */
        bool get_dependent_entities( std::set< gmme_t >& in_mesh_entities, std::set< gmge_t >& in_geological_entities ) const ;


        template< typename ENTITY >
        gmme_t create_mesh_entity( const MeshType mesh_type = "" )
        {
            const MeshEntityType entity_type = ENTITY::type_name_static() ;
            index_t nb_entities( geomodel_.nb_mesh_entities( entity_type ) ) ;
            index_t new_id( nb_entities ) ;
            ENTITY* new_entity = GeoModelMeshEntityAccess::create_entity< ENTITY >(
                geomodel_, new_id, mesh_type ) ;
            geomodel_access_.modifiable_mesh_entities( entity_type ).push_back(
                new_entity ) ;
            return new_entity->gmme_id() ;
        }

        bool create_mesh_entities(
            const MeshEntityType& type,
            index_t nb_additional_entities )
        {
            if( MeshEntityTypeManager::is_corner( type ) ) {
                return create_mesh_entities< Corner >( nb_additional_entities ) ;
            } else if( MeshEntityTypeManager::is_line( type ) ) {
                return create_mesh_entities< Line >( nb_additional_entities ) ;
            } else if( MeshEntityTypeManager::is_surface( type ) ) {
                return create_mesh_entities< Surface >( nb_additional_entities ) ;
            } else if( MeshEntityTypeManager::is_region( type ) ) {
                return create_mesh_entities< Region >( nb_additional_entities ) ;
            } else {
                ringmesh_assert_not_reached ;
                return false ;
            }
        }

        /*!
         * @brief Complete missing information in GeoModelEntities
         * boundaries - in_boundary - parent - children
         * @details For all 7 types of entities, check what information is available
         * for the first one and fill the entities of the same type accordingly
         * THIS MEANS that the all the entities of the same type have been initialized with
         * the same information
         */
        void complete_entity_connectivity() ;

        /*!
         * @brief Fill the boundaries of all entities of the given type
         * @details If the boundary entities do not have any in_boundary
         * information, nothing is done.
         */
        void fill_mesh_entities_boundaries( const MeshEntityType& type ) ;

        /*!
         * @brief Fill the in_boundary vector of all entities of the given type
         * @details If the in_boundary entities do not have any boundary
         * information, nothing is done, and geomodel construction will eventually fail.
         */
        void fill_mesh_entities_in_boundaries( const MeshEntityType& type ) ;

        void add_mesh_entity_boundary(
            const gmme_t& gme_id,
            index_t boundary_id,
            bool side = false ) ;

        void set_mesh_entity_boundary(
            const gmme_t& gme_id,
            index_t id,
            index_t boundary_id,
            bool side = false ) ;

        void add_universe_boundary( index_t boundary_id, bool side ) ;

        void set_universe_boundary( index_t id, index_t boundary_id, bool side ) ;

        void add_mesh_entity_in_boundary( const gmme_t& t, index_t in_boundary_id ) ;

        void set_mesh_entity_in_boundary(
            const gmme_t& gme_id,
            index_t id,
            index_t in_boundary_id ) ;

        void delete_mesh_entity( const MeshEntityType& type, index_t index ) ;

        /*!
         * @brief Finds or creates a corner at given coordinates.
         * @param[in] point Geometric location of the Corner
         * @return Index of the Corner
         */
        gmme_t find_or_create_corner( const vec3& point ) ;
        gmme_t find_or_create_corner( index_t geomodel_point_id ) ;
        gmme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;

        /*!
         * @brief Finds or creates a line knowing its topological adjacencies
         */
        gmme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gmme_t& first_corner,
            const gmme_t& second_corner ) ;

        void compute_universe() ;

    private:
        GeoModelBuilderTopology( GeoModelBuilder& builder, GeoModel& geomodel ) ;

        template< typename ENTITY >
        bool create_mesh_entities(
            index_t nb_additionnal_entities,
            const MeshType type = "" )
        {
            const MeshEntityType& entity_type = ENTITY::type_name_static() ;
            std::vector< GeoModelMeshEntity* >& store =
                geomodel_access_.modifiable_mesh_entities( entity_type ) ;
            index_t old_size = static_cast< index_t >( store.size() ) ;
            index_t new_size = old_size + nb_additionnal_entities ;
            store.resize( new_size, nil ) ;
            for( index_t i = old_size; i < new_size; i++ ) {
                store[i] = GeoModelMeshEntityAccess::create_entity< ENTITY >(
                    geomodel_, i, type ) ;
            }
            return true ;
        }

        void complete_mesh_entity_connectivity( const MeshEntityType& type ) ;

        template< typename ENTITY >
        void copy_mesh_entity_topology( const GeoModel& from )
        {
            const MeshEntityType& type = ENTITY::type_name_static() ;
            create_mesh_entities< ENTITY >( from.nb_mesh_entities( type ) ) ;

            RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < geomodel_.nb_mesh_entities( type ); ++e ) {
                gmme_t id( type, e ) ;
                GeoModelMeshEntityAccess gmme_access(
                    geomodel_access_.modifiable_mesh_entity( id ) ) ;
                gmme_access.copy( from.mesh_entity( id ) ) ;
            }
        }

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

}
