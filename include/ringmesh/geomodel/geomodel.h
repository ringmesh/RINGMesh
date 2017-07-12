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

#include <vector>

#include <geogram/basic/factory.h>

#include <ringmesh/basic/algorithm.h>

#include <ringmesh/geomodel/entity_type.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_mesh.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @file ringmesh/geomodel.h
 * @brief Class representing a geological structural geomodel: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {
    class WellGroup;
    class GeoModelGeologicalEntity;
    class GeoModelMeshEntity;
    class Corner;
    class Surface;
    class Line;
    class Region;
    class EntityTypeManager;
}

namespace RINGMesh {
    /*!
     * @brief The class to describe a geological geomodel represented 
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    class RINGMESH_API GeoModel {
    ringmesh_disable_copy( GeoModel );
        friend class GeoModelAccess;

    public:
        /*!
         * @brief Constructs an empty GeoModel
         */
        GeoModel();

        /*!
         * @brief Gets the name of the GeoModel
         */
        const std::string& name() const
        {
            return geomodel_name_;
        }

        /*!
         * @brief Gets the EntityTypeManager associated to the GeoModel
         */
        const EntityTypeManager& entity_type_manager() const
        {
            return entity_type_manager_;
        }

        /*!
         * @brief Returns the number of mesh entities of the given type
         * @details Default value is 0
         * @param[in] type the mesh entity type
         */
        index_t nb_mesh_entities( const MeshEntityType& type ) const;

        /*!
         * @brief Returns the number of geological entities of the given type
         * @details Default value is 0
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entities( const GeologicalEntityType& type ) const
        {
            return static_cast< index_t >( geological_entities( type ).size() );
        }
        /*!
         * @brief Returns the index of the geological entity type storage
         * @details Default value is NO_ID
         * @param[in] type the geological entity type
         */
        index_t nb_geological_entity_types() const
        {
            return entity_type_manager_.geological_entity_manager.nb_geological_entity_types();
        }

        const GeologicalEntityType& geological_entity_type( index_t index ) const
        {
            return entity_type_manager_.geological_entity_manager.geological_entity_type(
                index );
        }
        /*!
         * @brief Returns a const reference the identified GeoModelGeologicalEntity
         * @param[in] id Type and index of the entity. For the
         * pair (Region, NO_ID) universe region is returned.
         * @pre Entity identification is valid.
         */
        const GeoModelGeologicalEntity& geological_entity( gmge_id id ) const
        {
            return *geological_entities( id.type() )[id.index()];
        }
        /*!
         * Convenient overload of entity( gmge_id id )
         */
        const GeoModelGeologicalEntity& geological_entity(
            const GeologicalEntityType& entity_type,
            index_t entity_index ) const
        {
            return geological_entity( gmge_id( entity_type, entity_index ) );
        }
        /*!
         * @brief Generic access to a meshed entity
         * @pre Type of the entity is CORNER, LINE, SURFACE, or REGION
         */
        const GeoModelMeshEntity& mesh_entity( gmme_id id ) const;
        /*!
         * Convenient overload of mesh_entity( gmme_id id )
         */
        const GeoModelMeshEntity& mesh_entity(
            const MeshEntityType& entity_type,
            index_t entity_index ) const
        {
            return mesh_entity( gmme_id( entity_type, entity_index ) );
        }
        /*! @}
         * \name Specialized accessors.
         * @{
         */
        index_t nb_corners() const
        {
            return static_cast< index_t >( corners_.size() );
        }
        index_t nb_lines() const
        {
            return static_cast< index_t >( lines_.size() );
        }
        index_t nb_surfaces() const
        {
            return static_cast< index_t >( surfaces_.size() );
        }
        index_t nb_regions() const
        {
            return static_cast< index_t >( regions_.size() );
        }

        const Corner& corner( index_t index ) const;
        const Line& line( index_t index ) const;
        const Surface& surface( index_t index ) const;
        const Region& region( index_t index ) const;

        double epsilon() const;

        double epsilon2() const
        {
            return epsilon() * epsilon();
        }

        double epsilon3() const
        {
            return epsilon2() * epsilon();
        }

        /*!
         * @}
         */
        /*!
         * Associates a WellGroup to the GeoModel
         * @param[in] wells the WellGroup
         * @todo Review : What is this for ?
         * @todo Extend to other object types.
         */
        void set_wells( const WellGroup* wells );
        const WellGroup* wells() const
        {
            return wells_;
        }

        std::tuple< std::vector< index_t >, std::vector< bool > > get_voi_surfaces() const;

    public:
        GeoModelMesh mesh;

    private:
        /*!
         * Access to the position of the entity of that type in storage.
         */
        index_t geological_entity_type_index(
            const GeologicalEntityType& type ) const
        {
            return entity_type_manager_.geological_entity_manager.geological_entity_type_index(
                type );
        }
        /*!
         * @brief Generic accessor to the storage of mesh entities of the given type
         */
        const std::vector< std::unique_ptr< GeoModelMeshEntity > >& mesh_entities(
            const MeshEntityType& type ) const;

        /*!
         * @brief Generic accessor to the storage of geological entities of the given type
         */
        const std::vector< std::unique_ptr< GeoModelGeologicalEntity > >& geological_entities(
            const GeologicalEntityType& type ) const
        {
            index_t entity_index = geological_entity_type_index( type );
            return geological_entities( entity_index );
        }

        const std::vector< std::unique_ptr< GeoModelGeologicalEntity > >& geological_entities(
            index_t geological_entity_type_index ) const
        {
            ringmesh_assert( geological_entity_type_index != NO_ID );
            return geological_entities_[geological_entity_type_index];
        }

    private:
        std::string geomodel_name_;
        mutable double epsilon_;

        EntityTypeManager entity_type_manager_;

        /*!
         * \name Mandatory entities of the geomodel
         * @{
         */
        std::vector< std::unique_ptr< GeoModelMeshEntity > > corners_;
        std::vector< std::unique_ptr< GeoModelMeshEntity > > lines_;
        std::vector< std::unique_ptr< GeoModelMeshEntity > > surfaces_;
        std::vector< std::unique_ptr< GeoModelMeshEntity > > regions_;

        /*!
         * @brief Geological entities. They are optional.
         * The EntityTypes are managed by the EntityTypeManager of the class.
         */
        std::vector< std::vector< std::unique_ptr< GeoModelGeologicalEntity > > > geological_entities_;

        /*!
         * @}
         */

        /*! Optional WellGroup associated with the geomodel
         * @todo Move it out. It has nothing to do here. [JP]
         */
        const WellGroup* wells_;
    };

    class GeoModelAccess {
    ringmesh_disable_copy( GeoModelAccess );
        friend class GeoModelBuilder;
        friend class GeoModelBuilderGM;
        friend class GeoModelBuilderTopology;
        friend class GeoModelBuilderGeometry;
        friend class GeoModelBuilderGeology;
        friend class GeoModelBuilderRemoval;
        friend class GeoModelBuilderRepair;
        friend class GeoModelBuilderCopy;
        friend class GeoModelBuilderInfo;
        friend class GeoModelBuilderFromSurfaces;

    private:
        GeoModelAccess( GeoModel& geomodel )
            : geomodel_( geomodel )
        {
        }

        std::string& modifiable_name()
        {
            return geomodel_.geomodel_name_;
        }

        EntityTypeManager& modifiable_entity_type_manager()
        {
            return geomodel_.entity_type_manager_;
        }

        std::vector< std::unique_ptr< GeoModelMeshEntity > >& modifiable_mesh_entities(
            const MeshEntityType& type )
        {
            return const_cast< std::vector< std::unique_ptr< GeoModelMeshEntity > >& >( geomodel_.mesh_entities(
                type ) );
        }

        GeoModelMeshEntity& modifiable_mesh_entity( const gmme_id& id )
        {
            return *modifiable_mesh_entities( id.type() )[id.index()];
        }

        std::vector< std::vector< std::unique_ptr< GeoModelGeologicalEntity > > >& modifiable_geological_entities()
        {
            return geomodel_.geological_entities_;
        }

        std::vector< std::unique_ptr< GeoModelGeologicalEntity > >& modifiable_geological_entities(
            const GeologicalEntityType& type )
        {
            return const_cast< std::vector<
                std::unique_ptr< GeoModelGeologicalEntity > >& >( geomodel_.geological_entities(
                type ) );
        }

        GeoModelGeologicalEntity& modifiable_geological_entity( const gmge_id& id )
        {
            return *modifiable_geological_entities( id.type() )[id.index()];
        }

        double& modifiable_epsilon()
        {
            return geomodel_.epsilon_;
        }

    private:
        GeoModel& geomodel_;
    };

    using GeoModelGeologicalEntityFactory = GEO::Factory1< GeoModelGeologicalEntity, GeoModel >;

#define ringmesh_register_GeoModelGeologicalEntity_creator( type ) \
    geo_register_creator( GeoModelGeologicalEntityFactory, type, type::type_name_static() )

}
