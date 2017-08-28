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

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/entity_type.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/geomodel/geomodel_mesh.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @file ringmesh/geomodel.h
 * @brief Class representing a geological structural model: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {
    FORWARD_DECLARATION_DIMENSION_CLASS( WellGroup );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( Corner );
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    FORWARD_DECLARATION_DIMENSION_CLASS( Line );
    FORWARD_DECLARATION_DIMENSION_CLASS( Region );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelAccess );
    template< index_t DIMENSION > struct EntityTypeManager;
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometryBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometry );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemovalBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoval );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRepair );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderCopy );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderInfo );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGM );
} // namespace RINGMesh

namespace RINGMesh {
    /*!
     * @brief The class to describe a geological structural model represented
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    template< index_t DIMENSION >
    class GeoModelBase {
    ringmesh_disable_copy_and_move( GeoModelBase );
        ringmesh_template_assert_2d_or_3d (DIMENSION);
        friend class GeoModelAccess< DIMENSION > ;

    public:
        virtual ~GeoModelBase() = default;

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
        const EntityTypeManager< DIMENSION >& entity_type_manager() const
        {
            return entity_type_manager_;
        }

        /*!
         * @brief Returns the number of mesh entities of the given type
         * @details Default value is 0
         * @param[in] type the mesh entity type
         */
        virtual index_t nb_mesh_entities( const MeshEntityType& type ) const;

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
         * @param[in] id Type and index of the entity.
         * @pre Entity identification is valid.
         */
        const GeoModelGeologicalEntity< DIMENSION >& geological_entity(
            gmge_id id ) const
        {
            return *geological_entities( id.type() )[id.index()];
        }

        /*!
         * Convenient overload of entity( gmge_id id )
         */
        const GeoModelGeologicalEntity< DIMENSION >& geological_entity(
            const GeologicalEntityType& entity_type,
            index_t entity_index ) const
        {
            return geological_entity( gmge_id( entity_type, entity_index ) );
        }

        /*!
         * @brief Generic access to a meshed entity
         * @pre Type of the entity is CORNER, LINE, SURFACE, or REGION
         */
        virtual const GeoModelMeshEntity< DIMENSION >& mesh_entity(
            const gmme_id& id ) const;

        /*!
         * Convenient overload of mesh_entity( gmme_id id )
         */
        const GeoModelMeshEntity< DIMENSION >& mesh_entity(
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

        const Corner< DIMENSION >& corner( index_t index ) const;
        const Line< DIMENSION >& line( index_t index ) const;
        const Surface< DIMENSION >& surface( index_t index ) const;
        const Universe< DIMENSION >& universe() const
        {
            return universe_;
        }

        double epsilon() const;

        double epsilon2() const
        {
            return epsilon() * epsilon();
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
        void set_wells( const WellGroup< DIMENSION >* wells );
        const WellGroup< DIMENSION >* wells() const
        {
            return wells_;
        }

    public:
        GeoModelMesh< DIMENSION > mesh;

    protected:
        /*!
         * @brief Constructs an empty GeoModel
         */
        explicit GeoModelBase( GeoModel< DIMENSION >& geomodel );
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
        virtual const std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >& mesh_entities(
            const MeshEntityType& type ) const;

        /*!
         * @brief Generic accessor to the storage of geological entities of the given type
         */
        const std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& geological_entities(
            const GeologicalEntityType& type ) const
        {
            index_t entity_index = geological_entity_type_index( type );
            return geological_entities( entity_index );
        }

        const std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& geological_entities(
            index_t geological_entity_type_index ) const
        {
            ringmesh_assert( geological_entity_type_index != NO_ID );
            return geological_entities_[geological_entity_type_index];
        }

    protected:
        std::string geomodel_name_;
        mutable double epsilon_ { -1 };

        EntityTypeManager< DIMENSION > entity_type_manager_;

        /*!
         * \name Mandatory entities of the geomodel
         * @{
         */
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > > corners_;
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > > lines_;
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > > surfaces_;

        /*!
         * The Universe defines the extension of the GeoModel
         */
        Universe< DIMENSION > universe_;

        /*!
         * @brief Geological entities. They are optional.
         * The EntityTypes are managed by the EntityTypeManager of the class.
         */
        std::vector<
            std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > > > geological_entities_;

        /*!
         * @}
         */

        /*! Optional WellGroup associated with the geomodel
         * @todo Move it out. It has nothing to do here. [JP]
         */
        const WellGroup< DIMENSION >* wells_ { nullptr };
    };
    ALIAS_2D_AND_3D( GeoModelBase );

    template< index_t DIMENSION >
    class GeoModel final: public GeoModelBase< DIMENSION > {
        friend class GeoModelAccess< DIMENSION > ;
    public:
        GeoModel();

        corner_range< DIMENSION > corners() const
        {
            return corner_range< DIMENSION >( *this );
        }
        line_range< DIMENSION > lines() const
        {
            return line_range< DIMENSION >( *this );
        }
        surface_range< DIMENSION > surfaces() const
        {
            return surface_range< DIMENSION >( *this );
        }
    };

    template< >
    class RINGMESH_API GeoModel< 3 > final: public GeoModelBase< 3 > {
        friend class GeoModelAccess< 3 > ;
    public:
        GeoModel();

        corner_range< 3 > corners() const
        {
            return corner_range< 3 >( *this );
        }
        line_range< 3 > lines() const
        {
            return line_range< 3 >( *this );
        }
        surface_range< 3 > surfaces() const
        {
            return surface_range< 3 >( *this );
        }
        region_range< 3 > regions() const
        {
            return region_range< 3 >( *this );
        }

        index_t nb_regions() const
        {
            return static_cast< index_t >( regions_.size() );
        }

        const Region3D& region( index_t index ) const;

        const GeoModelMeshEntity3D& mesh_entity(
            const MeshEntityType& entity_type,
            index_t entity_index ) const
        {
            return GeoModelBase3D::mesh_entity( entity_type, entity_index );
        }

        const GeoModelMeshEntity3D& mesh_entity( const gmme_id& id ) const override;

        index_t nb_mesh_entities( const MeshEntityType& type ) const override;

        double epsilon3() const
        {
            return epsilon2() * epsilon();
        }
    private:
        const std::vector< std::unique_ptr< GeoModelMeshEntity3D > >& mesh_entities(
            const MeshEntityType& type ) const override;

    private:
        std::vector< std::unique_ptr< GeoModelMeshEntity3D > > regions_;
    };
    ALIAS_2D_AND_3D( GeoModel );

    template< index_t DIMENSION >
    class GeoModelAccess {
    ringmesh_disable_copy_and_move( GeoModelAccess );
        ringmesh_template_assert_2d_or_3d (DIMENSION);
        friend class GeoModelBuilderBase< DIMENSION > ;
        friend class GeoModelBuilder< DIMENSION > ;
        friend class GeoModelBuilderGM< DIMENSION > ;
        friend class GeoModelBuilderTopologyBase< DIMENSION > ;
        friend class GeoModelBuilderTopology< DIMENSION > ;
        friend class GeoModelBuilderGeometryBase< DIMENSION > ;
        friend class GeoModelBuilderGeometry< DIMENSION > ;
        friend class GeoModelBuilderGeology< DIMENSION > ;
        friend class GeoModelBuilderRemovalBase< DIMENSION > ;
        friend class GeoModelBuilderRemoval< DIMENSION > ;
        friend class GeoModelBuilderRepair< DIMENSION > ;
        friend class GeoModelBuilderCopy< DIMENSION > ;
        friend class GeoModelBuilderInfo< DIMENSION > ;

    private:
        explicit GeoModelAccess( GeoModel< DIMENSION >& geomodel )
            : geomodel_( geomodel )
        {
        }
        ~GeoModelAccess() = default;

        std::string& modifiable_name()
        {
            return geomodel_.geomodel_name_;
        }

        EntityTypeManager< DIMENSION >& modifiable_entity_type_manager()
        {
            return geomodel_.entity_type_manager_;
        }

        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >& modifiable_mesh_entities(
            const MeshEntityType& type )
        {
            return const_cast< std::vector<
                std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >& >( geomodel_.mesh_entities(
                type ) );
        }

        GeoModelMeshEntity< DIMENSION >& modifiable_mesh_entity( const gmme_id& id )
        {
            return *modifiable_mesh_entities( id.type() )[id.index()];
        }

        std::vector<
            std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > > >& modifiable_geological_entities()
        {
            return geomodel_.geological_entities_;
        }

        std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& modifiable_geological_entities(
            const GeologicalEntityType& type )
        {
            return const_cast< std::vector<
                std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& >( geomodel_.geological_entities(
                type ) );
        }

        GeoModelGeologicalEntity< DIMENSION >& modifiable_geological_entity(
            const gmge_id& id )
        {
            return *modifiable_geological_entities( id.type() )[id.index()];
        }

        Universe< DIMENSION >& modifiable_universe()
        {
            return geomodel_.universe_;
        }

        double& modifiable_epsilon()
        {
            return geomodel_.epsilon_;
        }

    private:
        GeoModel< DIMENSION >& geomodel_;
    };
} // namespace RINGMesh
