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

#include <ringmesh/geomodel/core/common.h>

#include <vector>

#include <ringmesh/basic/frame.h>

#include <ringmesh/geomodel/core/entity_type_manager.h>
#include <ringmesh/geomodel/core/geomodel_mesh.h>
#include <ringmesh/geomodel/core/geomodel_ranges.h>

/*!
 * @file ringmesh/geomodel.h
 * @brief Class representing a geological structural model: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( WellGroup );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( Corner );
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    FORWARD_DECLARATION_DIMENSION_CLASS( Line );
    FORWARD_DECLARATION_DIMENSION_CLASS( Region );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelAccess );
    FORWARD_DECLARATION_DIMENSION_STRUCT( EntityTypeManager );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometryBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeometry );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoveBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemove );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRepair );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderInfo );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGM );

    ALIAS_2D_AND_3D( GeoModelMeshEntity );
    ALIAS_2D_AND_3D( Region );
} // namespace RINGMesh

namespace RINGMesh
{
    struct LineSide
    {
        LineSide() = default;
        std::vector< index_t > lines_;
        std::vector< bool > sides_;
    };
    struct SurfaceSide
    {
        SurfaceSide() = default;
        std::vector< index_t > surfaces_;
        std::vector< bool > sides_;
    };

    template < index_t DIMENSION >
    class geomodel_core_api GeoModelBase
    {
        ringmesh_disable_copy_and_move( GeoModelBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelAccess< DIMENSION >;

    public:
        virtual ~GeoModelBase();

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
            return entity_type_manager_.geological_entity_manager
                .nb_geological_entity_types();
        }

        const GeologicalEntityType& geological_entity_type(
            index_t index ) const
        {
            return entity_type_manager_.geological_entity_manager
                .geological_entity_type( index );
        }

        /*!
         * @brief Returns a const reference the identified
         * GeoModelGeologicalEntity
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
            const MeshEntityType& entity_type, index_t entity_index ) const
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
        mutable GeoModelMesh< DIMENSION > mesh;

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
            return entity_type_manager_.geological_entity_manager
                .geological_entity_type_index( type );
        }
        /*!
         * @brief Generic accessor to the storage of mesh entities of the given
         * type
         */
        virtual const std::vector<
            std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >&
            mesh_entities( const MeshEntityType& type ) const;

        /*!
         * @brief Generic accessor to the storage of geological entities of the
         * given type
         */
        const std::vector<
            std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >&
            geological_entities( const GeologicalEntityType& type ) const
        {
            index_t entity_index = geological_entity_type_index( type );
            return geological_entities( entity_index );
        }

        const std::vector<
            std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >&
            geological_entities( index_t geological_entity_type_index ) const
        {
            ringmesh_assert( geological_entity_type_index != NO_ID );
            return geological_entities_[geological_entity_type_index];
        }

    protected:
        std::string geomodel_name_;
        mutable double epsilon_{ -1 };

        EntityTypeManager< DIMENSION > entity_type_manager_;

        /*!
         * \name Mandatory entities of the geomodel
         * @{
         */
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >
            corners_;
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >
            lines_;
        std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >
            surfaces_;

        /*!
         * @brief Geological entities. They are optional.
         * The EntityTypes are managed by the EntityTypeManager of the class.
         */
        std::vector< std::vector<
            std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > > >
            geological_entities_;

        /*!
         * @}
         */

        /*! Optional WellGroup associated with the geomodel
         * @todo Move it out. It has nothing to do here. [JP]
         */
        const WellGroup< DIMENSION >* wells_{ nullptr };
    };
    ALIAS_2D_AND_3D( GeoModelBase );

    template < index_t DIMENSION >
    class geomodel_core_api GeoModel final : public GeoModelBase< DIMENSION >
    {
        friend class GeoModelAccess< DIMENSION >;

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
        geol_entity_range< DIMENSION > geol_entities(
            const GeologicalEntityType& geol_type ) const
        {
            return geol_entity_range< DIMENSION >( *this, geol_type );
        }
    };

    template <>
    class geomodel_core_api GeoModel< 3 > final : public GeoModelBase< 3 >
    {
        friend class GeoModelAccess< 3 >;

    public:
        GeoModel();
        ~GeoModel() override;

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
        geol_entity_range< 3 > geol_entities(
            const GeologicalEntityType& geol_type ) const
        {
            return geol_entity_range< 3 >( *this, geol_type );
        }

        index_t nb_regions() const
        {
            return static_cast< index_t >( regions_.size() );
        }

        const Region3D& region( index_t index ) const;

        const GeoModelMeshEntity3D& mesh_entity(
            const MeshEntityType& entity_type, index_t entity_index ) const
        {
            return GeoModelBase3D::mesh_entity( entity_type, entity_index );
        }

        const GeoModelMeshEntity3D& mesh_entity(
            const gmme_id& id ) const override;

        index_t nb_mesh_entities( const MeshEntityType& type ) const override;

        double epsilon3() const
        {
            return epsilon2() * epsilon();
        }
        SurfaceSide voi_surfaces() const;

    private:
        const std::vector< std::unique_ptr< GeoModelMeshEntity3D > >&
            mesh_entities( const MeshEntityType& type ) const override;

    private:
        std::vector< std::unique_ptr< GeoModelMeshEntity3D > > regions_;
    };

    template <>
    class GeoModel< 2 > final : public GeoModelBase< 2 >
    {
        friend class GeoModelAccess< 2 >;

    public:
        GeoModel();

        explicit GeoModel( PlaneReferenceFrame3D plane_reference_frame );

        ~GeoModel() override;

        corner_range< 2 > corners() const
        {
            return corner_range< 2 >( *this );
        }
        line_range< 2 > lines() const
        {
            return line_range< 2 >( *this );
        }
        surface_range< 2 > surfaces() const
        {
            return surface_range< 2 >( *this );
        }
        geol_entity_range< 2 > geol_entities(
            const GeologicalEntityType& geol_type ) const
        {
            return geol_entity_range< 2 >( *this, geol_type );
        }
        LineSide voi_lines() const;

    private:
        PlaneReferenceFrame3D reference_frame_{};
    };

    ALIAS_2D_AND_3D( GeoModel );
} // namespace RINGMesh
