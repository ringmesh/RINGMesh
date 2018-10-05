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

/*!
 * @file Declaration of GeoModelGeologicalEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella
 */

#pragma once

#include <ringmesh/geomodel/core/common.h>

#include <memory>

#include <ringmesh/basic/factory.h>

#include <ringmesh/geomodel/core/geomodel_entity.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderGeology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoveBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemove );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderInfo );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntityAccess );

    class GeologicalEntityType;
    class MeshEntityType;
    struct gmge_id;
    struct gmme_id;
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class geomodel_core_api GeoModelGeologicalEntity
        : public GeoModelEntity< DIMENSION >
    {
    public:
        friend class GeoModelGeologicalEntityAccess< DIMENSION >;

        /*!
         * @brief Geological feature types for GeoModelEntity
         */
        enum struct GEOL_FEATURE
        {
            /// All geological features
            ALL_GEOL,
            /// Default value - No geological feature defined
            NO_GEOL,
            /// Stratigraphic surface - an horizon
            STRATI,
            /// Unconformity
            UNCONFORMITY,
            /// A normal fault
            NORMAL_FAULT,
            /// A reverse fault
            REVERSE_FAULT,
            /// An unspecified fault
            FAULT,
            /// Volume Of Interest
            VOI,
            // stratigraphic unit
            STRATI_UNIT,
            // geobody
            GEOBODY
        };

        /*!
         * @brief Map the name of a geological type with a value of GEOL_FEATURE
         *
         * @param[in] in Name of the feature. Can be
         * \li "reverse_fault"
         * \li "normal_fault"
         * \li "fault"
         * \li "top"
         * \li "none"
         * \li "topographic"
         * \li "unconformity"
         * \li "boundary"
         * Other strings will end up in \p NO_GEOL
         * @return The geological feature index
         * @todo Add other types of unconformity, see
         * RINGMesh::GeoModelEntity::TYPE. --GC
         */
        static GEOL_FEATURE determine_geological_type( const std::string& in );
        /*!
         * \return the (lowercase) string associated to a
         * GeoModelELement::GEOL_FEATURE
         */
        static std::string geol_name( GEOL_FEATURE feature );
        static bool is_fault( GEOL_FEATURE feature );
        static bool is_stratigraphic_limit( GEOL_FEATURE feature )
        {
            return feature == GEOL_FEATURE::STRATI
                   || feature == GEOL_FEATURE::UNCONFORMITY;
        }

        bool has_geological_feature() const
        {
            return geological_feature() != GEOL_FEATURE::NO_GEOL;
        }

        GEOL_FEATURE geological_feature() const
        {
            return geol_feature_;
        }

        static void initialize();

        gmge_id gmge() const;

        GeologicalEntityType entity_type() const;

        virtual MeshEntityType child_type_name() const = 0;
        virtual bool is_on_voi() const;
        virtual bool is_connectivity_valid() const;
        virtual bool is_valid() const;
        static GeologicalEntityType type_name_static();
        virtual GeologicalEntityType type_name() const;
        index_t nb_children() const
        {
            return static_cast< index_t >( children_.size() );
        }
        const gmme_id& child_gmme( index_t x ) const;
        const GeoModelMeshEntity< DIMENSION >& child( index_t x ) const;

        virtual bool is_identification_valid() const;

    protected:
        explicit GeoModelGeologicalEntity(
            const GeoModel< DIMENSION >& geomodel )
            : GeoModelEntity< DIMENSION >( geomodel, NO_ID )
        {
        }

        virtual bool is_index_valid() const;

    private:
        void copy_geological_entity(
            const GeoModelGeologicalEntity< DIMENSION >& from );

    protected:
        /// Children relations of this entity
        std::vector< index_t > children_{};

        /// Geological feature of this object - default is NO_GEOL
        GEOL_FEATURE geol_feature_{ GEOL_FEATURE::NO_GEOL };
    };

    ALIAS_2D_AND_3D( GeoModelGeologicalEntity );

    template < index_t DIMENSION >
    using GeoModelGeologicalEntityFactory = Factory< GeologicalEntityType,
        GeoModelGeologicalEntity< DIMENSION >,
        const GeoModel< DIMENSION >& >;

    ALIAS_2D_AND_3D( GeoModelGeologicalEntityFactory );

    template < index_t DIMENSION >
    class geomodel_core_api Contact
        : public GeoModelGeologicalEntity< DIMENSION >
    {
    public:
        explicit Contact( const GeoModel< DIMENSION >& geomodel )
            : GeoModelGeologicalEntity< DIMENSION >( geomodel )
        {
        }

        static GeologicalEntityType type_name_static();
        GeologicalEntityType type_name() const override;
        static MeshEntityType child_type_name_static();
        MeshEntityType child_type_name() const override;
    };

    ALIAS_2D_AND_3D( Contact );

    template < index_t DIMENSION >
    class geomodel_core_api Interface
        : public GeoModelGeologicalEntity< DIMENSION >
    {
    public:
        explicit Interface( const GeoModel< DIMENSION >& geomodel )
            : GeoModelGeologicalEntity< DIMENSION >( geomodel )
        {
        }

        static GeologicalEntityType type_name_static();
        GeologicalEntityType type_name() const override;
        static MeshEntityType child_type_name_static();
        MeshEntityType child_type_name() const override;
    };

    ALIAS_2D_AND_3D( Interface );

    template < index_t DIMENSION >
    class geomodel_core_api Layer : public GeoModelGeologicalEntity< DIMENSION >
    {
    public:
        explicit Layer( const GeoModel< DIMENSION >& geomodel )
            : GeoModelGeologicalEntity< DIMENSION >( geomodel )
        {
        }

        static GeologicalEntityType type_name_static();
        GeologicalEntityType type_name() const override;
        static MeshEntityType child_type_name_static();
        MeshEntityType child_type_name() const override;
    };

    ALIAS_2D_AND_3D( Layer );

} // namespace RINGMesh
