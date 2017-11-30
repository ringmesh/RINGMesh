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

/*!
 * @file Declaration of GeoModelEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella
 */

#pragma once

#include <ringmesh/geomodel/core/common.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopologyBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderTopology );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderRemoveBase );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Abstract base class describing one entity of a GeoModel
     */
    template < index_t DIMENSION >
    class geomodel_core_api GeoModelEntity
    {
        ringmesh_disable_copy_and_move( GeoModelEntity );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelEntity() = default;

        virtual bool is_on_voi() const = 0;
        virtual bool is_valid() const = 0;

        const GeoModel< DIMENSION >& geomodel() const
        {
            return geomodel_;
        }
        const std::string& name() const
        {
            return name_;
        }
        index_t index() const
        {
            return id_;
        }

    protected:
        /*!
         * @details Client code should only create GeoModelEntities through
         * GeoModelBuilderTopology class.
         *
         * @param[in] geomodel Geomodel owning the Entity to create
         * @param[in] id Index of the entity in the corresponding vector in the
         * geomodel
         * @param[in] name Name of the entity
         * @param[in] geological_feature Geological feature of the entity, none
         * by default.
         */
        GeoModelEntity( const GeoModel< DIMENSION >& geomodel, index_t id )
            : geomodel_( geomodel ), id_( id )
        {
        }

        void copy_name( const GeoModelEntity< DIMENSION >& from )
        {
            name_ = from.name_;
        }
        virtual bool is_index_valid() const = 0;

    protected:
        /// Reference to the GeoModel owning this entity
        const GeoModel< DIMENSION >& geomodel_;
        /// Name of the entity - default is "Unnamed"
        std::string name_ = std::string{ "Unnamed" };

        /// Index of the entity
        index_t id_{ NO_ID };
    };

    ALIAS_2D_AND_3D( GeoModelEntity );
} // namespace RINGMesh
