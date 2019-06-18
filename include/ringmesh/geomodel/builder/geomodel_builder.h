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

#include <ringmesh/basic/pimpl.h>
#include <ringmesh/geomodel/builder/common.h>

#include <ringmesh/geomodel/builder/geomodel_builder_geology.h>
#include <ringmesh/geomodel/builder/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder_remove.h>
#include <ringmesh/geomodel/builder/geomodel_builder_topology.h>

/*!
 * @brief Classes to build GeoModel from various inputs
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
    class geomodel_builder_api GeoModelBuilderInfo
    {
        ringmesh_disable_copy_and_move( GeoModelBuilderInfo );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;

    public:
        /*!
         *@brief Set the name of the geomodel
         */
        void set_geomodel_name( const std::string& name )
        {
            geomodel_access_.modifiable_name() = name;
        }

        /*!
         *@brief Set the name of a geomodel mesh entity
         */
        void set_mesh_entity_name(
            const gmme_id& gmme_id, const std::string& name )
        {
            GeoModelMeshEntityAccess< DIMENSION > gmme_access(
                geomodel_access_.modifiable_mesh_entity( gmme_id ) );
            gmme_access.modifiable_name() = name;
        }

        /*!
         *@brief Set the name of a geomodel geological entity
         */
        void set_geological_entity_name(
            const gmge_id& gmge_id, const std::string& name )
        {
            GeoModelGeologicalEntityAccess< DIMENSION > gmge_access(
                geomodel_access_.modifiable_geological_entity( gmge_id ) );
            gmge_access.modifiable_name() = name;
        }

    protected:
        GeoModelBuilderInfo( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );
        ~GeoModelBuilderInfo() = default;

    private:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };

    /*!
     * @brief Base class to build or edit a GeoModel
     * @details All needed functions are organized in several specific builder
     * in accordance with the kind of edition operation or
     * with the GeoModel part which is edited (topology, geometry, geology,
     * info)
     */
    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderBase
    {
        ringmesh_disable_copy_and_move( GeoModelBuilderBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelBuilderBase() = default;

        /*!
         * @brief Finish up geomodel building and complete missing information.
         */
        virtual void end_geomodel() = 0;

        void build_corners_from_lines();

        void build_lines_and_corners_from_surfaces();

    protected:
        GeoModelBuilderBase( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

        void cut_geomodel_on_internal_boundaries();

    public:
        GeoModelBuilderTopology< DIMENSION > topology;
        GeoModelBuilderGeometry< DIMENSION > geometry;
        GeoModelBuilderGeology< DIMENSION > geology;
        GeoModelBuilderRemove< DIMENSION > remove;
        GeoModelBuilderInfo< DIMENSION > info;

    protected:
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilder
        : public GeoModelBuilderBase< DIMENSION >
    {
    };

    template <>
    class geomodel_builder_api GeoModelBuilder< 2 >
        : public GeoModelBuilderBase< 2 >
    {
    public:
        explicit GeoModelBuilder( GeoModel2D& geomodel );

        ~GeoModelBuilder();

        void end_geomodel() final;

        void build_surfaces_from_corners_and_lines();

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };

    template <>
    class geomodel_builder_api GeoModelBuilder< 3 >
        : public GeoModelBuilderBase< 3 >
    {
    public:
        explicit GeoModelBuilder( GeoModel3D& geomodel );

        void end_geomodel() final;

        void build_regions_from_lines_and_surfaces();
    };

    ALIAS_2D_AND_3D( GeoModelBuilder );
} // namespace RINGMesh
