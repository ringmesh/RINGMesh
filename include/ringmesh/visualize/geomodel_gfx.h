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

#include <ringmesh/visualize/common.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <memory>

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

#include <ringmesh/visualize/geomodel_entity_gfx.h>
#include <ringmesh/visualize/mesh_entity_gfx.h>

/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class visualize_api GeoModelGfxBase
    {
        ringmesh_disable_copy_and_move( GeoModelGfxBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelGfxBase() = default;

        /*!
         * Sets the GeoModel associated to the graphics
         * @param[in] geomodel the GeoModel
         */
        void set_geomodel( const GeoModel< DIMENSION >& geomodel );
        /*!
         * Gets the GeoModel associated to the graphics
         * @return the GeoModel
         */
        const GeoModel< DIMENSION >* geomodel() const;

    protected:
        /*!
         * Initializes the database according the GeoModel dimensions
         */
        virtual void initialize();

        explicit GeoModelGfxBase( GeoModelGfx< DIMENSION >& gfx );

    private:
        /// The GeoModel associated to the graphics
        const GeoModel< DIMENSION >* geomodel_{ nullptr };

    public:
        CornerGfxEntity< DIMENSION > corners;
        LineGfxEntity< DIMENSION > lines;
        SurfaceGfxEntity< DIMENSION > surfaces;
        AttributeGfxManager< DIMENSION > attribute;
    };

    ALIAS_2D_AND_3D( GeoModelGfxBase );

    template < index_t DIMENSION >
    class visualize_api GeoModelGfx final : public GeoModelGfxBase< DIMENSION >
    {
    public:
        GeoModelGfx();
    };

    template <>
    class visualize_api GeoModelGfx< 3 > final : public GeoModelGfxBase< 3 >
    {
    public:
        GeoModelGfx();

        void initialize() final;

    public:
        RegionGfxEntity3D regions;
    };

    ALIAS_2D_AND_3D( GeoModelGfx );
} // namespace RINGMesh

#endif
