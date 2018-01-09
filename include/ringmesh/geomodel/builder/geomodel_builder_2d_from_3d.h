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

#include <ringmesh/geomodel/builder/common.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    ALIAS_2D_AND_3D( GeoModel );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Base class for GeoModel2D building from GeoModel3D.
     */
    class geomodel_builder_api GeoModelBuilder2DFrom3D
        : public GeoModelBuilder< 2 >
    {
    public:
        GeoModelBuilder2DFrom3D( GeoModel2D& geomodel2d,
            const GeoModel3D& geomodel3d_from,
            const Geometry::Plane& plane );

    protected:
        vec2 get_2d_coord( const vec3& coord3d );

    protected:
        const GeoModel3D& geomodel3d_from_;
        const Geometry::Plane& plane_;
        vec3 u_axis_{};
        vec3 v_axis_{};
    };

    /*!
     * @brief Builder of GeoModel2D which project a GeoModel3D onto a plane.
     * @note This builder is dedicated to planar or sub-planar
     * GeoModel3D without volume, i.e. cross-sections, map-view models
     * or GeoModel3D only made of 1 fault or 1 horizon.
     * @warning The result GeoModel2D is not guaranteed to be valid.
     * It depends of the projection.
     */
    class geomodel_builder_api GeoModelBuilder2DProjection
        : public GeoModelBuilder2DFrom3D
    {
    public:
        GeoModelBuilder2DProjection( GeoModel2D& geomodel2d,
            const GeoModel3D& geomodel3d_from,
            const Geometry::Plane& plane );

        void build_geomodel();

    private:
        void copy_geomodel_3d_topology();

        void copy_geomodel_3d_geological_informations();

        void project_geomodel_3d_mesh_entities();

        std::vector< vec2 > compute_projected_vertices(
            const GeoModelMeshEntity3D& entity );
    };
} // namespace RINGMesh
