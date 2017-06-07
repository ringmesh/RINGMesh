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
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/geomodel/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/geomodel_builder_remove.h>
#include <ringmesh/geomodel/geomodel_builder_repair.h>
#include <ringmesh/geomodel/geomodel_builder_topology.h>

/*!
 * @brief Classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModelBuilder;
}

namespace RINGMesh {
    /*!
     * @brief First draft of flags to build a GeoModel
     * @todo Implements functions to set, access the values, depending on what ?
     * To check the consistency of the options. What do we do about the other entities ? [JP]
     *
     * @todo We need to keep track of the status of the GeoModel when building it:
     * same flags or some others ?
     *
     * @todo To separate in two classes ? One providing the low level functions set, assign etc,
     * and the other one some high level functions. [JP]
     */
    class GeoModelBuildingFlags {
    public:
        GeoModelBuildingFlags()
        {
            compute_corners = false;
            compute_lines = false;
            compute_surfaces = false;
            compute_regions_brep = false;
            compute_regions_mesh = false;
        }
        bool compute_corners;
        bool compute_lines;
        bool compute_surfaces;
        bool compute_regions_brep;
        bool compute_regions_mesh;
    };

    // Implementation details
    class GeoModelRegionFromSurfaces;
}

namespace RINGMesh {

    template< index_t DIMENSION >
    class GeoModelBuilderFromSurfaces {
    ringmesh_disable_copy( GeoModelBuilderFromSurfaces );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilder< DIMENSION >;

    public:
        /*!
         * Delete all GeoModelRegionFromSurfaces owned by the builder
         */
        virtual ~GeoModelBuilderFromSurfaces();
        /*
         * @brief From a GeoModel in which only Surfaces are defined,
         * create Corners, Lines and Regions depending on the building flags
         * @note Validity is not checked
         */
        void build();

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces();

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces
         * must have been called before
         */
        bool build_brep_regions_from_surfaces();

    private:
        GeoModelBuilderFromSurfaces(
            GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

    public:
        /*! Options to toggle the building of entities from the available entities */
        GeoModelBuildingFlags options_;

    private:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_;
    };
}
