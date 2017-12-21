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

#include <ringmesh/tetrahedralize/common.h>

#include <memory>

#include <geogram/mesh/mesh.h>

#include <ringmesh/basic/factory.h>

#include <ringmesh/geomodel/builder/geomodel_builder.h>

/*!
 * @file ringmesh/tetragen.h
 * @brief API class interfacing GeoModel with external tetrahedral meshers
 * @author Arnaud Botella
 */

#ifdef USE_MG_TETRA
extern "C" {
#include <meshgems/meshgems.h>
#include <meshgems/tetra.h>
}
#endif

namespace RINGMesh
{
    class TetraGen;
    FORWARD_DECLARATION_DIMENSION_CLASS( Region );
    FORWARD_DECLARATION_DIMENSION_CLASS( WellGroup );

    ALIAS_3D( Region );
    ALIAS_3D( WellGroup );
} // namespace RINGMesh

namespace RINGMesh
{
    class tetrahedralize_api TetraGen
    {
        ringmesh_disable_copy_and_move( TetraGen );

    public:
        virtual ~TetraGen() = default;
        static std::unique_ptr< TetraGen > create(
            GeoModel3D& M, index_t region_id, const std::string& algo_name );
        static void initialize();

        /*!
         * Sets the boundaries of the domain
         * @param[in] region The Region of the GeoModel to mesh
         * @param[in] wells the wells to be conformal to
         */
        void set_boundaries(
            const Region3D& region, const WellGroup3D* wells = nullptr );

        /*!
         * Set additional points to be in the output tetrahedral mesh
         * @param[in] points the points to add
         */
        void set_internal_points( const std::vector< vec3 >& points );

        /*!
         * @brief Send the set of points/edges/triangles to MGTetra or TetGen
         * @details A set of points/edges/triangles are given to MGtetra or
         * Tetgen
         * The two mesh generators are configurated. Then check and repair
         * functions
         * are launched in order to control the outputs
         * @param[in] refine tells whether or not there are refined options to
         * set (true by defaults)
         */
        bool tetrahedralize( bool refine = true );

    protected:
        TetraGen( GeoModel3D& geomodel, index_t region_id )
            : builder_( geomodel ), output_region_( region_id )
        {
        }

        virtual bool do_tetrahedralize( bool refine ) = 0;

    protected:
        GeoModelBuilder3D builder_;
        index_t output_region_{ NO_ID };
        GEO::Mesh tetmesh_constraint_;
        const Region3D* region_{ nullptr };
        const WellGroup3D* wells_{ nullptr };
    };

    using TetraGenFactory =
        Factory< std::string, TetraGen, GeoModel3D&, index_t >;
} // namespace RINGMesh
