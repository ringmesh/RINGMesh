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

#include <ringmesh/geomodel/tools/common.h>

/*!
 * @author Benjamin Chauvin
 * This code is inspired from
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMesh );

    ALIAS_3D( GeoModel );
    ALIAS_3D( VolumeMesh );
} // namespace RINGMesh

namespace RINGMesh
{
    enum MeshQualityMode
    {
        INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS,
        INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH,
        VOLUME_BY_SUM_SQUARE_EDGE,
        MIN_SOLID_ANGLE
    };

    /*!
     * @brief Computes and stores mesh quality in the GeoModel.
     *
     * Enables to have a metrics on the quality of the 3D mesh. This is
     * important for instance for numerical processes. For now all the qualities
     * are based on the regular tetrahedron: a good element here is an element
     * near equilaterality (this definition may change in function of the
     * requirements of the numerical process).
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     * For more information about the mesh quality, see
     * <a
     * href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] mesh_qual_mode mesh quality to compute.
     * @param[in,out] geomodel GeoModel in which the mesh quality is performed.
     * The quality is stored on the cells of each Region.
     *
     * @warning The GeoModel must have at least one region. All the regions
     * must be meshed by simplexes (tetrahedra).
     */
    void geomodel_tools_api compute_prop_tet_mesh_quality(
        MeshQualityMode mesh_qual_mode, const GeoModel3D& geomodel );

    /*!
     * @brief Fill the /p output_mesh with cells of quality below \p min_quality
     * @param[in] mesh_qual_mode mesh quality of cells.
     * @param[in] min_quality Value of quality below which a cell is add in the
     * mesh.
     * @param[in] geomodel GeoModel in which the mesh quality is read.
     * @param[out] output_mesh VolumeMesh to fill with low quality cells.
     * @returns The minimum value of cell quality
     *
     * @warning The GeoModel must have at least one region. All the regions
     * must be meshed by simplexes (tetrahedra).
     */
    double geomodel_tools_api fill_mesh_with_low_quality_cells(
        MeshQualityMode mesh_qual_mode,
        double min_quality,
        const GeoModel3D& geomodel,
        VolumeMesh3D& output_mesh );
} // namespace RINGMesh
