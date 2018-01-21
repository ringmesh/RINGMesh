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
#include <ringmesh/geomodel/core/entity_type.h>

/*!
 * @file ringmesh/geomodel_api.h
 * @brief High level functions on GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 * @todo Encapsulate these functions in a namespace and TEST them.
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );

    ALIAS_3D( GeoModel );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Print in the console the geomodel statistics
     * @details Output number of polygons, vertices, and of the different entity
     * types
     * @todo Implement a test are_geomodels_equals to be able to check that
     * tests went well
     */
    template < index_t DIMENSION >
    void geomodel_core_api print_geomodel( const GeoModel< DIMENSION >& geomodel );

    /*!
     * Output the number of vertices, edges, polygons and cells.
     * Also output the number of triangles, quads and polygons if any.
     * Also output the number of tetra, prisms, pyramids, hex and polyhedra if
     * any.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    template < index_t DIMENSION >
    void geomodel_core_api print_geomodel_mesh_stats( const GeoModel< DIMENSION >& geomodel );

    /*!
     * Output the volume of the geomodel and the volume per cell type.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    void geomodel_core_api print_geomodel_mesh_cell_volumes(
        const GeoModel3D& geomodel );

    /*!
     * @return the index of the mesh entity \param gme_type named as \param name
     * in the GeoModel \param geomodel.
     * @note throw exception if no entities have this \param name or if two
     * entities
     * have the same \param name
     */
    template < index_t DIMENSION >
    index_t geomodel_core_api find_mesh_entity_id_from_name(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& gmme_type,
        const std::string& name );

    /*!
     * @return the index of the geological entity \param gme_type named as
     * \param name
     * in the GeoModel \param geomodel.
     * @note throw exception if no entities have this \param name or if two
     * entities
     * have the same \param name
     */
    template < index_t DIMENSION >
    index_t geomodel_core_api find_geological_entity_id_from_name(
        const GeoModel< DIMENSION >& geomodel,
        const GeologicalEntityType& gmge_type,
        const std::string& name );
} // namespace RINGMesh
