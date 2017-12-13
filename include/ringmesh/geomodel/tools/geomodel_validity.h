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

#include <ringmesh/geomodel/tools/common.h>

/*!
 * @file ringmesh/geomodel_validity.h
 * @brief Functions to check the validity of GeoModels
 * @author Jeanne Pellerin
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelEntity );
} // namespace RINGMesh

namespace RINGMesh
{
    /// Option to select what are checked.
    enum struct ValidityCheckMode
    {
        EMPTY = 0,
        FINITE_EXTENSION = 1,
        GEOMODEL_CONNECTIVITY = 1 << 1,
        GEOLOGICAL_ENTITIES = 1 << 2,
        SURFACE_LINE_MESH_CONFORMITY = 1 << 3,
        REGION_SURFACE_MESH_CONFORMITY = 1 << 4,
        MESH_ENTITIES = 1 << 5,
        NON_MANIFOLD_EDGES = 1 << 6,
        POLYGON_INTERSECTIONS = 1 << 7,
        GEOLOGY = GEOLOGICAL_ENTITIES,
        TOPOLOGY = FINITE_EXTENSION | GEOMODEL_CONNECTIVITY,
        GEOMETRY = SURFACE_LINE_MESH_CONFORMITY | REGION_SURFACE_MESH_CONFORMITY
                   | MESH_ENTITIES | NON_MANIFOLD_EDGES | POLYGON_INTERSECTIONS,
        ALL = TOPOLOGY | GEOMETRY | GEOLOGY
    };
    ENABLE_BITMASK_OPERATORS( ValidityCheckMode );

    /*!
     * @brief Set the directory where debugging information on
     * invalid entities shall be stored
     * @details If directory does not exist keep the previous value.
     */
    void geomodel_tools_api set_validity_errors_directory(
        const std::string& directory );

    /*!
     * @brief Get the directory where debugging information on
     * invalid entities shall be stored
     */
    std::string geomodel_tools_api get_validity_errors_directory();

    /*!
     * @brief Get validity mode from command line argument validity:do_not_check
     */
    ValidityCheckMode geomodel_tools_api get_validity_mode_from_arg();

    /*!
     * @brief Check global geomodel validity
     * @param[in] geomodel GeoModel to check
     * @param[in] validity_check_mode Mode to select what model feature should
     * be checked. Set by default to the most complete check option.
     */
    template < index_t DIMENSION >
    bool is_geomodel_valid( const GeoModel< DIMENSION >& geomodel,
        ValidityCheckMode validity_check_mode = get_validity_mode_from_arg() );

    /*!
     * @brief Check the validity of all individual entity meshes
     * @details Check that the entities belong to this geomodel,
     *          call the check validity for each entity
     */
    template < index_t DIMENSION >
    bool are_geomodel_mesh_entities_mesh_valid(
        const GeoModel< DIMENSION >& geomodel );

    /*!
     * @brief Check the connectivity of mesh entities
     */
    template < index_t DIMENSION >
    bool are_geomodel_mesh_entities_connectivity_valid(
        const GeoModel< DIMENSION >& geomodel );

    template < index_t DIMENSION >
    bool are_geomodel_mesh_entities_parent_valid(
        const GeoModel< DIMENSION >& geomodel );

    template < index_t DIMENSION >
    bool are_geomodel_geological_entities_valid(
        const GeoModel< DIMENSION >& geomodel );
} // namespace RINGMesh
