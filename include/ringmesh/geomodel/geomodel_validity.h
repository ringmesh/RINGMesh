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

#include <geogram/basic/file_system.h>

#include <ringmesh/geomodel/geomodel_indexing_types.h>

/*!
 * @file ringmesh/geomodel_validity.h
 * @brief Functions to check the validity of GeoModels
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModel;
    template< index_t DIMENSION > class GeoModelEntity;
}

namespace RINGMesh {
    /// Option to select what are checked.
    enum struct ValidityCheckMode {
        TOPOLOGY,
        GEOMETRY,
        GEOMETRY_EXCEPT_FACET_INTERSECTION,
        ALL,
        ALL_EXCEPT_FACET_INTERSECTION,
        UNDEFINED
    };

    /*! 
     * @brief Set the global default directory to store invalid entities of 
     *  geomodels to be the current working directory
     */
    static std::string validity_errors_directory =
        GEO::FileSystem::get_current_working_directory();

    /*!
     * @brief Set the directory where debugging information on 
     * invalid entities shall be stored
     * @details If directory does not exist keep the previous value.
     */
    void RINGMESH_API set_validity_errors_directory( const std::string& directory );

    /*!
     * @brief Check global geomodel validity
     * @param[in] geomodel GeoModel to check
     * @param[in] validity_check_mode Mode to select what model feature should
     * be checked. Set by default to the most complete check option.
     */
    template< index_t DIMENSION >
    bool is_geomodel_valid(
        const GeoModel< DIMENSION >& geomodel,
        ValidityCheckMode validity_check_mode = ValidityCheckMode::ALL );

    /*!
     * @brief Check the validity of all individual entity meshes
     * @details Check that the entities belong to this geomodel,
     *          call the check validity for each entity
     */
    bool RINGMESH_API are_geomodel_mesh_entities_mesh_valid(
        const GeoModel< 3 >& geomodel );

    /*!
     * @brief Check the connectivity of mesh entities
     */
    bool RINGMESH_API are_geomodel_mesh_entities_connectivity_valid(
        const GeoModel< 3 >& geomodel );

    bool RINGMESH_API are_geomodel_mesh_entities_parent_valid(
        const GeoModel< 3 >& geomodel );

    bool RINGMESH_API are_geomodel_geological_entities_valid(
        const GeoModel< 3 >& geomodel );
}
