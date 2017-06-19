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

/*!
 * @file ringmesh/geomodel_api.h
 * @brief High level functions on GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 * @todo Encapsulate these functions in a namespace and TEST them.
 */

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModel;
    class MeshEntityType;
    class GeologicalEntityType;
}

namespace RINGMesh {

    /*!
     * @brief Print in the console the geomodel statistics
     * @details Output number of polygons, vertices, and of the different entity types
     * @todo Implement a test are_geomodels_equals to be able to check that tests went well
     */
    template< index_t DIMENSION >
    void print_geomodel( const GeoModel< DIMENSION >& geomodel );

    /*!
     * Output the number of vertices, edges, polygons and cells.
     * Also output the number of triangles, quads and polygons if any.
     * Also output the number of tetra, prisms, pyramids, hex and polyhedra if any.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    template< index_t DIMENSION >
    void print_geomodel_mesh_stats( const GeoModel< DIMENSION >& geomodel );

    /*!
     * Output the volume of the geomodel and the volume per cell type.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    void RINGMESH_API print_geomodel_mesh_cell_volumes( const GeoModel< 3 >& geomodel );

    template< index_t DIMENSION >
    bool are_geomodel_surface_meshes_simplicial(
        const GeoModel< DIMENSION >& geomodel );

    bool RINGMESH_API are_geomodel_region_meshes_simplicial(
        const GeoModel< 3 >& geomodel );

    /*!
     * @return the index of the mesh entity \param gme_type named as \param name
     * in the GeoModel \param geomodel.
     * @note throw exception if no entities have this \param name or if two entities
     * have the same \param name
     */
    template< index_t DIMENSION >
    index_t find_mesh_entity_id_from_name(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& gmme_type,
        const std::string& name );

    /*!
     * @return the index of the geological entity \param gme_type named as \param name
     * in the GeoModel \param geomodel.
     * @note throw exception if no entities have this \param name or if two entities
     * have the same \param name
     */
    template< index_t DIMENSION >
    index_t find_geological_entity_id_from_name(
        const RINGMesh::GeoModel< DIMENSION >& geomodel,
        const RINGMesh::GeologicalEntityType& gmge_type,
        const std::string& name );

#ifdef RINGMESH_WITH_TETGEN

    /*!
     * Compute the tetrahedral mesh of the input structural geomodel
     * @param[in/out] geomodel GeoModel to tetrahedralize
     * @param[in] method External mesher used, Tetgen by default
     * @param[in] region_id Region to mesh. By default it set to NO_ID and all regions are meshed.
     * @param[in] add_steiner_points if true (default value), the mesher will add some points inside the region.
     */
    void RINGMESH_API tetrahedralize( GeoModel< 3 >& geomodel, const std::string& method =
        "TetGen", index_t region_id = NO_ID, bool add_steiner_points = true );

    /*!
     * Compute the tetrahedral mesh of the input structural geomodel
     * @param[in/out] geomodel GeoModel to tetrahedralize
     * @param[in] method External mesher used
     * @param[in] region_id Region to mesh. If set to NO_ID and all regions are meshed.
     * @param[in] add_steiner_points if true, the mesher will add some points inside the region.
     * @param[in] internal_vertices points inside the domain to constrain mesh generation.
     * There is one vector per region.
     */
    void RINGMESH_API tetrahedralize(
        GeoModel< 3 >& geomodel,
        const std::string& method,
        index_t region_id,
        bool add_steiner_points,
        const std::vector< std::vector< vec3 > >& internal_vertices );

#endif

    /*!
     * @brief Translates the boundary geomodel by a vector.
     *
     * Every single mesh of the boundary geomodel is translated:
     * corners, lines and surfaces.
     *
     * @param[in/out] geomodel GeoModel on which compute the translation
     * @param[in] translation_vector vector of translation.
     */
    template< index_t DIMENSION >
    void translate(
        GeoModel< DIMENSION >& geomodel,
        const vecn< DIMENSION >& translation_vector );

    /*!
     * \brief Rotate the boundary geomodel.
     *
     * Applies a rotation about the line defined by the point
     * \p origin and the vector \p axis. The rotation angle is
     * \p angle. If \p degrees is true the angle is in degrees,
     * else in radians. All the vertices of the boundary geomodel
     * undergo the rotation (each mesh inside the boundary geomodel:
     * corners, lines and surfaces).
     *
     * @param[in/out] geomodel GeoModel on which compute the rotation
     * @param[in] origin point in which passes the rotation axis.
     * @param[in] axis vector which defines the rotation axis.
     * @param[in] angle rotation angle (in radians or degrees).
     * @param[in] degrees true is \p angle is in degrees, false
     * if in radians.
     */
    void RINGMESH_API rotate(
        GeoModel< 3 >& geomodel,
        const vec3& origin,
        const vec3& axis,
        double angle,
        bool degrees = false );
}
