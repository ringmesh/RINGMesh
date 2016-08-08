/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#ifndef __RINGMESH_GEO_MODEL_API__
#define __RINGMESH_GEO_MODEL_API__

#include <ringmesh/basic/common.h>

#include <vector>

#include <geogram/basic/attributes.h>

#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_entity.h>
#include <ringmesh/geomodel/geo_model_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>

/*!
 * @file ringmesh/geo_model_api.h
 * @brief High level functions on GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 * @todo Encapsulate these functions in a namespace and TEST them.
 */

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    /*!
     * @brief Print in the console the model statistics
     * @details Output number of facets, vertices, and of the different entity types
     * @todo Implement a test are_geomodels_equals to be able to check that tests went well
     */
    void RINGMESH_API print_geomodel( const GeoModel& geomodel ) ;

    /*!
     * Output the number of vertices, edges, facets and cells.
     * Also output the number of triangles, quads and polygons if any.
     * Also output the number of tetra, prisms, pyramids, hex and polyhedra if any.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    void RINGMESH_API print_geomodel_mesh_stats( const GeoModel& geomodel ) ;

    /*!
     * Output the volume of the geomodel and the volume per cell type.
     * @param[in] geomodel the geomodel to compute the statistics on
     */
    void RINGMESH_API print_geomodel_mesh_cell_volumes( const GeoModel& geomodel ) ;

    bool RINGMESH_API are_geomodel_surface_meshes_simplicial(
        const GeoModel& geomodel ) ;

    bool RINGMESH_API are_geomodel_region_meshes_simplicial(
        const GeoModel& geomodel ) ;

    /*!
     * @brief Build a Mesh from the model non-duplicated vertices and its Surface facets.
     * @details Adjacencies are not set. Client should call mesh repair functions afterwards.
     * @todo Add flag options to specify which Mesh should be created, with what attributes.
     */
    void RINGMESH_API build_mesh_from_geomodel(
        const GeoModel& model,
        GEO::Mesh& M ) ;

    void RINGMESH_API build_mesh_from_geomodel(
        const GeoModel& model,
        GEO::Mesh& M,
        bool connect_facets ) ;

    void RINGMESH_API build_mesh_from_model_mesh_entities(
        const GeoModel& model,
        const std::vector< gme_t >& surface_entities,
        GEO::Mesh& M ) ;

    /*! 
     * @brief Bind named GEO::Attribute on the GeoModel entity facets
     * @pre Entities of geomodel_entity_type are GeoModelMeshEntity
     */
    template< typename T >
    void create_attributes_on_geomodel_surfaces_facets(
        const GeoModel& geomodel,
        const std::string& attribute_name,
        AttributeVector< T >& attributes )
    {
        index_t nb_entities = geomodel.nb_surfaces() ;
        attributes.resize( nb_entities ) ;
        for( index_t i = 0; i < nb_entities; ++i ) {
            const Surface& S = geomodel.surface( i ) ;
            GEO::AttributesManager& manager = S.facet_attribute_manager() ;
            attributes.bind_one_attribute( i, manager, attribute_name ) ;
        }
    }

    /*!
     * @brief Bind named GEO::Attribute on the GeoModel entities cells
     * @pre Entities of mesh_entity_type are GeoModelMeshEntity
     */
    template< typename T >
    void create_attributes_on_geomodel_regions_cells(
        const GeoModel& geomodel,
        const std::string& attribute_name,
        AttributeVector< T >& attributes )
    {
        index_t nb_entities = geomodel.nb_regions() ;
        attributes.resize( nb_entities ) ;
        for( index_t i = 0; i < nb_entities; ++i ) {
            const Region& R = geomodel.region( i ) ;
            GEO::AttributesManager& manager = R.cell_attribute_manager() ;
            attributes.bind_one_attribute( i, manager, attribute_name ) ;
        }
    }

#ifdef RINGMESH_WITH_TETGEN
    /*!
     * @brief Tetrahedralize the B-Rep defined regions of a GeoModel
     */
    void RINGMESH_API tetgen_tetrahedralize_geomodel_regions( GeoModel& geomodel ) ;
#endif

    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param[in] geomodel GeoModel to tetrahedralize
     * @param[in] method External mesher used, Tetgen by default
     * @param[in] region_id Region to mesh. By default it set to NO_ID and all regions are meshed.
     * @param[in] add_steiner_points if true (default value), the mesher will add some points inside the region.
     */
    void RINGMESH_API tetrahedralize( GeoModel& geomodel, const std::string& method =
        "TetGen", index_t region_id = NO_ID, bool add_steiner_points = true ) ;

    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param[in] geomodel GeoModel to tetrahedralize
     * @param[in] method External mesher used
     * @param[in] region_id Region to mesh. If set to NO_ID and all regions are meshed.
     * @param[in] add_steiner_points if true, the mesher will add some points inside the region.
     * @param[in] internal_vertices points inside the domain to constrain mesh generation.
     * There is one vector per region.
     */
    void RINGMESH_API tetrahedralize(
        GeoModel& geomodel,
        const std::string& method,
        index_t region_id,
        bool add_steiner_points,
        const std::vector< std::vector< vec3 > >& internal_vertices ) ;

    /*!
     * @brief Translates the boundary model by a vector.
     *
     * Every single mesh of the boundary model is translated:
     * corners, lines and surfaces.
     *
     * @param[in] geomodel GeoModel on which compute the translation
     * @param[in] translation_vector vector of translation.
     */
    void RINGMESH_API translate( GeoModel& geomodel, const vec3& translation_vector ) ;

    /*!
     * \brief Rotate the boundary model.
     *
     * Applies a rotation about the line defined by the point
     * \p origin and the vector \p axis. The rotation angle is
     * \p angle. If \p degrees is true the angle is in degrees,
     * else in radians. All the vertices of the boundary model
     * undergo the rotation (each mesh inside the boundary model:
     * corners, lines and surfaces).
     *
     * @param[in] geomodel GeoModel on which compute the rotation
     * @param[in] origin point in which passes the rotation axis.
     * @param[in] axis vector which defines the rotation axis.
     * @param[in] angle rotation angle (in radians or degrees).
     * @param[in] degrees true is \p angle is in degrees, false
     * if in radians.
     */
    void RINGMESH_API rotate(
        GeoModel& geomodel,
        const vec3& origin,
        const vec3& axis,
        float64 angle,
        bool degrees = false ) ;

    /*-----------------------------------------------------------------------*/

    /*!
     * @brief Compute the size (volume, area, length) of an Entity
     * @param[in] E Entity to evaluate
     */
    double RINGMESH_API model_entity_size( const GeoModelMeshEntity& E ) ;
    double RINGMESH_API model_entity_size( const GeoModelGeologicalEntity& E ) ;

    /*!
     * Compute the size (volume, area, length) of an Entity cell (cell, facet, edge)
     * @param[in] E Entity to evaluate
     * @param[in] c the cell index
     */
    double RINGMESH_API model_entity_cell_size(
        const GeoModelEntity& E,
        index_t c ) ;

    /*!
     * @brief Compute the barycenter of a GeoModelEntity
     * @param[in] E Entity to evaluate
     * @return The coordinates of the barycenter
     */
    vec3 RINGMESH_API model_entity_barycenter( const GeoModelEntity& E ) ;

    /*!
     * @brief Compute the centroid of a GeoModelMeshEntity cell (cell, facet, edge)
     * @param[in] E Entity to evaluate
     * @param[in] c the cell index
     * @return The coordinates of the center
     * @pre E has a valid mesh.
     */
    vec3 RINGMESH_API model_entity_cell_barycenter(
        const GeoModelMeshEntity& E,
        index_t c ) ;

    /*-----------------------------------------------------------------------*/

    /*!
     * @brief Gets the index of the Corner for a given point
     * @param[in] geomodel GeoModel to consider
     * @param[in] point Geometric location to look for
     * @return NO_ID or the index of the Corner
     */
    gme_t RINGMESH_API find_corner( const GeoModel& geomodel, const vec3& point ) ;

    /*!
     * @brief Gets the index of the Corner at a given model point
     * @param[in] geomodel GeoModel to consider
     * @param[in] model_point_id Index of the point in the GeoModel
     * @return NO_ID or the index of the Corner
     */
    gme_t RINGMESH_API find_corner(
        const GeoModel& geomodel,
        index_t model_point_id ) ;

}

#endif 
