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
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 *
 *
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */


#ifndef __RINGMESH_GEO_MODEL_BUILDER_SO__
#define __RINGMESH_GEO_MODEL_BUILDER_SO__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_builder.h>

namespace RINGMesh {
    /*!
     * @brief Builds a meshed GeoModel from a Gocad TSolid (file.so)
     */
    class RINGMESH_API GeoModelBuilderTSolid : public GeoModelBuilderFile {
    public:
        GeoModelBuilderTSolid( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename )
        {
            filename_ = filename ;
            z_sign_ = 1 ;
        }
        virtual ~GeoModelBuilderTSolid()
        {}
        bool load_file() ;

    private:

        /*!
         * \name Reads and sets Gocad Coordinate System information from .so file
         * @{
         */

        /*!
        * @brief Reads and sets the Gocad coordinate system information
        * from input .so file.
        */
        void read_and_set_gocad_coordinate_system() ;

        /*! @}
         * \name Volume mesh import
         * @{
         */

        /*!
         * @brief Reads vertex coordinates and adds it in the list
         * of region vertices
         * @param[in] region_id Index of the region
         * @param[in,out] region_vertices Vector of the coordinates of the
         * vertices of the region
         */
        void read_and_add_vertex_to_region_vertices(
            const index_t region_id,
            std::vector < vec3 >& region_vertices ) ;

        /*!
         * @brief Reads atom information and adds it in the list
         * of region vertices only if it refers to a vertex of another region
         * @param[in] region_id Index of the region
         * @param[in,out] region_vertices Vector of the coordinates of the
         * vertices of the region
         */
        void read_and_add_atom_to_region_vertices(
            const index_t region_id,
            std::vector < vec3 >& region_vertices ) ;

        /*! @}
         * \name Boundary model import
         * @{
         */

        /*!
         * Builds surface by setting the points and facets of the surface
         * @param[in] surface_id Index of the surface to build
         * @param[in,out] facet_corners Vector of the (gocad) indices of the
         * three corners of each facet (gocad) indices of the surface.
         * Re-initialized at the end of the function.
         * @param[in,out] facet_ptr Pointer to the beginning of a facet in
         * facets. Re-initialized at the end of the function.
         */
        void build_surface(
            const index_t surface_id,
            std::vector< index_t >& facet_corners,
            std::vector< index_t >& facet_ptr ) ;

        /*!
         * @brief Gets the points and the indices in the points vector to
         * build the facets
         * @param[in] facet_corners Vector of the (gocad) indices of the
         * three corners of each facet (gocad) indices
         * @param[out] cur_surf_points Vector of unique point coordinates
         * belonging to the surface
         * @param[out] cur_surf_facets Vector of each facet corner indices in
         * the cur_surf_points vector to build facets
         */
        void get_surface_points_and_facets_from_gocad_indices(
            const std::vector< index_t >& facet_corners,
            std::vector< vec3 >& cur_surf_points,
            std::vector< index_t >& cur_surf_facets ) const ;

        /*!
         * @brief Gets the point and the index in the points vector to
         * build the facets for one read gocad vertex
         * @param[in] vertex_gocad_id Gocad index of the vertex
         * @param[in] gocad_vertices2cur_surf_points Map between vertices with
         * gocad indices and position of the corresponding point
         * in the points vector
         * @param[out] cur_surf_points Vector of unique point coordinates
         * belonging to the surface
         * @param[out] cur_surf_facets Vector of each facet corner indices in
         * the cur_surf_points vector to build facets
         */
        void get_surface_point_and_facet_from_gocad_index(
            const index_t vertex_gocad_id,
            std::vector< index_t >& gocad_vertices2cur_surf_points,
            std::vector< vec3 >& cur_surf_points,
            std::vector< index_t >& cur_surf_facets ) const ;

        /*!
         * @brief Gets the coordinates of the point from gocad index
         * @param[in] point_gocad_id Gocad index of the point to get
         * @param[out] point Coordinates of the point
         */
        void get_point_from_gocad_id(
            const index_t point_gocad_id,
            vec3& point ) const ;

        /*! @}
         */

    private:
        std::string filename_ ;
        int z_sign_ ;
        std::string gocad_coordinates_system_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_unit_ ;
        /*!
         * Vector which maps the indices of vertices from Gocad .so file
         * to the local (in region) indices of vertices
         */
        std::vector< index_t > gocad_vertices2region_vertices_ ;
        /*!
         * Vector which maps the indices of vertices from Gocad .so file
         * to the index of the region they belong to
         */
        std::vector< index_t > gocad_vertices2region_id_ ;
    } ;
}

#endif
