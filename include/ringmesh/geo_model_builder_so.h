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
#include <ringmesh/utils.h>

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

        /*!
         * @brief Sets the Gocad coordinate system axis name
         */
        void set_gocad_coordinate_system_axis_name() ;

        /*!
         * @brief Sets the Gocad coordinate system axis units (meters, feet, ...)
         */
        void set_gocad_coordinate_system_axis_unit() ;

        /*!
         * @brief Sets the Gocad coordinate system Z sign
         * @details If the Z values increase upwards, z_sign_ is positive,
         * else (Z values increasing downwards), z_sign_ is negative.
         */
        void set_gocad_coordinate_system_z_sign() ;

        /*! @}
         * \name Information on the number of mesh elements from .so file
         * @{
         */

        /*!
         * @brief Counts number of vertices and tetras in each region
         * @details Reserve also space from the attributes which maps vertex
         * indices between gocad and region local indices.
         * @param[out] Vector built from number of vertices and tetras region
         * after region (i.e. [nb_v1, nb_t1, nb_v2, nb_t2, nb_v3, nb_t3, ...]
         */
        void count_nb_vertices_and_tetras_per_region(
                std::vector< index_t >& nb_elements_per_region ) ;

        /*!
         * @brief Shows number of vertices and tetras in each region
         * @param[in] Vector built from number of vertices and tetras region
         * after region (i.e. [nb_v1, nb_t1, nb_v2, nb_t2, nb_v3, nb_t3, ...]
         */
        void print_nb_vertices_and_tetras_per_region(
            const std::vector< index_t >& nb_elements_per_region ) const ;

        /*! @}
         * \name Properties import
         * @{
         */

        void add_new_property(
            std::vector< std::string >& property_names,
            GEO::AttributesManager& attribute_manager ) ;

        /*! @}
         * \name Volume mesh import
         * @{
         */

        /*!
         * @brief Creates an empty element of type GME::REGION and sets
         * its name from .so file
         * @return The index of the initialized region
         */
        index_t initialize_region() ;

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

        /*!
         * Builds region by setting the points and tetras of the region
         * @param[in] region_id Index of the surface to build
         * @param[in] nb_vertices_in_next_region Number of vertices in the
         * next region (to reverse space)
         * @param[in] nb_tetras_in_next_region Number of tetrahedra in the
         * next region (to reverse space)
         * @param[in,out] region_vertices Vector of the coordinates of the
         * vertices of the region. Re-initialized at the end of the function.
         * @param[in,out] tetra_corners Vector of the region tetrahedra corner
         * indices. Re-initialized at the end of the function.
         */
        void build_region(
            const index_t region_id,
            const index_t nb_vertices_in_next_region,
            const index_t nb_tetras_in_next_region,
            std::vector < vec3 >& region_vertices,
            std::vector < index_t >& tetra_corners ) ;

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
         * \name Read mesh elements (points, triangles, tetrehedra)
         * @{
         */

        /*!
         * Reads the coordinates of a vertex from file
         * @param[out] vertex Vertex
         */
        void read_vertex_coordinates( vec3& vertex ) const ;

        /*!
         * @brief Reads the four vertices index of a tetrahedron
         * @details Reads gocad indices (from .so file) and transforms
         * them to vertex local (region) indices
         * @param[out] corners_id Indices of the four vertices
         */
        void read_tetraedra( std::vector< index_t >& corners_id ) const ;

        /*!
         * @brief Reads the three vertices index of a triangle and adds
         * them to the facet corners of the currently built surface
         * @details Reads gocad indices
         * @param[out] cur_surf_facets Vector of each facet corner indices
         * to build facets
         */
        void read_triangle( std::vector< index_t >& cur_surf_facets ) const ;

        /*! @}
         * \name Linking surfaces and region boundaries
         * @{
         */

        /*!
         * @brief Sets the boundaries of the GeoModel regions
         */
        void compute_boundaries_of_geomodel_regions() ;

        /*!
         * @brief Computes the colocaters of the centers of cell facets for
         * each region
         * @param[out] region_anns Pointers to the ColocaterANNs of regions
         */
        void compute_cell_facet_centers_region_anns(
            std::vector< ColocaterANN* >& region_anns ) const ;

        /*!
         * @brief Builds a vector with the center of the cell
         * facets of a given region
         * @param[in] region_id Index of the region
         * @param[out] cell_facet_centers Vector of cell facet centers
         */
        void compute_region_cell_facet_centers(
            const index_t region_id,
            std::vector< vec3 >& cell_facet_centers ) const ;

        /*!
         * @brief Sets the given surface as regions boundaries
         * @details Based on ColocaterANN, retrieves the regions bounded by the
         * given surface. One side or the both sides of the surface
         * could bound model regions.
         * @param[in] surface_id Index of the surface
         * @param[in] region_anns Vector of ColocaterANN of the model regions
         */
        void add_surface_to_region_boundaries(
            const index_t surface_id,
            const std::vector< ColocaterANN* >& region_anns ) ;

        /*!
         * @brief Tests if a surface is a boundary of the region.
         * @details If it is the case, add the surface to the boundaries of
         * the region and the region to the in_boundaries of the surface
         * @param[in] surface_id Index of the surface
         * @param[in] region_id Index of the region
         * @param[in] region_anns Vector of ColocaterANN of the model regions
         * @param[out] colocated_cell_facet_centers Vector of colocated cell
         * facet centers
         * @return The number of surface sides bounding the region
         */
        index_t are_surface_sides_region_boundaries(
            const index_t surface_id,
            const index_t region_id,
            const ColocaterANN& region_ann,
            std::vector< index_t >& colocated_cell_facet_centers ) const ;

        /*!
         * @brief Adds the surface sides which bound the region to the
         * boundaries of the region (and add the region to in boundaries
         * of the surface)
         * @param[in] surface_id Index of the surface
         * @param[in] region_id Index of the region
         * @param[in] colocated_cell_facet_centers Vector of colocated cell
         * facet centers
         */
        void add_surface_sides_to_region_boundaries(
            const index_t surface_id,
            const index_t region_id,
            const std::vector< index_t >& colocated_cell_facet_centers ) ;

        /*!
         * @brief Adds one surface side in the boundaries of a region
         * and add the region to the in_boundaries of the surface
         * @details The index of the cell facet center is used for the
         * determination of the side to add.
         * @param[in] region_id Index of the region
         * @param[in] surface_id Index of the surface
         * @param[in] cell_facet_center_id Index of the cell facet center
         * (i.e., cell_id * 4 + local_facet_id)
         */
        void add_one_surface_side_to_region_boundaries(
            const index_t region_id,
            const index_t surface_id,
            const index_t cell_facet_center_id ) ;

        /*!
         * @brief Determines which side of the surface is to be added in the
         * region boundaries
         * @param[in] region_id Index of the region
         * @param[in] surface_id Index of the surface
         * @param[in] cell_facet_center_id Index of the cell facet center
         * (i.e., cell_id * 4 + local_facet_id)
         * @return The side of the surface to add in the region boundaries,
         * i.e. true for the '+' side (normal size) and false for the
         * '-' side (other side)
         */
        bool determine_surface_side_to_add(
            const index_t region_id,
            const index_t surface_id,
            const index_t cell_facet_center_id ) const ;

        /*!
         * @brief Adds the both surface sides in the boundaries of a region
         * (internal boundary) and add twice the region to the in_boundaries
         * of the surface
         * @param[in] region_id Index of the region
         * @param[in] surface_id Index of the surface
         */
        void add_both_surface_sides_to_region_boundaries(
            const index_t region_id,
            const index_t surface_id ) ;

        /*!
         * @brief Sets the boundaries of region Universe
         * @details A surface is set in the boundaries of region Universe if
         * only one of its sides belongs to the boundaries of other regions.
         */
        void compute_universe_boundaries() ;

        /*
         * @brief Determines if each side of the surfaces are
         * in the boundaries of model regions
         * @param[out] surf_side_minus Vector indicating if the '-' side of
         * surfaces are in the boundaries of model regions
         * @param[out] surf_side_plus Vector indicating if the '+' side of
         * surfaces are in the boundaries of model regions
         */
        void determine_if_surface_sides_bound_regions(
            std::vector< bool >& surf_minus_side,
            std::vector< bool >& surf_plus_side ) const ;

        /*
         * @brief Adds the right surface sides in universe boundaries
         * @param[in] surf_side_minus Vector indicating if the '-' side of
         * surfaces are in the boundaries of model regions
         * @param[in] surf_side_plus Vector indicating if the '+' side of
         * surfaces are in the boundaries of model regions
         */
        void add_surfaces_to_universe_boundaries(
            const std::vector< bool >& surf_minus_side,
            const std::vector< bool >& surf_plus_side ) ;

        /*!
         * @brief Both adds the surface in the boundaries of a region and
         * adds the region to the in_boundaries of the surface
         * @param[in] region_id Index of the region
         * @param[in] surface_id Index of the surface
         * @param[in] surf_side Side of the surface bounding the region
         */
        void fill_region_and_surface_boundaries_links(
            const index_t region_id,
            const index_t surface_id,
            const bool surf_side ) ;

        /*! @}
         * \name Surface internal borders determination
         * @{
         */

        /*!
         * @brief Computes internal borders of the model surfaces
         * @details An surface facet edge is an internal border if it is shared
         * by at least two surfaces. Adjacency of such a facet edge is set to
         * GEO::NO_FACET.
         */
        void compute_surfaces_internal_borders() ;

        /*!
         * @brief Computes internal borders of a given surface
         * @details A surface facet edge is an internal border if it is shared
         * by at least two surfaces. Adjacency of such a facet edge is set to
         * GEO::NO_FACET.
         * @param[in] surface_id Index of the surface
         * @param[in] surface_anns Pointers to the ColocaterANNs of surfaces
         */
        void compute_surface_internal_borders(
            const index_t surface_id,
            const std::vector< ColocaterANN* >& surface_anns,
            const std::vector< Box3d >& surface_boxes ) ;

        /*!
         * @brief Finds if a surface facet edge is an internal border
         * (i.e. shared by at least two surfaces)
         * @param[in] surface_id Index of the surface
         * @param[in] facet Index of the facet in the surface
         * @param[in] edge Index of the edge in the facet
         * @param[in] surface_anns Pointers to the ColocaterANNs of surfaces
         * @param[in] surface_boxes Bounding Box of surfaces
         * @return True is the edge is found in at least another surface
         */
        bool is_edge_in_several_surfaces(
            const index_t surface_id,
            const index_t facet,
            const index_t edge,
            const std::vector< ColocaterANN* >& surface_anns,
            const std::vector< Box3d >& surface_boxes ) const ;

        /*!
         * @brief Computes the colocaters of the centers of facet edges for
         * each surface and their Box3d
         * @param[out] surface_anns Pointers to the ColocaterANNs of surfaces
         * @param[out] surface_boxes Bounding Box of surfaces
         */
        void compute_facet_edge_centers_anns_and_surface_boxes(
            std::vector< ColocaterANN* >& surface_anns,
            std::vector< Box3d >& surface_boxes ) const ;

        /*!
         * @brief Gets the facet edge barycenters of a given surface
         * @param[in] surface_id Index of the surface
         * @param[out] border_edge_barycenters Vector of all the border
         * edge barycenters of the surface
         */
        void get_surface_border_edge_barycenters(
            const index_t surface_id,
            std::vector< vec3 >& border_edge_barycenters ) const ;

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
