/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
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
         * \name Reads and sets Gocad Coordinates System information from .so file
         * @{
         */

        /*!
        * Reads and sets the Gocad coordinates system information
        * from input .so file.
        */
        void read_and_set_gocad_coordinates_system() ;

        void set_gocad_coordinates_system_axis_name() ;

        void set_gocad_coordinates_system_axis_unit() ;

        void set_gocad_coordinates_system_z_sign() ;

        /*! @}
         * \name Information on the number of mesh elements from .so file
         * @{
         */

        ///@todo comment
        void read_number_of_mesh_elements(
                std::vector< index_t >& nb_elements_per_region ) const ;

        void print_number_of_mesh_elements(
            const std::vector< index_t >& nb_elements_per_region ) const ;

        /*! @}
         * \name Other helper functions
         * @{
         */

        void add_new_property(
            std::vector< std::string >& property_names,
            GEO::AttributesManager& attribute_manager ) ;

        GME::gme_t create_region() ;

        /*!
         * Reads the coordinates of a vertex from file
         * @param[out] vertex Vertex
         */
        void read_vertex_coordinates( vec3& vertex ) const ;

        /*!
         * @brief Reads the four vertices index
         * @details Reads gocad indices (from .so file) and transforms
         * them to vertex local (region) indices
         * @param[in] gocad_vertices2region_vertices Map from the gocad vertex
         * indices to vertex local indices (in region)
         * @param[out] corners_id Indices of the four vertices
         */
        void read_tetraedra(
            const std::vector< index_t >& gocad_vertices2region_vertices,
            std::vector< index_t >& corners_id ) const ;

        /*!
         * @brief Sets the boundaries of the GeoModel regions
         */
        void compute_boundaries_of_geomodel_regions() ;

        /*!
         * @brief Computes the colocaters of the centers of cell facets for
         * each region
         * @param[out] region_anns Pointers to the ColocaterANNs
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

        void build_surface(
            index_t surface_id,
            std::vector< index_t >& facet_corners,
            std::vector< index_t >& facet_ptr,
            const std::vector< index_t >& gocad_vertices2region_id,
            const std::vector< index_t >& gocad_vertices2region_vertices ) ;

        void compute_internal_borders() ;

        /*!
         * @brief Both adds the surface in the boundaries of a region and
         * add the region to the in_boundaries of the surface
         * @param[in] region_id Index of the region
         * @param[in] surface_id Index of the surface
         * @param[in] surf_side Side of the surface bounding the region
         */
        void fill_region_and_surface_boundaries_links(
            const index_t region_id,
            const index_t surface_id,
            const bool surf_side ) ;

        /*! @}
         */

    private:
        std::string filename_ ;
        int z_sign_ ;
        std::string gocad_coordinates_system_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_unit_ ;
    } ;
}

#endif
