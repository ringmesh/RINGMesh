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
                std::vector< index_t >& nb_elements_per_region ) ;

        void print_number_of_mesh_elements(
            const std::vector< index_t >& nb_elements_per_region ) const ;

        /*! @}
         * \name Other functions
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
        void read_vertex_coordinates( vec3& vertex ) ;

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
            std::vector< index_t >& corners_id ) ;

        /*!
         * @brief Sets the boundaries of the GeoModel regions
         */
        void compute_boundaries_of_geomodel_regions() ;

        /*!
         * @brief Sets the boundaries of region Universe
         * @details A surface is set in the boundaries of region Universe if
         * only one of its sides belongs to the boundaries of other regions.
         */
        void compute_universe_boundaries() ;

        void build_surface(
            index_t surface_id,
            std::vector< index_t >& facet_corners,
            std::vector< index_t >& facet_ptr,
            const std::vector< index_t >& gocad_vertices2region_id,
            const std::vector< index_t >& gocad_vertices2region_vertices ) ;

        void compute_internal_borders() ;

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
