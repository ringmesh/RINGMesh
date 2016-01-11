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
     * @brief Build a meshed GeoModel from a Gocad TSolid (file.so)
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
        * Read and set the Gocad coordinates system information
        * from input .so file.
        */
        void read_and_set_gocad_coordinates_system() ;
        ///@todo comment
        std::vector< index_t > read_number_of_mesh_elements() ;

        void print_number_of_mesh_elements(
            const std::vector< index_t >& nb_elements_per_region) const ;

        void add_new_property(
            std::vector < std::string >& property_names,
            GEO::AttributesManager& attribute_manager ) ;

        GME::gme_t create_region() ;

        void read_vertex_coordinates(
            GEO::Mesh* mesh,
            std::vector< index_t >& vertices_id_in_region ) ;

        /*
         * @brief Read the four vertices index
         * and regions depending on the building flags
         * @param[in] vertices_id_in_region Local indices (in region) of all the
         * model indices
         * @param[out] vertices Indices of the four vertices
         */
//        void read_tetraedra(
//            std::vector< index_t >& vertices_id_in_region,
//            std::vector< index_t >& vertices_id ) ;
        void read_tetraedra(
            std::vector< index_t >& map_gocad2gmm_vertices,
            std::vector< index_t >& map_gocad_vertices2region_id,
            std::vector< index_t >& vertices_id ) ;

    private:
        std::string filename_ ;
        int z_sign_ ;
        std::string gocad_coordinates_system_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_name_ ;
        std::vector< std::string > gocad_coordinates_system_axis_unit_ ;
    } ;
}

#endif
