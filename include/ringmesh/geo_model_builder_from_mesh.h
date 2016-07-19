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

#ifndef __RINGMESH_GEO_MODEL_BUILDER_FROM_MESH__
#define __RINGMESH_GEO_MODEL_BUILDER_FROM_MESH__

#include <ringmesh/common.h>

#include <vector>
#include <string>
#include <stack>

#include <geogram/basic/line_stream.h>
#include <ringmesh/geo_model_editor.h>
#include <third_party/zlib/unzip.h>

#define MAX_FILENAME 512
#define READ_SIZE 8192

/*!
 * @file ringmesh/geo_model_builder.h
 * @brief Classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    class GeoModelRegionFromSurfaces ;
    class GeoModelEntityFromMesh ;
}

namespace RINGMesh {

    /*!
     * @brief To build a GeoModel from a set of disconnected polygonal surfaces
     */
    class RINGMESH_API GeoModelBuilderSurfaceMesh: public GeoModelBuilder {
    public:
        GeoModelBuilderSurfaceMesh( GeoModel& model, const GEO::Mesh& mesh )
            : GeoModelBuilder( model ), mesh_( mesh )
        {
            options_.compute_lines = true ;
            options_.compute_corners = true ;
            options_.compute_regions_brep = true ;
        }
        void build_polygonal_surfaces_from_connected_components() ;

    private:
        const GEO::Mesh& mesh_ ;
    } ;

    /*!
     * @brief Builder of a GeoModel from a simplicial surface or volumetric Mesh 
     * @details Regions and Surfaces are identified with an attribute of type index_t
     * on the mesh cells or facet 
     */
    class RINGMESH_API GeoModelBuilderMesh: public GeoModelBuilder {
    public:
        GeoModelBuilderMesh(
            GeoModel& model,
            const GEO::Mesh& mesh,
            const std::string& surface_attribute_name,
            const std::string& region_attribute_name )
            :
                GeoModelBuilder( model ),
                mesh_( mesh ),
                surface_builder_( nil ),
                region_builder_( nil ),
                surface_attribute_name_( surface_attribute_name ),
                region_attribute_name_( region_attribute_name )
        {
            initialize_surface_builder() ;
            initialize_region_builder() ;
            add_mesh_vertices_to_model() ;
        }

        virtual ~GeoModelBuilderMesh() ;

        /*!
         * @brief Prepare a Mesh so that it can be used to build one GeoModel Surfaces
         * @details Repairs the mesh, triangulates it, computes a connected component 
         * attribute of type index_t on the mesh facets and removes colocated vertices. 
         */
        static void prepare_surface_mesh_from_connected_components(
            GEO::Mesh& mesh,
            const std::string& created_surface_attribute ) ;

        void create_and_build_surfaces() ;
        void build_surfaces() ;

        void create_and_build_regions() ;
        void build_regions() ;

        void copy_facet_attribute_from_mesh( const std::string& attribute_name ) ;
        void copy_cell_attribute_from_mesh( const std::string& attribute_name ) ;

    protected:
        /*!
         * @brief Set the unique vertices used to build the GeoModel
         * @details They are cleared when end_model() is called
         */
        void add_mesh_vertices_to_model() ;

        void initialize_surface_builder() ;
        void initialize_region_builder() ;

    private:
        bool is_mesh_valid_for_surface_building() const ;
        bool is_mesh_valid_for_region_building() const ;

    protected:
        const GEO::Mesh& mesh_ ;
        GeoModelEntityFromMesh* surface_builder_ ;
        GeoModelEntityFromMesh* region_builder_ ;
        std::string surface_attribute_name_ ;
        std::string region_attribute_name_ ;
        index_t nb_surface_attribute_values_ ;
        index_t nb_region_attribute_values_ ;
    } ;
}

#endif
