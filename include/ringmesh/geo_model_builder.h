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

#ifndef __RINGMESH_GEO_MODEL_BUILDER__
#define __RINGMESH_GEO_MODEL_BUILDER__

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
}

namespace RINGMesh {
    /*!
     * @brief First draft of flags to build a GeoModel
     * @todo Implements functions to set, access the values, depending on what ?
     * To check the consistency of the options. What do we do about the other entities ? [JP] 
     * @todo We need to keep track of the status of the GeoModel when building it:
     * same flags or some others ?    
     */
    class GeoModelBuildingFlags {
    public:
        GeoModelBuildingFlags()
        {
            compute_corners = false ;
            compute_lines = false ;
            compute_surfaces = false ;
            compute_regions_brep = false ;
            compute_regions_mesh = false ;
        }
        bool compute_corners ;
        bool compute_lines ;
        bool compute_surfaces ;
        bool compute_regions_brep ;
        bool compute_regions_mesh ;
    } ;

    /*!
     * @brief Base class for all classes building a GeoModel.
     * @details Derive from this class to build or modify a GeoModel. 
     * @note NON Geometry related modifications are in GeoModelEditor class.
     * @todo To refactor and rename. We need a GeoModelTopologyEditor 
     * and a GeoModelGeometryEditor
     */
    class RINGMESH_API GeoModelBuilder: public GeoModelEditor {
    public:
        GeoModelBuilder( GeoModel& model )
            : GeoModelEditor( model ), options_()
        {
        }
        virtual ~GeoModelBuilder() ;

        /*!
         * @todo Implements sot that it returns true if the input options are consistent
         */
        void set_options( const GeoModelBuildingFlags& options )
        {
            options_ = options ;
        }

        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The model under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;

        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;

        void assign_mesh_to_entity( const Mesh& mesh, GME::gme_t to ) ;

        /*!
         * \name Set entity geometry from geometrical positions
         * @{
         */
        void set_mesh_entity_vertex(
            const GME::gme_t& entity_id,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_mesh_entity_vertices(
            const GME::gme_t& entity_id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        void set_corner( index_t corner_id, const vec3& point ) ;

        void set_line( index_t id, const std::vector< vec3 >& vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */
        void set_mesh_entity_vertex(
            const GME::gme_t& id,
            index_t v,
            index_t model_vertex ) ;

        void set_mesh_entity_vertices(
            const GME::gme_t& entity_id,
            const std::vector< index_t >& model_vertices,
            bool clear ) ;

        index_t add_unique_vertex( const vec3& p ) ;

        void set_corner( index_t corner_id, index_t unique_vertex ) ;

        void set_line( index_t id, const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners ) ;

        void set_surface_geometry_with_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners,
            const std::vector< index_t >& adjacent_triangles ) ;

        void set_surface_element_geometry(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& corners ) ;

        void set_surface_element_adjacency(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& adjacents ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< index_t >& tet_corners ) ;

        void set_region_element_geometry(
            index_t region_id,
            index_t cell_id,
            const std::vector< index_t >& corners ) ;

        /*! @}
         * \name Create entity element
         * @{
         */

        index_t create_mesh_entity_vertices(
            const GME::gme_t& entity_id,
            index_t nb_vertices ) ;

        index_t create_surface_facet(
            index_t surface_id,
            const GEO::vector< index_t >& vertex_indices ) ;

        index_t create_region_cell(
            index_t region_id,
            GEO::MeshCellType type,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cells(
            index_t region_id,
            GEO::MeshCellType type,
            index_t nb_cells ) ;

        /*! @}
         * \name Delete mesh element entities
         * @{
         */

        void delete_mesh_entity_mesh( GME::gme_t E_id ) ;
        void delete_mesh_entity_vertices( GME::gme_t E_id, GEO::vector< index_t >& to_delete ) ;
        void delete_corner_vertex( index_t corner_id ) ;
        void delete_line_edges( index_t line_id, GEO::vector< index_t >& to_delete ) ;
        void delete_surface_facets( index_t surface_id, GEO::vector< index_t >& to_delete ) ;
        void delete_region_cells( index_t region_id, GEO::vector< index_t >& to_delete ) ;

        /*! @}
         * \name Misc
         * @{
         */
        index_t find_or_create_duplicate_vertex(
            const GME::gme_t& E_id,
            index_t model_vertex_id,
            index_t surface_vertex_id ) ;

        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;

        void compute_surface_adjacencies( index_t surface_id ) ;
        void compute_region_adjacencies( index_t region_id ) ;
        void triangulate_surface(
            const RINGMesh::Surface& surface_in,
            index_t surface_out ) ;

        GME::gme_t find_or_create_corner( const vec3& point ) ;
        GME::gme_t find_or_create_corner( index_t model_point_id ) ;
        GME::gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
        GME::gme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            GME::gme_t first_corner,
            GME::gme_t second_corner ) ;

        void recompute_geomodel_mesh() ;

        /*!
         * @}
         * \name Model building functions
         */

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces() ;

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces must have been called before
         */
        bool build_brep_regions_from_surfaces() ;

        /*
         * @brief From a GeoModel in which only Surface are defined, create corners, contacts
         * and regions depending on the building flags
         * @note Valdity is not checked
         */
        void build_model_from_surfaces() ;

        /*!
         * @brief Finish up model building and complete missing information.
         */
        void end_model() ;

    protected:
        void build_contacts() ;
        void set_surface_facet_adjacencies(
                index_t surface_id,
                const std::vector< index_t >& facets_id,
                const std::vector< index_t >& edges_id,
                const std::vector< index_t >& adjacent_triangles ) ;

    protected:
        /*! Entities to compute from the available entities */
        GeoModelBuildingFlags options_ ;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;

    private:
        void assign_surface_mesh_facets(
            index_t surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices ) ;
        void update_facet_corner(
            Surface& S,
            const std::vector< index_t >& facets,
            index_t old,
            index_t neu ) ;
        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) ;

        void duplicate_surface_vertices_along_line( index_t surface_id, index_t line_id ) ;
        void disconnect_surface_facets_along_line_edges(
            index_t surface_id, index_t line_id ) ;
    } ;

    /*!
     * @brief Abstract class to load and build GeoModels from files 
     */
    class RINGMESH_API GeoModelBuilderFile: public GeoModelBuilder {
    public:
        GeoModelBuilderFile( GeoModel& model, const std::string& filename ) ;

        virtual ~GeoModelBuilderFile()
        {
        }
        void build_model()
        {
            load_file() ;
            end_model() ;
        }

    private:
        virtual void load_file() = 0 ;
        /*! @todo Implement function to read the lines of the 
         *        file and wrap the GEO::LineInput which is not that easy to use 
         */
    protected:
        std::string filename_ ;
    } ;
}

#endif
