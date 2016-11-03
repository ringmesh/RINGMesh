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

#include <ringmesh/basic/common.h>

#include <vector>
#include <string>
#include <stack>

#include <ringmesh/geomodel/geo_model_editor.h>

/*!
 * @file ringmesh/geo_model_builder.h
 * @brief Classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
     * @brief First draft of flags to build a GeoModel
     * @todo Implements functions to set, access the values, depending on what ?
     * To check the consistency of the options. What do we do about the other entities ? [JP] 
     * 
     * @todo We need to keep track of the status of the GeoModel when building it:
     * same flags or some others ?    
     *
     * @todo To separate in two classes ? One providing the low level functions set, assign etc,
     * and the other one some high level functions. [JP]
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

    // Implementation details
    class GeoModelRegionFromSurfaces ;


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
         * @todo Implement it so that it returns true if the input options are consistent
         */
        void set_options( const GeoModelBuildingFlags& options )
        {
            options_ = options ;
        }

        void copy( const GeoModel& from )
        {
            copy_macro_topology( from ) ;
            copy_meshes( from ) ;
        }
        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The model under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;
        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;
        void copy_mesh( const GeoModel& from, const gme_t& mesh_entity ) ;

        void assign_mesh_to_entity( const Mesh& mesh, const gme_t& to ) ;

        /*!
         * \name Set entity geometry from geometrical positions
         * @{
         */
        /*!
         * @brief Sets a vertex coordinates of a GeoModelMeshEntity
         * @param[in] entity_id the entity to edit
         * @param[in] v the index of the vertex in the entity
         * @param[in] point the coordinates to set
         * @param[in] update if true, updates all the colocated vertices
         * to the new coordinates (ie if edit a Corner coordinates, it will updates
         * its Lines, Surfaces...)
         */
        void set_mesh_entity_vertex(
            const gme_t& entity_id,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        /*!
         * @brief Sets the coordinates of a given existing Corner
         * @param[in] corner_id the index of the corner in the GeoModel
         * @param[in] point the coordinates to set
         */
        void set_corner( index_t corner_id, const vec3& point ) ;
        /*!
         * @brief Sets the mesh of a given existing Line
         * @param[in] line_id the index of the line in the GeoModel
         * @param[in] vertices the coordinates to set
         * @warning the vertices should be ordered from the first boundary
         * corner to the second one
         */
        void set_line( index_t line_id, const std::vector< vec3 >& vertices ) ;
        /*!
         * @brief Sets the mesh of a given existing Surface
         * @param[in] surface_id the index of the surface in the GeoModel
         * @param[in] surface_vertices the coordinates to set
         * @param[in] surface_facets the vertex indices of the facets
         * corresponding to \p surface_vertices
         * @param[in] surface_facet_ptr the index of each new facet start in \p surface_facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;
        /*!
         * @brief Sets the tetrahedral mesh of a given existing Region
         * @param[in] region_id the index of the region in the GeoModel
         * @param[in] points the coordinates to set
         * @param[in] tetras the vertex indices of the cells (to read 4 by 4)
         * corresponding to \p points
         */
        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */
        void set_mesh_entity_vertex(
            const gme_t& id,
            index_t v,
            index_t model_vertex ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< index_t >& model_vertices,
            bool clear ) ;

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
            const gme_t& entity_id,
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

        void delete_mesh_entity_mesh( const gme_t& E_id ) ;
        void delete_mesh_entity_isolated_vertices( const gme_t& E_id ) ;
        void delete_mesh_entity_vertices(
            const gme_t& E_id,
            GEO::vector< index_t >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_corner_vertex( index_t corner_id ) ;
        void delete_line_edges(
            index_t line_id,
            GEO::vector< index_t >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_surface_facets(
            index_t surface_id,
            GEO::vector< index_t >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_region_cells(
            index_t region_id,
            GEO::vector< index_t >& to_delete,
            bool remove_isolated_vertices ) ;

        /*! @}
         * \name Misc
         * @{
         */

        void compute_surface_adjacencies( index_t surface_id ) ;
        void compute_region_adjacencies( index_t region_id ) ;
        void triangulate_surface(
            const RINGMesh::Surface& surface_in,
            index_t surface_out ) ;

        gme_t find_or_create_corner( const vec3& point ) ;
        gme_t find_or_create_corner( index_t model_point_id ) ;
        gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
        gme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gme_t& first_corner,
            const gme_t& second_corner ) ;

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
        /*!
         * @brief Build the Contacts
         * @details One contact is a group of lines shared by the same Interfaces
         */
        void build_contacts() ;
        void invert_surface_normals( index_t surface_id ) ;
        index_t mesh_nb_connected_components( const gme_t& gmme_id ) const ;
        void remove_isolated_vertices( const gme_t& gmme_id ) ;

        void set_surface_facet_adjacencies(
                index_t surface_id,
                const std::vector< index_t >& facets_id,
                const std::vector< index_t >& edges_id,
                const std::vector< index_t >& adjacent_triangles ) ;
        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;
        index_t disconnect_surface_facets_along_line_edges(
            index_t surface_id,
            index_t line_id ) ;

    protected:
        /*! Options to toggle the building of entities from the available entities */
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
        void update_facet_vertices_around_facet_vertex(
            Surface& surface,
            index_t facet,
            index_t old_vertex,
            index_t new_vertex ) ;
        void update_facet_vertex(
            Surface& surface,
            const std::vector< index_t >& facets,
            index_t old_vertex,
            index_t new_vertex ) ;
        void update_cell_vertex(
            Region& region,
            const std::vector< index_t >& cells,
            index_t old_vertex,
            index_t new_vertex ) ;
        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) ;

        void compute_universe() ;

        void cut_surfaces_by_internal_lines() ;
        void cut_regions_by_internal_surfaces() ;

        void cut_region_by_surface( index_t region_id, index_t surface_id ) ;
        void duplicate_surface_vertices_along_line(
            index_t surface_id,
            index_t line_id ) ;
        void duplicate_region_vertices_along_surface(
            index_t region_id,
            index_t surface_id ) ;
        index_t disconnect_region_cells_along_surface_facets(
            index_t region_id, index_t surface_id ) ;
    } ;

    /*!
     * @brief Abstract interface class to load and build GeoModels from files 
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

    protected:
        std::string filename_ ;
    } ;
}

#endif
