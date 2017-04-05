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

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @brief Builder tools to edit and build GeoModel geometry
 * (into meshes)
 * @author Pierre Anquez
 */

namespace RINGMesh {
    class GeoModelBuilder ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderGeometry {
    ringmesh_disable_copy( GeoModelBuilderGeometry ) ;
        friend class GeoModelBuilder ;

    public:
        void recompute_geomodel_mesh() ;
        /*!
         * @brief Transfer general mesh information from one mesh
         * data structure to another one
         * @param[in] id the GeoModelMeshEntity id to operate on
         * @param[in] type the new mesh data structure type
         */
        void change_mesh_data_structure( const gmme_t& id, const MeshType type )
        {
            GeoModelMeshEntityAccess gmme_access(
                geomodel_access_.modifiable_mesh_entity( id ) ) ;
            gmme_access.change_mesh_data_structure( type ) ;
        }

        /*!
         * @brief Create a Mesh0DBuilder for a given corner
         * @param[in] corner_id the corner index
         * @return The created Mesh0DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh0DBuilder_var.
         */
        Mesh0DBuilder* create_corner_builder( index_t corner_id )
        {
            gmme_t id( Corner::type_name_static(), corner_id ) ;
            GeoModelMeshEntity& corner = geomodel_access_.modifiable_mesh_entity(
                id ) ;
            GeoModelMeshEntityAccess corner_access( corner ) ;
            Mesh0D& corner_mesh =
                dynamic_cast< Mesh0D& >( *corner_access.modifiable_mesh() ) ;
            return Mesh0DBuilder::create_builder( corner_mesh ) ;
        }

        /*!
         * @brief Create a Mesh1DBuilder for a given line
         * @param[in] line_id the line index
         * @return The created Mesh1DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh1DBuilder_var.
         */
        Mesh1DBuilder* create_line_builder( index_t line_id )
        {
            gmme_t id( Line::type_name_static(), line_id ) ;
            GeoModelMeshEntity& line = geomodel_access_.modifiable_mesh_entity(
                id ) ;
            GeoModelMeshEntityAccess line_access( line ) ;
            Mesh1D& line_mesh =
                dynamic_cast< Mesh1D& >( *line_access.modifiable_mesh() ) ;
            return Mesh1DBuilder::create_builder( line_mesh ) ;
        }

        /*!
         * @brief Create a Mesh2DBuilder for a given surface
         * @param[in] surface_id the surface index
         * @return The created Mesh2DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh2DBuilder_var.
         */
        Mesh2DBuilder* create_surface_builder( index_t surface_id )
        {
            gmme_t id( Surface::type_name_static(), surface_id ) ;
            GeoModelMeshEntity& surface = geomodel_access_.modifiable_mesh_entity(
                id ) ;
            GeoModelMeshEntityAccess surface_access( surface ) ;
            Mesh2D& surface_mesh =
                dynamic_cast< Mesh2D& >( *surface_access.modifiable_mesh() ) ;
            return Mesh2DBuilder::create_builder( surface_mesh ) ;
        }

        /*!
         * @brief Create a Mesh3DBuilder for a given region
         * @param[in] region_id the region index
         * @return The created Mesh3DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh3DBuilder_var.
         */
        Mesh3DBuilder* create_region_builder( index_t region_id )
        {
            gmme_t id( Region::type_name_static(), region_id ) ;
            GeoModelMeshEntity& region = geomodel_access_.modifiable_mesh_entity(
                id ) ;
            GeoModelMeshEntityAccess region_access( region ) ;
            Mesh3D& region_mesh =
                dynamic_cast< Mesh3D& >( *region_access.modifiable_mesh() ) ;
            return Mesh3DBuilder::create_builder( region_mesh ) ;
        }

        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The geomodel under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;
        void copy_meshes( const GeoModel& from, const MeshEntityType& entity_type ) ;
        void copy_mesh( const GeoModel& from, const gmme_t& mesh_entity ) ;

        void assign_mesh_to_entity( const MeshBase& mesh, const gmme_t& to ) ;

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
         * to the new coordinates (i.e. if edit a Corner coordinates, it will updates
         * its Lines, Surfaces...)
         */
        void set_mesh_entity_vertex(
            const gmme_t& entity_id,
            index_t v,
            const vec3& point,
            bool update ) ;

        /*!
         * @brief Adds vertices to the mesh
         * @details No update of the geomodel vertices is done
         * @param[in] id Entity index
         * @param[in] points Geometric positions of the vertices to add
         * @param[in] clear If true the mesh is cleared, keeping its attributes
         */
        void set_mesh_entity_vertices(
            const gmme_t& entity_id,
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
         * @brief Sets the points and facets for a surface
         * @details If facet_adjacencies are not given they are computed.
         *
         * @param[in] surface_id Index of the surface
         * @param[in] surface_vertices Coordinates of the vertices
         * @param[in] surface_facets Indices in the vertices vector to build facets
         * @param[in] surface_facet_ptr Pointer to the beginning of a facet in facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;
        /*!
         * @brief Set the points and tetras for a region
         *
         * @param[in] region_id Index of the regions
         * @param[in] points Coordinates of the vertices
         * @param[in] tetras Indices in the vertices vector to build tetras
         */
        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */

        /*!
         * @brief Sets the geometrical position of a vertex from a geomodel vertex
         * @param[in] entity_id Entity index
         * @param[in] v Index of the vertex to modify
         * @param[in] geomodel_vertex Index in GeoModelMeshVertices of the vertex giving
         *                     the new position
         */
        void set_mesh_entity_vertex(
            const gmme_t& entity_id,
            index_t v,
            index_t geomodel_vertex ) ;

        /*!
         * @brief Adds vertices to the mesh
         * @details No update of the geomodel vertices is done
         * @param[in] entity_id Entity index
         * @param[in] geomodel_vertices Geometric positions of the vertices to add
         * @param[in] clear If true the mesh if cleared, keeping its attributes
         */
        void set_mesh_entity_vertices(
            const gmme_t& entity_id,
            const std::vector< index_t >& geomodel_vertices,
            bool clear ) ;

        /*!
         * @brief Sets the vertex for a Corner. Store the info in the geomodel vertices
         *
         * @param[in] corner_id Index of the corner
         * @param[in] geomodel_vertex_id Index of the vertex in the geomodel
         */
        void set_corner( index_t corner_id, index_t geomodel_vertex_id ) ;

        /*!
         * @brief Sets one Line vertices. Store the info in the geomodel vertices
         *
         * @param[in] id Line index
         * @param[in] unique_vertices Indices in the geomodel of the unique vertices with which to build the Line
         */
        void set_line(
            index_t line_id,
            const std::vector< index_t >& unique_vertices ) ;

        /*!
         * @brief Sets the vertices and facets for a surface
         *
         * @param[in] surface_id Index of the surface
         * @param[in] geomodel_vertex_ids Indices of unique vertices in the GeoModel
         * @param[in] facets Indices in the vertices vector to build facets
         * @param[in] facet_ptr Pointer to the beginning of a facet in facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& geomodel_vertex_ids,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        /*!
         * @brief Sets the facets of a surface
         * @param[in] surface_id Index of the surface
         * @param[in] facets Indices of the geomodel vertices defining the facets
         * @param[in] facet_ptr Pointer to the beginning of a facet in facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& facets,
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

        /*!
         * @brief Creates new vertices to the mesh
         * @param[in] entity_id Entity index
         * @param[in] nb_vertices Number of vertices to create
         * @return the first vertex index created
         */
        index_t create_mesh_entity_vertices(
            const gmme_t& entity_id,
            index_t nb_vertices ) ;

        index_t create_surface_facet(
            index_t surface_id,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cell(
            index_t region_id,
            GEO::MeshCellType type,
            const std::vector< index_t >& vertex_indices ) ;

        /*!
         * @brief Creates new cells in the mesh
         * @param[in] region_id Entity index
         * @param[in] type Type of cell
         * @param[in] nb_cells Number of cells to creates
         * @return the index of the first created cell
         */
        index_t create_region_cells(
            index_t region_id,
            GEO::MeshCellType type,
            index_t nb_cells ) ;

        /*! @}
         * \name Delete mesh element entities
         * @{
         */

        void delete_mesh_entity_mesh( const gmme_t& E_id ) ;
        void delete_mesh_entity_isolated_vertices( const gmme_t& E_id ) ;
        void delete_mesh_entity_vertices(
            const gmme_t& E_id,
            const std::vector< bool >& to_delete ) ;
        void delete_corner_vertex( index_t corner_id ) ;
        void delete_line_edges(
            index_t line_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_surface_facets(
            index_t surface_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_region_cells(
            index_t region_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;

        /*! @}
         * \name Misc
         * @{
         */

        void set_surface_facet_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& facets_ids,
            const std::vector< index_t >& edges_ids,
            const std::vector< index_t >& adjacent_triangles ) ;

        /*!
         * @brief Computes and sets the adjacencies between the facets
         * @details The adjacent facet is given for each vertex of each facet for the edge
         * starting at this vertex.
         * If there is no neighbor inside the same Surface adjacent is set to NO_ID
         *
         * @param[in] surface_id Index of the surface
         * @param[in] recompute_adjacency If true, recompute the existing adjacencies
         */
        void compute_surface_adjacencies(
            index_t surface_id,
            bool recompute_adjacency = true ) ;

        /*!
         * @brief Computes and sets the adjacencies between the cells
         * @details The adjacent cell is given for each facet of each cell
         * If there is no neighbor inside the same Region adjacent is set to NO_ID
         *
         * @param[in] region_id Index of the region
         * @param[in] recompute_adjacency If true, recompute the existing adjacencies
         */
        void compute_region_adjacencies(
            index_t region_id,
            bool recompute_adjacency = true ) ;

        void cut_surfaces_by_internal_lines() ;

        void cut_regions_by_internal_surfaces() ;

        /*!
         * @brief Cuts a Surface along a Line assuming that the edges of the Line are edges of the Surface
         * @pre Surface is not already cut. Line L does not cut the Surface S into 2 connected components.
         * @todo Add a test for this function.
         */
        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;

        void cut_region_by_surface( index_t region_id, index_t surface_id ) ;

    protected:
        GeoModelBuilderGeometry( GeoModelBuilder& builder, GeoModel& geomodel ) ;

    private:
        void assign_surface_mesh_facets(
            index_t surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) ;

        /*!
         * @brief Duplicates the surface vertices along the fake boundary
         * (NO_ID adjacencies but shared vertices) and duplicate the vertices
         */
        void duplicate_surface_vertices_along_line(
            index_t surface_id,
            index_t line_id ) ;
        void duplicate_region_vertices_along_surface(
            index_t region_id,
            index_t surface_id ) ;
        /*
         * @brief Resets the adjacencies for all Surface facets adjacent to the Line
         * @return The number of disconnection done
         * @pre All the edges of the Line are edges of at least one facet of the Surface
         */
        index_t disconnect_surface_facets_along_line_edges(
            index_t surface_id,
            index_t line_id ) ;
        index_t disconnect_region_cells_along_surface_facets(
            index_t region_id,
            index_t surface_id ) ;

        void update_facet_vertex(
            index_t surface_id,
            const std::vector< index_t >& facets,
            index_t old_vertex,
            index_t new_vertex ) ;

        void update_cell_vertex(
            index_t region_id,
            const std::vector< index_t >& cells,
            index_t old_vertex,
            index_t new_vertex ) ;

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

}
