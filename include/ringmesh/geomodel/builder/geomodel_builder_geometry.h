/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <ringmesh/geomodel/builder/common.h>

#include <ringmesh/geomodel/builder/geomodel_builder_access.h>

/*!
 * @brief Builder tools to edit and build GeoModel geometry
 * (into meshes)
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshBuilder );

    ALIAS_3D( GeoModel );
    ALIAS_3D( GeoModelBuilder );
    ALIAS_3D( VolumeMeshBuilder );
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderGeometryBase
    {
        ringmesh_disable_copy_and_move( GeoModelBuilderGeometryBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;

    public:
        virtual ~GeoModelBuilderGeometryBase() = default;

        void clear_geomodel_mesh();
        /*!
         * @brief Transfer general mesh information from one mesh
         * data structure to another one
         * @param[in] id the GeoModelMeshEntity id to operate on
         * @param[in] type the new mesh data structure type
         */
        void change_mesh_data_structure(
            const gmme_id& id, const MeshType type )
        {
            GeoModelMeshEntityAccess< DIMENSION > gmme_access(
                geomodel_access_.modifiable_mesh_entity( id ) );
            gmme_access.change_mesh_data_structure( type );
        }

        /*!
         * @brief Create a PointMeshBuilder for a given corner
         * @param[in] corner_id the corner index
         * @return The created Mesh0DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh0DBuilder_var.
         */
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > >
            create_corner_builder( index_t corner_id );

        /*!
         * @brief Create a LineMeshBuilder for a given line
         * @param[in] line_id the line index
         * @return The created LineMeshBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer LineMeshBuilder_var.
         */
        std::unique_ptr< LineMeshBuilder< DIMENSION > > create_line_builder(
            index_t line_id );

        /*!
         * @brief Create a Mesh2DBuilder for a given surface
         * @param[in] surface_id the surface index
         * @return The created Mesh2DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh2DBuilder_var.
         */
        std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >
            create_surface_builder( index_t surface_id );

        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The geomodel under construction has exaclty the same number of
         * entities
         * than the input geomodel.
         */
        virtual void copy_meshes( const GeoModel< DIMENSION >& from );
        void copy_meshes( const GeoModel< DIMENSION >& from,
            const MeshEntityType& entity_type );
        void copy_mesh(
            const GeoModel< DIMENSION >& from, const gmme_id& mesh_entity );

        void assign_mesh_to_entity(
            const MeshBase< DIMENSION >& mesh, const gmme_id& to );

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
         * to the new coordinates (i.e. if edit a Corner coordinates, it will
         * updates its Lines, Surfaces...)
         */
        void set_mesh_entity_vertex( const gmme_id& entity_id,
            index_t v,
            const vecn< DIMENSION >& point,
            bool update );

        /*!
         * @brief Sets a vertex coordinates of all GeoModelMeshEntity
         * from the GeoModelMesh vertex relationships
         * @param[in] geomodel_vertex Index in GeoModelMeshVertices of the
         * vertex giving
         * @param[in] point the coordinates to set
         */
        void set_mesh_entity_vertex(
            index_t geomodel_vertex_id, const vecn< DIMENSION >& point );

        /*!
         * @brief Adds vertices to the mesh
         * @details No update of the geomodel vertices is done
         * @param[in] entity_id Entity index
         * @param[in] points Geometric positions of the vertices to add
         * @param[in] clear If true the mesh is cleared, keeping its attributes
         */
        void set_mesh_entity_vertices( const gmme_id& entity_id,
            const std::vector< vecn< DIMENSION > >& points,
            bool clear );

        /*!
         * @brief Sets the coordinates of a given existing Corner
         * @param[in] corner_id the index of the corner in the GeoModel
         * @param[in] point the coordinates to set
         */
        void set_corner( index_t corner_id, const vecn< DIMENSION >& point );
        /*!
         * @brief Sets the mesh of a given existing Line
         * @param[in] line_id the index of the line in the GeoModel
         * @param[in] vertices the coordinates to set
         * @warning the vertices should be ordered from the first boundary
         * corner to the second one
         */
        void set_line(
            index_t line_id, const std::vector< vecn< DIMENSION > >& vertices );
        /*!
         * @brief Sets the points and polygons for a surface
         * @details If polygon_adjacencies are not given they are computed.
         *
         * @param[in] surface_id Index of the surface
         * @param[in] surface_vertices Coordinates of the vertices
         * @param[in] surface_polygons Indices in the vertices vector to build
         * polygons
         * @param[in] surface_polygon_ptr Pointer to the beginning of a polygon
         * in \p surface_polygons
         */
        void set_surface_geometry( index_t surface_id,
            const std::vector< vecn< DIMENSION > >& surface_vertices,
            const std::vector< index_t >& surface_polygons,
            const std::vector< index_t >& surface_polygon_ptr );

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */

        /*!
         * @brief Sets the geometrical position of a vertex from a geomodel
         * vertex
         * @param[in] entity_id Entity index
         * @param[in] v Index of the vertex to modify
         * @param[in] geomodel_vertex Index in GeoModelMeshVertices of the
         * vertex giving
         *                     the new position
         */
        void set_mesh_entity_vertex(
            const gmme_id& entity_id, index_t v, index_t geomodel_vertex );

        /*!
         * @brief Adds vertices to the mesh
         * @details No update of the geomodel vertices is done
         * @param[in] entity_id Entity index
         * @param[in] geomodel_vertices Geometric positions of the vertices to
         * add
         * @param[in] clear If true the mesh if cleared, keeping its attributes
         */
        void set_mesh_entity_vertices( const gmme_id& entity_id,
            const std::vector< index_t >& geomodel_vertices,
            bool clear );

        /*!
         * @brief Sets the vertex for a Corner. Store the info in the geomodel
         * vertices
         *
         * @param[in] corner_id Index of the corner
         * @param[in] geomodel_vertex_id Index of the vertex in the geomodel
         */
        void set_corner( index_t corner_id, index_t geomodel_vertex_id );

        /*!
         * @brief Sets one Line vertices. Store the info in the geomodel
         * vertices
         *
         * @param[in] id Line index
         * @param[in] unique_vertices Indices in the geomodel of the unique
         * vertices with which to build the Line
         */
        void set_line(
            index_t line_id, const std::vector< index_t >& unique_vertices );

        /*! @}
         * \name Set entity geometry using GeoModelMeshEntity vertices
         * @{
         */

        /*!
         * @brief Sets the polygons of a surface
         * @param[in] surface_id Index of the surface
         * @param[in] polygons Indices of the mesh vertices defining the
         * polygons
         * @param[in] polygon_ptr Pointer to the beginning of a polygon in \p
         * polygons
         */
        void set_surface_geometry( index_t surface_id,
            const std::vector< index_t >& polygons,
            const std::vector< index_t >& polygon_ptr );

        void set_surface_element_geometry( index_t surface_id,
            index_t polygon_id,
            const std::vector< index_t >& corners );

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
            const gmme_id& entity_id, index_t nb_vertices );

        index_t create_surface_polygon(
            index_t surface_id, const std::vector< index_t >& vertex_indices );

        /*! @}
         * \name Delete mesh element entities
         * @{
         */

        void delete_mesh_entity_mesh( const gmme_id& E_id );
        virtual void delete_mesh_entity_isolated_vertices(
            const gmme_id& E_id );
        void delete_mesh_entity_vertices(
            const gmme_id& E_id, const std::vector< bool >& to_delete );
        void delete_corner_vertex( index_t corner_id );
        void delete_line_edges( index_t line_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices );
        void delete_surface_polygons( index_t surface_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices );

        /*! @}
         * \name Misc
         * @{
         */

        void set_surface_element_adjacency( index_t surface_id,
            index_t polygon_id,
            const std::vector< index_t >& adjacents );

        /*!
         * @brief Computes and sets the adjacencies between the polygons
         * @details The adjacent polygon is given for each vertex of each
         * polygon for the edge
         * starting at this vertex.
         * If there is no neighbor inside the same Surface adjacent is set to
         * NO_ID
         *
         * @param[in] surface_id Index of the surface
         * @param[in] recompute_adjacency If true, recompute the existing
         * adjacencies
         */
        void compute_surface_adjacencies(
            index_t surface_id, bool recompute_adjacency = true );

        void cut_surfaces_by_internal_lines();

        /*!
         * @brief Cuts a Surface along a Line assuming that the edges of the
         * Line are edges of the Surface
         * @pre Surface is not already cut. Line L does not cut the Surface S
         * into 2 connected components.
         * @todo Add a test for this function.
         */
        void cut_surface_by_line( index_t surface_id, index_t line_id );

    protected:
        GeoModelBuilderGeometryBase( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

        /*!
         * @brief Duplicates the surface vertices along the fake boundary
         * (NO_ID adjacencies but shared vertices) and duplicate the vertices
         */
        void duplicate_surface_vertices_along_line(
            index_t surface_id, index_t line_id );
        /*
         * @brief Resets the adjacencies for all Surface polygons adjacent to
         * the Line
         * @return The number of disconnection done
         * @pre All the edges of the Line are edges of at least one polygon of
         * the Surface
         */
        index_t disconnect_surface_polygons_along_line_edges(
            index_t surface_id, index_t line_id );

        void update_polygon_vertex( index_t surface_id,
            const std::vector< index_t >& polygons,
            index_t old_vertex,
            index_t new_vertex );

    protected:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };

    ALIAS_2D_AND_3D( GeoModelBuilderGeometryBase );

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderGeometry final
        : public GeoModelBuilderGeometryBase< DIMENSION >
    {
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;

    protected:
        GeoModelBuilderGeometry( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGeometryBase< DIMENSION >( builder, geomodel )
        {
        }
    };

    template <>
    class geomodel_builder_api GeoModelBuilderGeometry< 3 > final
        : public GeoModelBuilderGeometryBase< 3 >
    {
        friend class GeoModelBuilderBase< 3 >;
        friend class GeoModelBuilder< 3 >;

    public:
        /*!
         * @brief Create a Mesh3DBuilder for a given region
         * @param[in] region_id the region index
         * @return The created Mesh3DBuilder
         * @warn The client code is responsible for the memory unallocation.
         * You can use the smartpointer Mesh3DBuilder_var.
         */
        std::unique_ptr< VolumeMeshBuilder3D > create_region_builder(
            index_t region_id );

        void copy_meshes( const GeoModel3D& geomodel ) override;

        void set_region_element_geometry( index_t region_id,
            index_t cell_id,
            const std::vector< index_t >& corners );

        /*!
         * @brief Set the points and tetras for a region
         *
         * @param[in] region_id Index of the regions
         * @param[in] points Coordinates of the vertices
         * @param[in] tetras Indices in the vertices vector to build tetras
         */
        void set_region_geometry( index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras );

        index_t create_region_cell( index_t region_id,
            CellType type,
            const std::vector< index_t >& vertex_indices );

        /*!
         * @brief Creates new cells in the mesh
         * @param[in] region_id Entity index
         * @param[in] type Type of cell
         * @param[in] nb_cells Number of cells to creates
         * @return the index of the first created cell
         */
        index_t create_region_cells(
            index_t region_id, CellType type, index_t nb_cells );

        void delete_region_cells( index_t region_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices );

        /*!
         * @brief Computes and sets the adjacencies between the cells
         * @details The adjacent cell is given for each facet of each cell
         * If there is no neighbor inside the same Region adjacent is set to
         * NO_ID
         *
         * @param[in] region_id Index of the region
         * @param[in] recompute_adjacency If true, recompute the existing
         * adjacencies
         */
        void compute_region_adjacencies(
            index_t region_id, bool recompute_adjacency = true );

        void cut_regions_by_internal_surfaces();

        void cut_region_by_surface( index_t region_id, index_t surface_id );

        void delete_mesh_entity_isolated_vertices(
            const gmme_id& E_id ) override;

    protected:
        GeoModelBuilderGeometry(
            GeoModelBuilder3D& builder, GeoModel3D& geomodel )
            : GeoModelBuilderGeometryBase< 3 >( builder, geomodel )
        {
        }

        void assign_region_tet_mesh(
            index_t region_id, const std::vector< index_t >& tet_vertices );

        void update_cell_vertex( index_t region_id,
            const std::vector< index_t >& cells,
            index_t old_vertex,
            index_t new_vertex );

        void duplicate_region_vertices_along_surface(
            index_t region_id, index_t surface_id );

        index_t disconnect_region_cells_along_surface_polygons(
            index_t region_id, index_t surface_id );
    };
} // namespace RINGMesh
