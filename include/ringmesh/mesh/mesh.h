/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <algorithm>
#include <memory>
#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/nn_search.h>
#include <ringmesh/mesh/common.h>
#include <ringmesh/mesh/mesh_base.h>

#include <ringmesh/mesh/mesh_aabb.h>



namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBaseBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );

    struct EdgeLocalVertex;
    struct ElementLocalVertex;
    struct PolygonLocalEdge;
    struct CellLocalFacet;
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * class base class for encapsulating Mesh structure
     * @brief encapsulate adimensional mesh functionalities in order to provide
     * an API
     * on which we base the RINGMesh algorithms
     * @note For now, we encapsulate the GEO::Mesh class.
     */
    /*!
     * class for encapsulating mesh composed of points
     */
    template < index_t DIMENSION >
    class PointSetMesh : public MeshBase< DIMENSION >
    {
        friend class PointSetMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< PointSetMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );
        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;
        bool is_mesh_valid() const override
        {
            return true;
        }

    protected:
        PointSetMesh() = default;
    };
    ALIAS_2D_AND_3D( PointSetMesh );

    template < index_t DIMENSION >
    using PointSetMeshFactory = Factory< MeshType, PointSetMesh< DIMENSION > >;
    ALIAS_2D_AND_3D( PointSetMeshFactory );

    /*!
     * class for encapsulating line mesh (composed of edges)
     */
    template < index_t DIMENSION >
    class LineMesh : public MeshBase< DIMENSION >
    {
        friend class LineMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< LineMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );

        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_local_vertex index of the edge and of the
         * local index of the vertex, in {0,1}
         * @return the global index of vertex in \p edge_local_vertex.
         */
        virtual index_t edge_vertex(
            const ElementLocalVertex& edge_local_vertex ) const = 0;

        /*!
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        virtual index_t nb_edges() const = 0;

        /*!
         * @brief Gets the length of the edge \param edge_id
         */
        double edge_length( index_t edge_id ) const;

        vecn< DIMENSION > edge_barycenter( index_t edge_id ) const;

        /*!
         * @brief return the NNSearch at edges
         * @warning the NNSearch is destroyed when calling the
         * Mesh::polygons_aabb()
         * and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& edge_nn_search() const;

        /*!
         * @brief Creates an AABB tree for a Mesh edges
         */
        const LineAABBTree< DIMENSION >& edge_aabb() const;

        virtual GEO::AttributesManager& edge_attribute_manager() const = 0;

        bool is_mesh_valid() const override;

        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;

    protected:
        LineMesh() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > edge_nn_search_{};
        mutable std::unique_ptr< LineAABBTree< DIMENSION > > edge_aabb_{};
    };
    ALIAS_2D_AND_3D( LineMesh );

    template < index_t DIMENSION >
    using LineMeshFactory = Factory< MeshType, LineMesh< DIMENSION > >;
    ALIAS_2D_AND_3D( LineMeshFactory );

    /*!
     * class for encapsulating surface mesh component
     */
    template < index_t DIMENSION >
    class SurfaceMeshBase : public MeshBase< DIMENSION >
    {
        friend class SurfaceMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< SurfaceMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );

        /*!
         * @brief Gets the vertex index by polygon index and local vertex index.
         * @param[in] polygon_local_vertex the polygon index and
         * the local edge index in the polygon.
         */
        virtual index_t polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const = 0;

        /*!
         * @brief Gets the number of all polygons in the whole Mesh.
         */
        virtual index_t nb_polygons() const = 0;

        /*!
         * @brief Gets the number of vertices in the polygon \param polygon_id.
         * @param[in] polygon_id polygon index
         */
        virtual index_t nb_polygon_vertices( index_t polygon_id ) const = 0;

        /*!
         * @brief Gets the number of edges in the polygon \param polygon_id.
         * @param[in] polygon_id polygon index
         */
        index_t nb_polygon_edges( index_t polygon_id ) const
        {
            return nb_polygon_vertices( polygon_id );
        }

        /*!
         * @brief Gets the next vertex index in the polygon vertex
         * \param polygon_local_vertex.
         * @param[in] polygon_local_vertex polygon index and the current
         * local vertex index.
         */
        ElementLocalVertex next_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const;

        /*!
         * @brief Get the next edge on the border
         * @warning the edge index is in fact the index of the vertex where the
         * edge starts.
         * @details The returned border edge is the next in the way of polygon
         * edges
         * orientation.
         * @param[in] polygon_local_edge input polygon index and its local
         * edge index in the polygon
         * @return the next polygon index
         * and the next edge index in the polygon
         *
         * @pre the given polygon edge must be on border
         */
        PolygonLocalEdge next_on_border(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*!
         * @brief Gets the previous vertex index in the polygon \param
         * polygon_id.
         * @param[in] polygon_local_vertex polygon index and its current
         * local vertex index
         */
        ElementLocalVertex prev_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const;

        /*!
         * @brief Get the previous edge on the border
         * @details The returned border edge is the previous in the way of
         * polygon edges
         * orientation.
         * @param[in] p Input polygon index
         * @param[in] e Edge index in the polygon
         * @return the previous polygon index
         * and the previous edge index in the polygon (tuple).
         *
         * @pre the surface must be correctly oriented and
         * the given polygon edge must be on border
         * @warning the edge index is in fact the index of the vertex where the
         * edge starts.
         */
        PolygonLocalEdge prev_on_border(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*!
         * @brief Get the vertex index in a polygon @param polygon_index from
         * its
         * global index in the SurfaceMesh @param vertex_id
         * @return NO_ID or index of the vertex in the polygon
         */
        index_t vertex_index_in_polygon(
            index_t polygon_index, index_t vertex_id ) const;

        /*!
         * @brief Compute closest vertex in a polygon to a point
         * @param[in] polygon_index Polygon index
         * @param[in] query_point Coordinates of the point to which distance is
         * measured
         * @return Index of the vertex of @param polygon_index closest to @param
         * query_point
         */
        index_t closest_vertex_in_polygon(
            index_t polygon_index, const vecn< DIMENSION >& query_point ) const;

        /*!
         * @brief Get the first polygon of the surface that has an edge linking
         * the two vertices (ids in the surface)
         *
         * @param[in] in0 Index of the first vertex in the surface
         * @param[in] in1 Index of the second vertex in the surface
         * @return NO_ID or the index of the polygon
         */
        index_t polygon_from_vertex_ids( index_t in0, index_t in1 ) const;

        /*!
         * @brief Determines the polygons around a vertex
         * @param[in] vertex_id Index of the vertex in the surface
         * @param[in] border_only If true only polygons on the border are
         * considered
         * @param[in] first_polygon (Optional) Index of one polygon containing
         * the vertex @param P
         * @return Indices of the polygons containing @param P
         * @note If a polygon containing the vertex is given, polygons around
         * this
         * vertex is search by propagation. Else, a first polygon is found by
         * brute
         * force algorithm, and then the other by propagation
         * @todo Try to use a AABB tree to remove @param first_polygon. [PA]
         */
        std::vector< index_t > polygons_around_vertex(
            index_t vertex_id, bool border_only, index_t first_polygon ) const;

        /*!
         * @brief Gets an adjacent polygon index by polygon index and local edge
         * index.
         * @param[in] polygon_id the polygon index.
         * @param[in] edge_id the local edge index in \param polygon_id.
         * @return the global polygon index adjacent to the \param edge_id of
         * the polygon \param polygon_id.
         * @precondition  \param edge_id < number of edge of the polygon \param
         * polygon_id .
         */
        virtual index_t polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge ) const = 0;

        virtual GEO::AttributesManager& polygon_attribute_manager() const = 0;

        /*!
         * @brief Tests whether all the polygons are triangles. when all the
         * polygons are triangles, storage and access is optimized.
         * @return True if all polygons are triangles and False otherwise.
         */
        virtual bool polygons_are_simplicies() const = 0;

        /*!
         * return true if the polygon \param polygon_id is a triangle
         */
        bool is_triangle( index_t polygon_id ) const
        {
            return nb_polygon_vertices( polygon_id ) == 3;
        }

        PolygonType polygone_type( index_t polygon_id ) const
        {
            if( is_triangle( polygon_id ) )
            {
                return PolygonType::TRIANGLE;
            }
            if( nb_polygon_vertices( polygon_id ) == 4 )
            {
                return PolygonType::QUAD;
            }
            return PolygonType::UNDEFINED;
        }

        /*!
         * Is the edge starting with the given vertex of the polygon on a border
         * of the Surface?
         */
        bool is_edge_on_border(
            const PolygonLocalEdge& polygon_local_edge ) const
        {
            return polygon_adjacent( polygon_local_edge ) == NO_ID;
        }

        /*!
         * Is one of the edges of the polygon on the border of the surface?
         */
        bool is_polygon_on_border( index_t polygon_index ) const;

        /*!
         * @brief Gets the length of the edge starting at a given vertex
         * @param[in] polygon_local_edge index of the polygon and its
         * local edge starting vertex index
         */
        double polygon_edge_length(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*!
         * @brief Gets the barycenter of the edge starting at a given vertex
         * @param[in] polygon_local_edge index of the polygon and
         * the local edge starting vertex index
         */
        vecn< DIMENSION > polygon_edge_barycenter(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*!
         * @brief Gets the vertex index on the polygon edge
         * @param[in] polygon_local_edge index of the polygon and
         * the local index of the edge in the polygon
         * @param[in] vertex_id index of the local vertex in the edge \param
         * edge_id (0 or 1)
         * @return the vertex index
         */
        index_t polygon_edge_vertex( const PolygonLocalEdge& polygon_local_edge,
            index_t vertex_id ) const;

        /*!
         * Computes the Mesh polygon barycenter
         * @param[in] polygon_id the polygon index
         * @return the polygon center
         */
        vecn< DIMENSION > polygon_barycenter( index_t polygon_id ) const;

        /*!
         * Computes the Mesh polygon area
         * @param[in] polygon_id the polygon index
         * @return the polygon area
         */
        virtual double polygon_area( index_t polygon_id ) const = 0;

        /*!
         * @brief return the NNSearch at polygons
         */
        const NNSearch< DIMENSION >& polygon_nn_search() const
        {
            if( !nn_search_ )
            {
                std::vector< vecn< DIMENSION > > polygon_centers(
                    nb_polygons() );
                for( auto p : range( nb_polygons() ) )
                {
                    polygon_centers[p] = polygon_barycenter( p );
                }
                nn_search_.reset(
                    new NNSearch< DIMENSION >( polygon_centers, true ) );
            }
            return *nn_search_.get();
        }
        /*!
         * @brief Creates an AABB tree for a Mesh polygons
         */
        const SurfaceAABBTree< DIMENSION >& polygon_aabb() const
        {
            if( !polygon_aabb_ )
            {
                polygon_aabb_.reset(
                    new SurfaceAABBTree< DIMENSION >( *this ) );
            }
            return *polygon_aabb_;
        }

        bool is_mesh_valid() const override;

        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;

    protected:
        SurfaceMeshBase() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > nn_search_{};
        mutable std::unique_ptr< SurfaceAABBTree< DIMENSION > > polygon_aabb_{};
    };
    ALIAS_2D_AND_3D( SurfaceMeshBase );

    template < index_t DIMENSION >
    class SurfaceMesh : public SurfaceMeshBase< DIMENSION >
    {
    };

    template < index_t DIMENSION >
    using SurfaceMeshFactory = Factory< MeshType, SurfaceMesh< DIMENSION > >;
    ALIAS_2D_AND_3D( SurfaceMeshFactory );

    template <>
    class mesh_api SurfaceMesh< 3 > : public SurfaceMeshBase< 3 >
    {
    public:
        /*!
         * Computes the Mesh polygon area
         * @param[in] polygon_id the polygon index
         * @return the polygon area
         */
        double polygon_area( index_t polygon_id ) const;

        /*!
         * Computes the Mesh polygon normal
         * @param[in] polygon_id the polygon index
         * @return the polygon normal
         */
        vec3 polygon_normal( index_t polygon_id ) const;

        /*!
         * @brief Computes the normal of the Mesh2D at the vertex location
         * it computes the average value of polygon normal neighbors
         * @param[in] vertex_id the vertex index
         * @param[in] p0 index of a polygon that contain the vertex \param
         * vertex_id
         * @return the normal at the given vertex
         */
        vec3 normal_at_vertex( index_t vertex_id, index_t p0 = NO_ID ) const;
    };

    template <>
    class mesh_api SurfaceMesh< 2 > : public SurfaceMeshBase< 2 >
    {
    public:
        /*!
         * Computes the Mesh polygon area
         * @param[in] polygon_id the polygon index
         * @return the polygon area
         */
        double polygon_area( index_t polygon_id ) const override;
    };

    ALIAS_2D_AND_3D( SurfaceMesh );

    /*!
     * class for encapsulating volume mesh component
     */
    template < index_t DIMENSION >
    class VolumeMesh : public MeshBase< DIMENSION >
    {
        ringmesh_template_assert_3d( DIMENSION );
        friend class VolumeMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< VolumeMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );

        /*!
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_vertex(
            const ElementLocalVertex& cell_local_vertex ) const = 0;

        /*!
         * @brief Gets a vertex index by cell and local edge and local vertex
         * index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_edge_vertex(
            index_t cell_id, index_t edge_id, index_t vertex_id ) const = 0;

        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param
         * facet_id
         * and facet_id number of facet in th cell \param cell_id
         */
        virtual index_t cell_facet_vertex(
            const CellLocalFacet& cell_local_facet,
            index_t vertex_id ) const = 0;

        /*!
         * @brief Gets a facet index by cell and local facet index.
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @return the global facet index.
         */
        virtual index_t cell_facet(
            const CellLocalFacet& cell_local_facet ) const = 0;

        /*!
         * Computes the Mesh cell edge length
         * @param[in] cell_id the facet index
         * @param[in] edge_id the edge index
         * @return the cell edge length
         */
        double cell_edge_length( index_t cell_id, index_t edge_id ) const;

        /*!
         * Computes the Mesh cell edge barycenter
         * @param[in] cell_id the facet index
         * @param[in] edge_id the edge index
         * @return the cell edge center
         */
        vecn< DIMENSION > cell_edge_barycenter(
            index_t cell_id, index_t edge_id ) const;

        /*!
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_facets( index_t cell_id ) const = 0;
        /*!
         * @brief Gets the total number of facet in a all cells
         */
        virtual index_t nb_cell_facets() const = 0;

        /*!
         * @brief Gets the number of edges in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_edges( index_t cell_id ) const = 0;

        /*!
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @return the number of vertices in the facet \param facet_id in the
         * cell \param cell_id
         */
        virtual index_t nb_cell_facet_vertices(
            const CellLocalFacet& cell_local_facet ) const = 0;

        /*!
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        virtual index_t nb_cell_vertices( index_t cell_id ) const = 0;

        /*!
         * @brief Gets the number of cells in the Mesh.
         */
        virtual index_t nb_cells() const = 0;

        virtual index_t cell_begin( index_t cell_id ) const = 0;

        virtual index_t cell_end( index_t cell_id ) const = 0;

        /*!
         * @return the index of the adjacent cell of \param cell_local_facet
         */
        virtual index_t cell_adjacent(
            const CellLocalFacet& cell_local_facet ) const = 0;

        virtual GEO::AttributesManager& cell_attribute_manager() const = 0;

        virtual GEO::AttributesManager&
            cell_facet_attribute_manager() const = 0;

        /*!
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        virtual CellType cell_type( index_t cell_id ) const = 0;

        /*!
         * @brief Tests whether all the cells are tetrahedra.
         * When all the cells are tetrahedra, storage and access is optimized.
         * @return True if all cells are tetrahedra and False otherwise.
         */
        virtual bool cells_are_simplicies() const = 0;

        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] cell_local_facet the cell index and
         * the local facet index in the cell
         * @return the cell facet center
         */
        vecn< DIMENSION > cell_facet_barycenter(
            const CellLocalFacet& cell_local_facet ) const;

        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        vecn< DIMENSION > cell_barycenter( index_t cell_id ) const;

        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_local_facet the cell index and
         * the local facet index in the cell
         * @return the cell facet normal
         */
        vecn< DIMENSION > cell_facet_normal(
            const CellLocalFacet& cell_local_facet ) const;

        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        virtual double cell_volume( index_t cell_id ) const = 0;

        std::vector< index_t > cells_around_vertex(
            index_t vertex_id, index_t cell_hint ) const;

        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const;

        bool find_cell_from_colocated_vertex_within_distance_if_any(
            const vecn< DIMENSION >& vertex_vec,
            double distance,
            index_t& cell_id,
            index_t& cell_vertex_id ) const;

        /*!
         * @brief return the NNSearch at cell facets
         * @warning the NNSearch is destroyed when calling the
         * Mesh::facets_aabb()
         *  and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& cell_facet_nn_search() const;

        /*!
         * @brief return the NNSearch at cells
         */
        const NNSearch< DIMENSION >& cell_nn_search() const
        {
            if( !cell_nn_search_ )
            {
                std::vector< vecn< DIMENSION > > cell_centers( nb_cells() );
                for( auto c : range( nb_cells() ) )
                {
                    cell_centers[c] = cell_barycenter( c );
                }
                cell_nn_search_.reset(
                    new NNSearch< DIMENSION >( cell_centers, true ) );
            }
            return *cell_nn_search_.get();
        }
        /*!
         * @brief Creates an AABB tree for a Mesh cells
         */
        const VolumeAABBTree< DIMENSION >& cell_aabb() const
        {
            if( !cell_aabb_ )
            {
                cell_aabb_.reset( new VolumeAABBTree< DIMENSION >( *this ) );
            }
            return *cell_aabb_.get();
        }

        bool is_mesh_valid() const override;
        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;

    protected:
        VolumeMesh() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > >
            cell_facet_nn_search_{};
        mutable std::unique_ptr< NNSearch< DIMENSION > > cell_nn_search_{};
        mutable std::unique_ptr< VolumeAABBTree< DIMENSION > > cell_aabb_{};
    };

    using VolumeMesh3D = VolumeMesh< 3 >;

    template < index_t DIMENSION >
    using VolumeMeshFactory = Factory< MeshType, VolumeMesh< DIMENSION > >;
    using VolumeMeshFactory3D = VolumeMeshFactory< 3 >;

    /*!
     * class composed of meshes from all the dimensions
     */
    template < index_t DIMENSION >
    class MeshSetBase
    {
        ringmesh_disable_copy_and_move( MeshSetBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        void create_point_set_mesh( MeshType type );
        void create_line_mesh( MeshType type );
        void create_well_mesh( MeshType type );
        void create_surface_mesh( MeshType type );

    protected:
        MeshSetBase();
        virtual ~MeshSetBase() = default;

    public:
        std::unique_ptr< PointSetMesh< DIMENSION > > point_set_mesh{};
        std::unique_ptr< LineMesh< DIMENSION > > well_mesh{};
        std::unique_ptr< LineMesh< DIMENSION > > line_mesh{};
        std::unique_ptr< SurfaceMesh< DIMENSION > > surface_mesh{};
    };

    template < index_t DIMENSION >
    class mesh_api MeshSet : public MeshSetBase< DIMENSION >
    {
    public:
        MeshSet() = default;
    };

    template <>
    class mesh_api MeshSet< 3 > : public MeshSetBase< 3 >
    {
    public:
        MeshSet();

        void create_volume_mesh( MeshType type );

    public:
        std::unique_ptr< VolumeMesh3D > volume_mesh{};
    };
} // namespace RINGMesh
