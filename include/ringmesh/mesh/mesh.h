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

#include <memory>
#include <stack>

#include <geogram/basic/attributes.h>

#include <geogram/mesh/mesh.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/nn_search.h>
#include <ringmesh/mesh/aabb.h>

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModel;
    template< index_t DIMENSION > class MeshBaseBuilder;
    template< index_t DIMENSION > class PointSetMeshBuilder;
    template< index_t DIMENSION > class LineMeshBuilder;
    template< index_t DIMENSION > class SurfaceMeshBuilder;
    template< index_t DIMENSION > class VolumeMeshBuilder;
    template< index_t DIMENSION > class SurfaceMesh;
    struct EdgeLocalVertex;
    struct PolygonLocalEdge;
    struct CellLocalFacet;
}

namespace RINGMesh {

    using MeshType = std::string;

    struct ElementLocalVertex {
        ElementLocalVertex() = default;
        ElementLocalVertex( index_t element_id, index_t local_vertex_id )
            : element_id_( element_id ), local_vertex_id_( local_vertex_id )
        {
        }
        ElementLocalVertex( EdgeLocalVertex edge_local_vertex );
        ElementLocalVertex( PolygonLocalEdge polygon_local_edge );
        ElementLocalVertex( CellLocalFacet cell_local_facet );
        index_t element_id_{ NO_ID };
        index_t local_vertex_id_{ NO_ID };
    };

    struct EdgeLocalVertex {
        EdgeLocalVertex() = default;
        EdgeLocalVertex( index_t edge_id, index_t local_vertex_id )
            : edge_id_( edge_id ), local_vertex_id_( local_vertex_id )
        {
        }
        EdgeLocalVertex( ElementLocalVertex edge_local_vertex )
            :
                edge_id_( std::move( edge_local_vertex.element_id_ ) ),
                local_vertex_id_( std::move( edge_local_vertex.local_vertex_id_ ) )
        {
        }
        index_t edge_id_{ NO_ID };
        index_t local_vertex_id_{ NO_ID };
    };

    struct PolygonLocalEdge {
        PolygonLocalEdge() = default;
        PolygonLocalEdge( index_t polygon_id, index_t local_edge_id )
            : polygon_id_( polygon_id ), local_edge_id_( local_edge_id )
        {
        }
        PolygonLocalEdge( ElementLocalVertex polygon_local_vertex )
            :
                polygon_id_( std::move( polygon_local_vertex.element_id_ ) ),
                local_edge_id_( std::move( polygon_local_vertex.local_vertex_id_ ) )
        {
        }
        index_t polygon_id_{ NO_ID };
        index_t local_edge_id_{ NO_ID };
    };

    struct CellLocalFacet {
        CellLocalFacet() = default;
        CellLocalFacet( index_t cell_id, index_t local_facet_id )
            : cell_id_( cell_id ), local_facet_id_( local_facet_id )
        {
        }
        index_t cell_id_{ NO_ID };
        index_t local_facet_id_{ NO_ID };
    };

    /*!
     * class base class for encapsulating Mesh structure
     * @brief encapsulate adimensional mesh functionalities in order to provide an API
     * on which we base the RINGMesh algorithms
     * @note For now, we encapsulate the GEO::Mesh class.
     */
    template< index_t DIMENSION >
    class MeshBase: public GEO::Counted {
    ringmesh_disable_copy( MeshBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class MeshBaseBuilder< DIMENSION > ;

    public:
        virtual ~MeshBase() = default;

        virtual void save_mesh( const std::string& filename ) const = 0;

        //TODO maybe reimplement the function with a RINGMesh::Mesh??
        virtual void print_mesh_bounded_attributes() const = 0;
        /*!
         * \name Vertex methods
         * @{
         */
        /*!
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return const reference to the point that corresponds to the vertex.
         */
        virtual const vecn< DIMENSION >& vertex( index_t v_id ) const = 0;
        /*
         * @brief Gets the number of vertices in the Mesh.
         */
        virtual index_t nb_vertices() const = 0;

        virtual GEO::AttributesManager& vertex_attribute_manager() const = 0;

        /*!
         * @brief return the NNSearch at vertices
         * @warning the NNSearch is destroyed when calling the Mesh::polygons_aabb()
         * and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& vertex_nn_search() const
        {
            if( !vertex_nn_search_ ) {
                std::vector< vecn< DIMENSION > > vec_vertices( nb_vertices() );
                for( index_t v : range( nb_vertices() ) ) {
                    vec_vertices[v] = vertex( v );
                }
                vertex_nn_search_.reset(
                    new NNSearch< DIMENSION >( vec_vertices, true ) );
            }
            return *vertex_nn_search_.get();
        }

        virtual MeshType type_name() const = 0;

        virtual std::string default_extension() const = 0;

        /*!
         * @}
         */
    protected:
        MeshBase() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > vertex_nn_search_;
    };

    CLASS_DIMENSION_ALIASES( MeshBase );

    /*!
     * class for encapsulating mesh composed of points
     */
    template< index_t DIMENSION >
    class PointSetMesh: public MeshBase< DIMENSION > {
    ringmesh_disable_copy( PointSetMesh );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class PointSetMeshBuilder< DIMENSION > ;

    public:
        virtual ~PointSetMesh() = default;

        static std::unique_ptr< PointSetMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );
    protected:
        PointSetMesh() = default;
    };

    CLASS_DIMENSION_ALIASES( PointSetMesh );

    template< index_t DIMENSION >
    using PointMeshFactory = GEO::Factory0< PointSetMesh< DIMENSION > >;

    using PointMeshFactory3D = PointMeshFactory< 3 >;
#define ringmesh_register_point_mesh_3d(type) \
    geo_register_creator(RINGMesh::PointMeshFactory3D, type, type::type_name_static())

    using PointMeshFactory2D = PointMeshFactory< 2 >;
#define ringmesh_register_point_mesh_2d(type) \
    geo_register_creator(RINGMesh::PointMeshFactory2D, type, type::type_name_static())

    /*!
     * class for encapsulating line mesh (composed of edges)
     */
    template< index_t DIMENSION >
    class LineMesh: public MeshBase< DIMENSION > {
    ringmesh_disable_copy( LineMesh );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class LineMeshBuilder< DIMENSION > ;
    public:
        virtual ~LineMesh() = default;

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
        double edge_length( index_t edge_id ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                edge_vertex( ElementLocalVertex( edge_id, 0 ) ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                edge_vertex( ElementLocalVertex( edge_id, 1 ) ) );
            return ( e1 - e0 ).length();
        }

        vecn< DIMENSION > edge_barycenter( index_t edge_id ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                edge_vertex( ElementLocalVertex( edge_id, 0 ) ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                edge_vertex( ElementLocalVertex( edge_id, 1 ) ) );
            return ( e1 + e0 ) / 2.;
        }

        /*!
         * @brief return the NNSearch at edges
         * @warning the NNSearch is destroyed when calling the Mesh::polygons_aabb()
         * and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& edge_nn_search() const
        {
            if( !edge_nn_search_ ) {
                std::vector< vecn< DIMENSION > > edge_centers( nb_edges() );
                for( index_t e : range( nb_edges() ) ) {
                    edge_centers[e] = edge_barycenter( e );
                }
                edge_nn_search_.reset(
                    new NNSearch< DIMENSION >( edge_centers, true ) );
            }
            return *edge_nn_search_.get();
        }
        /*!
         * @brief Creates an AABB tree for a Mesh edges
         */
        const LineAABBTree< DIMENSION >& edge_aabb() const
        {
            if( !edge_aabb_ ) {
                edge_aabb_.reset( new LineAABBTree< DIMENSION >( *this ) );
            }
            return *edge_aabb_.get();
        }

        virtual GEO::AttributesManager& edge_attribute_manager() const = 0;
    protected:
        LineMesh() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > edge_nn_search_;
        mutable std::unique_ptr< LineAABBTree< DIMENSION > > edge_aabb_;
    };

    CLASS_DIMENSION_ALIASES( LineMesh );

    template< index_t DIMENSION >
    using LineMeshFactory = GEO::Factory0< LineMesh< DIMENSION > >;

    using LineMeshFactory3D = LineMeshFactory< 3 >;
#define ringmesh_register_line_mesh_3d(type) \
    geo_register_creator(RINGMesh::LineMeshFactory3D, type, type::type_name_static())

    using LineMeshFactory2D = LineMeshFactory< 2 >;
#define ringmesh_register_line_mesh_2d(type) \
    geo_register_creator(RINGMesh::LineMeshFactory2D, type, type::type_name_static())

    /*!
     * class for encapsulating surface mesh component
     */
    template< index_t DIMENSION >
    class SurfaceMeshBase: public MeshBase< DIMENSION > {
    ringmesh_disable_copy( SurfaceMeshBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class SurfaceMeshBuilder< DIMENSION > ;

    public:
        virtual ~SurfaceMeshBase() = default;

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
         * @brief Gets the next vertex index in the polygon vertex
         * \param polygon_local_vertex.
         * @param[in] polygon_local_vertex polygon index and the current
         * local vertex index.
         */
        ElementLocalVertex next_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const
        {
            const index_t local_vertex_id = polygon_local_vertex.local_vertex_id_;
            ringmesh_assert(
                local_vertex_id
                    < nb_polygon_vertices( polygon_local_vertex.element_id_ ) );
            if( local_vertex_id
                != nb_polygon_vertices( polygon_local_vertex.element_id_ ) - 1 ) {
                return ElementLocalVertex( polygon_local_vertex.element_id_,
                    local_vertex_id + 1 );
            } else {
                return ElementLocalVertex( polygon_local_vertex.element_id_, 0 );
            }
        }
        /*!
         * @brief Get the next edge on the border
         * @warning the edge index is in fact the index of the vertex where the edge starts.
         * @details The returned border edge is the next in the way of polygon edges
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
         * @brief Gets the previous vertex index in the polygon \param polygon_id.
         * @param[in] polygon_local_vertex polygon index and its current
         * local vertex index
         */
        ElementLocalVertex prev_polygon_vertex(
            const ElementLocalVertex& polygon_local_vertex ) const
        {
            ringmesh_assert(
                polygon_local_vertex.local_vertex_id_
                    < nb_polygon_vertices( polygon_local_vertex.element_id_ ) );
            if( polygon_local_vertex.local_vertex_id_ > 0 ) {
                return ElementLocalVertex( polygon_local_vertex.element_id_,
                    polygon_local_vertex.local_vertex_id_ - 1 );
            } else {
                return ElementLocalVertex( polygon_local_vertex.element_id_,
                    nb_polygon_vertices( polygon_local_vertex.element_id_ ) - 1 );
            }
        }

        /*!
         * @brief Get the previous edge on the border
         * @details The returned border edge is the previous in the way of polygon edges
         * orientation.
         * @param[in] p Input polygon index
         * @param[in] e Edge index in the polygon
         * @return the previous polygon index
         * and the previous edge index in the polygon (tuple).
         *
         * @pre the surface must be correctly oriented and
         * the given polygon edge must be on border
         * @warning the edge index is in fact the index of the vertex where the edge starts.
         */
        PolygonLocalEdge prev_on_border(
            const PolygonLocalEdge& polygon_local_edge ) const;

        /*!
         * @brief Get the vertex index in a polygon @param polygon_index from its
         * global index in the SurfaceMesh @param vertex_id
         * @return NO_ID or index of the vertex in the polygon
         */
        index_t vertex_index_in_polygon(
            index_t polygon_index,
            index_t vertex_id ) const;

        /*!
         * @brief Compute closest vertex in a polygon to a point
         * @param[in] polygon_index Polygon index
         * @param[in] query_point Coordinates of the point to which distance is measured
         * @return Index of the vertex of @param polygon_index closest to @param query_point
         */
        index_t closest_vertex_in_polygon(
            index_t polygon_index,
            const vecn< DIMENSION >& query_point ) const;

        /*!
         * @brief Get the first polygon of the surface that has an edge linking the two vertices (ids in the surface)
         *
         * @param[in] in0 Index of the first vertex in the surface
         * @param[in] in1 Index of the second vertex in the surface
         * @return NO_ID or the index of the polygon
         */
        index_t polygon_from_vertex_ids( index_t in0, index_t in1 ) const;

        /*!
         * @brief Determines the polygons around a vertex
         * @param[in] vertex_id Index of the vertex in the surface
         * @param[in] border_only If true only polygons on the border are considered
         * @param[in] first_polygon (Optional) Index of one polygon containing the vertex @param P
         * @return Indices of the polygons containing @param P
         * @note If a polygon containing the vertex is given, polygons around this
         * vertex is search by propagation. Else, a first polygon is found by brute
         * force algorithm, and then the other by propagation
         * @todo Try to use a AABB tree to remove @param first_polygon. [PA]
         */
        std::vector< index_t > polygons_around_vertex(
            index_t vertex_id,
            bool border_only,
            index_t first_polygon ) const;

        /*!
         * @brief Gets an adjacent polygon index by polygon index and local edge index.
         * @param[in] polygon_id the polygon index.
         * @param[in] edge_id the local edge index in \param polygon_id.
         * @return the global polygon index adjacent to the \param edge_id of the polygon \param polygon_id.
         * @precondition  \param edge_id < number of edge of the polygon \param polygon_id .
         */
        virtual index_t polygon_adjacent(
            const PolygonLocalEdge& polygon_local_edge ) const = 0;

        virtual GEO::AttributesManager& polygon_attribute_manager() const = 0;
        /*!
         * @brief Tests whether all the polygons are triangles. when all the polygons are triangles, storage and access is optimized.
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
            if( is_triangle( polygon_id ) ) {
                return PolygonType::TRIANGLE;
            } else if( nb_polygon_vertices( polygon_id ) == 4 ) {
                return PolygonType::QUAD;
            } else {
                return PolygonType::UNDEFINED;
            }
        }

        /*!
         * Is the edge starting with the given vertex of the polygon on a border of the Surface?
         */
        bool is_edge_on_border( const PolygonLocalEdge& polygon_local_edge ) const
        {
            return polygon_adjacent( polygon_local_edge ) == NO_ID;
        }

        /*!
         * Is one of the edges of the polygon on the border of the surface?
         */
        bool is_polygon_on_border( index_t polygon_index ) const
        {
            for( index_t v : range( nb_polygon_vertices( polygon_index ) ) ) {
                if( is_edge_on_border( PolygonLocalEdge( polygon_index, v ) ) ) {
                    return true;
                }
            }
            return false;
        }

        /*!
         * @brief Gets the length of the edge starting at a given vertex
         * @param[in] polygon_local_edge index of the polygon and its
         * local edge starting vertex index
         */
        double polygon_edge_length(
            const PolygonLocalEdge& polygon_local_edge ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                polygon_edge_vertex( polygon_local_edge, 0 ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                polygon_edge_vertex( polygon_local_edge, 1 ) );
            return ( e1 - e0 ).length();
        }
        /*!
         * @brief Gets the barycenter of the edge starting at a given vertex
         * @param[in] polygon_local_edge index of the polygon and
         * the local edge starting vertex index
         */
        vecn< DIMENSION > polygon_edge_barycenter(
            const PolygonLocalEdge& polygon_local_edge ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                polygon_edge_vertex( polygon_local_edge, 0 ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                polygon_edge_vertex( polygon_local_edge, 1 ) );
            return ( e1 + e0 ) / 2.;
        }
        /*!
         * @brief Gets the vertex index on the polygon edge
         * @param[in] polygon_local_edge index of the polygon and
         * the local index of the edge in the polygon
         * @param[in] vertex_id index of the local vertex in the edge \param edge_id (0 or 1)
         * @return the vertex index
         */
        index_t polygon_edge_vertex(
            const PolygonLocalEdge& polygon_local_edge,
            index_t vertex_id ) const
        {
            ringmesh_assert( vertex_id < 2 );
            if( vertex_id == 0 ) {
                return polygon_vertex( polygon_local_edge );
            } else {
                return polygon_vertex(
                    ElementLocalVertex( polygon_local_edge.polygon_id_,
                        ( polygon_local_edge.local_edge_id_ + vertex_id )
                            % nb_polygon_vertices(
                                polygon_local_edge.polygon_id_ ) ) );
            }
        }

        /*!
         * Computes the Mesh polygon barycenter
         * @param[in] polygon_id the polygon index
         * @return the polygon center
         */
        vecn< DIMENSION > polygon_barycenter( index_t polygon_id ) const
        {
            vecn< DIMENSION > result;
            ringmesh_assert( nb_polygon_vertices( polygon_id ) >= 1 );
            for( index_t v : range( nb_polygon_vertices( polygon_id ) ) ) {
                result += this->vertex(
                    polygon_vertex( ElementLocalVertex( polygon_id, v ) ) );
            }
            return ( 1.0 / nb_polygon_vertices( polygon_id ) ) * result;
        }
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
            if( !nn_search_ ) {
                std::vector< vecn< DIMENSION > > polygon_centers( nb_polygons() );
                for( index_t p : range( nb_polygons() ) ) {
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
            if( !polygon_aabb_ ) {
                polygon_aabb_.reset( new SurfaceAABBTree< DIMENSION >( *this ) );
            }
            return *polygon_aabb_;
        }
    protected:
        SurfaceMeshBase() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > nn_search_;
        mutable std::unique_ptr< SurfaceAABBTree< DIMENSION > > polygon_aabb_;
    };

    CLASS_DIMENSION_ALIASES( SurfaceMeshBase );

    template< index_t DIMENSION >
    class SurfaceMesh: public SurfaceMeshBase< DIMENSION > {
    };

    template< index_t DIMENSION >
    using SurfaceMeshFactory = GEO::Factory0< SurfaceMesh< DIMENSION > >;

    using SurfaceMeshFactory3D = SurfaceMeshFactory< 3 >;
#define ringmesh_register_surface_mesh_3d(type) \
    geo_register_creator(RINGMesh::SurfaceMeshFactory3D, type, type::type_name_static())

    using SurfaceMeshFactory2D = SurfaceMeshFactory< 2 >;
#define ringmesh_register_surface_mesh_2d(type) \
    geo_register_creator(RINGMesh::SurfaceMeshFactory2D, type, type::type_name_static())

    template< >
    class SurfaceMesh< 3 > : public SurfaceMeshBase< 3 > {
    public:

        /*!
         * Computes the Mesh polygon area
         * @param[in] polygon_id the polygon index
         * @return the polygon area
         */
        double polygon_area( index_t polygon_id ) const override
        {
            double result = 0.0;
            if( nb_polygon_vertices( polygon_id ) == 0 ) {
                return result;
            }
            const vec3& p1 = vertex(
                polygon_vertex( ElementLocalVertex( polygon_id, 0 ) ) );
            for( index_t i : range( 1, nb_polygon_vertices( polygon_id ) - 1 ) ) {
                const vec3& p2 = vertex(
                    polygon_vertex( ElementLocalVertex( polygon_id, i ) ) );
                const vec3& p3 = vertex(
                    polygon_vertex( ElementLocalVertex( polygon_id, i + 1 ) ) );
                result += triangle_signed_area( p1, p2, p3,
                    polygon_normal( polygon_id ) );
            }
            return std::fabs( result );
        }

        /*!
         * Computes the Mesh polygon normal
         * @param[in] polygon_id the polygon index
         * @return the polygon normal
         */
        vec3 polygon_normal( index_t polygon_id ) const
        {
            const vec3& p1 = this->vertex(
                this->polygon_vertex( ElementLocalVertex( polygon_id, 0 ) ) );
            const vec3& p2 = this->vertex(
                this->polygon_vertex( ElementLocalVertex( polygon_id, 1 ) ) );
            const vec3& p3 = this->vertex(
                this->polygon_vertex( ElementLocalVertex( polygon_id, 2 ) ) );
            vec3 norm = cross( p2 - p1, p3 - p1 );
            return normalize( norm );
        }

        /*!
         * @brief Computes the normal of the Mesh2D at the vertex location
         * it computes the average value of polygon normal neighbors
         * @param[in] vertex_id the vertex index
         * @param[in] p0 index of a polygon that contain the vertex \param vertex_id
         * @return the normal at the given vertex
         */
        vec3 normal_at_vertex( index_t vertex_id, index_t p0 = NO_ID ) const
        {
            ringmesh_assert( vertex_id < nb_vertices() );
            index_t p = 0;
            while( p0 == NO_ID && p < nb_polygons() ) {
                for( index_t lv : range( nb_polygon_vertices( p ) ) ) {
                    if( polygon_vertex( ElementLocalVertex( p, lv ) )
                        == vertex_id ) {
                        p0 = p;
                        break;
                    }
                }
                p++;
            }

            std::vector< index_t > polygon_ids = polygons_around_vertex( vertex_id,
                false, p0 );
            vec3 norm;
            for( index_t polygon_id : polygon_ids ) {
                norm += polygon_normal( polygon_id );
            }
            return normalize( norm );
        }
    };

    template< >
    class SurfaceMesh< 2 > : public SurfaceMeshBase< 2 > {
    public:
        /*!
         * Computes the Mesh polygon area
         * @param[in] polygon_id the polygon index
         * @return the polygon area
         */
        double polygon_area( index_t polygon_id ) const override
        {
            double result = 0.0;
            if( nb_polygon_vertices( polygon_id ) == 0 ) {
                return result;
            }
            const vec2& p1 = vertex(
                polygon_vertex( ElementLocalVertex( polygon_id, 0 ) ) );
            for( index_t i : range( 1, nb_polygon_vertices( polygon_id ) - 1 ) ) {
                const vec2& p2 = vertex(
                    polygon_vertex( ElementLocalVertex( polygon_id, i ) ) );
                const vec2& p3 = vertex(
                    polygon_vertex( ElementLocalVertex( polygon_id, i + 1 ) ) );
                result += GEO::Geom::triangle_signed_area( p1, p2, p3 );
            }
            return std::fabs( result );
        }
    };

    CLASS_DIMENSION_ALIASES( SurfaceMesh );

    /*!
     * class for encapsulating volume mesh component
     */
    template< index_t DIMENSION >
    class VolumeMesh: public MeshBase< DIMENSION > {
    ringmesh_disable_copy( VolumeMesh );
        static_assert( DIMENSION == 3, "DIMENSION template should be 3" );
        friend class VolumeMeshBuilder< DIMENSION > ;
    public:
        virtual ~VolumeMesh() = default;

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
         * @brief Gets a vertex index by cell and local edge and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_edge_vertex(
            index_t cell_id,
            index_t edge_id,
            index_t vertex_id ) const = 0;

        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param facet_id
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
        double cell_edge_length( index_t cell_id, index_t edge_id ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                cell_edge_vertex( cell_id, edge_id, 0 ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                cell_edge_vertex( cell_id, edge_id, 1 ) );
            return ( e1 - e0 ).length();
        }

        index_t cells_around_vertex(
            index_t vertex_id,
            std::vector< index_t >& result,
            index_t cell_hint ) const
        {
            result.resize( 0 );

            if( cell_hint == NO_ID ) {
                return 0;
            }

            // Flag the visited cells
            std::vector< index_t > visited;
            visited.reserve( 10 );

            // Stack of the adjacent cells
            std::stack< index_t > S;
            S.push( cell_hint );
            visited.push_back( cell_hint );

            do {
                index_t c = S.top();
                S.pop();

                bool cell_includes_vertex = false;
                for( index_t v : range( nb_cell_vertices( c ) ) ) {
                    if( cell_vertex( ElementLocalVertex( c, v ) ) == vertex_id ) {
                        result.push_back( c );
                        cell_includes_vertex = true;
                        break;
                    }
                }
                if( !cell_includes_vertex ) {
                    continue;
                }

                for( index_t f : range( nb_cell_facets( c ) ) ) {
                    for( index_t v : range(
                        nb_cell_facet_vertices( CellLocalFacet( c, f ) ) ) ) {
                        index_t vertex = cell_facet_vertex( CellLocalFacet( c, f ),
                            v );
                        if( vertex == vertex_id ) {
                            index_t adj_P = cell_adjacent( CellLocalFacet( c, f ) );

                            if( adj_P != NO_ID ) {
                                if( !contains( visited, adj_P ) ) {
                                    S.push( adj_P );
                                    visited.push_back( adj_P );
                                }
                            }
                            break;
                        }
                    }
                }
            } while( !S.empty() );

            return static_cast< index_t >( result.size() );
        }

        /*!
         * Computes the Mesh cell edge barycenter
         * @param[in] cell_id the facet index
         * @param[in] edge_id the edge index
         * @return the cell edge center
         */
        vecn< DIMENSION > cell_edge_barycenter(
            index_t cell_id,
            index_t edge_id ) const
        {
            const vecn< DIMENSION >& e0 = this->vertex(
                cell_edge_vertex( cell_id, edge_id, 0 ) );
            const vecn< DIMENSION >& e1 = this->vertex(
                cell_edge_vertex( cell_id, edge_id, 1 ) );
            return ( e1 + e0 ) / 2.;
        }

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
         * @return the number of vertices in the facet \param facet_id in the cell \param cell_id
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

        virtual GEO::AttributesManager& cell_facet_attribute_manager() const = 0;

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
            const CellLocalFacet& cell_local_facet ) const
        {
            vecn< DIMENSION > result;
            index_t nb_vertices = nb_cell_facet_vertices( cell_local_facet );
            for( index_t v : range( nb_vertices ) ) {
                result += this->vertex( cell_facet_vertex( cell_local_facet, v ) );
            }
            ringmesh_assert( nb_vertices > 0 );

            return result / static_cast< double >( nb_vertices );
        }
        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        vecn< DIMENSION > cell_barycenter( index_t cell_id ) const
        {
            vecn< DIMENSION > result;
            ringmesh_assert( nb_cell_vertices( cell_id ) >= 1 );
            for( index_t v : range( nb_cell_vertices( cell_id ) ) ) {
                result += this->vertex(
                    cell_vertex( ElementLocalVertex( cell_id, v ) ) );
            }
            return ( 1.0 / nb_cell_vertices( cell_id ) ) * result;
        }
        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_local_facet the cell index and
         * the local facet index in the cell
         * @return the cell facet normal
         */
        vecn< DIMENSION > cell_facet_normal(
            const CellLocalFacet& cell_local_facet ) const
        {
            ringmesh_assert( cell_local_facet.cell_id_ < nb_cells() );
            ringmesh_assert(
                cell_local_facet.local_facet_id_
                    < nb_cell_facets( cell_local_facet.cell_id_ ) );

            const vecn< DIMENSION >& p1 = this->vertex(
                cell_facet_vertex( cell_local_facet, 0 ) );
            const vecn< DIMENSION >& p2 = this->vertex(
                cell_facet_vertex( cell_local_facet, 1 ) );
            const vecn< DIMENSION >& p3 = this->vertex(
                cell_facet_vertex( cell_local_facet, 2 ) );

            return cross( p2 - p1, p3 - p1 );
        }

        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        virtual double cell_volume( index_t cell_id ) const = 0;

        std::vector< index_t > cells_around_vertex(
            index_t vertex_id,
            index_t cell_hint ) const;

        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const
        {
            for( index_t v : range( nb_cell_vertices( cell_id ) ) ) {
                if( cell_vertex( ElementLocalVertex( cell_id, v ) ) == vertex_id ) {
                    return v;
                }
            }
            return NO_ID;
        }

        /*!
         * @brief return the NNSearch at cell facets
         * @warning the NNSearch is destroyed when calling the Mesh::facets_aabb()
         *  and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& cell_facet_nn_search() const
        {
            if( !cell_facet_nn_search_ ) {
                std::vector< vecn< DIMENSION > > cell_facet_centers(
                    nb_cell_facets() );
                index_t cf = 0;
                for( index_t c : range( nb_cells() ) ) {
                    for( index_t f : range( nb_cell_facets( c ) ) ) {
                        cell_facet_centers[cf] = cell_facet_barycenter(
                            CellLocalFacet( c, f ) );
                        ++cf;
                    }
                }
                cell_facet_nn_search_.reset(
                    new NNSearch< DIMENSION >( cell_facet_centers, true ) );
            }
            return *cell_facet_nn_search_.get();
        }
        /*!
         * @brief return the NNSearch at cells
         */
        const NNSearch< DIMENSION >& cell_nn_search() const
        {
            if( !cell_nn_search_ ) {
                std::vector< vecn< DIMENSION > > cell_centers( nb_cells() );
                for( index_t c : range( nb_cells() ) ) {
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
            if( !cell_aabb_ ) {
                cell_aabb_.reset( new VolumeAABBTree< DIMENSION >( *this ) );
            }
            return *cell_aabb_.get();
        }
    protected:
        VolumeMesh() = default;
    private:
        /// @TODO use find_cell_from_vertex function in geomodel_builder_geometry
        /// (in the unnamed namespace, geomodel2d branch) which uses C++11 a lambda
        /// function. BC.
        index_t find_first_cell_owing_vertex( index_t vertex_id_in_mesh ) const;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > cell_facet_nn_search_;
        mutable std::unique_ptr< NNSearch< DIMENSION > > cell_nn_search_;
        mutable std::unique_ptr< VolumeAABBTree< DIMENSION > > cell_aabb_;
    };

    using VolumeMesh3D = VolumeMesh< 3 >;

    template< index_t DIMENSION >
    using VolumeMeshFactory = GEO::Factory0< VolumeMesh< DIMENSION > >;

    using VolumeMeshFactory3D = VolumeMeshFactory< 3 >;
#define ringmesh_register_volume_mesh_3d(type) \
    geo_register_creator(RINGMesh::VolumeMeshFactory3D, type, type::type_name_static())

    /*!
     * class composed of meshes from all the dimensions
     */
    template< index_t DIMENSION >
    class MeshSetBase {
    ringmesh_disable_copy( MeshSetBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    public:
        void create_point_set_mesh( const MeshType type );
        void create_line_mesh( const MeshType type );
        void create_surface_mesh( const MeshType type );

    protected:
        MeshSetBase();

    public:
        std::unique_ptr< PointSetMesh< DIMENSION > > point_set_mesh;
        std::unique_ptr< LineMesh< DIMENSION > > line_mesh;
        std::unique_ptr< SurfaceMesh< DIMENSION > > surface_mesh;
    };

    template< index_t DIMENSION >
    class MeshSet: public MeshSetBase< DIMENSION > {
    public:
        MeshSet() = default;
    };

    template< >
    class MeshSet< 3 > : public MeshSetBase< 3 > {
    public:
        MeshSet();

        void create_volume_mesh( const MeshType type );

    public:
        std::unique_ptr< VolumeMesh3D > volume_mesh;
    };
}
