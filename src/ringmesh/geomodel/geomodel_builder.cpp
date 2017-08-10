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

#include <ringmesh/geomodel/geomodel_builder.h>

#include <stack>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/geomodel_api.h>

/*!
 * @file ringmesh/geomodel/geomodel_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace {
    using namespace RINGMesh;

    template< index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel,
        index_t geomodel_point_id )
    {
        const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
            geomodel.mesh.vertices;
        const std::vector< GMEVertex >& vertices = geomodel_vertices.gme_vertices(
            geomodel_point_id );
        for( const GMEVertex& vertex : vertices ) {
            if( vertex.gmme.type() == Corner< DIMENSION >::type_name_static() ) {
                return vertex.gmme;
            }
        }
        return gmme_id();
    }

    /*!
     * @brief Returns true if the Line< DIMENSION > has exactly the given vertices
     * @todo Reimplement using std::iterators
     */
    template< index_t DIMENSION >
    bool line_equal(
        const Line< DIMENSION >& line,
        const std::vector< index_t >& rhs_vertices )
    {
        if( line.nb_vertices() != rhs_vertices.size() ) {
            return false;
        }
        const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
            line.geomodel().mesh.vertices;
        bool equal = true;
        for( index_t i : range( line.nb_vertices() ) ) {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( line.gmme(), i ) ) {
                equal = false;
                break;
            }
        }
        if( equal ) {
            return true;
        }
        // If the order is the other one
        equal = true;
        for( index_t i : range( line.nb_vertices() ) ) {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( line.gmme(),
                    line.nb_vertices() - i - 1 ) ) {
                equal = false;
                break;
            }
        }
        return equal;
    }

    /*!
     * @brief Reorders the line so that front() is a corner.
     * @note Closed line has front()==back().
     */
    template< index_t DIMENSION >
    void reorder_closed_line_vertices_to_start_at_corner(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< index_t >& line_vertices )
    {
        if( geomodel.nb_corners() == 0 || line_vertices.empty() ) {
            return;
        }
        for( index_t i : range( 1, line_vertices.size() - 1 ) ) {
            gmme_id corner = find_corner( geomodel, line_vertices[i] );
            if( corner.is_defined() ) {
                line_vertices.pop_back();
                std::rotate( line_vertices.begin(), line_vertices.begin() + i,
                    line_vertices.end() );
                line_vertices.push_back( line_vertices.front() );
                break;
            }
        }
    }

    /*!
     * @brief Stores the vertices of a polygon which is on the boundary of a Surface
     */
    struct BorderPolygon {
        BorderPolygon( index_t surface, index_t polygon, index_t v0, index_t v1 )
            : v0_( v0 ), v1_( v1 ), surface_( surface ), polygon_( polygon )
        {
        }

        /*!
         * @brief Compares two polygons so those sharing an edge follow one another
         */
        bool operator<( const BorderPolygon& rhs ) const
        {
            if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) {
                return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ );
            } else if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) {
                return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ );
            } else if( surface_ != rhs.surface_ ) {
                return surface_ < rhs.surface_;
            } else {
                return polygon_ < rhs.polygon_;
            }
        }

        bool same_edge( const BorderPolygon& rhs ) const
        {
            return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ )
                && std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ );
        }

        /// Indices of the points in the geomodel.
        /// The edge v0v1 is on the one on the boundary.
        index_t v0_;
        index_t v1_;
        // Index of the surface containing the polygon
        index_t surface_;
        // Index of the polygon in the surface
        index_t polygon_;
    };

    template< index_t DIMENSION >
    class CommonDataFromGeoModelSurfaces {
    ringmesh_disable_copy( CommonDataFromGeoModelSurfaces );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
    protected:
        CommonDataFromGeoModelSurfaces( const GeoModel< DIMENSION >& geomodel )
            : geomodel_( geomodel )
        {
            const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
                geomodel_.mesh.vertices;
            for( const auto& surface : geomodel_.surfaces() ) {
                const auto& mesh = surface.low_level_mesh_storage();
                gmme_id S_id = surface.gmme();
                for( index_t p : range( surface.nb_mesh_elements() ) ) {
                    for( index_t v : range( surface.nb_mesh_element_vertices( p ) ) ) {
                        if( mesh.is_edge_on_border( PolygonLocalEdge( p, v ) ) ) {
                            index_t vertex = geomodel_vertices.geomodel_vertex_id(
                                S_id, ElementLocalVertex( p, v ) );
                            index_t next_vertex =
                                geomodel_vertices.geomodel_vertex_id( S_id,
                                    mesh.next_polygon_vertex(
                                        ElementLocalVertex( p, v ) ) );
                            border_polygons_.emplace_back( surface.index(), p,
                                vertex, next_vertex );
                        }
                    }
                }
            }
            std::sort( border_polygons_.begin(), border_polygons_.end() );
        }

        bool have_border_polygons_same_boundary_edge( index_t t0, index_t t1 ) const
        {
            return this->border_polygons_[t0].same_edge( this->border_polygons_[t1] );
        }

    protected:
        const GeoModel< DIMENSION >& geomodel_;
        // All the polygons on a boundary of all the Surfaces of the GeoModel
        std::vector< BorderPolygon > border_polygons_;
    };

    ALIAS_2D_AND_3D( CommonDataFromGeoModelSurfaces );

    /*!
     * @brief Utility class to sort a set of oriented polygons around a common edge
     */
    class GeoModelRegionFromSurfaces {
    public:
        /*!
         * @brief A polygon to sort around an edge, see GeoModelRegionFromSurfaces
         * @details This polygon belongs to a mesh connected component identified by its index.
         */
        struct PolygonToSort {
            /*!
             * @param index Index of this PolygonToSort in GeoModelRegionFromSurfaces
             * @param surface_index Index of the Surface
             * @param p0 point of the polygon
             * @param p1 point of the polygon
             */
            PolygonToSort(
                index_t index,
                index_t surface_index,
                const vec3& normal,
                const vec3& p0,
                const vec3& p1 )
                :
                    index_( index ),
                    surface_index_( surface_index ),
                    N_( normal ),
                    angle_( -max_float64() ),
                    side_( false )
            {
                ringmesh_assert( p0 != p1 );

                vec3 e1 = normalize( p1 - p0 );
                B_A_ = normalize( cross( e1, N_ ) );
            }

            bool operator<( const PolygonToSort& r ) const
            {
                return angle_ < r.angle_;
            }

            /// Index in GeoModelRegionFromSurfaces
            index_t index_;

            /// Global index of the surface owning this polygon
            index_t surface_index_;

            /// Normal to the polygon - normalized vector
            vec3 N_;

            /// Normal to the edge p0p1 in the plane defined by the polygon - normalized
            vec3 B_A_;

            // Values filled by sorting function in GeoModelRegionFromSurfaces
            double angle_;
            bool side_;
        };

        void add_polygon_edge(
            index_t surface_index,
            const vec3& normal,
            const vec3& p0,
            const vec3& p1 )
        {
            index_t polygon_id = static_cast< index_t >( polygons_.size() );
            polygons_.emplace_back( polygon_id, surface_index, normal, p0, p1 );
        }

        /*!
         * @details Rotation around axis of an angle A of the vector V
         *
         * @param axis Oriented axis
         * @param angle Angle in RADIANS
         * @param V the vector to rotate
         */
        static vec3 rotate( const vec3& axis, double angle, const vec3& V )
        {
            vec3 q = normalize( axis ) * std::sin( 0.5 * angle );
            vecn< 4 > quat( q[0], q[1], q[2], std::cos( 0.5 * angle ) );

            GEO::Matrix< 4, double > mat;
            mat( 0, 0 ) = 1 - 2.0 * ( quat[1] * quat[1] + quat[2] * quat[2] );
            mat( 1, 0 ) = 2.0 * ( quat[0] * quat[1] + quat[2] * quat[3] );
            mat( 2, 0 ) = 2.0 * ( quat[2] * quat[0] - quat[1] * quat[3] );
            mat( 3, 0 ) = 0.0;

            mat( 0, 1 ) = 2.0 * ( quat[0] * quat[1] - quat[2] * quat[3] );
            mat( 1, 1 ) = 1 - 2.0 * ( quat[2] * quat[2] + quat[0] * quat[0] );
            mat( 2, 1 ) = 2.0 * ( quat[1] * quat[2] + quat[0] * quat[3] );
            mat( 3, 1 ) = 0.0;

            mat( 0, 2 ) = 2.0 * ( quat[2] * quat[0] + quat[1] * quat[3] );
            mat( 1, 2 ) = 2.0 * ( quat[1] * quat[2] - quat[0] * quat[3] );
            mat( 2, 2 ) = 1 - 2.0 * ( quat[1] * quat[1] + quat[0] * quat[0] );
            mat( 3, 2 ) = 0.0;

            mat( 0, 3 ) = 0.0;
            mat( 1, 3 ) = 0.0;
            mat( 2, 3 ) = 0.0;
            mat( 3, 3 ) = 1.0;

            vecn< 4 > normalizedV( V.x, V.y, V.z, 1.0 );
            vecn< 4 > result;
            GEO::mult( mat, normalizedV.data(), result.data() );

            ringmesh_assert( std::fabs( result.w ) > global_epsilon );
            double inv_w = 1.0 / result.w;
            return vec3( result.x * inv_w, result.y * inv_w, result.z * inv_w );
        }

        void sort()
        {
            ringmesh_assert( polygons_.size() > 0 );

            std::pair< index_t, bool > default_pair( index_t( -1 ), false );
            sorted_polygons_.resize( 2 * polygons_.size(), default_pair );

            // If there is only one Polygon to sort - nothing to do
            if( polygons_.size() == 1 ) {
                sorted_polygons_[0] = std::pair< index_t, bool >(
                    polygons_[0].surface_index_, true );
                sorted_polygons_[1] = std::pair< index_t, bool >(
                    polygons_[0].surface_index_, false );
                return;
            }

            // Initialization
            // We start on the plus (true) side of the first Polygon
            sorted_polygons_[0] = std::pair< index_t, bool >(
                polygons_[0].surface_index_, true );

            // Reference vectors with wich angles will be computed
            vec3 N_ref = polygons_[0].N_;
            vec3 B_A_ref = polygons_[0].B_A_;
            vec3 Ax_ref = normalize( cross( B_A_ref, N_ref ) );

            // The minus (false) side of the start polygon will the last one encountered
            polygons_[0].angle_ = 2 * M_PI;
            polygons_[0].side_ = false;

            for( index_t i : range( 1, polygons_.size() ) ) {
                PolygonToSort& cur = polygons_[i];
                // Computes the angle RADIANS between the reference and the current
                // polygon
                double cos = dot( B_A_ref, cur.B_A_ );
                // Remove invalid values
                if( cos < -1 )
                    cos = -1;
                else if( cos > 1 ) cos = 1;
                cur.angle_ = std::acos( cos );
                // Put the angle between PI and 2PI if necessary
                if( dot( cross( B_A_ref, cur.B_A_ ), Ax_ref ) < 0. ) {
                    cur.angle_ = 2 * M_PI - cur.angle_;
                }

                // Get the side of the surface first encountered
                // when rotating in the N_ref direction
                vec3 N_rotate = rotate( Ax_ref, -cur.angle_, cur.N_ );
                cur.side_ = dot( N_rotate, N_ref ) > 0 ? false : true;
            }

            // Sorts the Surfaces according to the angle
            std::sort( polygons_.begin(), polygons_.end() );

            // Fills the sorted surfaces adding the side
            index_t it = 1;
            for( PolygonToSort& cur : polygons_ ) {
                if( cur.index_ == 0 ) { // The last to add
                    sorted_polygons_[it].first = cur.surface_index_;
                    sorted_polygons_[it].second = cur.side_;
                } else {
                    sorted_polygons_[it].first = cur.surface_index_;
                    sorted_polygons_[it].second = cur.side_;
                    sorted_polygons_[it + 1].first = cur.surface_index_;
                    sorted_polygons_[it + 1].second = !cur.side_;
                    it += 2;
                }
            }
            // All the surfaces must have been sorted
            ringmesh_assert(
                std::count( sorted_polygons_.begin(), sorted_polygons_.end(),
                    default_pair ) == 0 );
        }

        void clear()
        {
            polygons_.clear();
            sorted_polygons_.clear();
        }

        /*! Returns the next pair Polygon index (surface) + side
         */
        const std::pair< index_t, bool >& next(
            const std::pair< index_t, bool >& in ) const
        {
            for( index_t i : range( sorted_polygons_.size() ) ) {
                if( sorted_polygons_[i] == in ) {
                    if( i == sorted_polygons_.size() - 1 )
                        return sorted_polygons_[sorted_polygons_.size() - 2];
                    if( i == 0 ) return sorted_polygons_[1];

                    if( sorted_polygons_[i + 1].first
                        == sorted_polygons_[i].first ) {
                        // The next has the same surface id, check its sign
                        if( sorted_polygons_[i + 1].second
                            != sorted_polygons_[i].second ) {
                            return sorted_polygons_[i - 1];
                        } else {
                            // Sign is the same
                            return sorted_polygons_[i + 1];
                        }
                    } else {
                        ringmesh_assert(
                            sorted_polygons_[i - 1].first
                                == sorted_polygons_[i].first );
                        if( sorted_polygons_[i - 1].second
                            != sorted_polygons_[i].second ) {
                            return sorted_polygons_[i + 1];
                        } else {
                            return sorted_polygons_[i - 1];
                        }
                    }
                }
            }
            ringmesh_assert_not_reached;
            return sorted_polygons_.front();
        }

    private:
        std::vector< PolygonToSort > polygons_;
        // Pairs global polygon identifier (Surface index) and side reached
        std::vector< std::pair< index_t, bool > > sorted_polygons_;
    };

    class RegionTopologyFromGeoModelSurfaces: public CommonDataFromGeoModelSurfaces<
        3 > {
    public:
        RegionTopologyFromGeoModelSurfaces( const GeoModel3D& geomodel )
            :
                CommonDataFromGeoModelSurfaces3D( geomodel ),
                region_info_( geomodel.nb_lines() )
        {
        }

        void compute_region_info()
        {
            const GeoModelMeshVertices3D& vertices = this->geomodel_.mesh.vertices;
            for( const auto& line : geomodel_.lines() ) {
                BorderPolygon line_border( NO_ID, NO_ID,
                    vertices.geomodel_vertex_id( line.gmme(), 0 ),
                    vertices.geomodel_vertex_id( line.gmme(), 1 ) );
                for( const BorderPolygon& border : this->border_polygons_ ) {
                    if( line_border.same_edge( border ) ) {
                        index_t surface_id = border.surface_;
                        region_info_[line.index()].add_polygon_edge( surface_id,
                            this->geomodel_.surface( surface_id ).low_level_mesh_storage().polygon_normal(
                                border.polygon_ ), vertices.vertex( border.v0_ ),
                            vertices.vertex( border.v1_ ) );
                    }
                }
                region_info_[line.index()].sort();
            }
        }

        const std::vector< GeoModelRegionFromSurfaces >& region_info() const
        {
            return region_info_;
        }

    private:
        std::vector< GeoModelRegionFromSurfaces > region_info_;
    };

    struct LineDefinition {

        bool is_closed()
        {
            return vertices_.front() == vertices_.back();
        }

        template< index_t DIMENSION >
        void reoder_closed_line_vertices( const GeoModel< DIMENSION >& geomodel )
        {
            // Vertices can begin and end at any vertex
            reorder_closed_line_vertices_to_start_at_corner( geomodel, vertices_ );
        }

        void clear()
        {
            vertices_.clear();
            adjacent_surfaces_.clear();
        }

        std::vector< index_t > vertices_;
        std::vector< index_t > adjacent_surfaces_;
    };

    /*!
     * @brief Determines the geometry of the Lines of a GeoModel in which
     * the geometry of the Surfaces is given
     * @details All the polygons on the boundaries are classified as belonging to a Line< DIMENSION >
     * Two neighboring edges on the boundary belong to the same Line< DIMENSION > if their incident Surfaces
     * are the same.
     */
    template< index_t DIMENSION >
    class LineGeometryFromGeoModelSurfaces: public CommonDataFromGeoModelSurfaces<
        DIMENSION > {
    public:
        /*!
         * @param geomodel GeoModel providing the Surfaces
         * @param collect_region_info If true, information needed to determine closed Regions
         *  from a GeoModel Surfaces are collected for each Line< DIMENSION >.
         */
        LineGeometryFromGeoModelSurfaces( const GeoModel< DIMENSION >& geomodel )
            :
                CommonDataFromGeoModelSurfaces< DIMENSION >( geomodel ),
                cur_border_polygon_( 0 )
        {
            visited_.resize( this->border_polygons_.size(), false );
        }

        /*!
         * @brief Tries to compute a new line and returns true if one was.
         * @details To use in a while conditional loop, since the number of lines
         * is considered unknown.
         */
        bool compute_next_line_geometry()
        {
            for( ; cur_border_polygon_ < this->border_polygons_.size();
                cur_border_polygon_++ ) {
                if( is_visited( cur_border_polygon_ ) ) {
                    continue;
                }
                cur_line_.clear();
                compute_line_geometry();
                return true;
            }
            return false;
        }

        LineDefinition& current_line()
        {
            return cur_line_;
        }

    private:

        void compute_line_geometry()
        {
            visit_border_polygons_on_same_edge( cur_border_polygon_ );
            index_t p0 = this->border_polygons_[cur_border_polygon_].v0_;
            index_t p1 = this->border_polygons_[cur_border_polygon_].v1_;
            cur_line_.vertices_.push_back( p0 );
            cur_line_.vertices_.push_back( p1 );

            cur_line_.adjacent_surfaces_ = get_adjacent_surfaces(
                cur_border_polygon_ );

            bool backward = false;
            get_one_line_vertices( backward );
            backward = true;
            get_one_line_vertices( backward );
        }
        /*!
         * @brief Gets the geometry of a line propagating in a given direction
         * @details As long as the adjacent surfaces are the same, the vertices are added
         * belong to the line under construction
         */
        void get_one_line_vertices( bool backward )
        {
            ringmesh_assert( cur_border_polygon_ < this->border_polygons_.size() );
            ringmesh_assert( cur_border_polygon_ != NO_ID );

            index_t t = get_next_border_polygon( cur_border_polygon_, backward );
            ringmesh_assert( t != NO_ID );
            while( propagate( t ) ) {
                visit_border_polygons_on_same_edge( t );
                add_border_polygon_vertices_to_line( t, backward );
                t = get_next_border_polygon( t, backward );
                ringmesh_assert( t != NO_ID );
            }
        }

        bool propagate( index_t t ) const
        {
            return t != cur_border_polygon_ && !is_visited( t )
                && equal_to_line_adjacent_surfaces( t );
        }

        bool is_visited( index_t i ) const
        {
            return visited_[i];
        }

        bool equal_to_line_adjacent_surfaces( index_t t ) const
        {
            std::vector< index_t > input = get_adjacent_surfaces( t );
            if( input.size() != cur_line_.adjacent_surfaces_.size() ) {
                return false;
            } else {
                return std::equal( input.begin(), input.end(),
                    cur_line_.adjacent_surfaces_.begin() );
            }
        }

        void add_border_polygon_vertices_to_line(
            index_t polygon_index,
            bool backward )
        {
            const BorderPolygon& border_polygon =
                this->border_polygons_[polygon_index];
            add_vertices_to_line( border_polygon.v0_, border_polygon.v1_,
                !backward );
        }

        void add_vertices_to_line( index_t v0, index_t v1, bool at_the_end )
        {
            index_t to_add = v1;
            if( at_the_end ) {
                index_t end_vertex = cur_line_.vertices_.back();
                if( v1 == end_vertex ) {
                    to_add = v0;
                } else {
                    ringmesh_assert( v0 == end_vertex );
                }
                cur_line_.vertices_.push_back( to_add );
            } else {
                index_t first_vertex = cur_line_.vertices_.front();
                if( v1 == first_vertex ) {
                    to_add = v0;
                } else {
                    ringmesh_assert( v0 == first_vertex );
                }
                cur_line_.vertices_.insert( cur_line_.vertices_.begin(), to_add );
            }
        }

        /*!
         * @brief Gets the next BorderPolygon in the same surface
         */
        index_t get_next_border_polygon( index_t from, bool backward ) const
        {
            const BorderPolygon& border_polygon = this->border_polygons_[from];
            const Surface< DIMENSION >& S = this->geomodel_.surface(
                border_polygon.surface_ );
            const SurfaceMesh< DIMENSION >& mesh = S.low_level_mesh_storage();
            gmme_id surface_id = S.gmme();
            const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
                this->geomodel_.mesh.vertices;

            // Gets the next edge on border in the Surface
            index_t p = border_polygon.polygon_;
            std::vector< index_t > possible_v0_id =
                geomodel_vertices.mesh_entity_vertex_id( surface_id,
                    border_polygon.v0_ );
            ringmesh_assert( !possible_v0_id.empty() );
            index_t v0_id = NO_ID;
            for( index_t id : possible_v0_id ) {
                if( mesh.vertex_index_in_polygon( p, id ) != NO_ID ) {
                    v0_id = id;
                }
            }
            ringmesh_assert( v0_id != NO_ID );
            index_t v0_id_in_polygon = mesh.vertex_index_in_polygon( p, v0_id );
            ringmesh_assert( v0_id_in_polygon != NO_ID );

            const PolygonLocalEdge cur_polygon_local_edge( p, v0_id_in_polygon );
            PolygonLocalEdge next_polygon_local_edge0_on_border =
                backward ?
                    mesh.prev_on_border( cur_polygon_local_edge ) :
                    mesh.next_on_border( cur_polygon_local_edge );
            ringmesh_assert(
                next_polygon_local_edge0_on_border.polygon_id_ != NO_ID );
            ringmesh_assert(
                next_polygon_local_edge0_on_border.local_edge_id_ != NO_ID );

            PolygonLocalEdge next_polygon_local_edge1_on_border =
                mesh.next_polygon_vertex( next_polygon_local_edge0_on_border );
            ringmesh_assert(
                next_polygon_local_edge1_on_border.polygon_id_ != NO_ID );
            ringmesh_assert(
                next_polygon_local_edge1_on_border.local_edge_id_ != NO_ID );

            // Finds the BorderPolygon that is corresponding to this
            // It must exist and there is only one
            BorderPolygon bait( border_polygon.surface_,
                next_polygon_local_edge0_on_border.polygon_id_,
                geomodel_vertices.geomodel_vertex_id( surface_id,
                    next_polygon_local_edge0_on_border ),
                geomodel_vertices.geomodel_vertex_id( surface_id,
                    next_polygon_local_edge1_on_border ) );
            index_t result = find_sorted( this->border_polygons_, bait );

            ringmesh_assert( result != NO_ID );
            ringmesh_assert( this->border_polygons_[result].same_edge( bait ) );
            return result;
        }

        /*!
         * @brief Marks as visited all BorderPolygons whose first edge is equal to i's first edge
         */
        void visit_border_polygons_on_same_edge( index_t border_id )
        {
            visited_[border_id] = true;
            for( index_t next_border_id = border_id + 1;
                next_border_id < this->border_polygons_.size()
                    && this->have_border_polygons_same_boundary_edge( border_id,
                        next_border_id ); next_border_id++ ) {
                visited_[next_border_id] = true;
            }
            for( index_t prev_border_id = border_id - 1;
                prev_border_id != NO_ID
                    && this->have_border_polygons_same_boundary_edge( border_id,
                        prev_border_id ); prev_border_id-- ) {
                visited_[prev_border_id] = true;
            }
        }

        /*!
         * @brief Gets the sorted indices of the Surfaces incident to the first edge of the i-th BorderPolygon
         * @note When the surface appears twice (the line is an internal border)
         * both occurrences are kept.
         */
        std::vector< index_t > get_adjacent_surfaces( index_t border_id ) const
        {
            std::vector< index_t > adjacent_surfaces;
            adjacent_surfaces.reserve( 10 );
            adjacent_surfaces.push_back(
                this->border_polygons_[border_id].surface_ );

            for( index_t next_border_id = border_id + 1;
                next_border_id < this->border_polygons_.size()
                    && this->have_border_polygons_same_boundary_edge( border_id,
                        next_border_id ); next_border_id++ ) {
                adjacent_surfaces.push_back(
                    this->border_polygons_[next_border_id].surface_ );
            }

            for( index_t prev_border_id = border_id - 1;
                prev_border_id != NO_ID
                    && this->have_border_polygons_same_boundary_edge( border_id,
                        prev_border_id ); prev_border_id-- ) {
                adjacent_surfaces.push_back(
                    this->border_polygons_[prev_border_id].surface_ );
            }
            std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() );
            return adjacent_surfaces;
        }

    private:
        // Internal use to flag the visited border_polygons when computing the Lines
        std::vector< bool > visited_;

        // Currently computed line information
        index_t cur_border_polygon_;
        LineDefinition cur_line_;
    };

}

namespace RINGMesh {

    template< index_t DIMENSION >
    GeoModelBuilderCopy< DIMENSION >::GeoModelBuilderCopy(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderCopy< DIMENSION >::copy_geomodel(
        const GeoModel< DIMENSION >& from )
    {
        builder_.topology.copy_topology( from );
        builder_.geometry.copy_meshes( from );
        builder_.geology.copy_geology( from );
    }

    template< index_t DIMENSION >
    GeoModelBuilderInfo< DIMENSION >::GeoModelBuilderInfo(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    GeoModelBuilderBase< DIMENSION >::GeoModelBuilderBase(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        :
            topology( builder, geomodel ),
            geometry( builder, geomodel ),
            geology( builder, geomodel ),
            removal( builder, geomodel ),
            repair( builder, geomodel ),
            copy( builder, geomodel ),
            info( builder, geomodel ),
            geomodel_( geomodel ),
            geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderBase< DIMENSION >::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces< DIMENSION > line_computer( geomodel_ );

        while( line_computer.compute_next_line_geometry() ) {
            LineDefinition line = line_computer.current_line();

            if( line.is_closed() ) {
                line.reoder_closed_line_vertices( geomodel_ );
            }

            std::vector< index_t >& vertices = line.vertices_;
            gmme_id first_corner = topology.find_or_create_corner(
                vertices.front() );
            gmme_id second_corner = topology.find_or_create_corner(
                vertices.back() );

            const std::vector< index_t >& adjacent_surfaces = line.adjacent_surfaces_;
            index_t backup_nb_lines = geomodel_.nb_lines();
            gmme_id line_index = topology.find_or_create_line( adjacent_surfaces,
                first_corner, second_corner );

            bool created_line = geomodel_.nb_lines() != backup_nb_lines;
            if( created_line ) {
                geometry.set_line( line_index.index(), vertices );

                for( index_t j : adjacent_surfaces ) {
                    gmme_id surface_index( Surface< DIMENSION >::type_name_static(),
                        j );
                    topology.add_mesh_entity_boundary_relation( surface_index,
                        line_index );
                }
                topology.add_mesh_entity_boundary_relation( line_index,
                    first_corner );
                topology.add_mesh_entity_boundary_relation( line_index,
                    second_corner );
            } else {
                bool same_geometry = line_equal(
                    geomodel_.line( line_index.index() ), vertices );
                if( !same_geometry ) {
                    geometry.set_line( line_index.index(), vertices );
                }
            }
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderBase< DIMENSION >::end_geomodel()
    {
        if( geomodel_.name().empty() ) {
            info.set_geomodel_name( "model_default_name" );
        }

        cut_geomodel_on_internal_boundaries();
        topology.compute_universe();

        // Deliberate clear of the geomodel vertices used for geomodel building
        geometry.clear_geomodel_mesh();

        print_geomodel( geomodel_ );
    }

    GeoModelBuilder< 2 >::GeoModelBuilder( GeoModel2D& geomodel )
        : GeoModelBuilderBase< 2 >( *this, geomodel )
    {
    }

    template< >
    void GeoModelBuilderBase< 2 >::cut_geomodel_on_internal_boundaries()
    {
        geometry.cut_surfaces_by_internal_lines();
    }

    GeoModelBuilder< 3 >::GeoModelBuilder( GeoModel3D& geomodel )
        : GeoModelBuilderBase< 3 >( *this, geomodel )
    {
    }

    template< >
    void GeoModelBuilderBase< 3 >::cut_geomodel_on_internal_boundaries()
    {
        geometry.cut_surfaces_by_internal_lines();
        geometry.cut_regions_by_internal_surfaces();
    }

    void GeoModelBuilder< 3 >::build_regions_from_lines_and_surfaces()
    {
        RegionTopologyFromGeoModelSurfaces region_computer( geomodel_ );
        region_computer.compute_region_info();

        const std::vector< GeoModelRegionFromSurfaces >& region_info =
            region_computer.region_info();

        if( geomodel_.nb_surfaces() < 2 || geomodel_.nb_lines() == 0 ) {
            throw RINGMeshException( "GeoModel",
                "You need at least 1 line and 2 surfaces to use GeoModelBuilder::build_regions_from_lines_and_surfaces" );
        }

        // Each side of each Surface is in one Region( +side is first )
        std::vector< index_t > surf_2_region( 2 * geomodel_.nb_surfaces(), NO_ID );

        // Start with the first Surface on its + side
        std::stack< std::pair< index_t, bool > > S;
        S.emplace( 0, true );

        while( !S.empty() ) {
            std::pair< index_t, bool > cur = S.top();
            S.pop();
            // This side is already assigned
            if( surf_2_region[cur.second == true ? 2 * cur.first : 2 * cur.first + 1]
                != NO_ID ) {
                continue;
            }
            // Create a new region
            gmme_id cur_region_id( Region3D::type_name_static(),
                geomodel_.nb_regions() );
            topology.create_mesh_entities( Region3D::type_name_static(), 1 );
            // Get all oriented surfaces defining this region
            std::stack< std::pair< index_t, bool > > SR;
            SR.push( cur );
            while( !SR.empty() ) {
                std::pair< index_t, bool > s = SR.top();
                SR.pop();
                index_t s_id = s.second == true ? 2 * s.first : 2 * s.first + 1;
                // This oriented surface has already been visited
                if( surf_2_region[s_id] != NO_ID ) {
                    continue;
                }
                // Add the surface to the current region
                topology.add_mesh_entity_boundary_relation( cur_region_id,
                    gmme_id( Surface3D::type_name_static(), s.first ), s.second );
                surf_2_region[s_id] = cur_region_id.index();

                // Check the other side of the surface and push it in S
                index_t s_id_opp = !s.second == true ? 2 * s.first : 2 * s.first + 1;
                if( surf_2_region[s_id_opp] == NO_ID ) {
                    S.emplace( s.first, !s.second );
                }
                // For each contact, push the next oriented surface that is in the same region
                const Surface3D& surface = geomodel_.surface( s.first );
                for( index_t i : range( surface.nb_boundaries() ) ) {
                    const std::pair< index_t, bool >& n =
                        region_info[surface.boundary_gmme( i ).index()].next( s );
                    index_t n_id = n.second == true ? 2 * n.first : 2 * n.first + 1;

                    if( surf_2_region[n_id] == NO_ID ) {
                        SR.push( n );
                    }
                }
            }
        }

        // Check if all the surfaces were visited
        // If not, this means that there are additionnal regions included in those built
        if( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID ) != 0 ) {
            Logger::err( "GeoModel",
                "Small bubble regions were skipped at geomodel building " );
            // Or, most probably, we have a problem before
            ringmesh_assert( false );
            /// @todo handle the region building of small bubble regions
        }

        topology.compute_universe();
        // We need to remove from the regions_ the one corresponding
        // to the universe_, the one with the biggest volume
        double max_volume = -1.;
        index_t universe_id = NO_ID;
        for( const auto& region : geomodel_.regions() ) {
            double cur_volume = region.size();
            if( cur_volume > max_volume ) {
                max_volume = cur_volume;
                universe_id = region.index();
            }
        }
        const Region3D& cur_region = geomodel_.region( universe_id );
        for( index_t i : range( cur_region.nb_boundaries() ) ) {
            // Fill the Universe region boundaries
            // They are supposed to be empty
            topology.add_universe_boundary( cur_region.boundary( i ).index(),
                cur_region.side( i ) );
        }
        std::set< gmme_id > to_erase;
        to_erase.insert( cur_region.gmme() );
        removal.remove_mesh_entities( to_erase );
    }

    template class RINGMESH_API GeoModelBuilderBase< 2 > ;
    template class RINGMESH_API GeoModelBuilderInfo< 2 > ;
    template class RINGMESH_API GeoModelBuilderCopy< 2 > ;

    template class RINGMESH_API GeoModelBuilderBase< 3 > ;
    template class RINGMESH_API GeoModelBuilderInfo< 3 > ;
    template class RINGMESH_API GeoModelBuilderCopy< 3 > ;

}
