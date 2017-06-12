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
 * @file ringmesh/geomodel/geomodel_builder_from_surfaces.cpp
 * @brief Implementation of the classes to build GeoModel from surfaces
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
        const Line< DIMENSION >& L,
        const std::vector< index_t >& rhs_vertices )
    {
        if( L.nb_vertices() != rhs_vertices.size() ) {
            return false;
        }
        const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
            L.geomodel().mesh.vertices;
        bool equal = true;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( L.gmme(), i ) ) {
                equal = false;
                break;
            }
        }
        if( equal ) {
            return true;
        }
        // If the order is the other one
        equal = true;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( L.gmme(),
                    L.nb_vertices() - i - 1 ) ) {
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
        for( index_t i = 1; i + 1 < line_vertices.size(); ++i ) {
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

} // anonymous namespace

namespace RINGMesh {
    /*!
     * @brief Utility class to sort a set of oriented polygons around a common edge
     * Used in GeoModelBuilderSurface.
     *
     * @todo This code could certainly be improved [JP]
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
                    B_A_(),
                    angle_( -99999 ),
                    side_( false )
            {
                ringmesh_assert( p0 != p1 );
                vec3 e1 = normalize( p1 - p0 );
                ringmesh_assert( dot( N_, e1 ) < global_epsilon );
                B_A_ = normalize( cross( e1, N_ ) );
                ringmesh_assert( dot( B_A_, e1 ) < global_epsilon );
                ringmesh_assert( B_A_.length() > global_epsilon );
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
            vec3 q = axis;
            if( q.length() > 0 ) {
                double s = 1.0 / q.length();
                q[0] *= s;
                q[1] *= s;
                q[2] *= s;
            }
            q *= std::sin( 0.5 * angle );

            double quat[4] = { q[0], q[1], q[2], std::cos( 0.5 * angle ) };

            double m[4][4];

            m[0][0] = 1 - 2.0 * ( quat[1] * quat[1] + quat[2] * quat[2] );
            m[0][1] = 2.0 * ( quat[0] * quat[1] + quat[2] * quat[3] );
            m[0][2] = 2.0 * ( quat[2] * quat[0] - quat[1] * quat[3] );
            m[0][3] = 0.0;

            m[1][0] = 2.0 * ( quat[0] * quat[1] - quat[2] * quat[3] );
            m[1][1] = 1 - 2.0 * ( quat[2] * quat[2] + quat[0] * quat[0] );
            m[1][2] = 2.0 * ( quat[1] * quat[2] + quat[0] * quat[3] );
            m[1][3] = 0.0;

            m[2][0] = 2.0 * ( quat[2] * quat[0] + quat[1] * quat[3] );
            m[2][1] = 2.0 * ( quat[1] * quat[2] - quat[0] * quat[3] );
            m[2][2] = 1 - 2.0 * ( quat[1] * quat[1] + quat[0] * quat[0] );
            m[2][3] = 0.0;

            m[3][0] = 0.0;
            m[3][1] = 0.0;
            m[3][2] = 0.0;
            m[3][3] = 1.0;

            double x = V[0] * m[0][0] + V[1] * m[1][0] + V[2] * m[2][0] + m[3][0];
            double y = V[0] * m[0][1] + V[1] * m[1][1] + V[2] * m[2][1] + m[3][1];
            double z = V[0] * m[0][2] + V[1] * m[1][2] + V[2] * m[2][2] + m[3][2];
            double w = V[0] * m[0][3] + V[1] * m[1][3] + V[2] * m[2][3] + m[3][3];
            return vec3( x / w, y / w, z / w );
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

            for( index_t i = 1; i < polygons_.size(); ++i ) {
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
            for( index_t i = 0; i < sorted_polygons_.size(); ++i ) {
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
            region_information_.clear();
        }
        std::vector< index_t > vertices_;
        std::vector< index_t > adjacent_surfaces_;
        RINGMesh::GeoModelRegionFromSurfaces region_information_;
    };

    /*! 
     * @brief Determines the geometry of the Lines of a GeoModel in which
     * the geometry of the Surfaces is given
     * @details All the polygons on the boundaries are classified as belonging to a Line< 3 >
     * Two neighboring edges on the boundary belong to the same Line< 3 > if their incident Surfaces
     * are the same.
     */
    template< index_t DIMENSION >
    class LineGeometryFromGeoModelSurfaces {
    public:
        /*!
         * @param geomodel GeoModel providing the Surfaces
         * @param collect_region_info If true, information needed to determine closed Regions
         *  from a GeoModel Surfaces are collected for each Line< DIMENSION >.
         */
        LineGeometryFromGeoModelSurfaces(
            const GeoModel< DIMENSION >& geomodel,
            bool collect_region_info )
            :
                geomodel_( geomodel ),
                collect_region_information_( collect_region_info ),
                cur_border_polygon_( 0 )
        {
            initialize_border_polygons_from_model_surfaces();
            visited_.resize( border_polygons_.size(), false );
        }

        /*!
         * @brief Tries to compute a new line and returns true if one was.
         * @details To use in a while conditional loop, since the number of lines
         * is considered unknown. 
         */
        bool compute_next_line_geometry()
        {
            for( ; cur_border_polygon_ < border_polygons_.size();
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

        void compute_line_geometry()
        {
            visit_border_polygons_on_same_edge( cur_border_polygon_ );
            index_t p0 = border_polygons_[cur_border_polygon_].v0_;
            index_t p1 = border_polygons_[cur_border_polygon_].v1_;
            cur_line_.vertices_.push_back( p0 );
            cur_line_.vertices_.push_back( p1 );

            cur_line_.adjacent_surfaces_ = get_adjacent_surfaces(
                cur_border_polygon_ );

            bool backward = false;
            get_one_line_vertices( backward );
            backward = true;
            get_one_line_vertices( backward );

            if( collect_region_information_ ) {
                collect_region_information();
            }
        }
        /*!
         * @brief Gets the geometry of a line propagating in a given direction
         * @details As long as the adjacent surfaces are the same, the vertices are added
         * belong to the line under construction
         */
        void get_one_line_vertices( bool backward )
        {
            ringmesh_assert( cur_border_polygon_ < border_polygons_.size() );
            ringmesh_assert( cur_border_polygon_ != NO_ID );

            index_t t = get_next_border_polygon( cur_border_polygon_, backward );
            ringmesh_assert( t != NO_ID );
            while( t != cur_border_polygon_ && !is_visited( t )
                && equal_to_line_adjacent_surfaces( get_adjacent_surfaces( t ) ) ) {
                visit_border_polygons_on_same_edge( t );
                add_border_polygon_vertices_to_line( t, backward );
                t = get_next_border_polygon( t, backward );
                ringmesh_assert( t != NO_ID );
            }
        }

        bool is_visited( index_t i ) const
        {
            return visited_[i];
        }

        bool have_border_polygons_same_boundary_edge( index_t t0, index_t t1 ) const
        {
            return border_polygons_[t0].same_edge( border_polygons_[t1] );
        }

        /*!
         * @brief Gets polygons sharing the border edge of the current border polygon
         *        and adds them to current line region information 
         */
        void collect_region_information()
        {
            const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
                geomodel_.mesh.vertices;
            for( index_t j = cur_border_polygon_; j < border_polygons_.size();
                j++ ) {
                if( have_border_polygons_same_boundary_edge( cur_border_polygon_,
                    j ) ) {
                    const BorderPolygon& border = border_polygons_[j];
                    index_t surface_id = border.surface_;
                    cur_line_.region_information_.add_polygon_edge( surface_id,
                        geomodel_.surface( surface_id ).low_level_mesh_storage().polygon_normal(
                            border.polygon_ ),
                        geomodel_vertices.vertex( border.v0_ ),
                        geomodel_vertices.vertex( border.v1_ ) );
                }
            }
        }

        bool equal_to_line_adjacent_surfaces(
            const std::vector< index_t >& input ) const
        {
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
            const BorderPolygon& border_polygon = border_polygons_[polygon_index];
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

        void initialize_border_polygons_from_model_surfaces()
        {
            const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
                geomodel_.mesh.vertices;
            for( index_t s = 0; s < geomodel_.nb_surfaces(); s++ ) {
                const Surface< DIMENSION >& S = geomodel_.surface( s );
                const SurfaceMesh< DIMENSION >& mesh = S.low_level_mesh_storage();
                gmme_id S_id = S.gmme();
                for( index_t p = 0; p < S.nb_mesh_elements(); p++ ) {
                    for( index_t v = 0; v < S.nb_mesh_element_vertices( p ); ++v ) {
                        if( mesh.is_edge_on_border( p, v ) ) {
                            index_t vertex = geomodel_vertices.geomodel_vertex_id(
                                S_id, p, v );
                            index_t next_vertex =
                                geomodel_vertices.geomodel_vertex_id( S_id, p,
                                    mesh.next_polygon_vertex( p, v ) );
                            border_polygons_.emplace_back( s, p, vertex,
                                next_vertex );
                        }
                    }
                }
            }
            std::sort( border_polygons_.begin(), border_polygons_.end() );
        }

        /*!
         * @brief Gets the next BorderPolygon in the same surface
         */
        index_t get_next_border_polygon( index_t from, bool backward ) const
        {
            const BorderPolygon& border_polygon = border_polygons_[from];
            const Surface< DIMENSION >& S = geomodel_.surface(
                border_polygon.surface_ );
            const SurfaceMesh< DIMENSION >& mesh = S.low_level_mesh_storage();
            gmme_id surface_id = S.gmme();
            const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
                geomodel_.mesh.vertices;

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
                    break;
                }
            }
            ringmesh_assert( v0_id != NO_ID );
            index_t v0_id_in_polygon = mesh.vertex_index_in_polygon( p, v0_id );
            ringmesh_assert( v0_id_in_polygon != NO_ID );

            index_t next_f = NO_ID;
            index_t next_f_v0 = NO_ID;
            index_t next_f_v1 = NO_ID;

            if( !backward ) {
                mesh.next_on_border( p, v0_id_in_polygon, next_f, next_f_v0 );
                ringmesh_assert( next_f_v0 != NO_ID );
                next_f_v1 = mesh.next_polygon_vertex( next_f, next_f_v0 );
            } else {
                mesh.prev_on_border( p, v0_id_in_polygon, next_f, next_f_v0 );
                ringmesh_assert( next_f_v0 != NO_ID );
                next_f_v1 = mesh.next_polygon_vertex( next_f, next_f_v0 );
            }

            // Finds the BorderPolygon that is corresponding to this
            // It must exist and there is only one
            BorderPolygon bait( border_polygon.surface_, next_f,
                geomodel_vertices.geomodel_vertex_id( surface_id, next_f,
                    next_f_v0 ),
                geomodel_vertices.geomodel_vertex_id( surface_id, next_f,
                    next_f_v1 ) );
            index_t result = find_sorted( border_polygons_, bait );

            ringmesh_assert( result != NO_ID );
            ringmesh_assert( border_polygons_[result].same_edge( bait ) );
            return result;
        }

        /*!
         * @brief Marks as visited all BorderPolygons whose first edge is equal to i's first edge
         */
        void visit_border_polygons_on_same_edge( index_t border_id )
        {
            for( index_t i = 0; i < border_polygons_.size(); i++ ) {
                if( have_border_polygons_same_boundary_edge( border_id, i ) ) {
                    visited_[i] = true;
                }
            }
        }

        /*!
         * @brief Gets the sorted indices of the Surfaces incident to the first edge of the i-th BorderPolygon
         * @note When the surface appears twice (the line is an internal border)
         * both occurrences are kept.
         */
        std::vector< index_t > get_adjacent_surfaces( index_t border_id )
        {
            std::vector< index_t > adjacent_surfaces;
            for( index_t i = 0; i < border_polygons_.size(); i++ ) {
                if( have_border_polygons_same_boundary_edge( border_id, i ) ) {
                    adjacent_surfaces.push_back( border_polygons_[i].surface_ );
                }
            }
            std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() );
            return adjacent_surfaces;
        }

    private:
        const GeoModel< DIMENSION >& geomodel_;
        bool collect_region_information_;
        // All the polygons on a boundary of all the Surfaces of the GeoModel
        std::vector< BorderPolygon > border_polygons_;
        // Internal use to flag the visited border_polygons when computing the Lines
        std::vector< bool > visited_;

        // Currently computed line information
        index_t cur_border_polygon_;
        LineDefinition cur_line_;
    };

    /*************************************************************************/

    template< index_t DIMENSION >
    GeoModelBuilderFromSurfaces< DIMENSION >::GeoModelBuilderFromSurfaces(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        :
            options_(),
            builder_( builder ),
            geomodel_( geomodel ),
            geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    GeoModelBuilderFromSurfaces< DIMENSION >::~GeoModelBuilderFromSurfaces()
    {
        for( GeoModelRegionFromSurfaces*& info : regions_info_ ) {
            delete info;
        }
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderFromSurfaces< DIMENSION >::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces< DIMENSION > line_computer( geomodel_,
            options_.compute_regions_brep );

        while( line_computer.compute_next_line_geometry() ) {
            LineDefinition line = line_computer.current_line();

            if( line.is_closed() ) {
                line.reoder_closed_line_vertices( geomodel_ );
            }

            std::vector< index_t >& vertices = line.vertices_;
            gmme_id first_corner = builder_.topology.find_or_create_corner(
                vertices.front() );
            gmme_id second_corner = builder_.topology.find_or_create_corner(
                vertices.back() );

            const std::vector< index_t >& adjacent_surfaces = line.adjacent_surfaces_;
            index_t backup_nb_lines = geomodel_.nb_lines();
            gmme_id line_index = builder_.topology.find_or_create_line(
                adjacent_surfaces, first_corner, second_corner );

            bool created_line = geomodel_.nb_lines() != backup_nb_lines;
            if( created_line ) {
                builder_.geometry.set_line( line_index.index(), vertices );

                for( index_t j : adjacent_surfaces ) {
                    gmme_id surface_index( Surface< DIMENSION >::type_name_static(),
                        j );
                    builder_.topology.add_mesh_entity_boundary_relation(
                        surface_index, line_index );
                }
                builder_.topology.add_mesh_entity_boundary_relation( line_index,
                    first_corner );
                builder_.topology.add_mesh_entity_boundary_relation( line_index,
                    second_corner );

                // If the plan is to then build_regions, get the information
                if( options_.compute_regions_brep ) {
                    regions_info_.push_back(
                        new GeoModelRegionFromSurfaces( line.region_information_ ) );
                }
            } else {
                bool same_geometry = line_equal(
                    geomodel_.line( line_index.index() ), vertices );
                if( !same_geometry ) {
                    builder_.geometry.set_line( line_index.index(), vertices );
                }
            }
        }
        return true;
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderFromSurfaces< DIMENSION >::build_brep_regions_from_surfaces()
    {
        ringmesh_assert( geomodel_.nb_lines() == regions_info_.size() );

        // Sort surfaces around the contacts
        for( GeoModelRegionFromSurfaces*& info : regions_info_ ) {
            info->sort();
        }

        if( geomodel_.nb_surfaces() == 1 ) {
            if( geomodel_.nb_lines() != 0 ) {
                Logger::err( "GeoModel",
                    "The unique surface provided to build the geomodel has boundaries " );
                return false;
            } else {
                /// If there is only one surface, its inside is set to be
                /// the + side. No further check.
                bool inside = true;
                // Create the region - set the surface on its boundaries
                gmme_id region_id( Region< DIMENSION >::type_name_static(),
                    geomodel_.nb_regions() );
                builder_.topology.create_mesh_entities(
                    Region< DIMENSION >::type_name_static(), 1 );
                gmme_id surface_id( Surface< DIMENSION >::type_name_static(), 0 );
                builder_.topology.add_mesh_entity_boundary_relation( region_id,
                    surface_id, inside );

                // Set universe boundary
                builder_.topology.add_universe_boundary( 0, !inside );
            }
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2 * geomodel_.nb_surfaces(),
                NO_ID );

            // Start with the first Surface on its + side
            std::stack< std::pair< index_t, bool > > S;
            S.emplace( 0, true );

            while( !S.empty() ) {
                std::pair< index_t, bool > cur = S.top();
                S.pop();
                // This side is already assigned
                if( surf_2_region[
                    cur.second == true ? 2 * cur.first : 2 * cur.first + 1]
                    != NO_ID ) {
                    continue;
                }
                // Create a new region
                gmme_id cur_region_id( Region< DIMENSION >::type_name_static(),
                    geomodel_.nb_regions() );
                builder_.topology.create_mesh_entities(
                    Region< DIMENSION >::type_name_static(), 1 );
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
                    builder_.topology.add_mesh_entity_boundary_relation(
                        cur_region_id,
                        gmme_id( Surface< DIMENSION >::type_name_static(), s.first ),
                        s.second );
                    surf_2_region[s_id] = cur_region_id.index();

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp =
                        !s.second == true ? 2 * s.first : 2 * s.first + 1;
                    if( surf_2_region[s_id_opp] == NO_ID ) {
                        S.emplace( s.first, !s.second );
                    }
                    // For each contact, push the next oriented surface that is in the same region
                    const Surface< DIMENSION >& surface = geomodel_.surface(
                        s.first );
                    for( index_t i = 0; i < surface.nb_boundaries(); ++i ) {
                        const std::pair< index_t, bool >& n =
                            regions_info_[surface.boundary_gmme( i ).index()]->next(
                                s );
                        index_t n_id =
                            n.second == true ? 2 * n.first : 2 * n.first + 1;

                        if( surf_2_region[n_id] == NO_ID ) {
                            SR.push( n );
                        }
                    }
                }
            }

            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            if( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID )
                != 0 ) {
                Logger::err( "GeoModel",
                    "Small bubble regions were skipped at geomodel building " );
                // Or, most probably, we have a problem before
                ringmesh_assert( false );
            }

            builder_.topology.compute_universe();
            // We need to remove from the regions_ the one corresponding
            // to the universe_, the one with the biggest volume
            double max_volume = -1.;
            index_t universe_id = NO_ID;
            for( index_t i = 0; i < geomodel_.nb_regions(); ++i ) {
                double cur_volume = geomodel_.region( i ).size();
                if( cur_volume > max_volume ) {
                    max_volume = cur_volume;
                    universe_id = i;
                }
            }
            const Region< DIMENSION >& cur_region = geomodel_.region( universe_id );
            for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ) {
                // Fill the Universe region boundaries
                // They are supposed to be empty
                builder_.topology.add_universe_boundary(
                    cur_region.boundary( i ).index(), cur_region.side( i ) );
            }
            std::set< gmme_id > to_erase;
            to_erase.insert( cur_region.gmme() );
            builder_.removal.remove_mesh_entities( to_erase );
        }
        return true;
    }

    template class RINGMESH_API GeoModelBuilderFromSurfaces< 3 > ;

} // namespace
