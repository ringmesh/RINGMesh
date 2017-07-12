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

#include <ringmesh/geomodel/geomodel_api.h>

/*!
 * @file ringmesh/geomodel/geomodel_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace {
    using namespace RINGMesh;

    gmme_id find_corner( const GeoModel& geomodel, index_t geomodel_point_id )
    {
        const GeoModelMeshVertices& geomodel_vertices = geomodel.mesh.vertices;
        const std::vector< GMEVertex >& vertices = geomodel_vertices.gme_vertices(
            geomodel_point_id );
        for( const GMEVertex& vertex : vertices ) {
            if( vertex.gmme.type() == Corner::type_name_static() ) {
                return vertex.gmme;
            }
        }
        return gmme_id();
    }

    /*!
     * @brief Returns true if the Line has exactly the given vertices
     * @todo Reimplement using std::iterators
     */
    bool line_equal( const Line& L, const std::vector< index_t >& rhs_vertices )
    {
        if( L.nb_vertices() != rhs_vertices.size() ) {
            return false;
        }
        const GeoModelMeshVertices& geomodel_vertices = L.geomodel().mesh.vertices;
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
    void reorder_closed_line_vertices_to_start_at_corner(
        const GeoModel& geomodel,
        std::vector< index_t >& line_vertices )
    {
        if( geomodel.nb_corners() == 0 ) {
            // Maybe should throw an assertion, but I am not sure [JP]
            // this really may happen for sure, so no throw [RM]
            return;
        }
        if( line_vertices.empty() ) {
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
     * @brief Utility class to sort a set of oriented triangles around a common edge
     * Used in GeoModelBuilderSurface.
     *
     * @todo This code could certainly be improved [JP]
     */
    class GeoModelRegionFromSurfaces {
    public:
        /*!
         * @brief A triangle to sort around an edge, see GeoModelRegionFromSurfaces
         * @details This triangle belongs to a mesh connected component identified by its index.
         */
        struct TriangleToSort {
            /*!
             * @param index Index of this TriangleToSort in GeoModelRegionFromSurfaces
             * @param surface_index Index of the Surface
             * @param p0 point of the triangle
             * @param p1 point of the triangle
             * @param p2 point of the triangle
             */
            TriangleToSort(
                index_t index,
                index_t surface_index,
                const vec3& p0,
                const vec3& p1,
                const vec3& p2 )
                :
                    index_( index ),
                    surface_index_( surface_index ),
                    N_(),
                    B_A_(),
                    angle_( -99999 ),
                    side_( false )
            {
                ringmesh_assert( p0 != p1 );
                ringmesh_assert( p0 != p2 );
                ringmesh_assert( p1 != p2 );

                vec3 e1 = normalize( p1 - p0 );
                vec3 e2 = normalize( p2 - p0 );

                N_ = normalize( cross( e1, e2 ) );
                ringmesh_assert( dot( N_, e1 ) < global_epsilon );

                vec3 B = 0.5 * p1 + 0.5 * p0;
                vec3 p2B = p2 - B;
                B_A_ = normalize( p2B - dot( p2B, e1 ) * e1 );

                ringmesh_assert( dot( B_A_, e1 ) < global_epsilon );
                ringmesh_assert( B_A_.length() > global_epsilon );
            }

            bool operator<( const TriangleToSort& r ) const
            {
                return angle_ < r.angle_;
            }

            /// Index in GeoModelRegionFromSurfaces
            index_t index_;

            /// Global index of the surface owning this triangle
            index_t surface_index_;

            /// Normal to the triangle - normalized vector
            vec3 N_;

            /// Normal to the edge p0p1 in the plane defined by the triangle - normalized
            vec3 B_A_;

            // Values filled by sorting function in GeoModelRegionFromSurfaces
            double angle_;
            bool side_;
        };

        void add_triangle(
            index_t surface_index,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 )
        {
            index_t triangle_id = static_cast< index_t >( triangles_.size() );
            triangles_.push_back(
                TriangleToSort( triangle_id, surface_index, p0, p1, p2 ) );
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
            q *= sin( 0.5 * angle );

            double quat[4] = { q[0], q[1], q[2], cos( 0.5 * angle ) };

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
            ringmesh_assert( triangles_.size() > 0 );

            std::pair< index_t, bool > default_pair( index_t( -1 ), false );
            sorted_triangles_.resize( 2 * triangles_.size(), default_pair );

            // If there is only one Triangle to sort - nothing to do
            if( triangles_.size() == 1 ) {
                sorted_triangles_[0] = std::pair< index_t, bool >(
                    triangles_[0].surface_index_, true );
                sorted_triangles_[1] = std::pair< index_t, bool >(
                    triangles_[0].surface_index_, false );
                return;
            }

            // Initialization
            // We start on the plus (true) side of the first Triangle
            sorted_triangles_[0] = std::pair< index_t, bool >(
                triangles_[0].surface_index_, true );

            // Reference vectors with wich angles will be computed
            vec3 N_ref = triangles_[0].N_;
            vec3 B_A_ref = triangles_[0].B_A_;
            vec3 Ax_ref = normalize( cross( B_A_ref, N_ref ) );

            // The minus (false) side of the start triangle will the last one encountered
            triangles_[0].angle_ = 2 * M_PI;
            triangles_[0].side_ = false;

            for( index_t i = 1; i < triangles_.size(); ++i ) {
                TriangleToSort& cur = triangles_[i];
                // Computes the angle RADIANS between the reference and the current
                // triangle
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
            std::sort( triangles_.begin(), triangles_.end() );

            // Fills the sorted surfaces adding the side
            index_t it = 1;
            for( TriangleToSort& cur : triangles_ ) {
                if( cur.index_ == 0 ) { // The last to add
                    sorted_triangles_[it].first = cur.surface_index_;
                    sorted_triangles_[it].second = cur.side_;
                } else {
                    sorted_triangles_[it].first = cur.surface_index_;
                    sorted_triangles_[it].second = cur.side_;
                    sorted_triangles_[it + 1].first = cur.surface_index_;
                    sorted_triangles_[it + 1].second = !cur.side_;
                    it += 2;
                }
            }
            // All the surfaces must have been sorted
            ringmesh_assert(
                std::count( sorted_triangles_.begin(), sorted_triangles_.end(),
                    default_pair ) == 0 );
        }

        void clear()
        {
            triangles_.clear();
            sorted_triangles_.clear();
        }

        /*! Returns the next pair Triangle index (surface) + side
         */
        const std::pair< index_t, bool >& next(
            const std::pair< index_t, bool >& in ) const
        {
            for( index_t i = 0; i < sorted_triangles_.size(); ++i ) {
                if( sorted_triangles_[i] == in ) {
                    if( i == sorted_triangles_.size() - 1 )
                        return sorted_triangles_[sorted_triangles_.size() - 2];
                    if( i == 0 ) return sorted_triangles_[1];

                    if( sorted_triangles_[i + 1].first
                        == sorted_triangles_[i].first ) {
                        // The next has the same surface id, check its sign
                        if( sorted_triangles_[i + 1].second
                            != sorted_triangles_[i].second ) {
                            return sorted_triangles_[i - 1];
                        } else {
                            // Sign is the same
                            return sorted_triangles_[i + 1];
                        }
                    } else {
                        ringmesh_assert(
                            sorted_triangles_[i - 1].first
                                == sorted_triangles_[i].first );
                        if( sorted_triangles_[i - 1].second
                            != sorted_triangles_[i].second ) {
                            return sorted_triangles_[i + 1];
                        } else {
                            return sorted_triangles_[i - 1];
                        }
                    }
                }
            }
            ringmesh_assert_not_reached;
            return sorted_triangles_.front();
        }

    private:
        std::vector< TriangleToSort > triangles_;
        // Pairs global triangle identifier (Surface index) and side reached
        std::vector< std::pair< index_t, bool > > sorted_triangles_;
    };

    /*! 
     * @brief Determines the geometry of the Lines of a GeoModel in which 
     * the geometry of the Surfaces is given
     * @details All the triangles on the boundaries are classified as belonging to a Line
     * Two neighboring edges on the boundary belong to the same Line if their incident Surfaces
     * are the same.
     */
    class LineGeometryFromGeoModelSurfaces {
    public:
        /*!
         * @param geomodel GeoModel providing the Surfaces
         * @param collect_region_info If true, information needed to determine closed Regions
         *  from a GeoModel Surfaces are collected for each Line.
         */
        LineGeometryFromGeoModelSurfaces(
            const GeoModel& geomodel,
            bool collect_region_info )
            :
                geomodel_( geomodel ),
                collect_region_information_( collect_region_info ),
                cur_border_triangle_( 0 )
        {
            initialize_border_triangles_from_model_surfaces();
            visited_.resize( border_triangles_.size(), false );
        }

        /*!
         * @brief Tries to compute a new line and returns true if one was.
         * @details To use in a while conditional loop, since the number of lines
         * is considered unknown. 
         */
        bool compute_next_line_geometry()
        {
            go_to_next_non_visited_border_triangle();

            if( is_not_the_end() ) {
                init_next_line_computation();
                compute_line_geometry();
                return true;
            } else {
                return false;
            }
        }

        /*! 
         * @brief Returns the vertices of the last line built 
         */
        const std::vector< index_t >& vertices() const
        {
            return cur_line_vertices_;
        }

        /*!
         * @brief Returns the adjacent surfaces to the last line built
         */
        const std::vector< index_t >& adjacent_surfaces() const
        {
            return cur_line_adjacent_surfaces_;
        }

        /*!
         * @brief Returns region building information for the last line built.
         * @details Only computed if the collect_region_info was true at construction.
         */
        const GeoModelRegionFromSurfaces& region_information() const
        {
            return cur_line_region_information_;
        }

    private:
        /*!
         * @brief Stores the vertices of a triangle which is on the boundary of a Surface
         */
        struct BorderTriangle {
            BorderTriangle(
                index_t surface,
                index_t polygon,
                index_t v0,
                index_t v1,
                index_t v2 )
                :
                    v0_( v0 ),
                    v1_( v1 ),
                    v2_( v2 ),
                    surface_( surface ),
                    polygon_( polygon )
            {
            }

            /*!
             * @brief Compares two triangles so those sharing an edge follow one another
             */
            bool operator<( const BorderTriangle& rhs ) const
            {
                if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) {
                    return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ );
                }
                if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) {
                    return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ );
                }
                if( surface_ != rhs.surface_ ) {
                    return surface_ < rhs.surface_;
                }
                if( polygon_ != rhs.polygon_ ) {
                    return polygon_ < rhs.polygon_;
                }
                return rhs.v2_ == NO_ID ? false : v2_ < rhs.v2_;
            }

            bool same_edge( const BorderTriangle& rhs ) const
            {
                return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ )
                    && std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ );
            }

            /// Indices of the points in the geomodel. 
            /// The edge v0v1 is on the one on the boundary.
            index_t v0_;
            index_t v1_;
            index_t v2_;
            // Index of the surface containing the triangle
            index_t surface_;
            // Index of the polygon in the surface
            index_t polygon_;
        };

        void compute_line_geometry()
        {
            visit_border_triangles_on_same_edge( cur_border_triangle_ );
            index_t p0 = border_triangles_[cur_border_triangle_].v0_;
            index_t p1 = border_triangles_[cur_border_triangle_].v1_;
            cur_line_vertices_.push_back( p0 );
            cur_line_vertices_.push_back( p1 );

            get_adjacent_surfaces( cur_border_triangle_,
                cur_line_adjacent_surfaces_ );

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
            ringmesh_assert( is_not_the_end() );

            index_t start( cur_border_triangle_ );
            ringmesh_assert( start != NO_ID );

            index_t t = get_next_border_triangle( start, backward );

            while( t != start ) {
                ringmesh_assert( t != NO_ID );
                if( !is_visited( t ) ) {
                    std::vector< index_t > cur_adjacent_surfaces;
                    get_adjacent_surfaces( t, cur_adjacent_surfaces );

                    if( equal_to_line_adjacent_surfaces( cur_adjacent_surfaces ) ) {
                        visit_border_triangles_on_same_edge( t );
                        add_border_triangle_vertices_to_line( t, backward );
                    } else {
                        // Adjacent surface changed
                        break;
                    }
                } else {
                    // Adjacent surface changed
                    break;
                }
                t = get_next_border_triangle( t, backward );
            }
        }

        void init_next_line_computation()
        {
            cur_line_vertices_.resize( 0 );
            cur_line_adjacent_surfaces_.resize( 0 );
            cur_line_region_information_.clear();
        }

        bool is_visited( index_t i ) const
        {
            return visited_[i];
        }

        void go_to_next_non_visited_border_triangle()
        {
            while( is_not_the_end() && is_visited( cur_border_triangle_ ) ) {
                ++cur_border_triangle_;
            }
        }

        bool is_not_the_end() const
        {
            return cur_border_triangle_ < border_triangles_.size();
        }

        bool have_border_triangles_same_boundary_edge( index_t t0, index_t t1 ) const
        {
            return border_triangles_[t0].same_edge( border_triangles_[t1] );
        }

        /*!
         * @brief Gets triangles sharing the border edge of the current border triangle
         *        and adds them to current line region information 
         */
        void collect_region_information()
        {
            index_t i( cur_border_triangle_ );
            index_t j = i;
            while( j < border_triangles_.size()
                && have_border_triangles_same_boundary_edge( i, j ) ) {
                cur_line_region_information_.add_triangle(
                    border_triangles_[j].surface_,
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v0_ ),
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v1_ ),
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v2_ ) );
                j++;
            }
        }

        bool equal_to_line_adjacent_surfaces(
            const std::vector< index_t >& input ) const
        {
            if( input.size() != cur_line_adjacent_surfaces_.size() ) {
                return false;
            } else {
                return std::equal( input.begin(), input.end(),
                    cur_line_adjacent_surfaces_.begin() );
            }
        }

        void add_border_triangle_vertices_to_line(
            index_t triangle_index,
            bool backward )
        {
            const BorderTriangle& border_triangle = border_triangles_[triangle_index];
            index_t v0 = border_triangle.v0_;
            index_t v1 = border_triangle.v1_;

            add_vertices_to_line( v0, v1, !backward );
        }

        void add_vertices_to_line( index_t v0, index_t v1, bool at_the_end )
        {
            index_t to_add = v1;
            if( at_the_end ) {
                index_t end_vertex = cur_line_vertices_.back();
                if( v1 == end_vertex ) {
                    to_add = v0;
                } else {
                    ringmesh_assert( v0 == end_vertex );
                }
                cur_line_vertices_.push_back( to_add );
            } else {
                index_t first_vertex = cur_line_vertices_.front();
                if( v1 == first_vertex ) {
                    to_add = v0;
                } else {
                    ringmesh_assert( v0 == first_vertex );
                }
                cur_line_vertices_.insert( cur_line_vertices_.begin(), to_add );
            }
        }

        void initialize_border_triangles_from_model_surfaces()
        {
            const GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices;
            for( index_t s = 0; s < geomodel_.nb_surfaces(); ++s ) {
                const Surface& S = geomodel_.surface( s );
                gmme_id S_id = S.gmme();
                for( index_t p = 0; p < S.nb_mesh_elements(); ++p ) {
                    for( index_t v = 0; v < S.nb_mesh_element_vertices( p ); ++v ) {
                        if( S.is_on_border( p, v ) ) {
                            index_t vertex = geomodel_vertices.geomodel_vertex_id(
                                S_id, p, v );
                            index_t next_vertex =
                                geomodel_vertices.geomodel_vertex_id( S_id, p,
                                    S.next_polygon_vertex_index( p, v ) );
                            index_t previous_vertex =
                                geomodel_vertices.geomodel_vertex_id( S_id, p,
                                    S.prev_polygon_vertex_index( p, v ) );
                            border_triangles_.push_back(
                                BorderTriangle( s, p, vertex, next_vertex,
                                    previous_vertex ) );
                        }
                    }
                }
            }
            std::sort( border_triangles_.begin(), border_triangles_.end() );
        }

        /*!
         * @brief Gets the next BorderTriangle in the same surface
         */
        index_t get_next_border_triangle( index_t from, bool backward ) const
        {
            const BorderTriangle& border_triangle = border_triangles_[from];
            const Surface& S = geomodel_.surface( border_triangle.surface_ );

            const GeoModelMeshVertices& geomodel_vertices = geomodel_.mesh.vertices;

            // Gets the next edge on border in the Surface
            index_t p = border_triangle.polygon_;
            std::vector< index_t > possible_v0_id =
                geomodel_vertices.mesh_entity_vertex_id( S.gmme(),
                    border_triangle.v0_ );
            ringmesh_assert( !possible_v0_id.empty() );
            index_t v0_id = NO_ID;
            for( index_t id : possible_v0_id ) {
                if( S.vertex_index_in_polygon( p, id ) != NO_ID ) {
                    v0_id = id;
                }
            }
            ringmesh_assert( v0_id != NO_ID );
            index_t v0_id_in_polygon = S.vertex_index_in_polygon( p, v0_id );
            ringmesh_assert( v0_id_in_polygon != NO_ID );

            index_t next_f = NO_ID;
            index_t next_f_v0 = NO_ID;
            index_t next_f_v1 = NO_ID;

            if( !backward ) {
                S.next_on_border( p, v0_id_in_polygon, next_f, next_f_v0 );
                ringmesh_assert( next_f_v0 != NO_ID );
                next_f_v1 = S.next_polygon_vertex_index( next_f, next_f_v0 );
            } else {
                S.prev_on_border( p, v0_id_in_polygon, next_f, next_f_v0 );
                ringmesh_assert( next_f_v0 != NO_ID );
                next_f_v1 = S.next_polygon_vertex_index( next_f, next_f_v0 );
            }

            // Finds the BorderTriangle that is corresponding to this
            // It must exist and there is only one
            BorderTriangle bait( border_triangle.surface_, next_f,
                geomodel_vertices.geomodel_vertex_id( S.gmme(), next_f, next_f_v0 ),
                geomodel_vertices.geomodel_vertex_id( S.gmme(), next_f, next_f_v1 ),
                NO_ID );
            index_t result = static_cast< index_t >( std::lower_bound(
                border_triangles_.begin(), border_triangles_.end(), bait )
                - border_triangles_.begin() );

            ringmesh_assert( border_triangles_[result].same_edge( bait ) );
            ringmesh_assert( result < border_triangles_.size() );
            return result;
        }

        /*!
         * @brief Marks as visited all BorderTriangles whose first edge is equal to i's first edge
         */
        void visit_border_triangles_on_same_edge( index_t i )
        {
            index_t j = i;
            while( j < border_triangles_.size()
                && border_triangles_[i].same_edge( border_triangles_[j] ) ) {
                visited_[j] = true;
                j++;
            }
            index_t k = i - 1;
            while( k != NO_ID
                && border_triangles_[i].same_edge( border_triangles_[k] ) ) {
                visited_[k] = true;
                k--;
            }
        }

        /*!
         * @brief Gets the sorted indices of the Surfaces incident to the first edge of the i-th BorderTriangle
         * @note When the surface appears twice (the line is an internal border)
         * both occurrences are kept.
         */
        void get_adjacent_surfaces(
            index_t i,
            std::vector< index_t >& adjacent_surfaces )
        {
            index_t j = i;
            while( j < border_triangles_.size()
                && border_triangles_[i].same_edge( border_triangles_[j] ) ) {
                adjacent_surfaces.push_back( border_triangles_[j].surface_ );
                j++;
            }

            index_t k = i - 1;
            while( k != NO_ID
                && border_triangles_[i].same_edge( border_triangles_[k] ) ) {
                adjacent_surfaces.push_back( border_triangles_[k].surface_ );
                k--;
            }
            std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() );
        }

    private:
        const GeoModel& geomodel_;
        bool collect_region_information_;
        // All the triangles on a boundary of all the Surfaces of the GeoModel
        std::vector< BorderTriangle > border_triangles_;
        // Internal use to flag the visited border_triangles when computing the Lines
        std::vector< bool > visited_;

        // Currently computed line information
        index_t cur_border_triangle_;
        std::vector< index_t > cur_line_vertices_;
        std::vector< index_t > cur_line_adjacent_surfaces_;
        RINGMesh::GeoModelRegionFromSurfaces cur_line_region_information_;
    };

    /*************************************************************************/
    GeoModelBuilderFile::GeoModelBuilderFile(
        GeoModel& geomodel,
        const std::string& filename )
        : GeoModelBuilder( geomodel ), filename_( filename )
    {

    }

    GeoModelBuilderFromSurfaces::GeoModelBuilderFromSurfaces(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        :
            options_(),
            builder_( builder ),
            geomodel_( geomodel ),
            geomodel_access_( geomodel )
    {
    }

    GeoModelBuilderFromSurfaces::~GeoModelBuilderFromSurfaces()
    {
        for( GeoModelRegionFromSurfaces*& info : regions_info_ ) {
            delete info;
        }
    }

    bool GeoModelBuilderFromSurfaces::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces line_computer( geomodel_,
            options_.compute_regions_brep );

        bool new_line_was_built = true;
        while( new_line_was_built ) {
            new_line_was_built = line_computer.compute_next_line_geometry();

            // I know this is a copy - but should'nt be too big [JP]
            std::vector< index_t > vertices = line_computer.vertices();

            bool is_line_closed = vertices.front() == vertices.back();
            if( is_line_closed ) {
                // Vertices can begin and end at any vertex
                reorder_closed_line_vertices_to_start_at_corner( geomodel_,
                    vertices );
            }

            gmme_id first_corner = builder_.topology.find_or_create_corner(
                vertices.front() );
            gmme_id second_corner = builder_.topology.find_or_create_corner(
                vertices.back() );
            const std::vector< index_t >& adjacent_surfaces =
                line_computer.adjacent_surfaces();

            index_t backup_nb_lines = geomodel_.nb_lines();
            gmme_id line_index = builder_.topology.find_or_create_line(
                adjacent_surfaces, first_corner, second_corner );

            bool created_line = geomodel_.nb_lines() != backup_nb_lines;
            if( created_line ) {
                builder_.geometry.set_line( line_index.index(), vertices );

                for( index_t j : adjacent_surfaces ) {
                    gmme_id surface_index( Surface::type_name_static(), j );
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
                        new GeoModelRegionFromSurfaces(
                            line_computer.region_information() ) );
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

    bool GeoModelBuilderFromSurfaces::build_brep_regions_from_surfaces()
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
                gmme_id region_id = builder_.topology.create_mesh_entity< Region >();
                gmme_id surface_id( Surface::type_name_static(), 0 );
                builder_.topology.add_mesh_entity_boundary_relation( region_id,
                    surface_id, inside );

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
                gmme_id cur_region_id =
                    builder_.topology.create_mesh_entity< Region >();
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
                        gmme_id( Surface::type_name_static(), s.first ), s.second );
                    surf_2_region[s_id] = cur_region_id.index();

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp =
                        !s.second == true ? 2 * s.first : 2 * s.first + 1;
                    if( surf_2_region[s_id_opp] == NO_ID ) {
                        S.emplace( s.first, !s.second );
                    }
                    // For each contact, push the next oriented surface that is in the same region
                    const Surface& surface = geomodel_.surface( s.first );
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
            const Region& cur_region = geomodel_.region( universe_id );
            std::set< gmme_id > to_erase;
            to_erase.insert( cur_region.gmme() );
            builder_.removal.remove_mesh_entities( to_erase );
        }
        return true;
    }

    void GeoModelBuilderFromSurfaces::build()
    {
        if( geomodel_.nb_surfaces() == 0 ) {
            throw RINGMeshException( "GeoModel",
                "No surface to build the geomodel " );
        }

        // Initialize geomodel() global vertices
        geomodel_.mesh.vertices.test_and_initialize();

        build_lines_and_corners_from_surfaces();

        if( options_.compute_regions_brep ) {
            build_brep_regions_from_surfaces();
        }

        // Finish up the geomodel
        builder_.end_geomodel();
    }

    GeoModelBuilderCopy::GeoModelBuilderCopy(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderCopy::copy_geomodel( const GeoModel& from )
    {
        builder_.topology.copy_topology( from );
        builder_.geometry.copy_meshes( from );
        builder_.geology.copy_geology( from );
    }

    GeoModelBuilderInfo::GeoModelBuilderInfo(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    GeoModelBuilder::GeoModelBuilder( GeoModel& geomodel )
        :
            topology( *this, geomodel ),
            geometry( *this, geomodel ),
            geology( *this, geomodel ),
            removal( *this, geomodel ),
            repair( *this, geomodel ),
            copy( *this, geomodel ),
            info( *this, geomodel ),
            from_surfaces( *this, geomodel ),
            geomodel_( geomodel ),
            geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilder::end_geomodel()
    {
        if( geomodel_.name().empty() ) {
            info.set_geomodel_name( "model_default_name" );
        }
        // Get out if the geomodel has no surface
        if( geomodel_.nb_surfaces() == 0 ) {
            print_geomodel( geomodel_ );
            throw RINGMeshException( "GeoModel", "The GeoModel has no surface" );
        }

        geometry.cut_surfaces_by_internal_lines();
        geometry.cut_regions_by_internal_surfaces();

        // Deliberate clear of the geomodel vertices used for geomodel building
        geomodel_.mesh.vertices.clear();

        print_geomodel( geomodel_ );
    }

    GeoModelBuilderGeology::GeoModelBuilderGeology(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderGeology::copy_geology( const GeoModel& from )
    {
        for( index_t t = 0; t < from.nb_geological_entity_types(); t++ ) {
            builder_.geology.copy_geological_entity_topology( from,
                from.geological_entity_type( t ) );
        }
    }

    bool GeoModelBuilderGeology::create_geological_entities(
        const GeologicalEntityType& type,
        index_t nb_additional_entities )
    {
        find_or_create_geological_entity_type( type );
        std::vector< std::unique_ptr< GeoModelGeologicalEntity > >& store =
            geomodel_access_.modifiable_geological_entities( type );
        index_t old_size = static_cast< index_t >( store.size() );
        index_t new_size = old_size + nb_additional_entities;
        store.reserve( new_size );
        for( index_t i = old_size; i < new_size; i++ ) {
            store.emplace_back(
                GeoModelGeologicalEntityAccess::create_geological_entity( type,
                    geomodel_, i ) );
        }
        return true;
    }

    bool GeoModelBuilderGeology::check_if_boundary_incident_entity_relation_already_exists(
        const gmge_id& parent,
        const gmme_id& children )
    {
        const GeoModelMeshEntity& children_mesh_entity =
            geomodel_.mesh_entity( children );
        if( children_mesh_entity.has_parent( parent.type() ) ) {
            ringmesh_assert(parent == children_mesh_entity.parent_gmge(parent.type()));
            return true;
        } else {
            return false;
        }
    }

    void GeoModelBuilderGeology::add_parent_children_relation(
        const gmge_id& parent,
        const gmme_id& children )
    {
        GeoModelGeologicalEntity& parent_entity =
            geomodel_access_.modifiable_geological_entity( parent );
        const std::vector< GeologicalEntityType >& parent_entity_types =
            geomodel_.entity_type_manager().relationship_manager.parent_types(
                children.type() );
        if( !contains( parent_entity_types, parent.type() ) ) {
            std::ostringstream message;
            message << "Wrong Parent type in the parent children relation between "
                << parent << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }

        if( check_if_boundary_incident_entity_relation_already_exists( parent,
            children ) ) {
            return;
        }

        GeoModelMeshEntity& children_entity =
            geomodel_access_.modifiable_mesh_entity( children );
        const MeshEntityType& children_type =
            geomodel_.entity_type_manager().relationship_manager.child_type(
                parent.type() );

        if( children_type != children.type() ) {
            std::ostringstream message;
            message << "Wrong children type in the parent children relation between "
                << parent << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }

        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id = manager.add_parent_child_relationship( parent,
            children );
        GeoModelGeologicalEntityAccess parent_access( parent_entity );
        parent_access.modifiable_children().push_back( relation_id );
        GeoModelMeshEntityAccess children_access( children_entity );
        children_access.modifiable_parents().push_back( relation_id );
    }

    void GeoModelBuilderGeology::remove_parent_children_relation(
        const gmge_id& parent,
        const gmme_id& children )
    {
        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id = manager.find_parent_child_relationship( parent,
            children );
        if( relation_id == NO_ID ) {
            std::ostringstream message;
            message << "No parent children relation found between " << parent
                << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }
        GeoModelGeologicalEntityAccess parent_access(
            geomodel_access_.modifiable_geological_entity( parent ) );
        std::vector< index_t >& childs = parent_access.modifiable_children();
        std::remove_if( childs.begin(), childs.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
        GeoModelMeshEntityAccess children_access(
            geomodel_access_.modifiable_mesh_entity( children ) );
        std::vector< index_t >& parents = children_access.modifiable_parents();
        std::remove_if( parents.begin(), parents.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
    }
    void GeoModelBuilderGeology::delete_geological_entity(
        const GeologicalEntityType& type,
        index_t index )
    {
        geomodel_access_.modifiable_geological_entities( type )[index].reset();
    }

    gmge_id GeoModelBuilderGeology::create_geological_entity(
        const GeologicalEntityType& type )
    {
        index_t index = find_or_create_geological_entity_type( type );
        index_t id =
            static_cast< index_t >( geomodel_.nb_geological_entities( type ) );
        geomodel_access_.modifiable_geological_entities()[index].emplace_back(
            GeoModelGeologicalEntityAccess::create_geological_entity( type,
                geomodel_, id ) );
        return geomodel_access_.modifiable_geological_entities()[index].back()->gmge();
    }

    index_t GeoModelBuilderGeology::find_or_create_geological_entity_type(
        const GeologicalEntityType& type )
    {
        index_t type_index =
            geomodel_.entity_type_manager().geological_entity_manager.geological_entity_type_index(
                type );
        if( type_index == NO_ID ) {
            type_index = create_geological_entity_type( type );
        }
        return type_index;
    }

    index_t GeoModelBuilderGeology::create_geological_entity_type(
        const GeologicalEntityType& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory::has_creator( type ) );

        geomodel_access_.modifiable_entity_type_manager().geological_entity_manager.geological_entity_types_.push_back(
            type );
        geomodel_access_.modifiable_geological_entities().push_back(
            std::vector< std::unique_ptr< GeoModelGeologicalEntity > >() );
        std::unique_ptr< GeoModelGeologicalEntity > E(
            GeoModelGeologicalEntityFactory::create_object( type, geomodel_ ) );

        const MeshEntityType child_type = E->child_type_name();
        RelationshipManager& parentage =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        parentage.register_geology_relationship( type, child_type );

        return geomodel_.entity_type_manager().geological_entity_manager.nb_geological_entity_types()
            - 1;
    }

    void GeoModelBuilderGeology::copy_geological_entity_topology(
        const GeoModel& from,
        const GeologicalEntityType& type )
    {
        create_geological_entities( type, from.nb_geological_entities( type ) );

        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < geomodel_.nb_geological_entities( type ); ++e ) {
            gmge_id id( type, e );
            GeoModelGeologicalEntityAccess gmge_access(
                geomodel_access_.modifiable_geological_entity( id ) );
            gmge_access.copy( from.geological_entity( id ) );
        }
    }

    void GeoModelBuilderGeology::build_contacts()
    {
        if( geomodel_.entity_type_manager().geological_entity_manager.is_valid_type(
            Contact::type_name_static() )
            && geomodel_.nb_geological_entities( Contact::type_name_static() )
                > 0 ) {
            return;
        }

        std::vector< std::set< gmge_id > > interfaces;
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            const Line& L = geomodel_.line( i );
            std::set< gmge_id > cur_interfaces;
            for( index_t j = 0; j < L.nb_incident_entities(); ++j ) {
                const GeoModelMeshEntity& S = L.incident_entity( j );
                gmge_id parent_interface = S.parent_gmge(
                    Interface::type_name_static() );
                cur_interfaces.insert( parent_interface );
            }
            gmge_id contact_id;
            for( index_t j = 0; j < interfaces.size(); ++j ) {
                if( cur_interfaces.size() == interfaces[j].size()
                    && std::equal( cur_interfaces.begin(), cur_interfaces.end(),
                        interfaces[j].begin() ) ) {
                    contact_id = gmge_id( Contact::type_name_static(), j );
                    break;
                }
            }
            if( !contact_id.is_defined() ) {
                contact_id = create_geological_entity( Contact::type_name_static() );
                ringmesh_assert( contact_id.index() == interfaces.size() );
                interfaces.push_back( cur_interfaces );
                // Create a name for this contact
                std::string name = "contact";
                for( const gmge_id& it : cur_interfaces ) {
                    name += "_";
                    name += geomodel_.geological_entity( it ).name();
                }
                builder_.info.set_geological_entity_name( contact_id, name );
            }
            add_parent_children_relation( contact_id, L.gmme() );
        }
    }

} // namespace
