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

#include <ringmesh/geo_model_builder.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <set>
#include <stack>

#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/colocate.h>

#include <ringmesh/io.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_mesh_repair.h>
#include <ringmesh/utils.h>

/*!
 * @file ringmesh/geo_model_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

/*! @todo Split All functions of geo_model_builder.cpp into smaller functions
 * Split this file into at least 4 files.
 */

namespace {
    using namespace RINGMesh ;
    using GEO::Logger ;

    typedef GeoModelEntity::gme_t gme_t ;
    typedef GeoModelMeshEntity GMME ;

    void get_entity_vertices_and_update_corners(
        std::vector< index_t >& corners,
        std::vector< index_t >& vertices )
    {
        std::vector< index_t > corners_to_vertices ;
        get_unique_input_values_and_mapping< index_t >( corners, vertices,
            corners_to_vertices ) ;
        corners = corners_to_vertices ;
    }

    /*************************************************************************/
    /*!
     * @brief Gets the index of an Interface from its name
     * @param[in] geomodel GeoModel to consider
     * @param[in] interface_name Name of the interface to find
     * @return Index of the interface in the model, NO_ID if not found.
     */
    gme_t find_interface(
        const GeoModel& geomodel,
        const std::string& interface_name )
    {
        for( index_t i = 0; i < geomodel.nb_interfaces(); ++i ) {
            if( geomodel.one_interface( i ).name() == interface_name ) {
                return geomodel.one_interface( i ).gme_id() ;
            }
        }
        return gme_t() ;
    }

    /*!
     * @brief Structure used to build Line by GeoModelBuilderGocad
     */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1 )
            : part_id_( part ), corner_id_( corner ), p0_( p0 ), p1_( p1 )
        {
        }

        // Id of the Surface owning this Border
        index_t part_id_ ;

        // Id of p0 in the GeoModel corner vector
        index_t corner_id_ ;

        // Ids of the starting corner and second vertex on the border in the Surface
        // to which this Border belong
        index_t p0_ ;
        index_t p1_ ;
    } ;

    /*!
     * @brief Gets the index of the Corner for a given point
     * @param[in] geomodel GeoModel to consider
     * @param[in] point Geometric location to look for
     * @return NO_ID or the index of the Corner
     */
    gme_t find_corner( const GeoModel& geomodel, const vec3& point )
    {
        for( index_t i = 0; i < geomodel.nb_corners(); ++i ) {
            if( geomodel.corner( i ).vertex( 0 ) == point ) {
                return gme_t( GME::CORNER, i ) ;
            }
        }
        return gme_t() ;
    }

    /*!
     * @brief Gets the index of the Corner at a given model point
     * @param[in] geomodel GeoModel to consider
     * @param[in] model_point_id Index of the point in the GeoModel
     * @return NO_ID or the index of the Corner
     */
    gme_t find_corner( const GeoModel& geomodel, index_t model_point_id )
    {
        for( index_t i = 0; i < geomodel.nb_corners(); ++i ) {
            if( geomodel.corner( i ).model_vertex_id() == model_point_id ) {
                return gme_t( GME::CORNER, i ) ;
            }
        }
        return gme_t() ;
    }

    /*!
     * @brief Returns true if the Line has exactly the given vertices
     *
     * @param[in] L the line to compare to
     * @param[in] rhs_vertices Vertices to compare to
     */
    bool line_equal( const Line& L, const std::vector< vec3 >& rhs_vertices )
    {
        if( L.nb_vertices() != rhs_vertices.size() ) {
            return false ;
        }
        bool equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( i ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) return true ;

        equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( L.nb_vertices() - i - 1 ) ) {
                equal = false ;
                break ;
            }
        }
        return equal ;
    }

    /*!
     * @brief Returns true if the Line has exactly the given vertices
     * @todo Reimplement using std::iterators
     */
    bool line_equal( const Line& L, const std::vector< index_t >& rhs_vertices )
    {
        if( L.nb_vertices() != rhs_vertices.size() ) {
            return false ;
        }
        bool equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.model_vertex_id( i ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) {
            return true ;
        }
        equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.model_vertex_id( L.nb_vertices() - i - 1 ) ) {
                equal = false ;
                break ;
            }
        }
        return equal ;
    }

    /*!
     * Finds a facet and its edge index that are colocalised with an edge
     * defined by its two model vertex indices
     * @param[in] ann a ColocatorANN of the Surface \p surface using the keyword FACETS
     * @param[in] surface the surface where to find the facet
     * @param[in] model_v0 the first model vertex index of the edge
     * @param[in] model_v1 the second model vertex index of the edge
     * @param[out] f the found facet index
     * @param[out] e the found edge index
     * @return True if the facet and the edge indices are found
     * @todo RENAME these parameters and split in smaller functions !! [JP]
     */
    bool find_facet_and_edge(
        const ColocaterANN& ann,
        const Surface& surface,
        index_t model_v0,
        index_t model_v1,
        index_t& f,
        index_t& e )
    {
        // This is bad ! One level of abstraction is far far away
        const vec3& v0 = surface.model().mesh.vertices.vertex( model_v0 ) ;
        const vec3& v1 = surface.model().mesh.vertices.vertex( model_v1 ) ;
        vec3 v_bary = 0.5 * ( v0 + v1 ) ;

        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, surface.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = ann.get_neighbors( v_bary, cur_neighbor, neighbors,
                dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                f = neighbors[i] ;
                for( index_t j = 0; j < surface.nb_mesh_element_vertices( f ); j++ ) {
                    if( surface.model_vertex_id( f, j ) == model_v0 ) {
                        index_t j_next = surface.next_facet_vertex_index( f, j ) ;
                        if( surface.model_vertex_id( f, j_next ) == model_v1 ) {
                            e = j ;
                            return true ;
                        }
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor ) ;

        f = NO_ID ;
        e = NO_ID ;
        return false ;
    }

    bool is_corner_to_duplicate( const GeoModel& geomodel, index_t corner_id )
    {
        if( geomodel.corner( corner_id ).nb_in_boundary() > 3 ) {
            return true ;
        } else {
            return false ;
        }
    }

    void get_sorted_incident_surfaces(
        const GeoModelEntity& E,
        std::vector< index_t >& incident_surfaces )
    {
        index_t nb = E.nb_in_boundary() ;
        incident_surfaces.resize( nb ) ;
        for( index_t i = 0; i < nb; ++i ) {
            incident_surfaces[i] = E.in_boundary_gme( i ).index ;
        }
        std::sort( incident_surfaces.begin(), incident_surfaces.end() ) ;
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
                ringmesh_assert( p0 != p1 ) ;
                ringmesh_assert( p0 != p2 ) ;
                ringmesh_assert( p1 != p2 ) ;

                vec3 e1 = normalize( p1 - p0 ) ;
                vec3 e2 = normalize( p2 - p0 ) ;

                N_ = normalize( cross( e1, e2 ) ) ;
                ringmesh_assert( dot( N_, e1 ) < epsilon ) ;

                vec3 B = 0.5 * p1 + 0.5 * p0 ;
                vec3 p2B = p2 - B ;
                B_A_ = normalize( p2B - dot( p2B, e1 ) * e1 ) ;

                ringmesh_assert( dot( B_A_, e1 ) < epsilon ) ;
                ringmesh_assert( B_A_.length() > epsilon ) ;
            }
            ;

            bool operator<( const TriangleToSort& r ) const
            {
                return angle_ < r.angle_ ;
            }

            /// Index in GeoModelRegionFromSurfaces
            index_t index_ ;

            /// Global index of the surface owning this triangle
            index_t surface_index_ ;

            /// Normal to the triangle - normalized vector
            vec3 N_ ;

            /// Normal to the edge p0p1 in the plane defined by the triangle - normalized
            vec3 B_A_ ;

            // Values filled by sorting function in GeoModelRegionFromSurfaces
            double angle_ ;
            bool side_ ;
        } ;

        void add_triangle(
            index_t surface_index,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 )
        {
            triangles_.push_back(
                TriangleToSort( triangles_.size(), surface_index, p0, p1, p2 ) ) ;
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
            vec3 q = axis ;
            if( q.length() > 0 ) {
                double s = 1.0 / q.length() ;
                q[0] *= s ;
                q[1] *= s ;
                q[2] *= s ;
            }
            q *= sin( 0.5 * angle ) ;

            double quat[4] = { q[0], q[1], q[2], cos( 0.5 * angle ) } ;

            double m[4][4] ;

            m[0][0] = 1 - 2.0 * ( quat[1] * quat[1] + quat[2] * quat[2] ) ;
            m[0][1] = 2.0 * ( quat[0] * quat[1] + quat[2] * quat[3] ) ;
            m[0][2] = 2.0 * ( quat[2] * quat[0] - quat[1] * quat[3] ) ;
            m[0][3] = 0.0 ;

            m[1][0] = 2.0 * ( quat[0] * quat[1] - quat[2] * quat[3] ) ;
            m[1][1] = 1 - 2.0 * ( quat[2] * quat[2] + quat[0] * quat[0] ) ;
            m[1][2] = 2.0 * ( quat[1] * quat[2] + quat[0] * quat[3] ) ;
            m[1][3] = 0.0 ;

            m[2][0] = 2.0 * ( quat[2] * quat[0] + quat[1] * quat[3] ) ;
            m[2][1] = 2.0 * ( quat[1] * quat[2] - quat[0] * quat[3] ) ;
            m[2][2] = 1 - 2.0 * ( quat[1] * quat[1] + quat[0] * quat[0] ) ;
            m[2][3] = 0.0 ;

            m[3][0] = 0.0 ;
            m[3][1] = 0.0 ;
            m[3][2] = 0.0 ;
            m[3][3] = 1.0 ;

            double x = V[0] * m[0][0] + V[1] * m[1][0] + V[2] * m[2][0] + m[3][0] ;
            double y = V[0] * m[0][1] + V[1] * m[1][1] + V[2] * m[2][1] + m[3][1] ;
            double z = V[0] * m[0][2] + V[1] * m[1][2] + V[2] * m[2][2] + m[3][2] ;
            double w = V[0] * m[0][3] + V[1] * m[1][3] + V[2] * m[2][3] + m[3][3] ;
            return vec3( x / w, y / w, z / w ) ;
        }

        void sort()
        {
            ringmesh_assert( triangles_.size() > 0 ) ;

            std::pair< index_t, bool > default_pair( index_t( -1 ), false ) ;
            sorted_triangles_.resize( 2 * triangles_.size(), default_pair ) ;

            // If there is only one Triangle to sort - nothing to do
            if( triangles_.size() == 1 ) {
                sorted_triangles_[0] = std::pair< index_t, bool >(
                    triangles_[0].surface_index_, true ) ;
                sorted_triangles_[1] = std::pair< index_t, bool >(
                    triangles_[0].surface_index_, false ) ;
                return ;
            }

            // Initialization
            // We start on the plus (true) side of the first Triangle            
            sorted_triangles_[0] = std::pair< index_t, bool >(
                triangles_[0].surface_index_, true ) ;

            // Reference vectors with wich angles will be computed
            vec3 N_ref = triangles_[0].N_ ;
            vec3 B_A_ref = triangles_[0].B_A_ ;
            vec3 Ax_ref = normalize( cross( B_A_ref, N_ref ) ) ;

            // The minus (false) side of the start triangle will the last one encountered
            triangles_[0].angle_ = 2 * M_PI ;
            triangles_[0].side_ = false ;

            for( index_t i = 1; i < triangles_.size(); ++i ) {
                TriangleToSort& cur = triangles_[i] ;
                // Computes the angle RADIANS between the reference and the current
                // triangle 
                double cos = dot( B_A_ref, cur.B_A_ ) ;
                // Remove invalid values
                if( cos < -1 )
                    cos = -1 ;
                else if( cos > 1 ) cos = 1 ;
                cur.angle_ = std::acos( cos ) ;
                // Put the angle between PI and 2PI if necessary
                if( dot( cross( B_A_ref, cur.B_A_ ), Ax_ref ) < 0. ) {
                    cur.angle_ = 2 * M_PI - cur.angle_ ;
                }

                // Get the side of the surface first encountered
                // when rotating in the N_ref direction
                vec3 N_rotate = rotate( Ax_ref, -cur.angle_, cur.N_ ) ;
                cur.side_ = dot( N_rotate, N_ref ) > 0 ? false : true ;
            }

            // Sorts the Surfaces according to the angle
            std::sort( triangles_.begin(), triangles_.end() ) ;

            // Fills the sorted surfaces adding the side
            index_t it = 1 ;
            for( index_t i = 0; i < triangles_.size(); ++i ) {
                TriangleToSort& cur = triangles_[i] ;
                if( triangles_[i].index_ == 0 ) { // The last to add
                    ringmesh_assert( i == triangles_.size() - 1 ) ;
                    sorted_triangles_[it].first = cur.surface_index_ ;
                    sorted_triangles_[it].second = cur.side_ ;
                } else {
                    sorted_triangles_[it].first = cur.surface_index_ ;
                    sorted_triangles_[it].second = cur.side_ ;
                    sorted_triangles_[it + 1].first = cur.surface_index_ ;
                    sorted_triangles_[it + 1].second = !cur.side_ ;
                    it += 2 ;
                }
            }
            // All the surfaces must have been sorted
            ringmesh_assert(
                std::count( sorted_triangles_.begin(), sorted_triangles_.end(),
                    default_pair ) == 0 ) ;
        }

        void clear()
        {
            triangles_.clear() ;
            sorted_triangles_.clear() ;
        }

        /*! Returns the next pair Triangle index (surface) + side
         */
        const std::pair< index_t, bool >& next(
            const std::pair< index_t, bool >& in ) const
        {
            for( index_t i = 0; i < sorted_triangles_.size(); ++i ) {
                if( sorted_triangles_[i] == in ) {
                    if( i == sorted_triangles_.size() - 1 )
                        return sorted_triangles_[sorted_triangles_.size() - 2] ;
                    if( i == 0 ) return sorted_triangles_[1] ;

                    if( sorted_triangles_[i + 1].first
                        == sorted_triangles_[i].first ) {
                        // The next has the same surface id, check its sign
                        if( sorted_triangles_[i + 1].second
                            != sorted_triangles_[i].second ) {
                            return sorted_triangles_[i - 1] ;
                        } else {
                            // Sign is the same
                            return sorted_triangles_[i + 1] ;
                        }
                    } else {
                        ringmesh_assert(
                            sorted_triangles_[i - 1].first
                            == sorted_triangles_[i].first ) ;
                        if( sorted_triangles_[i - 1].second
                            != sorted_triangles_[i].second ) {
                            return sorted_triangles_[i + 1] ;
                        } else {
                            return sorted_triangles_[i - 1] ;
                        }
                    }
                }
            }
            ringmesh_assert_not_reached ;
            return sorted_triangles_.front() ;
        }

    private:
        std::vector< TriangleToSort > triangles_ ;
        // Pairs global triangle identifier (Surface index) and side reached
        std::vector< std::pair< index_t, bool > > sorted_triangles_ ;
    } ;

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
            initialize_border_triangles_from_model_surfaces() ;
            visited_.resize( border_triangles_.size(), false ) ;
        }

        /*!
         * @brief Tries to compute a new line and returns true if one was.
         * @details To use in a while conditional loop, since the number of lines
         * is considered unknown. 
         */
        bool compute_next_line_geometry()
        {
            go_to_next_non_visited_border_triangle() ;

            if( is_not_the_end() ) {
                init_next_line_computation() ;
                compute_line_geometry() ;
                return true ;
            } else {
                return false ;
            }
        }

        /*! 
         * @brief Returns the vertices of the last line built 
         */
        const std::vector< index_t >& vertices() const
        {
            return cur_line_vertices_ ;
        }

        /*!
         * @brief Returns the adjacent surfaces to the last line built
         */
        const std::vector< index_t >& adjacent_surfaces() const
        {
            return cur_line_adjacent_surfaces_ ;
        }

        /*!
         * @brief Returns region building information for the last line built.
         * @details Only computed if the collect_region_info was true at construction.
         */
        const GeoModelRegionFromSurfaces& region_information() const
        {
            return cur_line_region_information_ ;
        }

    private:
        /*!
         * @brief Stores the vertices of a triangle which is on the boundary of a Surface
         */
        struct BorderTriangle {
            BorderTriangle(
                index_t surface,
                index_t facet,
                index_t v0,
                index_t v1,
                index_t v2 )
                :
                    v0_( v0 ),
                    v1_( v1 ),
                    v2_( v2 ),
                    surface_( surface ),
                    facet_( facet )
            {
            }

            /*!
             * @brief Compares two triangles so those sharing an edge follow one another
             */
            bool operator<( const BorderTriangle& rhs ) const
            {
                if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) {
                    return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ ) ;
                }
                if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) {
                    return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ ) ;
                }
                if( surface_ != rhs.surface_ ) {
                    return surface_ < rhs.surface_ ;
                }
                if( facet_ != rhs.facet_ ) {
                    return facet_ < rhs.facet_ ;
                }
                return rhs.v2_ == NO_ID ? false : v2_ < rhs.v2_ ;
            }

            bool same_edge( const BorderTriangle& rhs ) const
            {
                return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ )
                    && std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ ) ;
            }

            /// Indices of the points in the model. 
            /// The edge v0v1 is on the one on the boundary.
            index_t v0_ ;
            index_t v1_ ;
            index_t v2_ ;
            // Index of the surface containing the triangle
            index_t surface_ ;
            // Index of the facet in the surface
            index_t facet_ ;
        } ;

        void compute_line_geometry()
        {
            visit_border_triangles_on_same_edge( cur_border_triangle_ ) ;
            index_t p0 = border_triangles_[cur_border_triangle_].v0_ ;
            index_t p1 = border_triangles_[cur_border_triangle_].v1_ ;
            cur_line_vertices_.push_back( p0 ) ;
            cur_line_vertices_.push_back( p1 ) ;

            get_adjacent_surfaces( cur_border_triangle_,
                cur_line_adjacent_surfaces_ ) ;

            bool backward = false ;
            get_one_line_vertices( backward ) ;
            backward = true ;
            get_one_line_vertices( backward ) ;

            if( collect_region_information_ ) {
                collect_region_information() ;
            }
        }
        /*!
         * @brief Gets the geometry of a line propagating in a given direction
         * @details As long as the adjacent surfaces are the same, the vertices are added
         * belong to the line under construction
         */
        void get_one_line_vertices( bool backward )
        {
            ringmesh_assert( is_not_the_end() ) ;

            index_t start( cur_border_triangle_ ) ;
            ringmesh_assert( start != NO_ID ) ;

            index_t t = get_next_border_triangle( start, backward ) ;
            ringmesh_assert( t != start ) ;

            bool same_surfaces = true ;
            while( same_surfaces && t != start ) {
                ringmesh_assert( t != NO_ID ) ;
                if( !is_visited( t ) ) {
                    std::vector< index_t > cur_adjacent_surfaces ;
                    get_adjacent_surfaces( t, cur_adjacent_surfaces ) ;

                    if( equal_to_line_adjacent_surfaces( cur_adjacent_surfaces ) ) {
                        visit_border_triangles_on_same_edge( t ) ;
                        add_border_triangle_vertices_to_line( t, backward ) ;
                    } else {
                        same_surfaces = false ;
                    }
                } else {
                    same_surfaces = false ;
                }
                t = get_next_border_triangle( t, backward ) ;
            }
        }

        void init_next_line_computation()
        {
            cur_line_vertices_.resize( 0 ) ;
            cur_line_adjacent_surfaces_.resize( 0 ) ;
            cur_line_region_information_.clear() ;
        }

        bool is_visited( index_t i ) const
        {
            return visited_[i] ;
        }

        void go_to_next_non_visited_border_triangle()
        {
            while( is_not_the_end() && is_visited( cur_border_triangle_ ) ) {
                ++cur_border_triangle_ ;
            }
        }

        bool is_not_the_end() const
        {
            return cur_border_triangle_ < border_triangles_.size() ;
        }

        bool have_border_triangles_same_boundary_edge( index_t t0, index_t t1 ) const
        {
            return border_triangles_[t0].same_edge( border_triangles_[t1] ) ;
        }

        /*!
         * @brief Gets triangles sharing the border edge of the current border triangle
         *        and adds them to current line region information 
         */
        void collect_region_information()
        {
            index_t i( cur_border_triangle_ ) ;
            index_t j = i ;
            while( j < border_triangles_.size()
                && have_border_triangles_same_boundary_edge( i, j ) ) {
                cur_line_region_information_.add_triangle(
                    border_triangles_[j].surface_,
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v0_ ),
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v1_ ),
                    geomodel_.mesh.vertices.vertex( border_triangles_[j].v2_ ) ) ;
                j++ ;
            }
        }

        bool equal_to_line_adjacent_surfaces(
            const std::vector< index_t >& input ) const
        {
            if( input.size() != cur_line_adjacent_surfaces_.size() ) {
                return false ;
            } else {
                return std::equal( input.begin(), input.end(),
                    cur_line_adjacent_surfaces_.begin() ) ;
            }
        }

        void add_border_triangle_vertices_to_line(
            index_t triangle_index,
            bool backward )
        {
            const BorderTriangle& border_triangle = border_triangles_[triangle_index] ;
            index_t v0 = border_triangle.v0_ ;
            index_t v1 = border_triangle.v1_ ;

            add_vertices_to_line( v0, v1, !backward ) ;
        }

        void add_vertices_to_line( index_t v0, index_t v1, bool at_the_end )
        {
            index_t to_add = v1 ;
            if( at_the_end ) {
                index_t end_vertex = cur_line_vertices_.back() ;
                if( v1 == end_vertex ) {
                    to_add = v0 ;
                } else {
                    ringmesh_assert( v0 == end_vertex ) ;
                }
                cur_line_vertices_.push_back( to_add ) ;
            } else {
                index_t first_vertex = cur_line_vertices_.front() ;
                if( v1 == first_vertex ) {
                    to_add = v0 ;
                } else {
                    ringmesh_assert( v0 == first_vertex ) ;
                }
                cur_line_vertices_.insert( cur_line_vertices_.begin(), to_add ) ;
            }
        }

        void initialize_border_triangles_from_model_surfaces()
        {
            for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
                const Surface& S = geomodel_.surface( i ) ;
                for( index_t j = 0; j < S.nb_mesh_elements(); ++j ) {
                    for( index_t v = 0; v < S.nb_mesh_element_vertices( j ); ++v ) {
                        if( S.is_on_border( j, v ) ) {
                            index_t vertex = S.model_vertex_id( j, v ) ;
                            index_t next_vertex = S.model_vertex_id( j,
                                S.next_facet_vertex_index( j, v ) ) ;
                            index_t previous_vertex = S.model_vertex_id( j,
                                S.prev_facet_vertex_index( j, v ) ) ;
                            border_triangles_.push_back(
                                BorderTriangle( i, j, vertex, next_vertex,
                                    previous_vertex ) ) ;
                        }
                    }
                }
            }
            std::sort( border_triangles_.begin(), border_triangles_.end() ) ;
        }

        /*!
         * @brief Gets the next BorderTriangle in the same surface
         */
        index_t get_next_border_triangle( index_t from, bool backward ) const
        {
            const BorderTriangle& border_triangle = border_triangles_[from] ;
            const Surface& S = geomodel_.surface( border_triangle.surface_ ) ;

            // Gets the next edge on border in the Surface
            index_t f = border_triangle.facet_ ;
            index_t f_v0 = S.facet_id_from_model( f, border_triangle.v0_ ) ;
            index_t f_v1 = S.facet_id_from_model( f, border_triangle.v1_ ) ;
            ringmesh_assert( f_v0 != NO_ID && f_v1 != NO_ID ) ;

            index_t next_f = NO_ID ;
            index_t next_f_v0 = NO_ID ;
            index_t next_f_v1 = NO_ID ;

            if( !backward ) {
                S.next_on_border( f, f_v0, f_v1, next_f, next_f_v0, next_f_v1 ) ;
            } else {
                S.next_on_border( f, f_v1, f_v0, next_f, next_f_v0, next_f_v1 ) ;
            }

            // Finds the BorderTriangle that is corresponding to this
            // It must exist and there is only one
            BorderTriangle bait( border_triangle.surface_, next_f,
                S.model_vertex_id( next_f, next_f_v0 ),
                S.model_vertex_id( next_f, next_f_v1 ), NO_ID ) ;
            index_t result(
                std::lower_bound( border_triangles_.begin(), border_triangles_.end(),
                    bait ) - border_triangles_.begin() ) ;

            ringmesh_assert( result < border_triangles_.size() ) ;
            return result ;
        }

        /*!
         * @brief Marks as visited all BorderTriangles whose first edge is equal to i's first edge
         */
        void visit_border_triangles_on_same_edge( index_t i )
        {
            index_t j = i ;
            while( j < border_triangles_.size()
                && border_triangles_[i].same_edge( border_triangles_[j] ) ) {
                visited_[j] = true ;
                j++ ;
            }
            index_t k = i - 1 ;
            while( k != NO_ID
                && border_triangles_[i].same_edge( border_triangles_[k] ) ) {
                visited_[k] = true ;
                k-- ;
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
            index_t j = i ;
            while( j < border_triangles_.size()
                && border_triangles_[i].same_edge( border_triangles_[j] ) ) {
                adjacent_surfaces.push_back( border_triangles_[j].surface_ ) ;
                j++ ;
            }

            index_t k = i - 1 ;
            while( k != NO_ID
                && border_triangles_[i].same_edge( border_triangles_[k] ) ) {
                adjacent_surfaces.push_back( border_triangles_[k].surface_ ) ;
                k-- ;
            }
            std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() ) ;
        }

    private:
        const GeoModel& geomodel_ ;
        bool collect_region_information_ ;
        // All the triangles on a boundary of all the Surfaces of the GeoModel
        std::vector< BorderTriangle > border_triangles_ ;
        // Internal use to flag the visited border_triangles when computing the Lines
        std::vector< bool > visited_ ;

        // Currently computed line information
        index_t cur_border_triangle_ ;
        std::vector< index_t > cur_line_vertices_ ;
        std::vector< index_t > cur_line_adjacent_surfaces_ ;
        RINGMesh::GeoModelRegionFromSurfaces cur_line_region_information_ ;
    } ;

    void reorder_line_vertices_to_start_at_corner(
        const GeoModel& geomodel,
        std::vector< index_t >& line_vertices )
    {
        if( geomodel.nb_corners() == 0 ) {
            // Maybe should throw an assertion, but I am not sure [JP]
            return ;
        }
        for( index_t i = 1; i + 1 < line_vertices.size(); ++i ) {
            gme_t corner = find_corner( geomodel, line_vertices[i] ) ;
            if( corner.is_defined() ) {
                std::vector< index_t > shuffled_vertices( line_vertices.size() ) ;
                std::copy( line_vertices.begin() + i, line_vertices.end(),
                    shuffled_vertices.begin() ) ;
                index_t nb_copied(
                    line_vertices.end() - line_vertices.begin() - i ) ;
                std::copy( line_vertices.begin() + 1, line_vertices.begin() + i + 1,
                    shuffled_vertices.begin() + nb_copied ) ;
                line_vertices = shuffled_vertices ;
                break ;
            }
        }
    }

    bool is_surface_mesh( const GEO::Mesh& mesh )
    {
        return mesh.facets.nb() != 0 ;
    }

    bool is_volume_mesh( const GEO::Mesh& mesh )
    {
        return mesh.cells.nb() != 0 ;
    }

    /*! Delete all GeoModelRegionFromSurfaces owned by the builder
     */
    GeoModelBuilder::~GeoModelBuilder()
    {
        for( index_t i = 0; i < regions_info_.size(); ++i ) {
            delete regions_info_[i] ;
        }
    }

    void GeoModelBuilder::copy_meshes( const GeoModel& geomodel )
    {
        for( index_t t = GME::CORNER; t <= GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            copy_meshes( geomodel, T ) ;
        }
    }

    void GeoModelBuilder::copy_meshes( const GeoModel& from, GME::TYPE entity_type )
    {
        for( index_t i = 0; i < model().nb_entities( entity_type ); ++i ) {
            const GeoModelMeshEntity& from_E = from.mesh_entity( entity_type,
                i ) ;
            assign_mesh_to_entity( from_E.mesh_, gme_t( entity_type, i ) ) ;
        }
    }

    void GeoModelBuilder::assign_mesh_to_entity( const Mesh& mesh, GME::gme_t to )
    {
        GeoModelMeshEntity& E = mesh_entity( to ) ;
        E.unbind_attributes() ;
        E.mesh_.copy( mesh, true, GEO::MESH_ALL_ELEMENTS ) ;
        E.bind_attributes() ;
        // Nothing else to do ? To test [JP]
    }

    /*!
     * @brief Finds or creates a corner at given coordinates.
     * @param[in] point Geometric location of the Corner
     * @return Index of the Corner
     */
    gme_t GeoModelBuilder::find_or_create_corner( const vec3& point )
    {
        gme_t result = find_corner( model(), point ) ;
        if( !result.is_defined() ) {
            result = create_entity( GME::CORNER ) ;
            set_corner( result.index, point ) ;
        }
        return result ;
    }

    gme_t GeoModelBuilder::find_or_create_corner( index_t model_point_id )
    {
        gme_t result = find_corner( model(), model_point_id ) ;
        if( !result.is_defined() ) {
            result = create_entity( GME::CORNER ) ;
            set_corner( result.index, model_point_id ) ;
        }
        return result ;
    }

    /*!
     * @brief Finds or creates a line
     * @param[in] vertices Coordinates of the vertices of the line
     * @return Index of the Line
     */
    gme_t GeoModelBuilder::find_or_create_line( const std::vector< vec3 >& vertices )
    {
        gme_t result ;
        for( index_t i = 0; i < model().nb_lines(); ++i ) {
            if( line_equal( model().line( i ), vertices ) ) {
                result = model().line( i ).gme_id() ;
            }
        }
        if( !result.is_defined() ) {
            result = create_entity( GME::LINE ) ;
            set_line( result.index, vertices ) ;

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_entity_boundary( result,
                find_or_create_corner( vertices.front() ) ) ;
            add_entity_boundary( result,
                find_or_create_corner( vertices.back() ) ) ;
        }
        return result ;
    }

    /*!
     * @brief Finds or creates a line knowing its topological adjacencies
     */
    gme_t GeoModelBuilder::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        gme_t first_corner,
        gme_t second_corner )
    {
        for( index_t i = 0; i < model().nb_lines(); ++i ) {
            const Line& L = model().line( i ) ;
            gme_t c0 = L.boundary_gme( 0 ) ;
            gme_t c1 = L.boundary_gme( 1 ) ;

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) ) {
                std::vector< index_t > cur_adjacent_surfaces ;
                get_sorted_incident_surfaces( L, cur_adjacent_surfaces ) ;
                if( cur_adjacent_surfaces.size() == sorted_adjacent_surfaces.size()
                    && std::equal( cur_adjacent_surfaces.begin(),
                        cur_adjacent_surfaces.end(),
                        sorted_adjacent_surfaces.begin() ) ) {
                    return L.gme_id() ;
                }
            }
        }
        return create_entity( GME::LINE ) ;
    }

    /*!
     * @brief Sets the geometrical position of a vertex
     * @param[in] t Entity index
     * @param[in] v Index of the vertex to modify
     * @param[in] point New coordinates
     * @param[in] update If true, all the vertices sharing the same geometrical position
     *               in the GeoModel have their position updated, if false they
     *               are not.
     * @warning Be careful with this update parameter, it is a very nice source of nasty bugs
     */
    void GeoModelBuilder::set_entity_vertex(
        const GME::gme_t& t,
        index_t v,
        const vec3& point,
        bool update )
    {
        GeoModelMeshEntity& E = mesh_entity( t ) ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        if( update ) {
            model().mesh.vertices.update_point( E.model_vertex_id( v ), point ) ;
        } else {
            MeshBuilder builder( E.mesh_ ) ;
            builder.vertex( v ) = point ;
        }
    }

    /*!
     * @brief Sets the geometrical position of a vertex from a model vertex
     * @details Sets also both mapping from (GeoModelMeshVertices::unique2bme)
     *          and to (model_vertex_id_) the model vertex.
     * @param[in] entity_id Entity index
     * @param[in] v Index of the vertex to modify
     * @param[in] model_vertex Index in GeoModelMeshVertices of the vertex giving
     *                     the new position
     */
    void GeoModelBuilder::set_entity_vertex(
        const gme_t& entity_id,
        index_t v,
        index_t model_vertex )
    {
        set_entity_vertex( entity_id, v,
            model().mesh.vertices.vertex( model_vertex ), false ) ;

        GeoModelMeshEntity& E = mesh_entity( entity_id ) ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        E.model_vertex_id_[v] = model_vertex ;
        model().mesh.vertices.add_to_bme( model_vertex, GMEVertex( entity_id, v ) ) ;
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     * @param[in] id Entity index
     * @param[in] points Geometric positions of the vertices to add
     * @param[in] clear If true the mesh is cleared, keeping its attributes
     */
    void GeoModelBuilder::set_entity_vertices(
        const gme_t& id,
        const std::vector< vec3 >& points,
        bool clear )
    {
        GeoModelMeshEntity& E = mesh_entity( id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder.clear( true, true ) ;
        }
        if( !points.empty() ) {
            index_t start = builder.create_vertices( points.size() ) ;
            GEO::Memory::copy( builder.point_ptr( start ), points.data()->data(),
                3 * sizeof(double) * points.size() ) ;
        }
    }

    /*!
     * @brief Creates new vertices to the mesh
     * @param[in] id Entity index
     * @param[in] nb_vertices Number of vertices to create
     * @return the first vertex index created
     */
    index_t GeoModelBuilder::create_entity_vertices(
        const GME::gme_t& id,
        index_t nb_vertices )
    {
        GeoModelMeshEntity& E = mesh_entity( id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        builder.delete_vertex_colocater() ;
        return builder.create_vertices( nb_vertices ) ;
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     *
     * @param[in] id Entity index
     * @param[in] model_vertices Geometric positions of the vertices to add
     * @param[in] clear If true the mesh if cleared, keeping its attributes
     */
    void GeoModelBuilder::set_entity_vertices(
        const gme_t& entity_id,
        const std::vector< index_t >& model_vertices,
        bool clear )
    {
        GeoModelMeshEntity& E = mesh_entity( entity_id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder.clear( true, true ) ;
        }
        index_t start = builder.create_vertices( model_vertices.size() ) ;
        for( index_t v = 0; v < model_vertices.size(); v++ ) {
            set_entity_vertex( entity_id, start + v, model_vertices[v] ) ;
        }
    }

    /*!
     * @brief Sets the geometric location of a Corner
     *
     * @param[in] corner_id Index of the corner
     * @param[in] point Coordinates of the vertex
     */
    void GeoModelBuilder::set_corner( index_t corner_id, const vec3& point )
    {
        set_entity_vertex( gme_t( GME::CORNER, corner_id ), 0, point, false ) ;
    }

    /*!
     * @brief Sets one Line points
     *
     * @param[in] line_id Line index
     * @param[in] vertices Coordinates of the vertices on the line
     */
    void GeoModelBuilder::set_line(
        index_t line_id,
        const std::vector< vec3 >& vertices )
    {
        set_entity_vertices( gme_t( GME::LINE, line_id ), vertices, true ) ;

        GeoModelMeshEntity& E = mesh_entity( GME::LINE, line_id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            builder.create_edge( e - 1, e ) ;
        }
    }

    /*!
     * @brief Sets the points and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] points Coordinates of the vertices
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void GeoModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        set_entity_vertices( gme_t( GME::SURFACE, surface_id ), points, true ) ;
        assign_surface_mesh_facets( surface_id, facets, facet_ptr ) ;
    }

    /*!
     * @brief Set the points and tetras for a region
     *
     * @param[in] region_id Index of the regions
     * @param[in] points Coordinates of the vertices
     * @param[in] tetras Indices in the vertices vector to build tetras
     */
    void GeoModelBuilder::set_region_geometry(
        index_t region_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& tetras )
    {
        set_entity_vertices( gme_t( GME::REGION, region_id ), points, true ) ;
        assign_region_tet_mesh( region_id, tetras ) ;
    }

    /*!
     * @brief Add a point to the GeoModel and not to one of its entities
     * @details To use when adding the points to the model before building its entities
     */
    index_t GeoModelBuilder::add_unique_vertex( const vec3& p )
    {
        return model().mesh.vertices.add_vertex( p ) ;
    }

    /*!
     * @brief Sets the vertex for a Corner. Store the info in the geomodel vertices
     *
     * @param[in] corner_id Index of the corner
     * @param[in] unique_vertex Index of the vertex in the model
     */
    void GeoModelBuilder::set_corner( index_t corner_id, index_t model_vertex_id )
    {
        set_entity_vertex( gme_t( GME::CORNER, corner_id ), 0, model_vertex_id ) ;
    }

    /*!
     * @brief Sets one Line vertices. Store the info in the geomodel vertices
     *
     * @param[in] id Line index
     * @param[in] unique_vertices Indices in the model of the unique vertices with which to build the Line
     */
    void GeoModelBuilder::set_line(
        index_t line_id,
        const std::vector< index_t >& unique_vertices )
    {
        bool clear_vertices = false ;
        GeoModelMeshEntity& E = mesh_entity( GME::LINE, line_id ) ;
        ringmesh_assert( E.nb_vertices() == 0 ) ; // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_entity_vertices( gme_t( GME::LINE, line_id ), unique_vertices,
            clear_vertices ) ;

        MeshBuilder builder( E.mesh_ ) ;
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            builder.create_edge( e - 1, e ) ;
        }
    }

    /*!
     * @brief Sest the vertices and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] model_vertex_ids Indices of unique vertices in the GeoModel
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void GeoModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& model_vertex_ids,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        set_entity_vertices( gme_t( GME::SURFACE, surface_id ), model_vertex_ids,
            false ) ;
        assign_surface_mesh_facets( surface_id, facets, facet_ptr ) ;
    }

    /*!
     * @brief Sets the facets of a surface
     * @param[in] surface_id Index of the surface
     * @param[in] facets Indices of the model vertices defining the facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void GeoModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_facets( facets ) ;
        get_entity_vertices_and_update_corners( new_facets, vertices ) ;
        set_surface_geometry( surface_id, vertices, new_facets, facet_ptr ) ;
    }

    void GeoModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& triangle_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_triangle_corners( triangle_corners ) ;
        get_entity_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_entity_vertices( gme_t( GME::SURFACE, surface_id ), vertices, false ) ;
        assign_surface_triangle_mesh( surface_id, new_triangle_corners ) ;
    }

    void GeoModelBuilder::set_surface_geometry_with_adjacencies(
        index_t surface_id,
        const std::vector< index_t >& triangle_corners,
        const std::vector< index_t >& adjacent_triangles )
    {
        /// @todo Reorganize to remove duplicated code in the class
        std::vector< index_t > vertices ;
        std::vector< index_t > new_triangle_corners( triangle_corners ) ;
        get_entity_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_entity_vertices( gme_t( GME::SURFACE, surface_id ), vertices, false ) ;

        assign_surface_triangle_mesh( surface_id, new_triangle_corners,
            adjacent_triangles ) ;
    }

    void GeoModelBuilder::set_surface_element_geometry(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = mesh_entity( GME::SURFACE, surface_id ) ;
        MeshBuilder builder( E.mesh_ ) ;

        for( index_t facet_vertex = 0; facet_vertex < corners.size(); facet_vertex++ ) {
            builder.set_facet_vertex( facet_id, facet_vertex, corners[facet_vertex] ) ;
        }
    }

    void GeoModelBuilder::set_surface_element_adjacency(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& adjacents )
    {
        GeoModelMeshEntity& E = mesh_entity( GME::SURFACE, surface_id ) ;
        MeshBuilder builder( E.mesh_ ) ;

        for( index_t facet_edge = 0; facet_edge < adjacents.size(); facet_edge++ ) {
            builder.set_facet_adjacent( facet_id, facet_edge, adjacents[facet_edge] ) ;
        }
    }

    void GeoModelBuilder::set_region_geometry(
        index_t region_id,
        const std::vector< index_t >& tet_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_tet_corners( tet_corners ) ;
        get_entity_vertices_and_update_corners( new_tet_corners, vertices ) ;

        set_entity_vertices( gme_t( GME::REGION, region_id ), vertices, false ) ;
        assign_region_tet_mesh( region_id, new_tet_corners ) ;
    }

    void GeoModelBuilder::set_region_element_geometry(
        index_t region_id,
        index_t cell_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = mesh_entity( GME::REGION, region_id ) ;
        MeshBuilder builder( E.mesh_ ) ;

        for( index_t cell_vertex = 0; cell_vertex < corners.size(); cell_vertex++ ) {
            builder.set_cell_vertex( cell_id, cell_vertex, corners[cell_vertex] ) ;
        }
    }

    /*!
     * @brief Creates new cells in the mesh
     * @param[in] region_id Entity index
     * @param[in] type Type of cell
     * @param[in] nb_cells Number of cells to creates
     * @return the index of the first created cell
     */
    index_t GeoModelBuilder::create_region_cells(
        index_t region_id,
        GEO::MeshCellType type,
        index_t nb_cells )
    {
        GeoModelMeshEntity& E = mesh_entity( GME::REGION, region_id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        builder.delete_vertex_colocater() ;
        builder.delete_cell_colocater() ;
        builder.delete_cell_aabb() ;

        return builder.create_cells( nb_cells, type ) ;
    }

    index_t GeoModelBuilder::create_region_cell(
        index_t region_id,
        GEO::MeshCellType type,
        const std::vector< index_t >& vertex_indices )
    {
        index_t cell_id = create_region_cells( region_id, type, 1 ) ;
        set_region_element_geometry( region_id, cell_id, vertex_indices ) ;
        return cell_id ;
    }

    index_t GeoModelBuilder::create_surface_facet(
        index_t surface_id,
        const GEO::vector< index_t >& vertex_indices )
    {
        GeoModelMeshEntity& E = mesh_entity( GME::SURFACE, surface_id ) ;
        MeshBuilder builder( E.mesh_ ) ;
        builder.delete_vertex_colocater() ;
        builder.delete_facet_colocater() ;
        builder.delete_facet_aabb() ;

        return builder.create_facet_polygon( vertex_indices ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices )
    {
        Mesh& M = mesh_entity( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.nb_vertices() > 0 ) ;
        MeshBuilder builder( M ) ;
        builder.assign_facet_triangle_mesh( triangle_vertices, true ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices,
        const std::vector< index_t >& adjacent_triangles )
    {
        Mesh& M = mesh_entity( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.nb_vertices() > 0 ) ;
        MeshBuilder builder( M ) ;
        builder.assign_facet_triangle_mesh( triangle_vertices, true ) ;

        ringmesh_assert( adjacent_triangles.size() == M.nb_facet_corners() ) ;
        for( index_t i = 0; i < adjacent_triangles.size(); i++ ) {
            builder.set_facet_corners_adjacent( i, adjacent_triangles[i] ) ;
        }
    }

    void GeoModelBuilder::set_surface_facet_adjacencies(
        index_t surface_id,
        const std::vector< index_t >& facets_id,
        const std::vector< index_t >& edges_id,
        const std::vector< index_t >& adjacent_triangles )
    {
        Mesh& M = mesh_entity( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.nb_vertices() > 0 ) ;
        MeshBuilder builder( M ) ;
        ringmesh_assert( facets_id.size() == edges_id.size() &&
            facets_id.size() == adjacent_triangles.size() ) ;
        for( index_t i = 0; i < facets_id.size(); ++i ) {
            builder.set_facet_adjacent( facets_id[i], edges_id[i],
                adjacent_triangles[i] ) ;
        }
    }

    void GeoModelBuilder::assign_surface_mesh_facets(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        Mesh& M = mesh_entity( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.nb_vertices() > 0 ) ;
        MeshBuilder builder( M ) ;
        builder.create_facet_polygons( facets, facet_ptr ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_region_tet_mesh(
        index_t region_id,
        const std::vector< index_t >& tet_vertices )
    {
        Mesh& M = mesh_entity( GME::REGION, region_id ).mesh_ ;
        ringmesh_assert( M.nb_vertices() > 0 ) ;
        MeshBuilder builder( M ) ;
        builder.assign_cell_tet_mesh( tet_vertices, true ) ;
        builder.connect_cells() ;
    }

    /*!
     * @brief Computes and sets the adjacencies between the facets
     * @details The adjacent facet is given for each vertex of each facet for the edge
     * starting at this vertex.
     * If there is no neighbor inside the same Surface adjacent is set to NO_ADJACENT
     *
     * @param[in] surface_id Index of the surface
     */
    void GeoModelBuilder::compute_surface_adjacencies( index_t surface_id )
    {
        Surface& S = dynamic_cast< Surface& >( *model().surfaces_[surface_id] ) ;

        index_t nb_facets = S.nb_mesh_elements() ;
        ringmesh_assert( nb_facets > 0 ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( nb_facets - 1 ), NO_ID ) ;

        index_t nb_vertices = S.nb_vertices() ;

        ///@todo Change representation of vertex_to_facets (without vector of vectors)

        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > empty_vector ;
        empty_vector.reserve( 10 ) ;
        std::vector< std::vector< index_t > > vertex_to_facets( nb_vertices,
            empty_vector ) ;

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_mesh_element_vertices( f ); v++ ) {
                vertex_to_facets[S.mesh_element_vertex_index( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ) {
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_mesh_element_vertices( f ); v++ ) {
                index_t cur = S.mesh_element_vertex_index( f, v ) ;
                index_t prev = S.mesh_element_vertex_index( f,
                    S.prev_facet_vertex_index( f, v ) ) ;

                const std::vector< index_t >& f_prev = vertex_to_facets[prev] ;
                const std::vector< index_t >& f_cur = vertex_to_facets[cur] ;

                std::vector< index_t > inter(
                    std::min( f_prev.size(), f_cur.size() ) ) ;
                index_t end = narrow_cast< index_t >(
                    std::set_intersection( f_prev.begin(), f_prev.end(),
                        f_cur.begin(), f_cur.end(), inter.begin() )
                        - inter.begin() ) ;

                if( end == 2 ) {
                    // There is one neighbor
                    index_t f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[S.facet_begin( f )
                        + S.prev_facet_vertex_index( f, v )] = f2 ;
                } else {
                    ringmesh_assert( end == 1 ) ;
                }
            }
        }

        ringmesh_assert( adjacent.size() == S.mesh_.nb_facet_corners() ) ;
        for( index_t i = 0; i < adjacent.size(); i++ ) {
            MeshBuilder builder( S.mesh_ ) ;
            builder.set_facet_corners_adjacent( i, adjacent[i] ) ;
        }
    }

    void GeoModelBuilder::compute_region_adjacencies( index_t region_id )
    {
        Mesh& mesh = mesh_entity( GME::gme_t( GME::REGION, region_id ) ).mesh_ ;
        MeshBuilder builder( mesh ) ;
        builder.connect_cells() ;
    }

    /*!
     * @brief Triangulate a surface
     * @details Compute the triangulation of a set of vertices using
     * a Restricted Voronoi Diagram of the vertices of the surface \p surface_in
     * @param[in] surface_in the surface used to compte the RVD
     * @param[in] surface_out the index of the triangulated surface
     */
    void GeoModelBuilder::triangulate_surface(
        const RINGMesh::Surface& surface_in,
        index_t surface_out )
    {
        Mesh& mesh = mesh_entity( GME::gme_t( GME::SURFACE, surface_out ) ).mesh_ ;
        MeshBuilder builder( mesh ) ;
        builder.triangulate( surface_in.mesh_ ) ;
    }

    /*!
     * Finds duplicate vertex or creates it
     */
    index_t GeoModelBuilder::find_or_create_duplicate_vertex(
        const GME::gme_t& E_id,
        index_t model_vertex_id,
        index_t surface_vertex_id )
    {
        GeoModelMeshEntity& E = mesh_entity( E_id ) ;

        const std::vector< GMEVertex >& vbme = model().mesh.vertices.gme_vertices(
            model_vertex_id ) ;
        index_t duplicate = NO_ID ;
        for( index_t i = 0; i < vbme.size(); ++i ) {
            if( vbme[i].gme_id == E.gme_id() ) {
                if( vbme[i].v_id != surface_vertex_id ) {
                    duplicate = vbme[i].v_id ;
                }
            }
        }
        if( duplicate == NO_ID ) {
            // Duplicate the vertex in the surface
            MeshBuilder builder( E.mesh_ ) ;
            duplicate = builder.create_vertex(
                model().mesh.vertices.vertex( model_vertex_id ).data() ) ;

            // Set its model vertex index
            ringmesh_assert( duplicate < E.nb_vertices() ) ;
            E.model_vertex_id_[duplicate] = model_vertex_id ;

            // Add the mapping from in the model vertices. Should we do this one ?
            model().mesh.vertices.add_to_bme( model_vertex_id,
                GMEVertex( E.gme_id(), duplicate ) ) ;
        }
        return duplicate ;
    }

    /*
     * @brief Resets the adjacencies for all Surface facets adjacent to the Line
     * @pre All the edges of the Line are edges of at least one facet of the Surface
     */
    void GeoModelBuilder::disconnect_surface_facets_along_line_edges(
        index_t surface_id,
        index_t line_id )
    {
        Surface& S = dynamic_cast< Surface& >( entity( gme_t( GME::SURFACE, surface_id ) ) ) ;
        const GMME& L = mesh_entity( gme_t( GME::LINE, line_id ) ) ;
        const ColocaterANN& ann = S.facet_colocater_ann() ;
        MeshBuilder builder( S.mesh_ ) ;
        for( index_t i = 0; i + 1 < L.nb_vertices(); ++i ) {
            index_t p0 = L.model_vertex_id( i ) ;
            index_t p1 = L.model_vertex_id( i + 1 ) ;

            index_t f = NO_ID ;
            index_t v = NO_ID ;
            bool found = find_facet_and_edge( ann, S, p0, p1, f, v ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && f != NO_ID && v != NO_ID ) ;

            index_t f2 = S.facet_adjacent_index( f, v ) ;
            if( f2 != NO_ID ) {
                index_t v2 = NO_ID ;
                // Get the edge in facet f2 matching model indices p0-p1
                S.oriented_edge_from_model_vertex_ids( p0, p1, f2, v2 ) ;
                if( v2 == NO_ID ) {
                    S.oriented_edge_from_model_vertex_ids( p1, p0, f2, v2 ) ;
                    ringmesh_assert( v2 != NO_ID ) ;
                }
                builder.set_facet_adjacent( f, v, NO_ID ) ;
                builder.set_facet_adjacent( f2, v2, NO_ID ) ;
            }
        }
    }

    // Internal function - not very clean I know [JP]
    void find_surface_vertices_adjacent_to_line_first_edge(
        const Surface& S,
        const Line& L,
        index_t& facet_index,
        index_t& surface_vertex_0,
        index_t& surface_vertex_1 )
    {
        // Reset outputs
        facet_index = NO_ID ;
        surface_vertex_0 = NO_ID ;
        surface_vertex_1 = NO_ID ;

        const ColocaterANN& ann = S.facet_colocater_ann() ;
        index_t p0 = L.model_vertex_id( 0 ) ;
        index_t p1 = L.model_vertex_id( 1 ) ;

        index_t v( NO_ID ) ;
        bool found = find_facet_and_edge( ann, S, p0, p1, facet_index, v ) ;
        ringmesh_unused( found ) ;
        ringmesh_assert( found && facet_index != NO_ID && v != NO_ID ) ;

        surface_vertex_0 = S.mesh_element_vertex_index( facet_index, v ) ;
        surface_vertex_1 = S.mesh_element_vertex_index( facet_index,
            S.next_facet_vertex_index( facet_index, v ) ) ;
    }

    /*!
     * @brief Duplicates the surface vertices along the fake boundary
     * (NO_ID adjacencies but shared vertices) and duplicate  the vertices
     * @note Bad written code - error prone
     * @todo Rewrite 
     */
    void GeoModelBuilder::duplicate_surface_vertices_along_line(
        index_t surface_id, index_t line_id )
    {
        gme_t S_id( GME::SURFACE, surface_id ) ;
        Surface& S = dynamic_cast< Surface& >( entity( S_id ) ) ;
        const Line& L = dynamic_cast< Line& >( mesh_entity( gme_t( GME::LINE, line_id ) ) ) ;
        const Corner& line_first_corner = dynamic_cast< const Corner& >( L.boundary(
            0 ) ) ;
        const Corner& line_second_corner = dynamic_cast< const Corner& >( L.boundary(
            1 ) ) ;

        // Surface vertex indices of the points along the line
        index_t id0( NO_ID ) ;
        index_t id1( NO_ID ) ;
        index_t f( NO_ID ) ;
        find_surface_vertices_adjacent_to_line_first_edge( S, L, f, id0, id1 ) ;

        // Backup the starting vertex in the Surface
        index_t first_corner_surf_id = id0 ;

        while( S.model_vertex_id( id1 ) != line_second_corner.model_vertex_id() ) {
            // Get the next facet and next triangle on this boundary
            // Same algorithm than in determine_line_vertices function
            index_t next_f( NO_ID ) ;
            index_t id1_in_next( NO_ID ) ;
            index_t next_id1_in_next( NO_ID ) ;

            S.next_on_border( f, S.vertex_index_in_facet( f, id0 ),
                S.vertex_index_in_facet( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;
            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.mesh_element_vertex_index( next_f, next_id1_in_next ) ;

            // Duplicate the vertex at id1
            // After having determined the next 1 we can probably get both at the same time
            // but I am lazy, and we must be careful not to break next_on_border function (Jeanne)
            std::vector< index_t > facets_around_id1 ;
            S.facets_around_vertex( id1, facets_around_id1, false, f ) ;

            index_t new_id1 = find_or_create_duplicate_vertex( S_id,
                S.model_vertex_id( id1 ), id1 ) ;

            // Update vertex index in facets
            update_facet_corner( S, facets_around_id1, id1, new_id1 ) ;

            // Update
            f = next_f ;
            id0 = new_id1 ;
            id1 = next_id1 ;
        }
        // Point where the process ended on the Surface
        index_t second_corner_surf_id = id1 ;

        // Take care of the corners
        bool duplicate_first_corner = is_corner_to_duplicate( S.model(),
            line_first_corner.index() ) ;
        bool duplicate_second_corner = is_corner_to_duplicate( S.model(),
            line_second_corner.index() ) ;

        if( duplicate_first_corner || duplicate_second_corner ) {
            index_t to_duplicate_model_index ;
            index_t to_duplicate_surface_index ;

            if( duplicate_first_corner ) {
                to_duplicate_model_index = line_first_corner.model_vertex_id() ;
                to_duplicate_surface_index = first_corner_surf_id ;
            } else if( duplicate_second_corner ) {
                to_duplicate_model_index = line_second_corner.model_vertex_id() ;
                to_duplicate_surface_index = second_corner_surf_id ;
            }

            index_t duplicated_surface_vertex = find_or_create_duplicate_vertex( S_id,
                to_duplicate_model_index, to_duplicate_surface_index ) ;

            std::vector< index_t > facets_around_to_duplicate ;
            S.facets_around_vertex( to_duplicate_surface_index,
                facets_around_to_duplicate, false ) ;
            update_facet_corner( S, facets_around_to_duplicate,
                to_duplicate_surface_index, duplicated_surface_vertex ) ;
        } else {
            // If both are duplicated, precondition is broken: line L cut S into 2 parts
            ringmesh_assert( !duplicate_first_corner || !duplicate_second_corner ) ;
        }
    }

    /*!
     * @brief Cuts a Surface along a Line assuming that the edges of the Line are edges of the Surface
     * @pre Surface is not already cut. Line L does not cut the Surface S into 2 connected components.
     * @todo Add a test for this function.
     */
    void GeoModelBuilder::cut_surface_by_line( index_t surface_id, index_t line_id )
    {
        /// @todo Replace the use of the model vertices by only a colocater
        /// of the surface vertices and the line vertices
        bool model_vertices_initialized = model().mesh.vertices.is_initialized() ;
        if( !model_vertices_initialized ) {
            model().mesh.vertices.test_and_initialize() ;
        }

        disconnect_surface_facets_along_line_edges( surface_id, line_id ) ;
        duplicate_surface_vertices_along_line( surface_id, line_id ) ;

        if( !model_vertices_initialized ) {
             model().mesh.vertices.clear() ;
        }
    }

    void GeoModelBuilder::end_model()
    {
        if( model().name() == "" ) {
            set_model_name( "model_default_name" ) ;
        }
        // Get out if the model has no surface
        if( model().nb_surfaces() == 0 ) {
            print_geomodel( model() ) ;
            throw RINGMeshException( "GeoModel", "The GeoModel has no surface" ) ;
        }

        model().init_global_model_entity_access() ;
        complete_entity_connectivity() ;

        // Fill geological feature if they are missing
        for( index_t i = 0; i < model().nb_entities( GME::ALL_TYPES ); ++i ) {
            GME& E = entity( gme_t( GME::ALL_TYPES, i ) ) ;
            if( !E.has_geological_feature() ) {
                if( E.has_parent() && E.parent().has_geological_feature() ) {
                    set_entity_geol_feature( E.gme_id(),
                        E.parent().geological_feature() ) ;
                } else if( E.nb_children() > 0
                    && E.child( 0 ).has_geological_feature() ) {
                    set_entity_geol_feature( E.gme_id(),
                        E.child( 0 ).geological_feature() ) ;
                }
            }
        }

        // Deliberate clear of the model vertices used for model building
        model().mesh.vertices.clear() ;
    }

    void GeoModelBuilder::recompute_geomodel_mesh()
    {
        model().mesh.vertices.clear() ;
        model().mesh.vertices.test_and_initialize() ;
    }

    bool GeoModelBuilder::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces line_computer( model(),
            options_.compute_regions_brep ) ;

        bool new_line_was_built = true ;
        while( new_line_was_built ) {
            new_line_was_built = line_computer.compute_next_line_geometry() ;

            // I know this is a copy - but should'nt be too big [JP]
            std::vector< index_t > vertices = line_computer.vertices() ;
            bool is_line_closed = vertices.front() == vertices.back() ;
            if( is_line_closed ) {
                // Vertices can begin and end at any vertex
                reorder_line_vertices_to_start_at_corner( model(), vertices ) ;
            }

            GME::gme_t first_corner = find_or_create_corner( vertices.front() ) ;
            GME::gme_t second_corner = find_or_create_corner( vertices.back() ) ;
            const std::vector< index_t >& adjacent_surfaces =
                line_computer.adjacent_surfaces() ;

            index_t backup_nb_lines = model().nb_lines() ;
            gme_t line_index = find_or_create_line( adjacent_surfaces, first_corner,
                second_corner ) ;

            bool created_line = model().nb_lines() != backup_nb_lines ;
            if( created_line ) {
                set_line( line_index.index, vertices ) ;

                for( index_t j = 0; j < adjacent_surfaces.size(); ++j ) {
                    GME::gme_t surface_id( GME::SURFACE, adjacent_surfaces[j] ) ;
                    add_entity_in_boundary( line_index, surface_id ) ;
                }
                add_entity_boundary( line_index, first_corner ) ;
                add_entity_boundary( line_index, second_corner ) ;

                // If the plan is to then build_regions, get the information
                if( options_.compute_regions_brep ) {
                    regions_info_.push_back(
                        new GeoModelRegionFromSurfaces(
                            line_computer.region_information() ) ) ;
                }
            } else {
                bool same_geometry = line_equal( model().line( line_index.index ),
                    vertices ) ;
                if( !same_geometry ) {
                    set_line( line_index.index, vertices ) ;
                }
            }
        }
        return true ;
    }

    bool GeoModelBuilder::build_brep_regions_from_surfaces()
    {
        ringmesh_assert( model().nb_lines() == regions_info_.size() ) ;

        // Complete boundary information for surfaces
        // to compute volumetric regions
        fill_entities_boundaries( GME::SURFACE ) ;

        // Sort surfaces around the contacts
        for( index_t i = 0; i < regions_info_.size(); ++i ) {
            regions_info_[i]->sort() ;
        }

        if( model().nb_surfaces() == 1 ) {
            if( model().nb_lines() != 0 ) {
                GEO::Logger::err( "GeoModel" )
                    << "The unique surface provided to build the model has boundaries "
                    << std::endl ;
                return false ;
            } else {
                /// If there is only one surface, its inside is set to be 
                /// the + side. No further check.
                bool inside = true ;
                gme_t surface_id( GME::SURFACE, 0 ) ;
                // Create the region - set the surface on its boundaries
                gme_t region_id = create_entity( GME::REGION ) ;
                add_entity_boundary( region_id, surface_id, inside ) ;

                // Set universe boundary
                gme_t universe_id( GME::REGION, NO_ID ) ;
                add_entity_boundary( region_id, surface_id, !inside ) ;
            }
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2 * model().nb_surfaces(), NO_ID ) ;

            // Start with the first Surface on its + side
            std::stack< std::pair< index_t, bool > > S ;
            S.push( std::pair< index_t, bool >( 0, true ) ) ;

            while( !S.empty() ) {
                std::pair< index_t, bool > cur = S.top() ;
                S.pop() ;
                // This side is already assigned
                if( surf_2_region[
                    cur.second == true ? 2 * cur.first : 2 * cur.first + 1]
                    != NO_ID ) {
                    continue ;
                }
                // Create a new region
                gme_t cur_region_id = create_entity( GME::REGION ) ;
                // Get all oriented surfaces defining this region
                std::stack< std::pair< index_t, bool > > SR ;
                SR.push( cur ) ;
                while( !SR.empty() ) {
                    std::pair< index_t, bool > s = SR.top() ;
                    SR.pop() ;
                    index_t s_id = s.second == true ? 2 * s.first : 2 * s.first + 1 ;
                    // This oriented surface has already been visited
                    if( surf_2_region[s_id] != NO_ID ) {
                        continue ;
                    }
                    // Add the surface to the current region
                    add_entity_boundary( cur_region_id,
                        gme_t( GME::SURFACE, s.first ), s.second ) ;
                    surf_2_region[s_id] = cur_region_id.index ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp =
                        !s.second == true ? 2 * s.first : 2 * s.first + 1 ;
                    if( surf_2_region[s_id_opp] == NO_ID ) {
                        S.push( std::pair< index_t, bool >( s.first, !s.second ) ) ;
                    }
                    // For each contact, push the next oriented surface that is in the same region
                    const GeoModelEntity& surface = model().surface( s.first ) ;
                    for( index_t i = 0; i < surface.nb_boundaries(); ++i ) {
                        const std::pair< index_t, bool >& n =
                            regions_info_[surface.boundary_gme( i ).index]->next(
                                s ) ;
                        index_t n_id =
                            n.second == true ? 2 * n.first : 2 * n.first + 1 ;

                        if( surf_2_region[n_id] == NO_ID ) {
                            SR.push( n ) ;
                        }
                    }
                }
            }

            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            if( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID )
                != 0 ) {
                GEO::Logger::err( "GeoModel" )
                    << "Small bubble regions were skipped at model building "
                    << std::endl ;
                // Or, most probably, we have a problem before
                ringmesh_assert( false ) ;
            }

            // We need to remove from the regions_ the one corresponding
            // to the universe_, the one with the biggest volume
            double max_volume = -1. ;
            index_t universe_id = NO_ID ;
            for( index_t i = 0; i < model().nb_regions(); ++i ) {
                double cur_volume = model_entity_size( model().region( i ) ) ;
                if( cur_volume > max_volume ) {
                    max_volume = cur_volume ;
                    universe_id = i ;
                }
            }
            const Region& cur_region = model().region( universe_id ) ;
            for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ) {
                // Fill the Universe region boundaries
                // They are supposed to be empty
                add_entity_boundary( gme_t( GME::REGION, NO_ID ),
                    cur_region.boundary( i ).gme_id(), cur_region.side( i ) ) ;
            }
            std::set< gme_t > to_erase ;
            to_erase.insert( cur_region.gme_id() ) ;
            remove_entities( to_erase ) ;
        }
        return true ;
    }

    void GeoModelBuilder::build_model_from_surfaces()
    {
        if( model().nb_surfaces() == 0 ) {
            throw RINGMeshException( "GeoModel", "No surface to build the model " ) ;
        }

        // Initialize model() global vertices
        model().mesh.vertices.test_and_initialize() ;

        build_lines_and_corners_from_surfaces() ;

        if( options_.compute_regions_brep ) {
            build_brep_regions_from_surfaces() ;
        }

        // Finish up the model
        end_model() ;
    }

    /*!
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces
     */
    void GeoModelBuilder::build_contacts()
    {
        std::vector< std::set< gme_t > > interfaces ;
        for( index_t i = 0; i < model().nb_lines(); ++i ) {
            const Line& L = model().line( i ) ;
            std::set< gme_t > cur_interfaces ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                cur_interfaces.insert(
                    model().entity( L.in_boundary_gme( j ) ).parent().gme_id() ) ;
            }
            gme_t contact_id ;
            for( index_t j = 0; j < interfaces.size(); ++j ) {
                if( cur_interfaces.size() == interfaces[j].size()
                    && std::equal( cur_interfaces.begin(), cur_interfaces.end(),
                        interfaces[j].begin() ) ) {
                    contact_id = gme_t( GME::CONTACT, j ) ;
                    break ;
                }
            }
            if( !contact_id.is_defined() ) {
                contact_id = create_entity( GME::CONTACT ) ;
                ringmesh_assert( contact_id.index == interfaces.size() ) ;
                interfaces.push_back( cur_interfaces ) ;
                // Create a name for this contact
                std::string name = "contact_" ;
                for( std::set< gme_t >::const_iterator it( cur_interfaces.begin() );
                    it != cur_interfaces.end(); ++it ) {
                    name += model().entity( *it ).name() ;
                    name += "_" ;
                }
                set_entity_name( contact_id, name ) ;
            }
            add_entity_child( contact_id, gme_t( GME::LINE, i ) ) ;
        }
    }

    void GeoModelBuilder::update_facet_corner(
        Surface& S,
        const std::vector< index_t >& facets,
        index_t old,
        index_t neu )
    {
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t cur_f = facets[i] ;
            for( index_t cur_v = 0; cur_v < S.nb_mesh_element_vertices( cur_f );
                cur_v++ ) {
                if( S.mesh_element_vertex_index( cur_f, cur_v ) == old ) {
                    MeshBuilder builder( S.mesh_ ) ;
                    builder.set_facet_vertex( cur_f, cur_v, neu ) ;
                }
            }
        }
    }


    void GeoModelBuilder::delete_entity_mesh( GME::gme_t E_id )
    {
        Mesh& M = mesh_entity( E_id ).mesh_ ;
        MeshBuilder builder( M ) ;
        builder.clear( true, false ) ;
    }

    void GeoModelBuilder::delete_entity_vertices(
        GME::gme_t E_id,
        GEO::vector< index_t >& to_delete )
    {
        Mesh& M = mesh_entity( E_id ).mesh_ ;
        MeshBuilder builder( M ) ;
        builder.delete_vertices( to_delete, false ) ;
    }

    void GeoModelBuilder::delete_corner_vertex( index_t corner_id )
    {
        GME::gme_t corner( GME::CORNER, corner_id ) ;
        GEO::vector< index_t > to_delete ;
        to_delete.push_back( 1 ) ;
        delete_entity_vertices( corner, to_delete ) ;
    }
    void GeoModelBuilder::delete_line_edges(
        index_t line_id,
        GEO::vector< index_t >& to_delete )
    {
        Mesh& M = mesh_entity( GME::CORNER, line_id ).mesh_ ;
        MeshBuilder builder( M ) ;
        builder.delete_edges( to_delete, false ) ;
    }
    void GeoModelBuilder::delete_surface_facets(
        index_t surface_id,
        GEO::vector< index_t >& to_delete )
    {
        Mesh& M = mesh_entity( GME::SURFACE, surface_id ).mesh_ ;
        MeshBuilder builder( M ) ;
        builder.delete_facets( to_delete, false ) ;
        builder.connect_facets() ;
    }
    void GeoModelBuilder::delete_region_cells(
        index_t region_id,
        GEO::vector< index_t >& to_delete )
    {
        Mesh& M = mesh_entity( GME::CORNER, region_id ).mesh_ ;
        MeshBuilder builder( M ) ;
        builder.delete_cells( to_delete, false ) ;
        builder.connect_cells() ;
    }


    /*************************************************************************/
    GeoModelBuilderFile::GeoModelBuilderFile(
        GeoModel& model,
        const std::string& filename )
        : GeoModelBuilder( model ), filename_( filename )
    {

    }

    GME::TYPE GeoModelBuilderFile::match_nb_entities( const char* s )
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3 ) != 0 ) {
            return GME::NO_TYPE ;
        } else {
            for( index_t i = GME::CORNER; i < GME::NO_TYPE; i++ ) {
                GME::TYPE type = (GME::TYPE) i ;
                if( strstr( s, GME::type_name( type ).data() ) != NULL ) {
                    return type ;
                }
            }
            return GME::NO_TYPE ;
        }
    }

    GME::TYPE GeoModelBuilderFile::match_type( const char* s )
    {
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; i++ ) {
            GME::TYPE type = (GME::TYPE) i ;
            if( strcmp( s, GME::type_name( type ).data() ) == 0 ) {
                return type ;
            }
        }
        return GME::NO_TYPE ;
    }

} // namespace
