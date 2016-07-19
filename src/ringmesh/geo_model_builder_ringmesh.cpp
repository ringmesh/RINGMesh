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

#include <ringmesh/geo_model_builder_ringmesh.h>

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
    /*************************************************************************/

    void GeoModelBuilderGM::load_topology( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                // Name of the model
                if( file_line.field_matches( 0, "NAME" ) ) {
                    if( file_line.nb_fields() > 1 ) {
                        set_model_name( file_line.field( 1 ) ) ;
                    }
                }
                // Number of entities of a given type
                else if( match_nb_entities( file_line.field( 0 ) )
                    != GME::NO_TYPE ) {
                    // Allocate the space
                    if( file_line.nb_fields() > 1 ) {
                        create_entities( match_nb_entities( file_line.field( 0 ) ),
                            file_line.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level entities
                else if( match_high_level_type( file_line.field( 0 ) ) ) {
                    // Read this entity
                    // First line : type - id - name - geol_feature
                    if( file_line.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    GME::TYPE t = match_type( file_line.field( 0 ) ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( t, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        add_entity_child( entity,
                            gme_t( GME::child_type( t ),
                                file_line.field_as_uint( c ) ) ) ;
                    }
                }
                // Regions
                else if( match_type( file_line.field( 0 ) ) == GME::REGION ) {
                    // First line : type - id - name
                    if( file_line.nb_fields() < 3 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
                    }
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( GME::REGION, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    // Second line : signed indices of boundaries
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_entity_boundary( entity, gme_t( GME::SURFACE, s ),
                            side ) ;
                    }
                }

                // Universe
                else if( file_line.field_matches( 0, "UNIVERSE" ) ) {
                    // Second line: signed indices of boundaries
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_entity_boundary( gme_t( GME::REGION, NO_ID ),
                            gme_t( GME::SURFACE, s ), side ) ;
                    }
                }

            }
        }
    }
    void GeoModelBuilderGM::load_file()
    {
        unzFile uz = unzOpen( filename_.c_str() ) ;
        unz_global_info global_info ;
        if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not read file global info" ) ;
        }

        std::string topology = "topology.txt" ;
        unzip_one_file( uz, topology.c_str() ) ;

        GEO::LineInput line_topo( topology ) ;

        load_topology( line_topo ) ;
        GEO::FileSystem::delete_file( topology ) ;

        for( index_t t = GME::CORNER; t <= GME::REGION; t++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( t ) ;
            load_entities( type, uz ) ;
        }

        std::string connectivity = "connectivity.txt" ;
        unzip_one_file( uz, connectivity.c_str() ) ;

        GEO::LineInput line_connectivity( connectivity ) ;
        load_connectivities( line_connectivity ) ;
        GEO::FileSystem::delete_file( connectivity ) ;

        unzClose( uz ) ;
    }

    void GeoModelBuilderGM::load_connectivities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                if( file_line.field_matches( 0, "GME" ) ) {
                    GME::TYPE t = match_type( file_line.field( 1 ) ) ;
                    index_t id = file_line.field_as_uint( 2 ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    const GeoModelMeshEntity& cur_gme = model().mesh_entity( t,
                        id ) ;
                    gme_t cur_gme_type( t, id ) ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_entity_in_boundary( cur_gme_type,
                            gme_t( cur_gme.in_boundary_type( t ),
                                file_line.field_as_uint( in_b ) ) ) ;
                    }
                }

            }

        }

    }

    void GeoModelBuilderGM::load_entities( GME::TYPE gme_t, unzFile& uz )
    {
        for( index_t el = 0; el < model().nb_entities( gme_t ); el++ ) {
            std::string file_to_extract_and_load ;
            build_string_for_geo_model_entity_export( gme_t, el,
                file_to_extract_and_load ) ;
            std::string str_try = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                str_try = file_to_extract_and_load + ".meshb" ;
                if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                    if( gme_t != GME::REGION ) {
                        std::string message = "Invalid format of .gm file" ;
                        message += "\n.geogram file (defining mesh) is missing." ;
                        throw RINGMeshException( "I/O", message ) ;
                    }
                    return ; // a region is not necessary meshed.
                } else {
                    std::string message =
                        "Warning! you are using an old file (*.gm). \n" ;
                    message += "Please use ringmeshconvert to update this file. \n" ;
                    message +=
                        "ringmeshconvert in:geomodel=old_geomodel.gm out:geomodel=new_geomodel.gm" ;
                    GEO::Logger::warn( "I/O" ) << message << std::endl ;
                }
            }
            unzip_one_file( uz, str_try.c_str() ) ;
            Mesh cur_mesh( model(), 3, false ) ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            GEO::Logger::instance()->set_minimal( true ) ;
            MeshBuilder builder(cur_mesh);
            builder.load_mesh( str_try, flags ) ;
            assign_mesh_to_entity( cur_mesh,
                model().entity( gme_t, el ).gme_id() ) ;
            GEO::Logger::instance()->set_minimal( false ) ;

            unzip_one_file( uz, str_try.c_str() ) ;

//            set_connectivities
            GEO::FileSystem::delete_file( str_try ) ;
        }

    }

    void GeoModelBuilderGM::unzip_one_file(
        unzFile& uz,
        const char filename[MAX_FILENAME] )
    {
        unzLocateFile( uz, filename, 0 ) ;
        char read_buffer[ READ_SIZE] ;

        if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open file" ) ;
        }
        FILE *out = fopen( filename, "wb" ) ;
        if( out == NULL ) {
            unzCloseCurrentFile( uz ) ;
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open destination file" ) ;
        }
        int error = UNZ_OK ;
        do {
            error = unzReadCurrentFile( uz, read_buffer, READ_SIZE ) ;
            if( error < 0 ) {
                unzCloseCurrentFile( uz ) ;
                unzClose( uz ) ;
                fclose( out ) ;
                throw RINGMeshException( "ZLIB",
                    "Invalid error: " + GEO::String::to_string( error ) ) ;
            }
            if( error > 0 ) {
                fwrite( read_buffer, error, 1, out ) ;
            }
        } while( error > 0 ) ;
        fclose( out ) ;
        unzCloseCurrentFile( uz ) ;

    }

} // namespace
