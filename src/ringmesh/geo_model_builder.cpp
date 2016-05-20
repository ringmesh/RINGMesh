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
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/colocate.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/io.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geogram_extension.h>
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

    typedef GeoModelElement::gme_t gme_t ;
    typedef GeoModelMeshElement GMME ;

    void get_element_vertices_and_update_corners(
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
            if( geomodel.corner( i ).vertex() == point ) {
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

        index_t nb_neighbors = std::min( index_t( 5 ), surface.nb_cells() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, surface.nb_cells() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = ann.get_neighbors( v_bary, cur_neighbor, neighbors,
                dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                f = neighbors[i] ;
                for( index_t j = 0; j < surface.nb_vertices_in_facet( f ); j++ ) {
                    if( surface.model_vertex_id( f, j ) == model_v0 ) {
                        index_t j_next = surface.next_in_facet( f, j ) ;
                        if( surface.model_vertex_id( f, j_next ) == model_v1 ) {
                            e = j ;
                            return true ;
                        }
                    }
                }
            }
        } while( surface.nb_cells() != cur_neighbor ) ;

        f = Surface::NO_ID ;
        e = Surface::NO_ID ;
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

    void update_facet_corner(
        Surface& S,
        const std::vector< index_t >& facets,
        index_t old,
        index_t neu )
    {
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t cur_f = facets[i] ;
            for( index_t cur_v = 0; cur_v < S.nb_vertices_in_facet( cur_f );
                cur_v++ ) {
                if( S.surf_vertex_id( cur_f, cur_v ) == old ) {
                    S.mesh().facets.set_vertex( cur_f, cur_v, neu ) ;
                }
            }
        }
    }

    void get_sorted_incident_surfaces(
        const GeoModelElement& E,
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
                for( index_t j = 0; j < S.nb_cells(); ++j ) {
                    for( index_t v = 0; v < S.nb_vertices_in_facet( j ); ++v ) {
                        if( S.is_on_border( j, v ) ) {
                            index_t vertex = S.model_vertex_id( j, v ) ;
                            index_t next_vertex = S.model_vertex_id( j,
                                S.next_in_facet( j, v ) ) ;
                            index_t previous_vertex = S.model_vertex_id( j,
                                S.prev_in_facet( j, v ) ) ;
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

    void GeoModelBuilder::copy_meshes( const GeoModel& from, GME::TYPE element_type )
    {
        for( index_t i = 0; i < model_.elements( element_type ).size(); ++i ) {
            const GeoModelMeshElement& from_E = from.mesh_element( element_type,
                i ) ;
            assign_mesh_to_element( from_E.mesh(), gme_t( element_type, i ) ) ;
        }
    }

    void GeoModelBuilder::assign_mesh_to_element(
        const GEO::Mesh& mesh,
        GME::gme_t to )
    {
        GeoModelMeshElement& E = mesh_element( to ) ;
        E.unbind_attributes() ;
        E.mesh().copy( mesh ) ;
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
        gme_t result = find_corner( model_, point ) ;
        if( !result.is_defined() ) {
            result = create_element( GME::CORNER ) ;
            set_corner( result.index, point ) ;
        }
        return result ;
    }

    gme_t GeoModelBuilder::find_or_create_corner( index_t model_point_id )
    {
        gme_t result = find_corner( model_, model_point_id ) ;
        if( !result.is_defined() ) {
            result = create_element( GME::CORNER ) ;
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
            result = create_element( GME::LINE ) ;
            set_line( result.index, vertices ) ;

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_element_boundary( result,
                find_or_create_corner( vertices.front() ) ) ;
            add_element_boundary( result,
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
        return create_element( GME::LINE ) ;
    }

    /*!
     * @brief Sets the geometrical position of a vertex
     * @param[in] t Element index
     * @param[in] v Index of the vertex to modify
     * @param[in] point New coordinates
     * @param[in] update If true, all the vertices sharing the same geometrical position
     *               in the GeoModel have their position updated, if false they
     *               are not.
     * @warning Be careful with this update parameter, it is a very nice source of nasty bugs
     */
    void GeoModelBuilder::set_element_vertex(
        const GME::gme_t& t,
        index_t v,
        const vec3& point,
        bool update )
    {
        GeoModelMeshElement& E = mesh_element( t ) ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        if( update ) {
            model_.mesh.vertices.update_point( E.model_vertex_id( v ), point ) ;
        } else {
            E.mesh_.vertices.point( v ) = point ;
        }
    }

    /*!
     * @brief Sets the geometrical position of a vertex from a model vertex
     * @details Sets also both mapping from (GeoModelMeshVertices::unique2bme)
     *          and to (model_vertex_id_) the model vertex.
     * @param[in] element_id Element index
     * @param[in] v Index of the vertex to modify
     * @param[in] model_vertex Index in GeoModelMeshVertices of the vertex giving
     *                     the new position
     */
    void GeoModelBuilder::set_element_vertex(
        const gme_t& element_id,
        index_t v,
        index_t model_vertex )
    {
        set_element_vertex( element_id, v,
            model_.mesh.vertices.vertex( model_vertex ), false ) ;

        GeoModelMeshElement& E = mesh_element( element_id ) ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        E.model_vertex_id_[v] = model_vertex ;
        model_.mesh.vertices.add_to_bme( model_vertex, GMEVertex( element_id, v ) ) ;
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     * @param[in] id Element index
     * @param[in] points Geometric positions of the vertices to add
     * @param[in] clear If true the mesh is cleared, keeping its attributes
     */
    void GeoModelBuilder::set_element_vertices(
        const gme_t& id,
        const std::vector< vec3 >& points,
        bool clear )
    {
        GeoModelMeshElement& E = mesh_element( id ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            E.mesh_.clear( true, true ) ;
        }
        if( !points.empty() ) {
            index_t start = E.mesh_.vertices.create_vertices( points.size() ) ;
            GEO::Memory::copy( E.mesh_.vertices.point_ptr( start ),
                points.data()->data(), 3 * sizeof(double) * points.size() ) ;
        }
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     *
     * @param[in] id Element index
     * @param[in] model_vertices Geometric positions of the vertices to add
     * @param[in] clear If true the mesh if cleared, keeping its attributes
     */
    void GeoModelBuilder::set_element_vertices(
        const gme_t& element_id,
        const std::vector< index_t >& model_vertices,
        bool clear )
    {
        GeoModelMeshElement& E = mesh_element( element_id ) ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            E.mesh_.clear( true, true ) ;
        }
        index_t start = E.mesh_.vertices.create_vertices( model_vertices.size() ) ;
        for( index_t v = 0; v < model_vertices.size(); v++ ) {
            set_element_vertex( element_id, start + v, model_vertices[v] ) ;
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
        set_element_vertex( gme_t( GME::CORNER, corner_id ), 0, point, false ) ;
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
        set_element_vertices( gme_t( GME::LINE, line_id ), vertices, false ) ;

        GeoModelMeshElement& E = mesh_element( GME::LINE, line_id ) ;
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            E.mesh_.edges.create_edge( e - 1, e ) ;
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
        set_element_vertices( gme_t( GME::SURFACE, surface_id ), points, false ) ;
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
        set_element_vertices( gme_t( GME::REGION, region_id ), points, false ) ;
        assign_region_tet_mesh( region_id, tetras ) ;
    }

    /*!
     * @brief Add a point to the GeoModel and not to one of its elements
     * @details To use when adding the points to the model before building its elements
     */
    index_t GeoModelBuilder::add_unique_vertex( const vec3& p )
    {
        return model_.mesh.vertices.add_vertex( p ) ;
    }

    /*!
     * @brief Sets the vertex for a Corner. Store the info in the geomodel vertices
     *
     * @param[in] corner_id Index of the corner
     * @param[in] unique_vertex Index of the vertex in the model
     */
    void GeoModelBuilder::set_corner( index_t corner_id, index_t model_vertex_id )
    {
        set_element_vertex( gme_t( GME::CORNER, corner_id ), 0, model_vertex_id ) ;
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
        GeoModelMeshElement& E = mesh_element( GME::LINE, line_id ) ;
        ringmesh_assert( E.nb_vertices() == 0 ) ; // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_element_vertices( gme_t( GME::LINE, line_id ), unique_vertices,
            clear_vertices ) ;

        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            E.mesh_.edges.create_edge( e - 1, e ) ;
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
        set_element_vertices( gme_t( GME::SURFACE, surface_id ), model_vertex_ids,
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
        get_element_vertices_and_update_corners( new_facets, vertices ) ;
        set_surface_geometry( surface_id, vertices, new_facets, facet_ptr ) ;
    }

    void GeoModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& triangle_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_triangle_corners( triangle_corners ) ;
        get_element_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_element_vertices( gme_t( GME::SURFACE, surface_id ), vertices, false ) ;
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
        get_element_vertices_and_update_corners( new_triangle_corners, vertices ) ;

        set_element_vertices( gme_t( GME::SURFACE, surface_id ), vertices, false ) ;

        assign_surface_triangle_mesh( surface_id, new_triangle_corners,
            adjacent_triangles ) ;
    }

    void GeoModelBuilder::set_region_geometry(
        index_t region_id,
        const std::vector< index_t >& tet_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_tet_corners( tet_corners ) ;
        get_element_vertices_and_update_corners( new_tet_corners, vertices ) ;

        set_element_vertices( gme_t( GME::REGION, region_id ), vertices, false ) ;
        assign_region_tet_mesh( region_id, new_tet_corners ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices )
    {
        GEO::Mesh& M = mesh_element( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.vertices.nb() > 0 ) ;
        GEO::vector< index_t > copy ;
        copy_std_vector_to_geo_vector( triangle_vertices, copy ) ;
        M.facets.assign_triangle_mesh( copy, true ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices,
        const std::vector< index_t >& adjacent_triangles )
    {
        // COPY !! REMOVE !! // Access Mesh.... 
        GEO::Mesh& M = mesh_element( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.vertices.nb() > 0 ) ;
        GEO::vector< index_t > copy ;
        copy_std_vector_to_geo_vector( triangle_vertices, copy ) ;
        M.facets.assign_triangle_mesh( copy, true ) ;

        ringmesh_assert( adjacent_triangles.size() == M.facet_corners.nb() ) ;
        for( index_t i = 0; i < adjacent_triangles.size(); i++ ) {
            M.facet_corners.set_adjacent_facet( i, adjacent_triangles[i] ) ;
        }
    }

    void GeoModelBuilder::assign_surface_mesh_facets(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        GEO::Mesh& M = mesh_element( GME::SURFACE, surface_id ).mesh_ ;
        ringmesh_assert( M.vertices.nb() > 0 ) ;
        for( index_t f = 0; f + 1 < facet_ptr.size(); f++ ) {
            index_t start = facet_ptr[f] ;
            index_t end = facet_ptr[f + 1] ;
            GEO::vector< index_t > facet_vertices ;
            copy_std_vector_to_geo_vector( facets, start, end, facet_vertices ) ;

            M.facets.create_polygon( facet_vertices ) ;
        }
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_region_tet_mesh(
        index_t region_id,
        const std::vector< index_t >& tet_vertices ) const
    {
        GEO::Mesh& M = mesh_element( GME::REGION, region_id ).mesh_ ;
        ringmesh_assert( M.vertices.nb() > 0 ) ;
        GEO::vector< index_t > copy ;
        copy_std_vector_to_geo_vector( tet_vertices, copy ) ;
        M.cells.assign_tet_mesh( copy, true ) ;
        M.cells.connect() ;
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
        Surface& S = dynamic_cast< Surface& >( *model_.surfaces_[surface_id] ) ;
        ringmesh_assert( S.nb_cells() > 0 ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( S.nb_cells() - 1 ), Surface::NO_ADJACENT ) ;

        index_t nb_facets = S.nb_cells() ;
        index_t nb_vertices = S.nb_vertices() ;

        ///@todo Change representation of vertex_to_facets (without vector of vectors)

        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > empty_vector ;
        empty_vector.reserve( 10 ) ;
        std::vector< std::vector< index_t > > vertex_to_facets( nb_vertices,
            empty_vector ) ;

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                vertex_to_facets[S.surf_vertex_id( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ) {
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                index_t cur = S.surf_vertex_id( f, v ) ;
                index_t prev = S.surf_vertex_id( f, S.prev_in_facet( f, v ) ) ;

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
                    adjacent[S.facet_begin( f ) + S.prev_in_facet( f, v )] = f2 ;
                } else {
                    ringmesh_assert( end == 1 ) ;
                }
            }
        }

        ringmesh_assert( adjacent.size() == S.mesh_.facet_corners.nb() ) ;
        for( index_t i = 0; i < adjacent.size(); i++ ) {
            S.mesh_.facet_corners.set_adjacent_facet( i, adjacent[i] ) ;
        }
    }

    /*!
     * Finds duplicate vertex or creates it
     */
    index_t GeoModelBuilder::find_or_create_duplicate_vertex(
        GeoModelMeshElement& E,
        index_t model_vertex_id,
        index_t surface_vertex_id )
    {
        GeoModel& M = const_cast< GeoModel& >( E.model() ) ;

        const std::vector< GMEVertex >& vbme = M.mesh.vertices.gme_vertices(
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
            duplicate = E.mesh_.vertices.create_vertex(
                M.mesh.vertices.vertex( model_vertex_id ).data() ) ;

            // Set its model vertex index
            ringmesh_assert( duplicate < E.nb_vertices() ) ;
            E.model_vertex_id_[duplicate] = model_vertex_id ;

            // Add the mapping from in the model vertices. Should we do this one ?
            M.mesh.vertices.add_to_bme( model_vertex_id,
                GMEVertex( E.gme_id(), duplicate ) ) ;
        }
        return duplicate ;
    }

    /*
     * @brief Resets the adjacencies for all Surface facets adjacent to the Line
     * @pre All the edges of the Line are edges of at least one facet of the Surface
     */
    void GeoModelBuilder::disconnect_surface_facets_along_line_edges(
        Surface& S,
        const Line& L )
    {
        ColocaterANN ann( S.mesh(), ColocaterANN::FACETS ) ;
        for( index_t i = 0; i + 1 < L.nb_vertices(); ++i ) {
            index_t p0 = L.model_vertex_id( i ) ;
            index_t p1 = L.model_vertex_id( i + 1 ) ;

            index_t f = NO_ID ;
            index_t v = NO_ID ;
            bool found = find_facet_and_edge( ann, S, p0, p1, f, v ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && f != NO_ID && v != NO_ID ) ;

            index_t f2 = S.adjacent( f, v ) ;
            if( f2 != NO_ID ) {
                index_t v2 = NO_ID ;
                // Get the edge in facet f2 matching model indices p0-p1
                S.oriented_edge_from_model_vertex_ids( p0, p1, f2, v2 ) ;
                if( v2 == NO_ID ) {
                    S.oriented_edge_from_model_vertex_ids( p1, p0, f2, v2 ) ;
                    ringmesh_assert( v2 != NO_ID ) ;
                }
                S.mesh_.facets.set_adjacent( f, v, Surface::NO_ADJACENT ) ;
                S.mesh_.facets.set_adjacent( f2, v2, Surface::NO_ADJACENT ) ;
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

        ColocaterANN ann( S.mesh(), ColocaterANN::FACETS ) ;
        index_t p0 = L.model_vertex_id( 0 ) ;
        index_t p1 = L.model_vertex_id( 1 ) ;

        index_t v( NO_ID ) ;
        bool found = find_facet_and_edge( ann, S, p0, p1, facet_index, v ) ;
        ringmesh_unused( found ) ;
        ringmesh_assert( found && facet_index != NO_ID && v != NO_ID ) ;

        surface_vertex_0 = S.surf_vertex_id( facet_index, v ) ;
        surface_vertex_1 = S.surf_vertex_id( facet_index,
            S.next_in_facet( facet_index, v ) ) ;
    }

    /*!
     * @brief Duplicates the surface vertices along the fake boundary
     * (NO_ID adjacencies but shared vertices) and duplicate  the vertices
     * @note Bad written code - error prone
     * @todo Rewrite 
     */
    void GeoModelBuilder::duplicate_surface_vertices_along_line(
        Surface& S,
        const Line& L )
    {
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

            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;
            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Duplicate the vertex at id1
            // After having determined the next 1 we can probably get both at the same time
            // but I am lazy, and we must be careful not to break next_on_border function (Jeanne)
            std::vector< index_t > facets_around_id1 ;
            S.facets_around_vertex( id1, facets_around_id1, false, f ) ;

            index_t new_id1 = find_or_create_duplicate_vertex( S,
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

            index_t duplicated_surface_vertex = find_or_create_duplicate_vertex( S,
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
    void GeoModelBuilder::cut_surface_by_line( Surface& S, const Line& L )
    {
        /// @todo Replace the use of the model vertices by only a colocater
        /// of the surface vertices and the line vertices
        bool model_vertices_initialized = model().mesh.vertices.is_initialized() ;
        if( !model_vertices_initialized ) {
            model().mesh.vertices.test_and_initialize() ;
        }

        disconnect_surface_facets_along_line_edges( S, L ) ;
        duplicate_surface_vertices_along_line( S, L ) ;

        if( !model_vertices_initialized ) {
            const_cast< GeoModel& >( model() ).mesh.vertices.clear() ;
        }
    }

    void GeoModelBuilder::end_model()
    {
        if( model_.name() == "" ) {
            set_model_name( "model_default_name" ) ;
        }
        // Get out if the model has no surface
        if( model_.nb_surfaces() == 0 ) {
            print_geomodel( model_ ) ;
            throw RINGMeshException( "GeoModel", "The GeoModel has no surface" ) ;
        }

        model_.init_global_model_element_access() ;
        complete_element_connectivity() ;

        // Fill geological feature if they are missing
        for( index_t i = 0; i < model_.nb_elements( GME::ALL_TYPES ); ++i ) {
            GME& E = element( gme_t( GME::ALL_TYPES, i ) ) ;
            if( !E.has_geological_feature() ) {
                if( E.has_parent() && E.parent().has_geological_feature() ) {
                    set_element_geol_feature( E.gme_id(),
                        E.parent().geological_feature() ) ;
                } else if( E.nb_children() > 0
                    && E.child( 0 ).has_geological_feature() ) {
                    set_element_geol_feature( E.gme_id(),
                        E.child( 0 ).geological_feature() ) ;
                }
            }
        }

        // Deliberate clear of the model vertices used for model building
        model_.mesh.vertices.clear() ;
    }

    void GeoModelBuilder::recompute_geomodel_mesh()
    {
        model_.mesh.vertices.clear() ;
        model_.mesh.vertices.test_and_initialize() ;
    }

    bool GeoModelBuilder::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces line_computer( model_,
            options_.compute_regions_brep ) ;

        bool new_line_was_built = true ;
        while( new_line_was_built ) {
            new_line_was_built = line_computer.compute_next_line_geometry() ;

            // I know this is a copy - but should'nt be too big [JP]
            std::vector< index_t > vertices = line_computer.vertices() ;
            bool is_line_closed = vertices.front() == vertices.back() ;
            if( is_line_closed ) {
                // Vertices can begin and end at any vertex
                reorder_line_vertices_to_start_at_corner( model_, vertices ) ;
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
                    add_element_in_boundary( line_index, surface_id ) ;
                }
                add_element_boundary( line_index, first_corner ) ;
                add_element_boundary( line_index, second_corner ) ;

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
        ringmesh_assert( model_.nb_lines() == regions_info_.size() ) ;

        // Complete boundary information for surfaces
        // to compute volumetric regions
        fill_elements_boundaries( GME::SURFACE ) ;

        // Sort surfaces around the contacts
        for( index_t i = 0; i < regions_info_.size(); ++i ) {
            regions_info_[i]->sort() ;
        }

        if( model_.nb_surfaces() == 1 ) {
            if( model_.nb_lines() != 0 ) {
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
                gme_t region_id = create_element( GME::REGION ) ;
                add_element_boundary( region_id, surface_id, inside ) ;

                // Set universe boundary
                gme_t universe_id( GME::REGION, NO_ID ) ;
                add_element_boundary( region_id, surface_id, !inside ) ;
            }
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2 * model_.nb_surfaces(), NO_ID ) ;

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
                gme_t cur_region_id = create_element( GME::REGION ) ;
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
                    add_element_boundary( cur_region_id,
                        gme_t( GME::SURFACE, s.first ), s.second ) ;
                    surf_2_region[s_id] = cur_region_id.index ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp =
                        !s.second == true ? 2 * s.first : 2 * s.first + 1 ;
                    if( surf_2_region[s_id_opp] == NO_ID ) {
                        S.push( std::pair< index_t, bool >( s.first, !s.second ) ) ;
                    }
                    // For each contact, push the next oriented surface that is in the same region
                    const GeoModelElement& surface = model_.surface( s.first ) ;
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
            for( index_t i = 0; i < model_.nb_regions(); ++i ) {
                double cur_volume = model_element_size( model_.region( i ) ) ;
                if( cur_volume > max_volume ) {
                    max_volume = cur_volume ;
                    universe_id = i ;
                }
            }
            const Region& cur_region = model_.region( universe_id ) ;
            for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ) {
                // Fill the Universe region boundaries
                // They are supposed to be empty
                add_element_boundary( gme_t( GME::REGION, NO_ID ),
                    cur_region.boundary( i ).gme_id(), cur_region.side( i ) ) ;
            }
            std::set< gme_t > to_erase ;
            to_erase.insert( cur_region.gme_id() ) ;
            remove_elements( to_erase ) ;
        }
        return true ;
    }

    void GeoModelBuilder::build_model_from_surfaces()
    {
        if( model_.nb_surfaces() == 0 ) {
            throw RINGMeshException( "GeoModel", "No surface to build the model " ) ;
        }

        // Initialize model_ global vertices 
        model_.mesh.vertices.test_and_initialize() ;

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
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            const Line& L = model_.line( i ) ;
            std::set< gme_t > cur_interfaces ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                cur_interfaces.insert(
                    model_.element( L.in_boundary_gme( j ) ).parent().gme_id() ) ;
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
                contact_id = create_element( GME::CONTACT ) ;
                ringmesh_assert( contact_id.index == interfaces.size() ) ;
                interfaces.push_back( cur_interfaces ) ;
                // Create a name for this contact
                std::string name = "contact_" ;
                for( std::set< gme_t >::const_iterator it( cur_interfaces.begin() );
                    it != cur_interfaces.end(); ++it ) {
                    name += model_.element( *it ).name() ;
                    name += "_" ;
                }
                set_element_name( contact_id, name ) ;
            }
            add_element_child( contact_id, gme_t( GME::LINE, i ) ) ;
        }
    }

    /*************************************************************************/

    /*!
     * @brief Implementation detail: abstract base class to create a GeoModelElements 
     *        from SIMPLICIAL meshes
     * @details Manages the correspondence between a Mesh elements and
     *          GeoModelElements only known by indices.
     *
     * @warning Implemented only for TRIANGULATED Surface and TETRAHEDRALIZED Region.
     * @note Used by GeoModelBuilderMesh.
     */
    class GeoModelElementFromMesh {
    public:
        typedef std::pair< index_t, index_t > index_pair ;
        typedef std::map< index_t, index_t > index_map ;

        virtual ~GeoModelElementFromMesh()
        {
        }

        /*! Check that the attribute is defined.
         * If not, returns false otherwise bind it.
         */
        void initialize()
        {
            if( is_gme_attribute_defined() ) {
                bind_gme_attribute() ;
            }
        }

        bool is_valid()
        {
            bool attribute_is_bounded = gme_attribute_.is_bound() ;
            return attribute_is_bounded ;
        }

        index_t count_attribute_values_and_simplexes()
        {
            if( !is_valid() ) {
                return 0 ;
            }
            index_t nb = nb_mesh_simplexes() ;
            for( index_t i = 0; i != nb; ++i ) {
                index_t value = gme_attribute_[i] ;
                if( !is_attribute_value( value ) ) {
                    nb_simplexes_per_attribute_value_[value] = 0 ;
                }
                ++nb_simplexes_per_attribute_value_[value] ;
            }
            return nb_attribute_values() ;
        }

        /*! Number of different values for the attribute */
        index_t nb_attribute_values()
        {
            return index_t( nb_simplexes_per_attribute_value_.size() ) ;
        }

        /*! Sets a mapping from the attribute values on the Mesh and 
         * the indices of the GeoModelElements to fill
         */
        void set_gme_id_attribute_mapping(
            const std::vector< index_t >& gme_id_to_attribute_in )
        {
            gme_id_to_attribute_value_ = gme_id_to_attribute_in ;

            for( index_t i = 0; i != nb_gme(); ++i ) {
                index_t value = attribute_value_from_gme( i ) ;
                if( is_attribute_value( value ) ) {
                    attribute_value_to_gme_id_[value] = i ;
                } else {
                    GEO::Logger::err( "Debug" )
                        << "Invalid mapping between Mesh attribute and GeoModelElement"
                        << std::endl ;
                    gme_id_to_attribute_value_[i] = NO_ID ;
                }
            }
        }

        /*! Default mapping: first GeoModelElement index corresponds to first attribute
         *  value, second to second, etc.
         */
        void set_default_gme_id_attribute_mapping( index_t nb_geomodel_elements )
        {
            std::vector< index_t > default_mapping( nb_geomodel_elements, NO_ID ) ;
            index_t count = 0 ;
            index_map::const_iterator it(
                nb_simplexes_per_attribute_value_.begin() ) ;
            while( it != nb_simplexes_per_attribute_value_.end()
                && count < nb_geomodel_elements ) {
                default_mapping[count] = it->first ;
                ++count ;
                ++it ;
            }
            set_gme_id_attribute_mapping( default_mapping ) ;
        }

        /*!
         * Computes the simplex vertex indices for each GeoModelElement
         * as well as the mapping from the mesh simplices to the GME mesh simplices
         * for eventual attribute copying
         */
        void compute_gme_simplexes()
        {
            if( !is_valid() ) {
                return ;
            }
            allocate_mesh_simplex_to_gme() ;
            allocate_gme_vertices() ;

            index_t nb = nb_mesh_simplexes() ;
            std::vector< index_t > gme_simplex_counter_( nb_gme(), 0 ) ;
            for( index_t mesh_simplex = 0; mesh_simplex < nb; ++mesh_simplex ) {
                index_t attribute_value = gme_attribute_[mesh_simplex] ;
                if( attribute_value_has_gme_id( attribute_value ) ) {
                    index_t gme_id = attribute_value_to_gme_id_[attribute_value] ;
                    index_t gme_simplex_id = gme_simplex_counter_[gme_id] ;

                    assign_one_gme_simplex_vertices( mesh_simplex, gme_id,
                        gme_simplex_id ) ;
                    assign_mesh_simplex_to_gme_simplex( mesh_simplex, gme_id,
                        gme_simplex_id ) ;
                    ++gme_simplex_counter_[gme_id] ;
                } else {
                    assign_mesh_simplex_to_no_gme_simplex( mesh_simplex ) ;
                }
            }
        }

        /*!
         * Computes adjacencies between simplices for each GeoModelElement
         * compute_gme_simplexes SHOULD have been called before, nothing done if it wasn't
         */
        void compute_gme_adjacencies()
        {
            if( gme_simplex_vertices_.empty() ) {
                return ;
            }
            allocate_gme_corner_adjacent_gme_simplex() ;

            index_t nb = nb_mesh_simplexes() ;
            for( index_t mesh_simplex = 0; mesh_simplex < nb; ++mesh_simplex ) {
                GMESimplex gme = mesh_simplex_to_gme_simplex_[mesh_simplex] ;
                index_t gme_id = gme.gme_id ;
                index_t gme_simplex = gme.gme_simplex_id ;

                if( gme_id == NO_ID || gme_simplex == NO_ID ) {
                    continue ;
                }

                for( index_t v = 0; v < nb_vertices_per_simplex(); ++v ) {
                    index_t gme_vertex = nb_vertices_per_simplex() * gme_simplex
                        + v ;
                    index_t mesh_vertex = nb_vertices_per_simplex() * mesh_simplex
                        + v ;
                    // Number of vertices = number of edges / number of facets
                    index_t adjacent_simplex = adjacent_simplex_index(
                        mesh_vertex ) ;

                    GMESimplex adjacent_gme_simplex ;
                    if( adjacent_simplex != NO_ID ) {
                        adjacent_gme_simplex =
                            mesh_simplex_to_gme_simplex_[adjacent_simplex] ;
                    }
                    if( adjacent_gme_simplex.gme_id == gme_id ) {
                        gme_corner_adjacent_gme_simplex_[gme_id][gme_vertex] =
                            adjacent_gme_simplex.gme_simplex_id ;
                    }
                }
            }
        }

        /*! Simplex vertex indices in the Mesh for one GeoModelElement
         */
        const std::vector< index_t >& gme_simplices( index_t gme_id ) const
        {
            ringmesh_assert( gme_id < nb_gme() ) ;
            return gme_simplex_vertices_[gme_id] ;
        }

        const std::vector< index_t >& adjacent_gme_simplices( index_t gme_id ) const
        {
            ringmesh_assert( gme_id < nb_gme() ) ;
            return gme_corner_adjacent_gme_simplex_[gme_id] ;
        }

        template< typename T >
        void copy_simplex_attribute_from_mesh_to_geomodel(
            GEO::Attribute< T >& mesh_attribute,
            AttributeVector< T >& model_attributes ) const
        {
            for( index_t i = 0; i < nb_mesh_simplexes(); ++i ) {
                const GMESimplex& copy_to = mesh_simplex_to_gme_simplex_[i] ;
                model_attributes[copy_to.gme_id][copy_to.gme_simplex_id] =
                    mesh_attribute[i] ;
            }
        }

    protected:
        // A simplex in a GeoModelElement
        struct GMESimplex {
            GMESimplex()
                : gme_id( NO_ID ), gme_simplex_id( NO_ID )
            {
            }
            GMESimplex( index_t gme_in, index_t simplex_in )
                : gme_id( gme_in ), gme_simplex_id( simplex_in )
            {
            }
            index_t gme_id ;
            index_t gme_simplex_id ;
        } ;

    protected:
        GeoModelElementFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : mesh_( M ), gme_attribute_name_( attribute_name )
        {
        }

        /*! Number of GeoModelElements of the considered type */
        index_t nb_gme() const
        {
            return gme_id_to_attribute_value_.size() ;
        }

        /*! Is the value indeed taken by the attribute on the mesh */
        bool is_attribute_value( index_t value )
        {
            return nb_simplexes_per_attribute_value_.count( value ) == 1 ;
        }

        bool attribute_value_has_gme_id( index_t attribute_value ) const
        {
            return attribute_value_to_gme_id_.count( attribute_value ) == 1 ;
        }

        index_t attribute_value_from_gme( index_t gme_id ) const
        {
            return gme_id_to_attribute_value_[gme_id] ;
        }

        bool is_gme_attribute_defined()
        {
            GEO::AttributesManager& manager = mesh_simplex_attribute_manager() ;
            return manager.is_defined( gme_attribute_name_ ) ;
        }

        void bind_gme_attribute()
        {
            GEO::AttributesManager& manager = mesh_simplex_attribute_manager() ;
            gme_attribute_.bind( manager, gme_attribute_name_ ) ;
        }

        void allocate_gme_vertices()
        {
            gme_simplex_vertices_.resize( nb_gme() ) ;
            for( index_t i = 0; i < nb_gme(); ++i ) {
                index_t value = attribute_value_from_gme( i ) ;
                if( is_attribute_value( value ) ) {
                    index_t nb_simplexes = nb_simplexes_per_attribute_value_[value] ;
                    gme_simplex_vertices_[i].resize(
                        nb_vertices_per_simplex() * nb_simplexes ) ;
                }
            }
        }

        void allocate_gme_corner_adjacent_gme_simplex()
        {
            gme_corner_adjacent_gme_simplex_.resize( nb_gme() ) ;
            for( index_t i = 0; i < nb_gme(); ++i ) {
                gme_corner_adjacent_gme_simplex_[i].resize(
                    gme_simplex_vertices_[i].size(), NO_ID ) ;
            }
        }

        void allocate_mesh_simplex_to_gme()
        {
            mesh_simplex_to_gme_simplex_.resize( nb_mesh_simplexes() ) ;
        }

        virtual void assign_one_gme_simplex_vertices(
            index_t mesh_simplex_id,
            index_t gme_id,
            index_t gme_simplex_id )
        {
            index_t from = gme_simplex_id * nb_vertices_per_simplex() ;
            for( index_t v = 0; v != nb_vertices_per_simplex(); ++v ) {
                gme_simplex_vertices_[gme_id][from + v] = mesh_vertex_index(
                    mesh_simplex_id, v ) ;
            }
        }

        void assign_mesh_simplex_to_no_gme_simplex( index_t mesh_simplex_id )
        {
            assign_mesh_simplex_to_gme_simplex( mesh_simplex_id, NO_ID, NO_ID ) ;
        }

        void assign_mesh_simplex_to_gme_simplex(
            index_t mesh_simplex_id,
            index_t gme_id,
            index_t gme_simplex_id )
        {
            mesh_simplex_to_gme_simplex_[mesh_simplex_id] = GMESimplex( gme_id,
                gme_simplex_id ) ;
        }

        virtual index_t nb_mesh_simplexes() const = 0 ;
        virtual GEO::AttributesManager& mesh_simplex_attribute_manager() = 0 ;
        virtual index_t nb_vertices_per_simplex() const = 0 ;
        virtual index_t mesh_vertex_index(
            index_t simplex_id,
            index_t lv ) const = 0 ;
        virtual index_t adjacent_simplex_index(
            index_t facet_or_edge_id ) const = 0 ;

    protected:
        // THE Mesh
        const GEO::Mesh& mesh_ ;
        // Name of the attribute on the Mesh simplices identifying the GMEs
        std::string gme_attribute_name_ ;
        // Attribute giving the GeoModelElement index on Mesh simplices
        GEO::Attribute< index_t > gme_attribute_ ;
        // Number of simplices of the Mesh with a given attribute value
        std::map< index_t, index_t > nb_simplexes_per_attribute_value_ ;
        // Mapping from the attribute value on the Mesh to a GME index in a GeoModel
        std::map< index_t, index_t > attribute_value_to_gme_id_ ;
        // Mapping from a GME index in a GeoModel to the attribute value on the Mesh  
        std::vector< index_t > gme_id_to_attribute_value_ ;
        // Vertex indices (in the Mesh) of the corners of the simplices per GME
        std::vector< std::vector< index_t > > gme_simplex_vertices_ ;
        // Mapping from a mesh simplex index to a simplex index in a GME
        std::vector< GMESimplex > mesh_simplex_to_gme_simplex_ ;
        // For each GME, store the adjacent simplex (in the GME) for each simplex corner
        std::vector< std::vector< index_t > > gme_corner_adjacent_gme_simplex_ ;
    } ;

    class GeoModelSurfaceFromMesh: public GeoModelElementFromMesh {
    public:
        GeoModelSurfaceFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : GeoModelElementFromMesh( M, attribute_name )
        {
        }

        GEO::AttributesManager& mesh_simplex_attribute_manager()
        {
            return mesh_.facets.attributes() ;
        }

        index_t nb_mesh_simplexes() const
        {
            return mesh_.facets.nb() ;
        }

        index_t nb_vertices_per_simplex() const
        {
            return 3 ;
        }

        index_t mesh_vertex_index( index_t simplex_id, index_t vertex ) const
        {
            return mesh_.facets.vertex( simplex_id, vertex ) ;
        }
        index_t adjacent_simplex_index( index_t corner_id ) const
        {
            return mesh_.facet_corners.adjacent_facet( corner_id ) ;
        }
    } ;

    class GeoModelRegionFromMesh: public GeoModelElementFromMesh {
    public:
        GeoModelRegionFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : GeoModelElementFromMesh( M, attribute_name )
        {
        }

        GEO::AttributesManager& mesh_simplex_attribute_manager()
        {
            return mesh_.cells.attributes() ;
        }

        index_t nb_mesh_simplexes() const
        {
            return mesh_.cells.nb() ;
        }

        virtual index_t nb_vertices_per_simplex() const
        {
            return 4 ;
        }

        index_t mesh_vertex_index( index_t simplex_id, index_t vertex ) const
        {
            return mesh_.cells.vertex( simplex_id, vertex ) ;
        }

        index_t adjacent_simplex_index( index_t facet_id ) const
        {
            return mesh_.cell_facets.adjacent_cell( facet_id ) ;
        }
    } ;

    /*************************************************************************/

    GeoModelBuilderMesh::~GeoModelBuilderMesh()
    {
        delete surface_builder_ ;
        surface_builder_ = nil ;
        delete region_builder_ ;
        region_builder_ = nil ;
    }

    bool GeoModelBuilderMesh::is_mesh_valid_for_surface_building() const
    {
        bool valid = true ;
        if( !is_surface_mesh( mesh_ ) ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not a surface mesh "
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.facets.are_simplices() ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not triangulated"
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.facets.attributes().is_defined( surface_attribute_name_ ) ) {
            Logger::warn( "GMBuilder" ) << "The attribute "
                << surface_attribute_name_ << " is not defined on the Mesh facets"
                << std::endl ;
            valid = false ;
        }
        if( RINGMesh::has_mesh_colocate_vertices( mesh_, epsilon ) ) {
            Logger::warn( "GMBuilder" )
                << " The Mesh has colocated vertices. Repair it beforehand "
                << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    bool GeoModelBuilderMesh::is_mesh_valid_for_region_building() const
    {
        bool valid = true ;
        if( !is_volume_mesh( mesh_ ) ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not a volumetric mesh "
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.cells.are_simplices() ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not tetrahedralized"
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.cells.attributes().is_defined( region_attribute_name_ ) ) {
            Logger::warn( "GMBuilder" ) << "The attribute " << region_attribute_name_
                << " is not defined on the Mesh cells " << std::endl ;
            valid = false ;
        }
        if( RINGMesh::has_mesh_colocate_vertices( mesh_, epsilon ) ) {
            Logger::warn( "GMBuilder" )
                << " The Mesh has colocated vertices. Repair it beforehand. "
                << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    void create_and_fill_connected_component_attribute(
        GEO::Mesh& mesh,
        const std::string& connected_component_attribute )
    {
        GEO::Attribute< index_t > connected_component ;
        GEO::AttributesManager& manager = mesh.facets.attributes() ;
        connected_component.bind( manager, connected_component_attribute ) ;

        index_t nb_facets = mesh.facets.nb() ;
        std::vector< bool > visited( nb_facets, false ) ;

        ///@todo This algorithm is implemented over and over again in RINGMesh
        /// and Geogram. Couldn't we do better ?
        index_t nb_connected_components = 0 ;
        for( index_t f = 0; f < nb_facets; f++ ) {
            if( !visited[f] ) {
                nb_connected_components++ ;

                std::stack< index_t > facet_stack ;
                facet_stack.push( f ) ;

                while( !facet_stack.empty() ) {
                    index_t f_from_stack = facet_stack.top() ;
                    facet_stack.pop() ;
                    visited[f_from_stack] = true ;
                    connected_component[f_from_stack] = nb_connected_components ;

                    for( index_t v = 0; v < 3; ++v ) {
                        index_t neighbor_facet = mesh.facets.adjacent( f_from_stack,
                            v ) ;
                        if( neighbor_facet != NO_ID && !visited[neighbor_facet] ) {
                            visited[neighbor_facet] = true ;
                            facet_stack.push( neighbor_facet ) ;
                        }
                    }
                }
            }
        }
    }

    void GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components(
        GEO::Mesh& mesh,
        const std::string& created_facet_attribute )
    {
        // Remove duplicate facets, and triangulate the mesh
        // Side effects: fixes facet orientation and split non-manifold vertices
        // AND empty all attributes.
        GEO::mesh_repair( mesh,
            GEO::MeshRepairMode(
                GEO::MESH_REPAIR_DUP_F | GEO::MESH_REPAIR_TRIANGULATE ) ) ;

        create_and_fill_connected_component_attribute( mesh,
            created_facet_attribute ) ;

        // Remove colocated vertices
        repair_colocate_vertices( mesh, epsilon ) ;
    }

    /*! @details Adds separately each connected component of the mesh
     *          as a Surface of the model under construction.
     *          All the facets of the input mesh are visited and added to a
     *          Surface of the GeoModel.
     *          Connected components of the mesh are determined with a
     *          propagation (or "coloriage" algorithm) using the adjacent_facet
     *          information provided on the input GEO::Mesh.
     *
     * @todo Old code - old building - to delimit connected components
     * vertices are duplicated in the input mesh
     *
     */
    void GeoModelBuilderSurfaceMesh::build_polygonal_surfaces_from_connected_components()
    {
        std::vector< index_t > global_vertex_id_to_id_in_cc( mesh_.vertices.nb(),
            NO_ID ) ;

        std::vector< bool > visited( mesh_.facets.nb(), false ) ;
        for( index_t i = 0; i < mesh_.facets.nb(); i++ ) {
            if( !visited[i] ) {
                std::vector< index_t > cc_corners ;
                std::vector< index_t > cc_facets_ptr ;
                std::vector< vec3 > cc_vertices ;

                /// @todo Review : This should not be necessary as each vertex should
                /// be in one and only one connected component. To test. [JP]
                std::fill( global_vertex_id_to_id_in_cc.begin(),
                    global_vertex_id_to_id_in_cc.end(), NO_ID ) ;

                // First facet begin at corner 0
                cc_facets_ptr.push_back( 0 ) ;

                // Propagate from facet #i 
                std::stack< index_t > S ;
                S.push( i ) ;
                while( !S.empty() ) {
                    index_t f = S.top() ;
                    S.pop() ;
                    visited[f] = true ;

                    for( index_t c = mesh_.facets.corners_begin( f );
                        c < mesh_.facets.corners_end( f ); ++c ) {
                        index_t v = mesh_.facet_corners.vertex( c ) ;
                        if( global_vertex_id_to_id_in_cc[v] == NO_ID ) {
                            global_vertex_id_to_id_in_cc[v] = cc_vertices.size() ;
                            cc_vertices.push_back( mesh_.vertices.point( v ) ) ;
                        }
                        cc_corners.push_back( global_vertex_id_to_id_in_cc[v] ) ;

                        index_t n = mesh_.facet_corners.adjacent_facet( c ) ;
                        if( n != NO_ID && !visited[n] ) {
                            visited[n] = true ;
                            S.push( n ) ;
                        }
                    }
                    cc_facets_ptr.push_back( cc_corners.size() ) ;
                }

                gme_t surface_gme = create_element( GME::SURFACE ) ;
                set_surface_geometry( surface_gme.index, cc_vertices, cc_corners,
                    cc_facets_ptr ) ;
            }
        }
    }

    void GeoModelBuilderMesh::create_and_build_surfaces()
    {
        create_geomodel_elements( GME::SURFACE, nb_surface_attribute_values_ ) ;
        build_surfaces() ;
    }

    void GeoModelBuilderMesh::build_surfaces()
    {
        if( !is_mesh_valid_for_surface_building() ) {
            return ;
        }

        index_t nb_surfaces = model_.nb_surfaces() ;
        surface_builder_->set_default_gme_id_attribute_mapping( nb_surfaces ) ;
        surface_builder_->compute_gme_simplexes() ;
        surface_builder_->compute_gme_adjacencies() ;
        for( index_t i = 0; i != nb_surfaces; ++i ) {
            const std::vector< index_t >& triangle_vertices =
                surface_builder_->gme_simplices( i ) ;
            const std::vector< index_t >& adjacent_triangles =
                surface_builder_->adjacent_gme_simplices( i ) ;
            // Set the Surface facets
            // Set the adjacencies so that internal boundaries are not lost 
            // when connecting the surface facets.
            set_surface_geometry_with_adjacencies( i, triangle_vertices,
                adjacent_triangles ) ;
        }
    }

    void GeoModelBuilderMesh::create_and_build_regions()
    {
        create_geomodel_elements( GME::REGION, nb_region_attribute_values_ ) ;
        build_regions() ;
    }

    void GeoModelBuilderMesh::build_regions()
    {
        if( !is_mesh_valid_for_region_building() ) {
            return ;
        }

        index_t nb_regions = model_.nb_regions() ;
        region_builder_->set_default_gme_id_attribute_mapping( nb_regions ) ;

        region_builder_->compute_gme_simplexes() ;
        for( index_t i = 0; i != nb_regions; ++i ) {
            const std::vector< index_t >& tet_vertices =
                region_builder_->gme_simplices( i ) ;
            set_region_geometry( i, tet_vertices ) ;
        }
    }

    void GeoModelBuilderMesh::add_mesh_vertices_to_model()
    {
        index_t nb_vertices = mesh_.vertices.nb() ;
        for( index_t i = 0; i < nb_vertices; ++i ) {
            model_.mesh.vertices.add_vertex( mesh_.vertices.point( i ) ) ;
        }
    }

    void GeoModelBuilderMesh::initialize_surface_builder()
    {
        surface_builder_ = new GeoModelSurfaceFromMesh( mesh_,
            surface_attribute_name_ ) ;
        surface_builder_->initialize() ;
        nb_surface_attribute_values_ =
            surface_builder_->count_attribute_values_and_simplexes() ;
    }

    void GeoModelBuilderMesh::initialize_region_builder()
    {
        region_builder_ = new GeoModelRegionFromMesh( mesh_,
            region_attribute_name_ ) ;
        region_builder_->initialize() ;
        nb_region_attribute_values_ =
            region_builder_->count_attribute_values_and_simplexes() ;
    }

    void GeoModelBuilderMesh::copy_facet_attribute_from_mesh(
        const std::string& attribute_name )
    {
        if( !is_facet_attribute_defined< index_t >( mesh_, attribute_name ) ) {
            GEO::Logger::warn( "GMBuilder" ) << "No INDEX_T attribute named "
                << attribute_name << " on mesh facets to copy " << std::endl ;
            return ;
        }
        GEO::Attribute< index_t > attribute( mesh_.facets.attributes(),
            attribute_name ) ;
        AttributeVector< index_t > attributes ;
        create_attributes_on_geomodel_element_facets< index_t >( model_,
            GME::SURFACE, attribute_name, attributes ) ;
        surface_builder_->copy_simplex_attribute_from_mesh_to_geomodel< index_t >(
            attribute, attributes ) ;
    }

    void GeoModelBuilderMesh::copy_cell_attribute_from_mesh(
        const std::string& attribute_name )
    {
        if( !is_cell_attribute_defined< index_t >( mesh_, attribute_name ) ) {
            GEO::Logger::warn( "GMBuilder" ) << "No INDEX_T attribute named "
                << attribute_name << " on mesh cells to copy " << std::endl ;
            return ;
        }
        GEO::Attribute< index_t > attribute( mesh_.cells.attributes(),
            attribute_name ) ;
        AttributeVector< index_t > attributes ;
        create_attributes_on_geomodel_element_cells< index_t >( model_, GME::REGION,
            attribute_name, attributes ) ;
        region_builder_->copy_simplex_attribute_from_mesh_to_geomodel< index_t >(
            attribute, attributes ) ;
    }

    /*************************************************************************/
    GeoModelBuilderFile::GeoModelBuilderFile(
        GeoModel& model,
        const std::string& filename )
        : GeoModelBuilder( model ), filename_( filename )
    {

    }

    GME::TYPE GeoModelBuilderFile::match_nb_elements( const char* s )
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

    /*************************************************************************/

    /*!
     * @brief Loads and builds a GeoModel from a Gocad .ml file
     * @warning Pretty unstable. Crashes if the file is not exactly what is expected.
     * @details Correspondance between Gocad::Model3D elements
     * and GeoModel elements is :
     *  - Gocad TSurf  <-> GeoModel Interface
     *  - Gocad TFace  <-> GeoModel Surface
     *  - Gocad Region <-> GeoModel Region
     *  - Gocad Layer  <-> GeoModel Layer
     * @param[in] ml_file_name Input .ml file stream
     * @param[in] ignore_file_borders If true, BORDER and BSTONE entries in the files
     * are ignored and the Lines and Corners of the GeoModel are deduced from the
     * connectivity of its Surfaces. By default set to false.
     */
    void GeoModelBuilderGocad::load_file()
    {
        // Count the number of TSurf - Interface
        index_t nb_tsurf = 0 ;

        // Count the number of TFace - Surface
        index_t nb_tface = 0 ;

        // Counters identifying the currently read TSurf or TFace
        index_t tsurf_count = 0 ;
        index_t tface_count = 0 ;

        index_t nb_tface_in_prev_tsurf = 0 ;

        /// The file contains 2 parts and is read in 2 steps
        /// 1. Read global information on model elements
        /// 2. Read surface geometries and info to build corners and contacts
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and between TSurfs

        // Coordinates of the vertices of the currently built TSurf in the model
        std::vector< vec3 > tsurf_vertices ;

        // Where the vertices of a TFace start in the vertices of the TSurf (offset)
        std::vector< index_t > tface_vertex_start ;

        // Triangles of the currently built TFace
        std::vector< index_t > tface_facets ;

        // Starting and ending indices of each facet triangle in the tface_facets vector
        /// @todo This is useless. Facets are all triangles.
        /// Write functions to be able to use the Mesh.facets.assign_triangle_mesh function [JP]
        std::vector< index_t > tface_facets_ptr ;
        tface_facets_ptr.push_back( 0 ) ;

        // Intermediate information for contact parts building
        std::vector< Border > borders_to_build ;

        // Surfaces for which the KeyFacet orientation should be changed
        // because it does not match triangle orientations.
        std::vector< index_t > change_key_facet ;

        ///@todo Add assert everywhere when doing substraction operations on unsigned int [JP]

        while( !file_line_.eof() && file_line_.get_line() ) {
            file_line_.get_fields() ;
            if( file_line_.nb_fields() > 0 ) {
                if( read_model ) {
                    if( strncmp( file_line_.field( 0 ), "name:", 5 ) == 0 ) {
                        // Sometimes there is a space after name:
                        // Sometimes not
                        if( file_line_.nb_fields() > 1 ) {
                            set_model_name( file_line_.field( 1 ) ) ;
                        } else {
                            set_model_name( &file_line_.field( 0 )[5] ) ;
                        }
                    } else if( file_line_.field_matches( 0, "TSURF" ) ) {
                        /// 1.1 Create Interface from its name
                        index_t f = 1 ;
                        std::ostringstream oss ;
                        do {
                            oss << file_line_.field( f++ ) ;
                        } while( f < file_line_.nb_fields() ) ;
                        // Create an interface and set its name
                        set_element_name( create_element( GME::INTERFACE ),
                            oss.str() ) ;

                        nb_tsurf++ ;
                    } else if( file_line_.field_matches( 0, "TFACE" ) ) {
                        /// 1.2 Create Surface from the name of its parent Interface
                        /// and its geological feature
                        std::string geol = file_line_.field( 2 ) ;
                        index_t f = 3 ;
                        std::ostringstream oss ;
                        do {
                            oss << file_line_.field( f++ ) ;
                        } while( f < file_line_.nb_fields() ) ;
                        std::string interface_name = oss.str() ;

                        // And its key facet that give the orientation of the surface part
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        vec3 p0( file_line_.field_as_double( 0 ),
                            file_line_.field_as_double( 1 ),
                            file_line_.field_as_double( 2 ) ) ;
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        vec3 p1( file_line_.field_as_double( 0 ),
                            file_line_.field_as_double( 1 ),
                            file_line_.field_as_double( 2 ) ) ;
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        vec3 p2( file_line_.field_as_double( 0 ),
                            file_line_.field_as_double( 1 ),
                            file_line_.field_as_double( 2 ) ) ;

                        create_surface( interface_name, geol, p0, p1, p2 ) ;
                        nb_tface++ ;
                    } else if( file_line_.field_matches( 0, "REGION" ) ) {
                        /// 1.3 Read Region information and create them from their name,
                        /// and the surfaces on their boundary
                        std::string name = file_line_.field( 2 ) ;

                        std::vector< std::pair< index_t, bool > > region_boundaries ;
                        bool end_region = false ;
                        while( !end_region ) {
                            file_line_.get_line() ;
                            file_line_.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                signed_index_t signed_id = file_line_.field_as_int(
                                    i ) ;
                                if( signed_id == 0 ) {
                                    end_region = true ;
                                    break ;
                                }
                                bool side = signed_id > 0 ;
                                index_t id = static_cast< index_t >( std::abs(
                                    signed_id ) - 1 ) ;
                                region_boundaries.push_back(
                                    std::pair< index_t, bool >( id, side ) ) ;
                            }
                        }

                        // By default the region id is the universe id
                        gme_t region_id( GME::REGION, NO_ID ) ;
                        // Create the element if it is not the universe
                        if( name != "Universe" ) {
                            region_id = create_element( GME::REGION ) ;
                        }
                        // Set the region name and boundaries
                        set_element_name( region_id, name ) ;
                        for( index_t i = 0; i < region_boundaries.size(); ++i ) {
                            add_element_boundary( region_id,
                                gme_t( GME::SURFACE, region_boundaries[i].first ),
                                region_boundaries[i].second ) ;
                        }
                    } else if( file_line_.field_matches( 0, "LAYER" ) ) {
                        /// 1.4 Build the volumetric layers from their name and
                        /// the ids of the regions they contain
                        gme_t layer_id = create_element( GME::LAYER ) ;
                        set_element_name( layer_id, file_line_.field( 1 ) ) ;
                        bool end_layer = false ;
                        while( !end_layer ) {
                            file_line_.get_line() ;
                            file_line_.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                index_t region_id = file_line_.field_as_uint( i ) ;
                                if( region_id == 0 ) {
                                    end_layer = true ;
                                    break ;
                                } else {
                                    region_id -= nb_tface + 1 ; // Remove Universe region
                                    // Correction because ids begin at 1 in the file
                                    add_element_child( layer_id,
                                        gme_t( GME::REGION, region_id - 1 ) ) ;
                                }
                            }
                        }
                    } else if( file_line_.field_matches( 0, "END" ) ) {
                        // End of the high level information on the model
                        // Switch to reading the geometry of the model surfaces
                        read_model = false ;
                        continue ;
                    }
                } else {
                    if( file_line_.field_matches( 0, "GOCAD" ) ) {
                        // This is the beginning of a new TSurf = Interface
                        tsurf_count++ ;
                    }
                    if( file_line_.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( file_line_.field_matches( 1, "Elevation" ) ) {
                            z_sign = 1 ;
                        } else if( file_line_.field_matches( 1, "Depth" ) ) {
                            z_sign = -1 ;
                        } else {
                            ringmesh_assert_not_reached ;
                        }
                    } else if( file_line_.field_matches( 0, "END" ) ) {
                        // This the END of a TSurf
                        if( tsurf_count > 0 ) {
                            // End the last TFace - Surface of this TSurf
                            set_surface_geometry( tface_count - 1,
                                std::vector< vec3 >(
                                    tsurf_vertices.begin()
                                        + tface_vertex_start.back(),
                                    tsurf_vertices.end() ), tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;

                            // End this TSurf - Interface
                            nb_tface_in_prev_tsurf += tface_vertex_start.size() ;
                            tsurf_vertices.clear() ;
                            tface_vertex_start.clear() ;
                        }
                    } else if( file_line_.field_matches( 0, "TFACE" ) ) {
                        // Beginning of a new TFace - Surface
                        if( tface_vertex_start.size() > 0 ) {
                            // End the previous TFace - Surface  (copy from line 1180)
                            set_surface_geometry( tface_count - 1,
                                std::vector< vec3 >(
                                    tsurf_vertices.begin()
                                        + tface_vertex_start.back(),
                                    tsurf_vertices.end() ), tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;
                        }

                        // Register where begin the new TFace vertices
                        tface_vertex_start.push_back( tsurf_vertices.size() ) ;

                        tface_count++ ;
                    }

                    /// 2.1 Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    else if( file_line_.field_matches( 0, "VRTX" )
                        || file_line_.field_matches( 0, "PVRTX" ) ) {
                        vec3 p( file_line_.field_as_double( 2 ),
                            file_line_.field_as_double( 3 ),
                            z_sign * file_line_.field_as_double( 4 ) ) ;
                        tsurf_vertices.push_back( p ) ;
                    } else if( file_line_.field_matches( 0, "PATOM" )
                        || file_line_.field_matches( 0, "ATOM" ) ) {
                        tsurf_vertices.push_back(
                            tsurf_vertices[file_line_.field_as_uint( 2 ) - 1] ) ;
                    } else if( file_line_.field_matches( 0, "TRGL" ) ) {
                        // Read ids of the vertices of each triangle in the TSurf
                        // and switch to ids in the TFace
                        tface_facets.push_back(
                            (index_t) file_line_.field_as_uint( 1 )
                                - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back(
                            (index_t) file_line_.field_as_uint( 2 )
                                - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back(
                            (index_t) file_line_.field_as_uint( 3 )
                                - tface_vertex_start.back() - 1 ) ;
                        tface_facets_ptr.push_back( tface_facets.size() ) ;
                    }

                    // 2.2 Build the corners from their position and the surface parts
                    //    containing them
                    else if( file_line_.field_matches( 0, "BSTONE" )
                        && !options_.compute_corners ) {
                        index_t v_id = file_line_.field_as_uint( 1 ) - 1 ;
                        if( !find_corner( model_, tsurf_vertices[v_id] ).is_defined() ) {
                            // Create the corner
                            gme_t corner_gme = create_element( GME::CORNER ) ;
                            set_corner( corner_gme.index, tsurf_vertices[v_id] ) ;
                        }
                    }

                    /// 2.3 Read the Border information and store it
                    else if( file_line_.field_matches( 0, "BORDER" )
                        && !options_.compute_lines ) {
                        index_t p1 = file_line_.field_as_uint( 2 ) - 1 ;
                        index_t p2 = file_line_.field_as_uint( 3 ) - 1 ;

                        // Get the global corner id
                        gme_t corner_id = find_corner( model_, tsurf_vertices[p1] ) ;
                        ringmesh_assert( corner_id.is_defined() ) ;

                        // Get the surface
                        index_t part_id = NO_ID ;
                        for( index_t i = 0; i < tface_vertex_start.size(); ++i ) {
                            if( p1 < tface_vertex_start[i] ) {
                                ringmesh_assert( p2 < tface_vertex_start[i] ) ;

                                // Get vertices ids in the surface
                                p1 = p1 - tface_vertex_start[i - 1] ;
                                p2 = p2 - tface_vertex_start[i - 1] ;

                                // i-1 is the id of the TFace in this TSurf
                                part_id = i - 1 ;
                                break ;
                            }
                        }
                        if( part_id == NO_ID ) {
                            // It is in the last built Tface
                            p1 = p1
                                - tface_vertex_start[tface_vertex_start.size() - 1] ;
                            p2 = p2
                                - tface_vertex_start[tface_vertex_start.size() - 1] ;

                            part_id = tface_vertex_start.size() - 1 ;
                        }

                        // The number of tfaces in previous tsurf is also to add
                        part_id += nb_tface_in_prev_tsurf ;

                        borders_to_build.push_back(
                            Border( part_id, corner_id.index, p1, p2 ) ) ;
                    }
                }
            }
        }

        // I agree that we do not need to compute the GeoModelMeshVertices here
        // But perhaps the computation of Lines would be faster and safer [JP]

        /// 3. Build the Lines        
        if( !options_.compute_lines ) {
            // Use info of the .ml file to fill the Lines
            std::vector< vec3 > line_vertices ;
            for( index_t i = 0; i < borders_to_build.size(); ++i ) {
                const Border& b = borders_to_build[i] ;
                // 1- Build the boundary : construct the vector
                // of vertices on the border
                const Surface& S = model_.surface( b.part_id_ ) ;
                determine_line_vertices( S, b.p0_, b.p1_, line_vertices ) ;
                if( line_vertices.empty() ) {
                    GEO::Logger::out( "I/O" )
                        << "One Line vertices determination failed in SURFACE "
                        << S.index() << std::endl ;
                } else {
                    // 2 - Check if this border already exists
                    gme_t line_id = find_or_create_line( line_vertices ) ;
                    // Add the surface in which this line is
                    add_element_in_boundary( line_id, S.gme_id() ) ;
                }
            }
        } else {
            // Ignore BORDER and CORNER information of the file
            // Create them now from the topology of the Surfaces
            model_.mesh.vertices.test_and_initialize() ;
            build_lines_and_corners_from_surfaces() ;
        }

        /// 4. Build the Contacts
        build_contacts() ;

        // Modify in the Region the side of the Surface for which the key facet
        // orientation was not the same than their facet orientations
        for( index_t i = 0; i < change_key_facet.size(); i++ ) {
            const Surface& S = model_.surface( change_key_facet[i] ) ;
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {
                Region& R = dynamic_cast< Region& >( element(
                    S.in_boundary_gme( j ) ) ) ;
                for( index_t b = 0; b < R.nb_boundaries(); ++b ) {
                    if( R.boundary_gme( b ).index == change_key_facet[i] ) {
                        bool old_side = R.side( b ) ;
                        set_element_boundary( R.gme_id(), b, R.boundary_gme( b ),
                            !old_side ) ;
                    }
                }
            }
        }
    }

    /*!
     * @brief Find the facet which first 3 vertices are given
     * 
     * @param[in] surface_id Index of the surface
     * @param[in] p0 First point coordinates
     * @param[in] p1 Second point coordinates
     * @param[in] p2 Third point coordinates
     * @param[out] same_sign Is true if the found facet has the same orientation than triangle p0p1p2
     * @return Index of the found facet, NO_ID if none found
     */
    index_t GeoModelBuilderGocad::find_key_facet(
        index_t surface_id,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool& same_sign ) const
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;

        for( index_t t = 0; t < surface.nb_cells(); ++t ) {
            const vec3& pp0 = surface.vertex( t, 0 ) ;
            const vec3& pp1 = surface.vertex( t, 1 ) ;
            const vec3& pp2 = surface.vertex( t, 2 ) ;

            if( p0 == pp0 ) {
                if( p1 == pp1 && p2 == pp2 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp1 ) {
                    same_sign = false ;
                    return t ;
                }
            }
            if( p0 == pp1 ) {
                if( p1 == pp0 && p2 == pp2 ) {
                    same_sign = false ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp0 ) {
                    same_sign = true ;
                    return t ;
                }
            }
            if( p0 == pp2 ) {
                if( p1 == pp0 && p2 == pp1 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp1 && p2 == pp0 ) {
                    same_sign = false ;
                    return t ;
                }
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Verify that a surface key facet has an orientation consistent with the surface facets.
     * 
     * @param[in] surface_id Index of the surface
     * @return False if the key_facet orientation is not the same than the surface facets, else true.
     */
    bool GeoModelBuilderGocad::check_key_facet_orientation(
        index_t surface_id ) const
    {
        const KeyFacet& key_facet = key_facets_[surface_id] ;

        const vec3& p0 = key_facet.p0_ ;
        const vec3& p1 = key_facet.p1_ ;
        const vec3& p2 = key_facet.p2_ ;
        bool same_sign = false ;

        index_t t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ;
        if( t == NO_ID ) {
            vec3 p00( p0 ) ;
            p00.z *= -1 ;
            vec3 p10( p1 ) ;
            p10.z *= -1 ;
            vec3 p20( p2 ) ;
            p20.z *= -1 ;

            // It is because of the sign of Z that is not the same
            t = find_key_facet( surface_id, p00, p10, p20, same_sign ) ;
        }
        ringmesh_assert( t != NO_ID ) ;
        return same_sign ;
    }

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_vertices Coordinates of the vertices on the Line (emptied and filled again)
     * @return Index of the Corner at which the Line ends
     */
    gme_t GeoModelBuilderGocad::determine_line_vertices(
        const Surface& S,
        index_t id0,
        index_t id1,
        std::vector< vec3 >& border_vertex_model_vertices ) const
    {
        ringmesh_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_vertices.resize( 0 ) ;

        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
//        ringmesh_assert( f != Surface::NO_ID ) ;
        if( f == NO_ID ) {
            border_vertex_model_vertices.resize( 0 ) ;
            return gme_t() ;
        }

        vec3 p0 = S.vertex( id0 ) ;
        vec3 p1 = S.vertex( id1 ) ;

        border_vertex_model_vertices.push_back( p0 ) ;
        border_vertex_model_vertices.push_back( p1 ) ;

        gme_t p1_corner = find_corner( model(), p1 ) ;
        while( !p1_corner.is_defined() ) {
            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;

            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.vertex( next_id1 ) ;
            border_vertex_model_vertices.push_back( p1 ) ;
            p1_corner = find_corner( model(), p1 ) ;
        }
        return p1_corner ;
    }

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] p0 Coordinates of the 1 point of the TFace key facet 
     * @param[in] p1 Coordinates of the 2 point of the TFace key facet 
     * @param[in] p2 Coordinates of the 3 point of the TFace key facet 
     */
    void GeoModelBuilderGocad::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        gme_t parent = find_interface( model_, interface_name ) ;
        if( interface_name != "" ) {
            ringmesh_assert( parent.is_defined() ) ;
        }

        gme_t id = create_element( GME::SURFACE ) ;
        set_element_parent( id, parent ) ;
        set_element_geol_feature( parent, GME::determine_geological_type( type ) ) ;
        key_facets_.push_back( KeyFacet( p0, p1, p2 ) ) ;
    }

    int GeoModelBuilderGocad::read_gocad_coordinates_system( const std::string& in )
    {
        if( in == "Elevation" ) {
            return 1 ;
        } else if( in == "Depth" ) {
            return -1 ;
        } else {
            ringmesh_assert_not_reached ;
            return 0 ;
        }
    }

    /*************************************************************************/

    void GeoModelBuilderBM::load_file()
    {
        while( !file_line_.eof() && file_line_.get_line() ) {
            file_line_.get_fields() ;
            if( file_line_.nb_fields() > 0 ) {
                // Name of the model
                if( file_line_.field_matches( 0, "NAME" ) ) {
                    if( file_line_.nb_fields() > 1 ) {
                        set_model_name( file_line_.field( 1 ) ) ;
                    }
                }
                // Number of elements of a given type
                else if( match_nb_elements( file_line_.field( 0 ) )
                    != GME::NO_TYPE ) {
                    // Allocate the space
                    if( file_line_.nb_fields() > 1 ) {
                        create_elements( match_nb_elements( file_line_.field( 0 ) ),
                            file_line_.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level elements
                else if( match_high_level_type( file_line_.field( 0 ) ) ) {
                    // Read this element
                    // First line : type - id - name - geol_feature
                    if( file_line_.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line_.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    GME::TYPE t = match_type( file_line_.field( 0 ) ) ;
                    index_t id = file_line_.field_as_uint( 1 ) ;
                    gme_t element( t, id ) ;
                    set_element_name( element, file_line_.field( 2 ) ) ;
                    set_element_geol_feature( element,
                        GME::determine_geological_type( file_line_.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    for( index_t c = 0; c < file_line_.nb_fields(); c++ ) {
                        add_element_child( element,
                            gme_t( GME::child_type( t ),
                                file_line_.field_as_uint( c ) ) ) ;
                    }
                }
                // Regions
                else if( match_type( file_line_.field( 0 ) ) == GME::REGION ) {
                    // First line : type - id - name
                    if( file_line_.nb_fields() < 3 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line_.line_number() )
                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
                    }
                    index_t id = file_line_.field_as_uint( 1 ) ;
                    gme_t element( GME::REGION, id ) ;
                    set_element_name( element, file_line_.field( 2 ) ) ;
                    // Second line : signed indices of boundaries
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    for( index_t c = 0; c < file_line_.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line_.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line_.field( c )[1], s ) ;

                        add_element_boundary( element, gme_t( GME::SURFACE, s ),
                            side ) ;
                    }
                }

                // Universe
                else if( file_line_.field_matches( 0, "UNIVERSE" ) ) {
                    // Second line: signed indices of boundaries
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    for( index_t c = 0; c < file_line_.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line_.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line_.field( c )[1], s ) ;

                        add_element_boundary( gme_t( GME::REGION, NO_ID ),
                            gme_t( GME::SURFACE, s ), side ) ;
                    }
                }

                // Model vertices
//                else if( in_.field_matches( 0, "MODEL_VERTICES" ) ) {
//                    index_t nb_vertices = in_.field_as_uint( 1 ) ;
//
//                    // Attributes
//                    in_.get_line() ;
//                    in_.get_fields() ;
//                    ringmesh_assert( in_.field_matches( 0, "MODEL_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in_.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute< GeoModel::VERTEX > >
//                    vertex_attribs( nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            model_.vertex_attribute_manager(), in_.field(
//                                1 + 2 * i ), in_.field( 2 + 2 * i ), nb_vertices ) ;
//                    }
//                    for( index_t i = 0; i < nb_vertices; ++i ) {
//                        in_.get_line() ;
//                        in_.get_fields() ;
//                        add_vertex( vec3(
//                                read_double( in_,
//                                    0 ), read_double( in_, 1 ), read_double( in_, 2 ) ) ) ;
//                        serialize_read_attributes( in_, 3, i, vertex_attribs ) ;
//                    }
//                }

                // Corners
                else if( match_type( file_line_.field( 0 ) ) == GME::CORNER ) {
                    // First line: CORNER - id - vertex id
                    if( file_line_.nb_fields() < 5 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line_.line_number() )
                                + ", 5 fields are expected to describe a corner: "
                                + " CORNER, index, and X, Y, Z coordinates " ) ;
                    }
                    index_t id = file_line_.field_as_uint( 1 ) ;
                    vec3 point( file_line_.field_as_double( 2 ),
                        file_line_.field_as_double( 3 ),
                        file_line_.field_as_double( 4 ) ) ;
                    set_corner( id, point ) ;
                }

                // Lines
                else if( match_type( file_line_.field( 0 ) ) == GME::LINE ) {
                    index_t id = file_line_.field_as_uint( 1 ) ;
                    gme_t cur_element( GME::LINE, id ) ;

                    // Following information: vertices of the line
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    ringmesh_assert( file_line_.field_matches( 0, "LINE_VERTICES" ) ) ;
                    index_t nb_vertices = file_line_.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        vec3 point( file_line_.field_as_double( 0 ),
                            file_line_.field_as_double( 1 ),
                            file_line_.field_as_double( 2 ) ) ;
                        vertices[i] = point ;
                    }

                    // Set the line points
                    set_line( cur_element.index, vertices ) ;

                    // Attributes on line vertices
//                    in_.get_line() ;
//                    in_.get_fields() ;
//                    ringmesh_assert( in_.field_matches( 0, "LINE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in_.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute > vertex_attribs(
//                        nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            L.vertex_attribute_manager(), in_.field(
//                                1 + 2 * i ), in_.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices indices and attributes on vertices
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in_.get_line() ;
//                        in_.get_fields() ;
//                        serialize_read_attributes( in_, 1, i, vertex_attribs ) ;
//                    }

                    // Read attributes on line segments
//                    in_.get_line() ;
//                    in_.get_fields() ;
//                   ringmesh_assert( in_.field_matches( 0, "LINE_SEGMENT_ATTRIBUTES" ) ) ;
//                    index_t nb_segment_attribs = ( in_.nb_fields() - 1 ) / 2 ;
//                    if( nb_segment_attribs > 0 ) {
//                        std::vector< SerializedAttribute< GME::FACET > >
//                        segment_attribs( nb_segment_attribs ) ;
//                        for( index_t i = 0; i < nb_segment_attribs; i++ ) {
//                            segment_attribs[ i ].bind(
//                                L.facet_attribute_manager(), in_.field(
//                                    1 + 2 * i ), in_.field( 2 + 2 * i ), L.nb_cells() ) ;
//                        }
//                        for( index_t i = 0; i < L.nb_cells(); i++ ) {
//                            in_.get_line() ;
//                            in_.get_fields() ;
//                            serialize_read_attributes( in_, 1, in_.field_as_uint(
//                                    0 ), segment_attribs ) ;
//                        }
//                    }

                    // Set the corners - they can be the same
                    add_element_boundary( cur_element,
                        find_corner( model(), vertices.front() ) ) ;
                    add_element_boundary( cur_element,
                        find_corner( model(), vertices.back() ) ) ;

                    // Finally we have the in_boundary information
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    ringmesh_assert( file_line_.field_matches( 0, "IN_BOUNDARY" ) ) ;
                    for( index_t b = 1; b < file_line_.nb_fields(); b++ ) {
                        add_element_in_boundary( cur_element,
                            gme_t( GME::SURFACE, file_line_.field_as_uint( b ) ) ) ;
                    }
                }

                // Surfaces
                else if( match_type( file_line_.field( 0 ) ) == GME::SURFACE ) {
                    index_t id = file_line_.field_as_uint( 1 ) ;
                    gme_t cur_element( GME::SURFACE, id ) ;

                    // Read the surface vertices and their attributes
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    ringmesh_assert(
                        file_line_.field_matches( 0, "SURFACE_VERTICES" ) ) ;
                    index_t nb_vertices = file_line_.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        vec3 point( file_line_.field_as_double( 0 ),
                            file_line_.field_as_double( 1 ),
                            file_line_.field_as_double( 2 ) ) ;
                        vertices[i] = point ;
                    }

//                    in_.get_line() ;
//                    in_.get_fields() ;
//                    ringmesh_assert( in_.field_matches( 0,
//                            "SURFACE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_vertex_attribs = ( in_.nb_fields() - 1 ) / 2 ;
//
                    // Bind the vertex attributes
//                    std::vector< SerializedAttribute< GME::VERTEX > > vertex_attribs(
//                        nb_vertex_attribs ) ;
//                    for( index_t i = 0; i < nb_vertex_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            S.vertex_attribute_manager(), in_.field(
//                                1 + 2 * i ), in_.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices global ids and attributes
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in_.get_line() ;
//                        in_.get_fields() ;
//                        serialize_read_attributes( in_, 1, i, vertex_attribs ) ;
//                    }

                    // Read the surface facets
                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    ringmesh_assert(
                        file_line_.field_matches( 0, "SURFACE_CORNERS" ) ) ;
                    index_t nb_corners = file_line_.field_as_uint( 1 ) ;

                    file_line_.get_line() ;
                    file_line_.get_fields() ;
                    ringmesh_assert(
                        file_line_.field_matches( 0, "SURFACE_FACETS" ) ) ;
                    index_t nb_facets = file_line_.field_as_uint( 1 ) ;

//                    in_.get_line() ;
//                    in_.get_fields() ;
//                    ringmesh_assert( in_.field_matches( 0, "SURFACE_FACET_ATTRIBUTES" ) ) ;
//                    index_t nb_facet_attribs = ( in_.nb_fields() - 1 ) / 2 ;

                    // Bind the facet attributes
//                    std::vector< SerializedAttribute< GME::FACET > > facet_attribs(
//                        nb_facet_attribs ) ;
//                    for( index_t i = 0; i < nb_facet_attribs; i++ ) {
//                        facet_attribs[ i ].bind(
//                            S.facet_attribute_manager(), in_.field(
//                                1 + 2 * i ), in_.field( 2 + 2 * i ), nb_facets ) ;
//                    }

                    std::vector< index_t > corners( nb_corners ) ;
                    std::vector< index_t > facet_ptr( nb_facets + 1, 0 ) ;
                    index_t count_facets = 0 ;
                    for( index_t f = 0; f < nb_facets; f++ ) {
                        file_line_.get_line() ;
                        file_line_.get_fields() ;
                        index_t nb_v = file_line_.field_as_uint( 0 ) ;
                        for( index_t v = 0; v < nb_v; ++v ) {
                            corners[count_facets + v] = file_line_.field_as_uint(
                                v + 1 ) ;
                        }
                        count_facets += nb_v ;
                        facet_ptr[f + 1] = count_facets ;
//                        serialize_read_attributes( in_, nb_v + 1, f, facet_attribs ) ;
                    }

                    set_surface_geometry( cur_element.index, vertices, corners,
                        facet_ptr ) ;
                    compute_surface_adjacencies( cur_element.index ) ;
                }
            }
        }
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
                // Number of elements of a given type
                else if( match_nb_elements( file_line.field( 0 ) )
                    != GME::NO_TYPE ) {
                    // Allocate the space
                    if( file_line.nb_fields() > 1 ) {
                        create_elements( match_nb_elements( file_line.field( 0 ) ),
                            file_line.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level elements
                else if( match_high_level_type( file_line.field( 0 ) ) ) {
                    // Read this element
                    // First line : type - id - name - geol_feature
                    if( file_line.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    GME::TYPE t = match_type( file_line.field( 0 ) ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t element( t, id ) ;
                    set_element_name( element, file_line.field( 2 ) ) ;
                    set_element_geol_feature( element,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        add_element_child( element,
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
                    gme_t element( GME::REGION, id ) ;
                    set_element_name( element, file_line.field( 2 ) ) ;
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

                        add_element_boundary( element, gme_t( GME::SURFACE, s ),
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

                        add_element_boundary( gme_t( GME::REGION, NO_ID ),
                            gme_t( GME::SURFACE, s ), side ) ;
                    }
                }

            }
        }
    }
    void GeoModelBuilderGM::load_file()
    {

        // ZLib fails to load a relative path (case of a directory
        // inside the geomodel relative path).
        std::string normalized_path = GEO::FileSystem::normalized_path(
        filename_.c_str() ) ;
        unzFile uz = unzOpen( normalized_path.c_str() ) ;
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
            load_elements( type, uz ) ;
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
                    const GeoModelMeshElement& cur_gme = model_.mesh_element( t,
                        id ) ;
                    gme_t cur_gme_type( t, id ) ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_element_in_boundary( cur_gme_type,
                            gme_t( cur_gme.in_boundary_type( t ),
                                file_line.field_as_uint( in_b ) ) ) ;
                    }
                }

            }

        }

    }

    void GeoModelBuilderGM::load_elements( GME::TYPE gme_t, unzFile& uz )
    {
        for( index_t el = 0; el < model_.nb_elements( gme_t ); el++ ) {
            std::string file_to_extract_and_load ;
            build_string_for_geo_model_element_export( gme_t, el,
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
            GEO::Mesh cur_mesh ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            GEO::Logger::instance()->set_minimal( true ) ;
            GEO::mesh_load( str_try, cur_mesh, flags ) ;
            assign_mesh_to_element( cur_mesh,
                model_.element( gme_t, el ).gme_id() ) ;
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
