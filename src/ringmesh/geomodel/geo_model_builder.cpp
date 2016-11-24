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

#include <ringmesh/geomodel/geo_model_builder.h>

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

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/box3d.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/geogram_extension/geogram_mesh_repair.h>
#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_api.h>
#include <ringmesh/geomodel/geo_model_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file ringmesh/geomodel/geo_model_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

/*! @todo Split All functions of geo_model_builder.cpp into smaller functions
 * Split this file into at least 4 files.
 */

namespace {
    using namespace RINGMesh ;

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

    void get_internal_borders(
        const GeoModelMeshEntity& entity,
        std::set< index_t >& internal_borders )
    {
        for( index_t i = 0; i < entity.nb_boundaries(); ++i ) {
            const GeoModelMeshEntity& border = entity.boundary( i ) ;
            if( border.is_inside_border( entity ) ) {
                internal_borders.insert( border.index() ) ;
            }
        }
    }

    bool inexact_equal( const vec3& v1, const vec3& v2, double epsilon )
    {
        return length( v2 - v1 ) < epsilon ;
    }

    /*************************************************************************/

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
        const GeoModelMeshVertices& model_vertices = L.model().mesh.vertices ;
        bool equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i]
                != model_vertices.model_vertex_id( L.gme_id(), i ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) {
            return true ;
        }
        // If the order is the other one
        equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i]
                != model_vertices.model_vertex_id( L.gme_id(),
                    L.nb_vertices() - i - 1 ) ) {
                equal = false ;
                break ;
            }
        }
        return equal ;
    }

    /*!
     * Finds a facet and its edge index that are colocalised with an edge
     * defined by its two model vertex indices
     * @param[in] surface the surface where to find the facet
     * @param[in] v0 the first vertex of the edge
     * @param[in] v1 the second vertex of the edge
     * @param[out] f the found facet index
     * @param[out] e the found edge index
     * @return True if the facet and the edge indices are found
     * @todo RENAME these parameters and split in smaller functions !! [JP]
     */
    bool find_facet_from_edge_vertices(
        const Surface& surface,
        const vec3& v0,
        const vec3& v1,
        index_t& f,
        index_t& e )
    {
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
            nb_neighbors = surface.facet_colocater_ann().get_neighbors( v_bary,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                f = neighbors[i] ;
                for( index_t j = 0; j < surface.nb_mesh_element_vertices( f );
                    j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( f, j ), v0,
                        surface.model().epsilon() ) ) {
                        index_t j_next = surface.next_facet_vertex_index( f, j ) ;
                        if( inexact_equal( surface.mesh_element_vertex( f, j_next ),
                            v1, surface.model().epsilon() ) ) {
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

    bool are_cell_facet_and_facet_equal(
        const Region& region,
        index_t cell,
        index_t cell_facet,
        const Surface& surface,
        index_t facet )
    {
        index_t nb_cell_facet_vertices = region.nb_cell_facet_vertices( cell,
            cell_facet ) ;
        index_t nb_facet_vertices = surface.nb_mesh_element_vertices( facet ) ;
        if( nb_cell_facet_vertices != nb_facet_vertices ) {
            return false ;
        }
        vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell,
            cell_facet ) ;
        vec3 facet_barycenter = surface.mesh_element_barycenter( facet ) ;
        return inexact_equal( cell_facet_barycenter, facet_barycenter,
            region.model().epsilon() ) ;
    }

    bool find_cell_facet_from_facet(
        const Region& region,
        const Surface& surface,
        index_t facet,
        index_t& cell,
        index_t& cell_facet )
    {
        vec3 v_bary = surface.mesh_element_barycenter( facet ) ;
        index_t nb_neighbors = std::min( index_t( 5 ), region.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, region.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = region.cell_colocater_ann().get_neighbors( v_bary,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                cell = neighbors[i] ;
                for( cell_facet = 0 ; cell_facet < region.nb_cell_facets( cell );
                    cell_facet++ ) {
                    if( are_cell_facet_and_facet_equal( region, cell, cell_facet,
                        surface, facet ) ) {
                        return true ;
                    }
                }
            }
        } while( region.nb_mesh_elements() != cur_neighbor ) ;

        cell = NO_ID ;
        cell_facet = NO_ID ;
        return false ;
    }

    bool find_facet_from_vertex(
        const Surface& surface,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
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
            nb_neighbors = surface.facet_colocater_ann().get_neighbors( v,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i] ;
                for( index_t j = 0;
                    j < surface.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( surface.mesh_element_vertex( element_id, j ),
                        v, surface.model().epsilon() ) ) {
                        vertex_id = surface.mesh_element_vertex_index( element_id,
                            j ) ;
                        return true ;
                    }
                }
            }
        } while( surface.nb_mesh_elements() != cur_neighbor ) ;

        element_id = NO_ID ;
        vertex_id = NO_ID ;
        return false ;
    }

    bool find_cell_from_vertex(
        const Region& entity,
        const vec3& v,
        index_t& element_id,
        index_t& vertex_id )
    {
        index_t nb_neighbors = std::min( index_t( 5 ), entity.nb_mesh_elements() ) ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            cur_neighbor = std::min( cur_neighbor, entity.nb_mesh_elements() ) ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = entity.cell_colocater_ann().get_neighbors( v,
                cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                element_id = neighbors[i] ;
                for( index_t j = 0;
                    j < entity.nb_mesh_element_vertices( element_id ); j++ ) {
                    if( inexact_equal( entity.mesh_element_vertex( element_id, j ),
                        v, entity.model().epsilon() ) ) {
                        vertex_id = entity.mesh_element_vertex_index( element_id,
                            j ) ;
                        return true ;
                    }
                }
            }
        } while( entity.nb_mesh_elements() != cur_neighbor ) ;

        element_id = NO_ID ;
        vertex_id = NO_ID ;
        return false ;
    }

    index_t edge_index_from_facet_and_edge_vertex_indices(
        const Surface& surface,
        index_t f,
        const vec3& v0,
        const vec3& v1 )
    {
        for( index_t v = 0; v < surface.nb_mesh_element_vertices( f ); v++ ) {
            if( !inexact_equal( surface.mesh_element_vertex( f, v ), v0,
                surface.model().epsilon() ) ) continue ;
            index_t prev_v = surface.prev_facet_vertex_index( f, v ) ;
            index_t next_v = surface.next_facet_vertex_index( f, v ) ;
            if( inexact_equal( surface.mesh_element_vertex( f, prev_v ), v1,
                surface.model().epsilon() ) ) {
                return prev_v ;
            } else if( inexact_equal( surface.mesh_element_vertex( f, next_v ), v1,
                surface.model().epsilon() ) ) {
                return v ;
            }
        }
        return NO_ID ;
    }

    index_t cell_facet_index_from_cell_and_facet(
        const Region& region,
        index_t cell,
        const Surface& surface,
        index_t facet )
    {
        vec3 facet_barycenter = surface.mesh_element_barycenter( facet ) ;
        for( index_t f = 0; f < region.nb_cell_facets( cell ); f++ ) {
            vec3 cell_facet_barycenter = region.cell_facet_barycenter( cell, f ) ;
            if( inexact_equal( cell_facet_barycenter, facet_barycenter,
                surface.model().epsilon() ) ) {
                return f ;
            }
        }
        return NO_ID ;
    }

    bool is_corner_to_duplicate(
        const GeoModel& geomodel,
        index_t corner_id,
        index_t surface_id )
    {
        const Corner& corner = geomodel.corner( corner_id ) ;
        std::vector< index_t > surfaces ;
        for( index_t l = 0; l < corner.nb_in_boundary(); l++ ) {
            const GMME& line = corner.in_boundary( l ) ;
            for( index_t s = 0; s < line.nb_in_boundary(); s++ ) {
                surfaces.push_back( line.in_boundary_gme( s ).index ) ;
            }
        }
        return std::count( surfaces.begin(), surfaces.end(), surface_id ) != 2 ;
    }

    void get_sorted_incident_surfaces(
        const GeoModelMeshEntity& E,
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
                ringmesh_assert( dot( N_, e1 ) < global_epsilon ) ;

                vec3 B = 0.5 * p1 + 0.5 * p0 ;
                vec3 p2B = p2 - B ;
                B_A_ = normalize( p2B - dot( p2B, e1 ) * e1 ) ;

                ringmesh_assert( dot( B_A_, e1 ) < global_epsilon ) ;
                ringmesh_assert( B_A_.length() > global_epsilon ) ;
            }

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
            index_t triangle_id = static_cast< index_t >( triangles_.size() ) ;
            triangles_.push_back(
                TriangleToSort( triangle_id, surface_index, p0, p1, p2 ) ) ;
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

            while( t != start ) {
                ringmesh_assert( t != NO_ID ) ;
                if( !is_visited( t ) ) {
                    std::vector< index_t > cur_adjacent_surfaces ;
                    get_adjacent_surfaces( t, cur_adjacent_surfaces ) ;

                    if( equal_to_line_adjacent_surfaces( cur_adjacent_surfaces ) ) {
                        visit_border_triangles_on_same_edge( t ) ;
                        add_border_triangle_vertices_to_line( t, backward ) ;
                    } else {
                        // Adjacent surface changed
                        break ;
                    }
                } else {
                    // Adjacent surface changed
                    break ;
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
            const GeoModelMeshVertices& model_vertices = geomodel_.mesh.vertices ;
            for( index_t s = 0; s < geomodel_.nb_surfaces(); ++s ) {
                const Surface& S = geomodel_.surface( s ) ;
                for( index_t f = 0; f < S.nb_mesh_elements(); ++f ) {
                    for( index_t v = 0; v < S.nb_mesh_element_vertices( f ); ++v ) {
                        if( S.is_on_border( f, v ) ) {
                            index_t vertex = model_vertices.model_vertex_id(
                                S.gme_id(), f, v ) ;
                            index_t next_vertex = model_vertices.model_vertex_id(
                                S.gme_id(), f, S.next_facet_vertex_index( f, v ) ) ;
                            index_t previous_vertex = model_vertices.model_vertex_id(
                                S.gme_id(), f, S.prev_facet_vertex_index( f, v ) ) ;
                            border_triangles_.push_back(
                                BorderTriangle( s, f, vertex, next_vertex,
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

            const GeoModelMeshVertices& model_vertices = geomodel_.mesh.vertices ;

            // Gets the next edge on border in the Surface
            index_t f = border_triangle.facet_ ;
            std::vector< index_t > possible_v0_id ;
            model_vertices.mesh_entity_vertex_id( S.gme_id(), border_triangle.v0_,
                possible_v0_id ) ;
            ringmesh_assert( !possible_v0_id.empty() ) ;
            index_t v0_id = NO_ID ;
            for( index_t i = 0; i < possible_v0_id.size(); i++ ) {
                if( S.vertex_index_in_facet( f, possible_v0_id[i] ) != NO_ID ) {
                    v0_id = possible_v0_id[i] ;
                }
            }
            ringmesh_assert( v0_id != NO_ID ) ;
            index_t v0_id_in_facet = S.vertex_index_in_facet( f,
                 v0_id ) ;
            ringmesh_assert( v0_id_in_facet != NO_ID ) ;

            index_t next_f = NO_ID ;
            index_t next_f_v0 = NO_ID ;
            index_t next_f_v1 = NO_ID ;

            if( !backward ) {
                S.next_on_border( f, v0_id_in_facet, next_f, next_f_v0 ) ;
                ringmesh_assert( next_f_v0 != NO_ID ) ;
                next_f_v1 = S.next_facet_vertex_index( next_f, next_f_v0 ) ;
            } else {
                S.prev_on_border( f, v0_id_in_facet, next_f, next_f_v0 ) ;
                ringmesh_assert( next_f_v0 != NO_ID ) ;
                next_f_v1 = S.next_facet_vertex_index( next_f, next_f_v0 ) ;
            }

            // Finds the BorderTriangle that is corresponding to this
            // It must exist and there is only one
            BorderTriangle bait( border_triangle.surface_, next_f,
                model_vertices.model_vertex_id( S.gme_id(), next_f, next_f_v0 ),
                model_vertices.model_vertex_id( S.gme_id(), next_f, next_f_v1 ),
                NO_ID ) ;
            index_t result = static_cast< index_t >( std::lower_bound(
                border_triangles_.begin(), border_triangles_.end(), bait )
                - border_triangles_.begin() ) ;

            ringmesh_assert( border_triangles_[result].same_edge( bait ) ) ;
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
                index_t nb_copied = static_cast< index_t >( line_vertices.end()
                    - line_vertices.begin() - i ) ;
                std::copy( line_vertices.begin() + 1, line_vertices.begin() + i + 1,
                    shuffled_vertices.begin() + nb_copied ) ;
                line_vertices = shuffled_vertices ;
                break ;
            }
        }
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
        copy_meshes( geomodel, Corner::type_name_static() ) ;
        copy_meshes( geomodel, Line::type_name_static() ) ;
        copy_meshes( geomodel, Surface::type_name_static() ) ;
        copy_meshes( geomodel, Region::type_name_static() ) ;
    }

    void GeoModelBuilder::copy_meshes(
        const GeoModel& from,
        const std::string& entity_type )
    {
        for( index_t i = 0; i < model().nb_mesh_entities( entity_type ); ++i ) {
            copy_mesh( from, gme_t( entity_type, i ) ) ;
        }
    }

    void GeoModelBuilder::copy_mesh( const GeoModel& from, const gme_t& mesh_entity )
    {
        const GeoModelMeshEntity& from_E = from.mesh_entity( mesh_entity ) ;
        assign_mesh_to_entity( *from_E.mesh_, mesh_entity ) ;
    }

    void GeoModelBuilder::assign_mesh_to_entity(
        const MeshBase& mesh,
        const gme_t& to )
    {
        GeoModelMeshEntity& E = mesh_entity( to ) ;
        MeshBaseBuilder* builder = E.mesh_->get_mesh_base_builder() ;
        builder->copy( mesh, true, GEO::MESH_ALL_ELEMENTS ) ;
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
            result = create_mesh_entity< Corner >() ;
            set_corner( result.index, point ) ;
        }
        return result ;
    }

    gme_t GeoModelBuilder::find_or_create_corner( index_t model_point_id )
    {
        gme_t result = find_corner( model(), model_point_id ) ;
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >() ;
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
            result = create_mesh_entity< Line >() ;
            set_line( result.index, vertices ) ;

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_mesh_entity_boundary( result,
                find_or_create_corner( vertices.front() ).index ) ;
            add_mesh_entity_boundary( result,
                find_or_create_corner( vertices.back() ).index ) ;
        }
        return result ;
    }

    /*!
     * @brief Finds or creates a line knowing its topological adjacencies
     */
    gme_t GeoModelBuilder::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        const gme_t& first_corner,
        const gme_t& second_corner )
    {
        for( index_t i = 0; i < model().nb_lines(); ++i ) {
            const Line& line = model().line( i ) ;
            gme_t c0 = line.boundary_gme( 0 ) ;
            gme_t c1 = line.boundary_gme( 1 ) ;

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) ) {
                std::vector< index_t > cur_adjacent_surfaces ;
                get_sorted_incident_surfaces( line, cur_adjacent_surfaces ) ;
                if( cur_adjacent_surfaces.size() == sorted_adjacent_surfaces.size()
                    && std::equal( cur_adjacent_surfaces.begin(),
                        cur_adjacent_surfaces.end(),
                        sorted_adjacent_surfaces.begin() ) ) {
                    return line.gme_id() ;
                }
            }
        }
        return create_mesh_entity< Line >() ;
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
    void GeoModelBuilder::set_mesh_entity_vertex(
        const gme_t& t,
        index_t v,
        const vec3& point,
        bool update )
    {
        GeoModelMeshEntity& E = mesh_entity( t ) ;
        GeoModelMeshVertices& model_vertices = model().mesh.vertices ;
        ringmesh_assert( v < E.nb_vertices() ) ;
        if( update ) {
            model_vertices.update_point(
                model_vertices.model_vertex_id( E.gme_id(), v ), point ) ;
        } else {
            MeshBaseBuilder* builder = E.mesh_->get_mesh_base_builder() ;
            builder->set_vertex( v, point ) ;
        }
    }

    /*!
     * @brief Sets the geometrical position of a vertex from a model vertex
     * @param[in] entity_id Entity index
     * @param[in] v Index of the vertex to modify
     * @param[in] model_vertex Index in GeoModelMeshVertices of the vertex giving
     *                     the new position
     */
    void GeoModelBuilder::set_mesh_entity_vertex(
        const gme_t& entity_id,
        index_t v,
        index_t model_vertex )
    {
        GeoModelMeshVertices& model_vertices = model().mesh.vertices ;
        set_mesh_entity_vertex( entity_id, v,
            model_vertices.vertex( model_vertex ), false ) ;

        ringmesh_assert( v < mesh_entity( entity_id ).nb_vertices() ) ;
        model_vertices.update_vertex_mapping( entity_id, v, model_vertex ) ;
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     * @param[in] id Entity index
     * @param[in] points Geometric positions of the vertices to add
     * @param[in] clear If true the mesh is cleared, keeping its attributes
     */
    void GeoModelBuilder::set_mesh_entity_vertices(
        const gme_t& id,
        const std::vector< vec3 >& points,
        bool clear )
    {
        GeoModelMeshEntity& E = mesh_entity( id ) ;
        MeshBaseBuilder* builder = E.mesh_->get_mesh_base_builder() ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true ) ;
        }
        if( !points.empty() ) {
            index_t nb_points = static_cast< index_t >( points.size() ) ;
            index_t start = builder->create_vertices( nb_points ) ;
            for( index_t v = 0; v < nb_points; v++ ) {
                builder->set_vertex( start + v, points[v] ) ;
            }
        }
    }

    /*!
     * @brief Creates new vertices to the mesh
     * @param[in] id Entity index
     * @param[in] nb_vertices Number of vertices to create
     * @return the first vertex index created
     */
    index_t GeoModelBuilder::create_mesh_entity_vertices(
        const gme_t& id,
        index_t nb_vertices )
    {
        GeoModelMeshEntity& E = mesh_entity( id ) ;
        MeshBaseBuilder* builder = E.mesh_->get_mesh_base_builder() ;
        return builder->create_vertices( nb_vertices ) ;
    }

    /*!
     * @brief Adds vertices to the mesh
     * @details No update of the model vertices is done
     *
     * @param[in] id Entity index
     * @param[in] model_vertices Geometric positions of the vertices to add
     * @param[in] clear If true the mesh if cleared, keeping its attributes
     */
    void GeoModelBuilder::set_mesh_entity_vertices(
        const gme_t& entity_id,
        const std::vector< index_t >& model_vertices,
        bool clear )
    {
        GeoModelMeshEntity& E = mesh_entity( entity_id ) ;
        MeshBaseBuilder* builder = E.mesh_->get_mesh_base_builder() ;
        // Clear the mesh, but keep the attributes and the space
        if( clear ) {
            builder->clear( true, true ) ;
        }
        index_t nb_model_vertices = static_cast< index_t >( model_vertices.size() ) ;
        index_t start = builder->create_vertices( nb_model_vertices ) ;
        for( index_t v = 0; v < nb_model_vertices; v++ ) {
            set_mesh_entity_vertex( entity_id, start + v, model_vertices[v] ) ;
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
        set_mesh_entity_vertex( gme_t( Corner::type_name_static(), corner_id ), 0,
            point, false ) ;
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
        set_mesh_entity_vertices( gme_t( Line::type_name_static(), line_id ),
            vertices, true ) ;

        Line& line = dynamic_cast< Line& >( mesh_entity( Line::type_name_static(),
            line_id ) ) ;
        Mesh1DBuilder* builder = line.mesh1d_->get_mesh1d_builder() ;
        for( index_t e = 1; e < line.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e ) ;
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
        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            points, true ) ;
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
        set_mesh_entity_vertices( gme_t( Region::type_name_static(), region_id ),
            points, true ) ;
        assign_region_tet_mesh( region_id, tetras ) ;
    }

    /*!
     * @brief Sets the vertex for a Corner. Store the info in the geomodel vertices
     *
     * @param[in] corner_id Index of the corner
     * @param[in] unique_vertex Index of the vertex in the model
     */
    void GeoModelBuilder::set_corner( index_t corner_id, index_t model_vertex_id )
    {
        set_mesh_entity_vertex( gme_t( Corner::type_name_static(), corner_id ), 0,
            model_vertex_id ) ;
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
        GeoModelMeshEntity& E = mesh_entity( Line::type_name_static(), line_id ) ;

        ringmesh_assert( E.nb_vertices() == 0 ) ; // If there are already some vertices
        // we are doomed because they are not removed
        /// @todo Do this test for all others set_something
        set_mesh_entity_vertices( E.gme_id(), unique_vertices, clear_vertices ) ;

        Line& line = dynamic_cast< Line& >( E ) ;
        Mesh1DBuilder* builder = line.mesh1d_->get_mesh1d_builder() ;
        for( index_t e = 1; e < E.nb_vertices(); e++ ) {
            builder->create_edge( e - 1, e ) ;
        }
    }

    /*!
     * @brief Sets the vertices and facets for a surface
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
        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            model_vertex_ids, false ) ;
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

        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            vertices, false ) ;
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

        set_mesh_entity_vertices( gme_t( Surface::type_name_static(), surface_id ),
            vertices, false ) ;

        assign_surface_triangle_mesh( surface_id, new_triangle_corners,
            adjacent_triangles ) ;
    }

    void GeoModelBuilder::set_surface_element_geometry(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = mesh_entity( Surface::type_name_static(),
            surface_id ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;

        for( index_t facet_vertex = 0; facet_vertex < corners.size();
            facet_vertex++ ) {
            builder->set_facet_vertex( facet_id, facet_vertex,
                corners[facet_vertex] ) ;
        }
    }

    void GeoModelBuilder::set_surface_element_adjacency(
        index_t surface_id,
        index_t facet_id,
        const std::vector< index_t >& adjacents )
    {
        GeoModelMeshEntity& E = mesh_entity( Surface::type_name_static(),
            surface_id ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;

        for( index_t facet_edge = 0; facet_edge < adjacents.size(); facet_edge++ ) {
            builder->set_facet_adjacent( facet_id, facet_edge,
                adjacents[facet_edge] ) ;
        }
    }

    void GeoModelBuilder::set_region_geometry(
        index_t region_id,
        const std::vector< index_t >& tet_corners )
    {
        std::vector< index_t > vertices ;
        std::vector< index_t > new_tet_corners( tet_corners ) ;
        get_entity_vertices_and_update_corners( new_tet_corners, vertices ) ;

        set_mesh_entity_vertices( gme_t( Region::type_name_static(), region_id ),
            vertices, false ) ;
        assign_region_tet_mesh( region_id, new_tet_corners ) ;
    }

    void GeoModelBuilder::set_region_element_geometry(
        index_t region_id,
        index_t cell_id,
        const std::vector< index_t >& corners )
    {
        GeoModelMeshEntity& E = mesh_entity( Region::type_name_static(),
            region_id ) ;

        Region& region = dynamic_cast< Region& >( E ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;

        for( index_t cell_vertex = 0; cell_vertex < corners.size(); cell_vertex++ ) {
            builder->set_cell_vertex( cell_id, cell_vertex, corners[cell_vertex] ) ;
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
        GeoModelMeshEntity& E = mesh_entity( Region::type_name_static(),
            region_id ) ;
        Region& region = dynamic_cast< Region& >( E ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        return builder->create_cells( nb_cells, type ) ;
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
        GeoModelMeshEntity& E = mesh_entity( Surface::type_name_static(),
            surface_id ) ;
        Surface& surface = dynamic_cast< Surface& >( E ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        return builder->create_facet_polygon( vertex_indices ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices )
    {
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        builder->assign_facet_triangle_mesh( triangle_vertices, true ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_surface_triangle_mesh(
        index_t surface_id,
        const std::vector< index_t >& triangle_vertices,
        const std::vector< index_t >& adjacent_triangles )
    {
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        builder->assign_facet_triangle_mesh( triangle_vertices, true ) ;

        ringmesh_assert( adjacent_triangles.size() == surface.nb_mesh_elements() * 3 ) ;
        for( index_t f = 0; f < surface.nb_mesh_elements(); f++ ) {
            for( index_t v = 0; v < 3; v++ ) {
                builder->set_facet_adjacent( f, v, adjacent_triangles[3 * f + v] ) ;
            }
        }
    }

    void GeoModelBuilder::set_surface_facet_adjacencies(
        index_t surface_id,
        const std::vector< index_t >& facet_ids,
        const std::vector< index_t >& edge_ids,
        const std::vector< index_t >& adjacent_triangles )
    {
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        ringmesh_assert( facet_ids.size() == edge_ids.size() &&
            facet_ids.size() == adjacent_triangles.size() ) ;
        for( index_t i = 0; i < facet_ids.size(); ++i ) {
            builder->set_facet_adjacent( facet_ids[i], edge_ids[i],
                adjacent_triangles[i] ) ;
        }
    }

    void GeoModelBuilder::assign_surface_mesh_facets(
        index_t surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        ringmesh_assert( surface.nb_vertices() > 0 ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        builder->create_facet_polygons( facets, facet_ptr ) ;
        compute_surface_adjacencies( surface_id ) ;
    }

    void GeoModelBuilder::assign_region_tet_mesh(
        index_t region_id,
        const std::vector< index_t >& tet_vertices )
    {
        Region& region = dynamic_cast< Region& >( mesh_entity(
            Region::type_name_static(), region_id ) ) ;
        ringmesh_assert( region.nb_vertices() > 0 ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        builder->assign_cell_tet_mesh( tet_vertices, true ) ;
        builder->connect_cells() ;
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
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;

        for( index_t f = 0; f < surface.nb_mesh_elements(); f++ ) {
            for( index_t v = 0; v < surface.nb_mesh_element_vertices( f ); v++ ) {
                builder->set_facet_adjacent( f, v, NO_ID ) ;
            }
        }
        builder->connect_facets() ;
    }

    void GeoModelBuilder::compute_region_adjacencies( index_t region_id )
    {
        Region& region = dynamic_cast< Region& >( mesh_entity(
            Region::type_name_static(), region_id ) ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
            for( index_t f = 0; f < region.nb_cell_facets( c ); f++ ) {
                builder->set_cell_adjacent( c, f, NO_ID ) ;
            }
        }
        builder->connect_cells() ;
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
        index_t surface_out_id )
    {
        Surface& surface_out = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_out_id ) ) ;
        Mesh2DBuilder* builder = surface_out.mesh2d_->get_mesh2d_builder() ;
        builder->triangulate( *surface_in.mesh2d_ ) ;
    }

    /*
     * @brief Resets the adjacencies for all Surface facets adjacent to the Line
     * @return The number of disconnection done
     * @pre All the edges of the Line are edges of at least one facet of the Surface
     */
    index_t GeoModelBuilder::disconnect_surface_facets_along_line_edges(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < model().nb_surfaces() ) ;
        ringmesh_assert( line_id < model().nb_lines() ) ;

        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            gme_t( Surface::type_name_static(), surface_id ) ) ) ;
        const Line& line = model().line( line_id ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        index_t nb_disconnected_edges = 0 ;
        for( index_t i = 0; i < line.nb_mesh_elements(); ++i ) {
            const vec3& p0 = line.vertex( i ) ;
            const vec3& p1 = line.vertex( i + 1 ) ;

            index_t f = NO_ID ;
            index_t e = NO_ID ;
            bool found = find_facet_from_edge_vertices( surface, p0, p1, f, e ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && f != NO_ID && e != NO_ID ) ;

            index_t adj_f = surface.facet_adjacent_index( f, e ) ;
            if( adj_f != NO_ID ) {
                index_t adj_e = edge_index_from_facet_and_edge_vertex_indices(
                    surface, adj_f, p0, p1 ) ;
                ringmesh_assert( adj_e != NO_ID ) ;
                builder->set_facet_adjacent( f, e, NO_ID ) ;
                builder->set_facet_adjacent( adj_f, adj_e, NO_ID ) ;
                nb_disconnected_edges++ ;
            }
        }
        return nb_disconnected_edges ;
    }

    index_t GeoModelBuilder::disconnect_region_cells_along_surface_facets(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < model().nb_regions() ) ;
        ringmesh_assert( surface_id < model().nb_surfaces() ) ;

        Region& region = dynamic_cast< Region& >( mesh_entity(
            gme_t( Region::type_name_static(), region_id ) ) ) ;
        const Surface& surface = model().surface( surface_id ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        index_t nb_disconnected_facets = 0 ;
        for( index_t facet = 0; facet < surface.nb_mesh_elements(); ++facet ) {
            index_t cell = NO_ID ;
            index_t cell_facet = NO_ID ;
            bool found = find_cell_facet_from_facet( region, surface, facet, cell,
                cell_facet ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && cell != NO_ID && cell_facet != NO_ID ) ;

            index_t adj_cell = region.cell_adjacent_index( cell, cell_facet ) ;
            if( adj_cell != NO_ID ) {
                index_t adj_cell_facet = cell_facet_index_from_cell_and_facet(
                    region, adj_cell, surface, facet ) ;
                ringmesh_assert( adj_cell_facet != NO_ID ) ;
                builder->set_cell_adjacent( cell, cell_facet, NO_ID ) ;
                builder->set_cell_adjacent( adj_cell, adj_cell_facet, NO_ID ) ;
                nb_disconnected_facets++ ;
            }
        }
        return nb_disconnected_facets ;
    }

    struct ElementVertex {
        index_t element_ ;
        index_t vertex_ ;
    } ;

    /*!
     * @brief Duplicates the surface vertices along the fake boundary
     * (NO_ID adjacencies but shared vertices) and duplicate the vertices
     */
    void GeoModelBuilder::duplicate_surface_vertices_along_line(
        index_t surface_id,
        index_t line_id )
    {
        ringmesh_assert( surface_id < model().nb_surfaces() ) ;
        ringmesh_assert( line_id < model().nb_lines() ) ;

        gme_t surface_gme( Surface::type_name_static(), surface_id ) ;
        Surface& surface = dynamic_cast< Surface& >( mesh_entity( surface_gme ) ) ;
        const Line& line = model().line( line_id ) ;

        std::vector< ElementVertex > facet_vertices( line.nb_vertices() ) ;
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v ) ;

            index_t& facet_vertex = facet_vertices[v].vertex_ ;
            index_t& facet = facet_vertices[v].element_ ;
            bool found = find_facet_from_vertex( surface, p, facet, facet_vertex ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && facet != NO_ID && facet_vertex != NO_ID ) ;
        }

        index_t vertex_id = create_mesh_entity_vertices( surface_gme,
            line.nb_vertices() ) ;
        Mesh2DBuilder* surface_mesh_builder = surface.mesh2d_->get_mesh2d_builder() ;
        for( index_t v = 0; v < line.nb_vertices(); v++ ) {
            const vec3& p = line.vertex( v ) ;
            const index_t& facet_vertex = facet_vertices[v].vertex_ ;
            const index_t& facet = facet_vertices[v].element_ ;

            std::vector< index_t > facets ;
            surface.facets_around_vertex( facet_vertex, facets, false, facet ) ;
            update_facet_vertex( surface, facets, facet_vertex, vertex_id ) ;
            surface_mesh_builder->set_vertex( vertex_id, p ) ;
            vertex_id++ ;
        }
    }

    void GeoModelBuilder::duplicate_region_vertices_along_surface(
        index_t region_id,
        index_t surface_id )
    {
        ringmesh_assert( region_id < model().nb_regions() ) ;
        ringmesh_assert( surface_id < model().nb_surfaces() ) ;

        gme_t region_gme( Region::type_name_static(), region_id ) ;
        Region& region = dynamic_cast< Region& >( mesh_entity( region_gme ) ) ;
        const Surface& surface = model().surface( surface_id ) ;

        std::vector< ElementVertex > cell_vertices( surface.nb_vertices() ) ;
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v ) ;

            index_t& cell = cell_vertices[v].element_ ;
            index_t& cell_vertex = cell_vertices[v].vertex_ ;
            bool found = find_cell_from_vertex( region, p, cell, cell_vertex ) ;
            ringmesh_unused( found ) ;
            ringmesh_assert( found && cell != NO_ID && cell_vertex != NO_ID ) ;

        }

        index_t vertex_id = create_mesh_entity_vertices( region_gme,
            surface.nb_vertices() ) ;
        Mesh3DBuilder* region_mesh_builder = region.mesh3d_->get_mesh3d_builder() ;
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            const vec3& p = surface.vertex( v ) ;
            const index_t& cell = cell_vertices[v].element_ ;
            const index_t& cell_vertex = cell_vertices[v].vertex_ ;

            std::vector< index_t > cells ;
            region.cells_around_vertex( cell_vertex, cells, cell ) ;
            update_cell_vertex( region, cells, cell_vertex, vertex_id ) ;
            region_mesh_builder->set_vertex( vertex_id, p ) ;
            vertex_id++ ;
        }
    }

    /*!
     * @brief Cuts a Surface along a Line assuming that the edges of the Line are edges of the Surface
     * @pre Surface is not already cut. Line L does not cut the Surface S into 2 connected components.
     * @todo Add a test for this function.
     */
    void GeoModelBuilder::cut_surface_by_line( index_t surface_id, index_t line_id )
    {
        index_t nb_disconnected_edges = disconnect_surface_facets_along_line_edges(
            surface_id, line_id ) ;
        if( nb_disconnected_edges > 0 ) {
            duplicate_surface_vertices_along_line( surface_id, line_id ) ;
        }
    }

    void GeoModelBuilder::cut_region_by_surface(
        index_t region_id,
        index_t surface_id )
    {
        index_t nb_disconnected_facets =
            disconnect_region_cells_along_surface_facets( region_id, surface_id ) ;
        if( nb_disconnected_facets > 0 ) {
            duplicate_region_vertices_along_surface( region_id, surface_id ) ;
        }
    }

    void GeoModelBuilder::end_model()
    {
        if( model().name().empty() ) {
            set_model_name( "model_default_name" ) ;
        }
        // Get out if the model has no surface
        if( model().nb_surfaces() == 0 ) {
            print_geomodel( model() ) ;
            throw RINGMeshException( "GeoModel", "The GeoModel has no surface" ) ;
        }

        complete_entity_connectivity() ;

        // Fill geological feature if they are missing
        complete_mesh_entities_geol_feature_from_first_parent(
            Corner::type_name_static() ) ;
        complete_mesh_entities_geol_feature_from_first_parent(
            Line::type_name_static() ) ;
        complete_mesh_entities_geol_feature_from_first_parent(
            Surface::type_name_static() ) ;
        complete_mesh_entities_geol_feature_from_first_parent(
            Region::type_name_static() ) ;

        for( index_t i = 0; i < model().nb_geological_entity_types(); i++ ) {
            const std::string& type = model().geological_entity_type( i ) ;
            complete_geological_entities_geol_feature_from_first_child( type ) ;
        }

        cut_surfaces_by_internal_lines() ;
        cut_regions_by_internal_surfaces() ;
        compute_universe() ;

        // Deliberate clear of the model vertices used for model building
        model().mesh.vertices.clear() ;
    }

    void GeoModelBuilder::cut_surfaces_by_internal_lines()
    {
        for( index_t s = 0; s < model().nb_surfaces(); s++ ) {
            Surface& surface = dynamic_cast< Surface& >( mesh_entity(
                Surface::type_name_static(), s ) ) ;
            std::set< index_t > cutting_lines ;
            get_internal_borders( surface, cutting_lines ) ;
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                cut_surface_by_line( s, *it ) ;
            }
            if( !cutting_lines.empty() ) {
                Mesh2DBuilder* surface_mesh_builder =
                    surface.mesh2d_->get_mesh2d_builder() ;
                surface_mesh_builder->remove_isolated_vertices() ;
            }
        }
    }

    void GeoModelBuilder::cut_regions_by_internal_surfaces()
    {
        for( index_t r = 0; r < model().nb_regions(); r++ ) {
            Region& region = dynamic_cast< Region& >( mesh_entity(
                Region::type_name_static(), r ) ) ;
            if( region.nb_mesh_elements() == 0 ) continue ;
            std::set< index_t > cutting_surfaces ;
            get_internal_borders( region, cutting_surfaces ) ;
            for( std::set< index_t >::iterator it = cutting_surfaces.begin();
                it != cutting_surfaces.end(); ++it ) {
                cut_region_by_surface( r, *it ) ;
            }
            if( !cutting_surfaces.empty() ) {
                Mesh3DBuilder* region_mesh_builder =
                    region.mesh3d_->get_mesh3d_builder() ;
                region_mesh_builder->remove_isolated_vertices() ;
            }
        }
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

            gme_t first_corner = find_or_create_corner( vertices.front() ) ;
            gme_t second_corner = find_or_create_corner( vertices.back() ) ;
            const std::vector< index_t >& adjacent_surfaces =
                line_computer.adjacent_surfaces() ;

            index_t backup_nb_lines = model().nb_lines() ;
            gme_t line_index = find_or_create_line( adjacent_surfaces, first_corner,
                second_corner ) ;

            bool created_line = model().nb_lines() != backup_nb_lines ;
            if( created_line ) {
                set_line( line_index.index, vertices ) ;

                for( index_t j = 0; j < adjacent_surfaces.size(); ++j ) {
                    add_mesh_entity_in_boundary( line_index, adjacent_surfaces[j] ) ;
                }
                add_mesh_entity_boundary( line_index, first_corner.index ) ;
                add_mesh_entity_boundary( line_index, second_corner.index ) ;

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
        fill_mesh_entities_boundaries( Surface::type_name_static() ) ;

        // Sort surfaces around the contacts
        for( index_t i = 0; i < regions_info_.size(); ++i ) {
            regions_info_[i]->sort() ;
        }

        if( model().nb_surfaces() == 1 ) {
            if( model().nb_lines() != 0 ) {
                Logger::err( "GeoModel" )
                    << "The unique surface provided to build the model has boundaries "
                    << std::endl ;
                return false ;
            } else {
                /// If there is only one surface, its inside is set to be 
                /// the + side. No further check.
                bool inside = true ;
                // Create the region - set the surface on its boundaries
                gme_t region_id = create_mesh_entity< Region >() ;
                add_mesh_entity_boundary( region_id, 0, inside ) ;

                // Set universe boundary
                add_universe_boundary( 0, !inside ) ;
            }
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2 * model().nb_surfaces(),
                NO_ID ) ;

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
                gme_t cur_region_id = create_mesh_entity< Region >() ;
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
                    add_mesh_entity_boundary( cur_region_id, s.first, s.second ) ;
                    surf_2_region[s_id] = cur_region_id.index ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp =
                        !s.second == true ? 2 * s.first : 2 * s.first + 1 ;
                    if( surf_2_region[s_id_opp] == NO_ID ) {
                        S.push( std::pair< index_t, bool >( s.first, !s.second ) ) ;
                    }
                    // For each contact, push the next oriented surface that is in the same region
                    const Surface& surface = model().surface( s.first ) ;
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
                Logger::err( "GeoModel" )
                    << "Small bubble regions were skipped at model building "
                    << std::endl ;
                // Or, most probably, we have a problem before
                ringmesh_assert( false ) ;
            }

            compute_universe() ;
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
                add_universe_boundary( cur_region.boundary( i ).index(),
                    cur_region.side( i ) ) ;
            }
            std::set< gme_t > to_erase ;
            to_erase.insert( cur_region.gme_id() ) ;
            remove_mesh_entities( to_erase ) ;
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

    void GeoModelBuilder::build_contacts()
    {
        std::vector< std::set< gme_t > > interfaces ;
        for( index_t i = 0; i < model().nb_lines(); ++i ) {
            const Line& L = model().line( i ) ;
            std::set< gme_t > cur_interfaces ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                const GeoModelMeshEntity& S = L.in_boundary( j ) ;
                gme_t parent_interface = S.parent_gme(
                    Interface::type_name_static() ) ;
                cur_interfaces.insert( parent_interface ) ;
            }
            gme_t contact_id ;
            for( index_t j = 0; j < interfaces.size(); ++j ) {
                if( cur_interfaces.size() == interfaces[j].size()
                    && std::equal( cur_interfaces.begin(), cur_interfaces.end(),
                        interfaces[j].begin() ) ) {
                    contact_id = gme_t( Contact::type_name_static(), j ) ;
                    break ;
                }
            }
            if( !contact_id.is_defined() ) {
                contact_id = create_geological_entity(
                    Contact::type_name_static() ) ;
                ringmesh_assert( contact_id.index == interfaces.size() ) ;
                interfaces.push_back( cur_interfaces ) ;
                // Create a name for this contact
                std::string name = "contact" ;
                for( std::set< gme_t >::const_iterator it( cur_interfaces.begin() );
                    it != cur_interfaces.end(); ++it ) {
                    name += "_" ;
                    name += model().geological_entity( *it ).name() ;
                }
                set_entity_name( contact_id, name ) ;
            }
            add_geological_entity_child( contact_id, i ) ;
        }
    }

    void GeoModelBuilder::invert_surface_normals( index_t surface_id )
    {
        ringmesh_assert( surface_id < model().nb_surfaces() ) ;
        Surface& surface = dynamic_cast<Surface&> (
            mesh_entity( Surface::type_name_static(), surface_id ) ) ;
        Mesh2DBuilder* surf_mesh_builder = surface.mesh2d_->get_mesh2d_builder() ;
        surf_mesh_builder->invert_normals() ;
    }

    void GeoModelBuilder::update_facet_vertex(
        Surface& surface,
        const std::vector< index_t >& facets,
        index_t old_vertex,
        index_t new_vertex )
    {
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t cur_f = facets[i] ;
            for( index_t cur_v = 0;
                cur_v < surface.nb_mesh_element_vertices( cur_f ); cur_v++ ) {
                if( surface.mesh_element_vertex_index( cur_f, cur_v )
                    == old_vertex ) {
                    builder->set_facet_vertex( cur_f, cur_v, new_vertex ) ;
                }
            }
        }
    }

    index_t GeoModelBuilder::get_connected_commoponents(
        const gme_t& gmme_id,
        GEO::vector< index_t >& component ) const
    {
        const MeshBase* mesh = model().mesh_entity( gmme_id ).mesh_ ;
        ringmesh_assert( mesh != nil ) ;
        return mesh->get_connected_commoponents( component ) ;
    }

    void GeoModelBuilder::update_cell_vertex(
        Region& region,
        const std::vector< index_t >& cells,
        index_t old_vertex,
        index_t new_vertex )
    {
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        for( index_t i = 0; i < cells.size(); ++i ) {
            index_t cur_c = cells[i] ;
            for( index_t cur_v = 0; cur_v < region.nb_mesh_element_vertices( cur_c );
                cur_v++ ) {
                if( region.mesh_element_vertex_index( cur_c, cur_v )
                    == old_vertex ) {
                    builder->set_cell_vertex( cur_c, cur_v, new_vertex ) ;
                }
            }
        }
    }

    void GeoModelBuilder::delete_mesh_entity_mesh( const gme_t& E_id )
    {
        MeshBase& M = *mesh_entity( E_id ).mesh_ ;
        MeshBaseBuilder* builder = M.get_mesh_base_builder() ;
        builder->clear( true, false ) ;
    }

    void GeoModelBuilder::delete_mesh_entity_isolated_vertices( const gme_t& E_id )
    {
        MeshBase& M = *mesh_entity( E_id ).mesh_ ;
        MeshBaseBuilder* builder = M.get_mesh_base_builder() ;
        builder->remove_isolated_vertices() ;
    }

    void GeoModelBuilder::delete_mesh_entity_vertices(
        const gme_t& E_id,
        GEO::vector< index_t >& to_delete,
        bool remove_isolated_vertices )
    {
        MeshBase& M = *mesh_entity( E_id ).mesh_ ;
        MeshBaseBuilder* builder = M.get_mesh_base_builder() ;
        builder->delete_vertices( to_delete, remove_isolated_vertices ) ;
    }

    void GeoModelBuilder::delete_corner_vertex( index_t corner_id )
    {
        gme_t corner( Corner::type_name_static(), corner_id ) ;
        GEO::vector< index_t > to_delete ;
        to_delete.push_back( 1 ) ;
        delete_mesh_entity_vertices( corner, to_delete, false ) ;
    }
    void GeoModelBuilder::delete_line_edges(
        index_t line_id,
        GEO::vector< index_t >& to_delete,
        bool remove_isolated_vertices )
    {
        Line& line = dynamic_cast< Line& >( mesh_entity( Line::type_name_static(),
            line_id ) ) ;
        Mesh1DBuilder* builder = line.mesh1d_->get_mesh1d_builder() ;
        builder->delete_edges( to_delete, remove_isolated_vertices ) ;
    }
    void GeoModelBuilder::delete_surface_facets(
        index_t surface_id,
        GEO::vector< index_t >& to_delete,
        bool remove_isolated_vertices )
    {
        Surface& surface = dynamic_cast< Surface& >( mesh_entity(
            Surface::type_name_static(), surface_id ) ) ;
        Mesh2DBuilder* builder = surface.mesh2d_->get_mesh2d_builder() ;
        builder->delete_facets( to_delete, remove_isolated_vertices ) ;
    }
    void GeoModelBuilder::delete_region_cells(
        index_t region_id,
        GEO::vector< index_t >& to_delete,
        bool remove_isolated_vertices )
    {
        Region& region = dynamic_cast< Region& >( mesh_entity(
            Region::type_name_static(), region_id ) ) ;
        Mesh3DBuilder* builder = region.mesh3d_->get_mesh3d_builder() ;
        builder->delete_cells( to_delete, remove_isolated_vertices ) ;
    }

    void GeoModelBuilder::compute_universe()
    {
        if( model().universe().nb_boundaries() != 0 ) return ;
        std::vector< bool > is_surface_universe_boundary( model().nb_surfaces(),
            false ) ;
        std::vector< bool > surface_side( model().nb_surfaces() ) ;
        for( index_t r = 0; r < model().nb_regions(); r++ ) {
            const Region& region = model().region( r ) ;
            for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                index_t surface_id = region.boundary_gme( s ).index ;
                is_surface_universe_boundary[surface_id] =
                    !is_surface_universe_boundary[surface_id] ;
                surface_side[surface_id] = region.side( s ) ;
            }
        }

        for( index_t s = 0; s < model().nb_surfaces(); s++ ) {
            if( !is_surface_universe_boundary[s] ) continue ;
            add_universe_boundary( s, surface_side[s] ) ;
        }
    }

    /*************************************************************************/
    GeoModelBuilderFile::GeoModelBuilderFile(
        GeoModel& model,
        const std::string& filename )
        : GeoModelBuilder( model ), filename_( filename )
    {

    }

} // namespace
