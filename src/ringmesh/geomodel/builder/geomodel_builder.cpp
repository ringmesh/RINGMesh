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

#include <ringmesh/geomodel/builder/geomodel_builder.h>

#include <stack>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/pimpl_impl.h>
#include <ringmesh/geomodel/builder/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder_remove.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>

#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/surface_mesh.h>

/*!
 * @file ringmesh/geomodel/builder/geomodel_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel, index_t geomodel_point_id )
    {
        const auto& geomodel_vertices = geomodel.mesh.vertices;
        const auto& vertices =
            geomodel_vertices.gme_vertices( geomodel_point_id );
        for( const auto& vertex : vertices )
        {
            if( vertex.gmme.type() == Corner< DIMENSION >::type_name_static() )
            {
                return vertex.gmme;
            }
        }
        return gmme_id{};
    }

    /*!
     * @brief Returns true if the Line< DIMENSION > has exactly the given
     * vertices
     * @todo Reimplement using std::iterators
     */
    template < index_t DIMENSION >
    bool line_equal( const Line< DIMENSION >& line,
        const std::vector< index_t >& rhs_vertices )
    {
        if( line.nb_vertices() != rhs_vertices.size() )
        {
            return false;
        }
        const auto& geomodel_vertices = line.geomodel().mesh.vertices;
        bool equal{ true };
        for( auto i : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( line.gmme(), i ) )
            {
                equal = false;
                break;
            }
        }
        if( equal )
        {
            return true;
        }
        // If the order is the other one
        equal = true;
        for( auto i : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id(
                       line.gmme(), line.nb_vertices() - i - 1 ) )
            {
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
    template < index_t DIMENSION >
    void reorder_closed_line_vertices_to_start_at_corner(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< index_t >& line_vertices )
    {
        if( geomodel.nb_corners() == 0 || line_vertices.empty() )
        {
            return;
        }
        for( auto i : range( 1, line_vertices.size() - 1 ) )
        {
            auto corner = find_corner( geomodel, line_vertices[i] );
            if( corner.is_defined() )
            {
                line_vertices.pop_back();
                std::rotate( line_vertices.begin(), line_vertices.begin() + i,
                    line_vertices.end() );
                line_vertices.push_back( line_vertices.front() );
                break;
            }
        }
    }

    /*!
     * @brief Stores the vertices of a polygon which is on the boundary of a
     * Surface
     */
    struct BorderPolygon
    {
        BorderPolygon( index_t surface_in,
            index_t polygon_in,
            index_t v0_in,
            index_t v1_in )
            : v0( v0_in ),
              v1( v1_in ),
              surface( surface_in ),
              polygon( polygon_in )
        {
        }

        /*!
         * @brief Compares two polygons so those sharing an edge follow one
         * another
         */
        bool operator<( const BorderPolygon& rhs ) const
        {
            if( std::min( v0, v1 ) != std::min( rhs.v0, rhs.v1 ) )
            {
                return std::min( v0, v1 ) < std::min( rhs.v0, rhs.v1 );
            }
            if( std::max( v0, v1 ) != std::max( rhs.v0, rhs.v1 ) )
            {
                return std::max( v0, v1 ) < std::max( rhs.v0, rhs.v1 );
            }
            if( surface != rhs.surface )
            {
                return surface < rhs.surface;
            }
            return polygon < rhs.polygon;
        }

        bool same_edge( const BorderPolygon& rhs ) const
        {
            return std::min( v0, v1 ) == std::min( rhs.v0, rhs.v1 )
                   && std::max( v0, v1 ) == std::max( rhs.v0, rhs.v1 );
        }

        /// Indices of the points in the geomodel.
        /// The edge v0v1 is on the one on the boundary.
        index_t v0;
        index_t v1;
        /// Index of the surface containing the polygon
        index_t surface;
        /// Index of the polygon in the surface
        index_t polygon;
    };

    template < index_t DIMENSION >
    class CommonDataFromGeoModelSurfaces
    {
        ringmesh_disable_copy_and_move( CommonDataFromGeoModelSurfaces );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    protected:
        explicit CommonDataFromGeoModelSurfaces(
            const GeoModel< DIMENSION >& geomodel )
            : geomodel_( geomodel )
        {
            const auto& geomodel_vertices = geomodel_.mesh.vertices;
            for( const auto& surface : geomodel_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                auto S_id = surface.gmme();
                for( auto p : range( surface.nb_mesh_elements() ) )
                {
                    for( auto v :
                        range( surface.nb_mesh_element_vertices( p ) ) )
                    {
                        if( mesh.is_edge_on_border( { p, v } ) )
                        {
                            auto vertex = geomodel_vertices.geomodel_vertex_id(
                                S_id, { p, v } );
                            auto next_vertex =
                                geomodel_vertices.geomodel_vertex_id( S_id,
                                    mesh.next_polygon_vertex( { p, v } ) );
                            border_polygons_.emplace_back(
                                surface.index(), p, vertex, next_vertex );
                        }
                    }
                }
            }
            std::sort( border_polygons_.begin(), border_polygons_.end() );
        }
        ~CommonDataFromGeoModelSurfaces() = default;

        bool have_border_polygons_same_boundary_edge(
            index_t t0, index_t t1 ) const
        {
            return this->border_polygons_[t0].same_edge(
                this->border_polygons_[t1] );
        }

    protected:
        const GeoModel< DIMENSION >& geomodel_;
        /// All the polygons on a boundary of all the Surfaces of the GeoModel
        std::vector< BorderPolygon > border_polygons_;
    };

    ALIAS_2D_AND_3D( CommonDataFromGeoModelSurfaces );

    /*!
     * @brief Utility class to sort a set of oriented polygons around a common
     * edge
     */
    class GeoModelRegionFromSurfaces
    {
    public:
        /*!
         * @brief A polygon to sort around an edge, see
         * GeoModelRegionFromSurfaces
         * @details This polygon belongs to a mesh connected component
         * identified by its index.
         */
        struct PolygonToSort
        {
            /*!
             * @param index Index of this PolygonToSort in
             * GeoModelRegionFromSurfaces
             * @param surface_index Index of the Surface
             * @param normal normal to the polygon - normalized vector
             * @param p0 point of the polygon
             * @param p1 point of the polygon
             */
            PolygonToSort( index_t index,
                index_t surface_index,
                const vec3& normal,
                const vec3& p0,
                const vec3& p1 )
                : index_( index ), surface_index_( surface_index ), N_( normal )
            {
                ringmesh_assert( p0 != p1 );
                auto e1 = normalize( p1 - p0 );
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

            /// Normal to the edge p0p1 in the plane defined by the polygon -
            /// normalized
            vec3 B_A_;

            // Values filled by sorting function in GeoModelRegionFromSurfaces
            double angle_{ -max_float64() };
            bool side_{ false };
        };

        void add_polygon_edge( index_t surface_index,
            const vec3& normal,
            const vec3& p0,
            const vec3& p1 )
        {
            auto polygon_id = static_cast< index_t >( polygons_.size() );
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
            auto q = normalize( axis ) * std::sin( 0.5 * angle );
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
            double inv_w{ 1.0 / result.w };
            return { result.x * inv_w, result.y * inv_w, result.z * inv_w };
        }

        void sort()
        {
            ringmesh_assert( !polygons_.empty() );

            std::pair< index_t, bool > default_pair( index_t( -1 ), false );
            sorted_polygons_.resize( 2 * polygons_.size(), default_pair );

            // If there is only one Polygon to sort, nothing to do
            if( polygons_.size() == 1 )
            {
                sorted_polygons_[0] = std::pair< index_t, bool >(
                    polygons_[0].surface_index_, true );
                sorted_polygons_[1] = std::pair< index_t, bool >(
                    polygons_[0].surface_index_, false );
                return;
            }

            // Initialization
            // We start on the plus (true) side of the first Polygon
            sorted_polygons_[0] =
                std::pair< index_t, bool >( polygons_[0].surface_index_, true );

            // Reference vectors with wich angles will be computed
            vec3 N_ref{ polygons_[0].N_ };
            vec3 B_A_ref{ polygons_[0].B_A_ };
            vec3 Ax_ref{ normalize( cross( B_A_ref, N_ref ) ) };

            // The minus (false) side of the start polygon will the last one
            // encountered
            polygons_[0].angle_ = 2 * M_PI;
            polygons_[0].side_ = false;

            for( auto i : range( 1, polygons_.size() ) )
            {
                auto& cur = polygons_[i];
                // Computes the angle RADIANS between the reference and the
                // current
                // polygon
                double cos = dot( B_A_ref, cur.B_A_ );
                // Remove invalid values
                if( cos < -1 )
                {
                    cos = -1;
                }
                else if( cos > 1 )
                {
                    cos = 1;
                }
                cur.angle_ = std::acos( cos );
                // Put the angle between PI and 2PI if necessary
                if( dot( cross( B_A_ref, cur.B_A_ ), Ax_ref ) < 0. )
                {
                    cur.angle_ = 2 * M_PI - cur.angle_;
                }

                // Get the side of the surface first encountered
                // when rotating in the N_ref direction
                auto N_rotate = rotate( Ax_ref, -cur.angle_, cur.N_ );
                cur.side_ = dot( N_rotate, N_ref ) < 0;
            }

            // Sorts the Surfaces according to the angle
            std::sort( polygons_.begin(), polygons_.end() );

            // Fills the sorted surfaces adding the side
            index_t it{ 1 };
            for( auto& cur : polygons_ )
            {
                if( cur.index_ == 0 )
                { // The last to add
                    sorted_polygons_[it].first = cur.surface_index_;
                    sorted_polygons_[it].second = cur.side_;
                }
                else
                {
                    sorted_polygons_[it].first = cur.surface_index_;
                    sorted_polygons_[it].second = cur.side_;
                    sorted_polygons_[it + 1].first = cur.surface_index_;
                    sorted_polygons_[it + 1].second = !cur.side_;
                    it += 2;
                }
            }
            // All the surfaces must have been sorted
            ringmesh_assert( std::count( sorted_polygons_.begin(),
                                 sorted_polygons_.end(), default_pair )
                             == 0 );
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
            for( auto i : range( sorted_polygons_.size() ) )
            {
                if( sorted_polygons_[i] == in )
                {
                    if( i == sorted_polygons_.size() - 1 )
                    {
                        return sorted_polygons_[sorted_polygons_.size() - 2];
                    }
                    if( i == 0 )
                    {
                        return sorted_polygons_[1];
                    }

                    if( sorted_polygons_[i + 1].first
                        == sorted_polygons_[i].first )
                    {
                        // The next has the same surface id, check its sign
                        if( sorted_polygons_[i + 1].second
                            != sorted_polygons_[i].second )
                        {
                            return sorted_polygons_[i - 1];
                        }
                        // Sign is the same
                        return sorted_polygons_[i + 1];
                    }
                    ringmesh_assert( sorted_polygons_[i - 1].first
                                     == sorted_polygons_[i].first );
                    if( sorted_polygons_[i - 1].second
                        != sorted_polygons_[i].second )
                    {
                        return sorted_polygons_[i + 1];
                    }
                    return sorted_polygons_[i - 1];
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

    class RegionTopologyFromGeoModelSurfaces
        : public CommonDataFromGeoModelSurfaces< 3 >
    {
    public:
        explicit RegionTopologyFromGeoModelSurfaces(
            const GeoModel3D& geomodel )
            : CommonDataFromGeoModelSurfaces3D( geomodel ),
              region_info_( geomodel.nb_lines() )
        {
        }

        void compute_region_info()
        {
            const auto& vertices = this->geomodel_.mesh.vertices;
            for( const auto& line : geomodel_.lines() )
            {
                BorderPolygon line_border{ NO_ID, NO_ID,
                    vertices.geomodel_vertex_id( line.gmme(), 0 ),
                    vertices.geomodel_vertex_id( line.gmme(), 1 ) };
                for( const auto& border : this->border_polygons_ )
                {
                    if( line_border.same_edge( border ) )
                    {
                        auto surface_id = border.surface;
                        region_info_[line.index()].add_polygon_edge( surface_id,
                            this->geomodel_.surface( surface_id )
                                .mesh()
                                .polygon_normal( border.polygon ),
                            vertices.vertex( border.v0 ),
                            vertices.vertex( border.v1 ) );
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

    struct LineDefinition
    {
        bool is_closed()
        {
            return vertices_.front() == vertices_.back();
        }

        template < index_t DIMENSION >
        void reoder_closed_line_vertices(
            const GeoModel< DIMENSION >& geomodel )
        {
            // Vertices can begin and end at any vertex
            reorder_closed_line_vertices_to_start_at_corner(
                geomodel, vertices_ );
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
     * @details All the polygons on the boundaries are classified as belonging
     * to a Line< DIMENSION >
     * Two neighboring edges on the boundary belong to the same Line< DIMENSION
     * > if their incident Surfaces
     * are the same.
     */
    template < index_t DIMENSION >
    class LineGeometryFromGeoModelSurfaces
        : public CommonDataFromGeoModelSurfaces< DIMENSION >
    {
    public:
        /*!
         * @param geomodel GeoModel providing the Surfaces
         * @param collect_region_info If true, information needed to determine
         * closed Regions
         *  from a GeoModel Surfaces are collected for each Line< DIMENSION >.
         */
        explicit LineGeometryFromGeoModelSurfaces(
            const GeoModel< DIMENSION >& geomodel )
            : CommonDataFromGeoModelSurfaces< DIMENSION >( geomodel ),
              cur_border_polygon_( 0 )
        {
            visited_.resize( this->border_polygons_.size(), false );
            const auto& geomodel_vertices = this->geomodel_.mesh.vertices;
            surfaces_around_vertices_.resize( geomodel_vertices.nb() );
            for( const auto& border : this->border_polygons_ )
            {
                fill_surfaces_around_vertex( border.v0 );
                fill_surfaces_around_vertex( border.v1 );
            }
        }

        /*!
         * @brief Tries to compute a new line and returns true if one was.
         * @details To use in a while conditional loop, since the number of
         * lines
         * is considered unknown.
         */
        bool compute_next_line_geometry()
        {
            for( ; cur_border_polygon_ < this->border_polygons_.size();
                 cur_border_polygon_++ )
            {
                if( is_visited( cur_border_polygon_ ) )
                {
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
        void fill_surfaces_around_vertex( index_t vertex )
        {
            auto& surfaces = surfaces_around_vertices_[vertex];
            if( !surfaces.empty() )
            {
                return;
            }
            auto gme_vertices = this->geomodel_.mesh.vertices.gme_type_vertices(
                surface_type_name_static(), vertex );
            for( const auto& gme_vertex : gme_vertices )
            {
                surfaces.push_back( gme_vertex.gmme.index() );
            }
        }

        void compute_line_geometry()
        {
            visit_border_polygons_on_same_edge( cur_border_polygon_ );
            auto p0 = this->border_polygons_[cur_border_polygon_].v0;
            auto p1 = this->border_polygons_[cur_border_polygon_].v1;
            cur_line_.vertices_.push_back( p0 );
            cur_line_.vertices_.push_back( p1 );

            cur_line_.adjacent_surfaces_ =
                get_adjacent_surfaces( cur_border_polygon_ );

            bool backward{ false };
            get_one_line_vertices( backward );
            backward = true;
            get_one_line_vertices( backward );
        }
        /*!
         * @brief Gets the geometry of a line propagating in a given direction
         * @details As long as the adjacent surfaces are the same, the vertices
         * are added
         * belong to the line under construction
         */
        void get_one_line_vertices( bool backward )
        {
            ringmesh_assert(
                cur_border_polygon_ < this->border_polygons_.size() );
            ringmesh_assert( cur_border_polygon_ != NO_ID );

            auto t = get_next_border_polygon( cur_border_polygon_, backward );
            ringmesh_assert( t != NO_ID );
            while( propagate( t, backward ) )
            {
                visit_border_polygons_on_same_edge( t );
                add_border_polygon_vertices_to_line( t, backward );
                t = get_next_border_polygon( t, backward );
                ringmesh_assert( t != NO_ID );
            }
        }

        bool propagate( index_t t, bool backward ) const
        {
            return t != cur_border_polygon_ && !is_visited( t )
                   && equal_to_line_adjacent_surfaces( t )
                   && !vertex_is_on_corner( t, backward );
        }

        bool is_visited( index_t i ) const
        {
            return visited_[i];
        }

        bool vertex_is_on_corner( index_t t, bool backward ) const
        {
            const auto& surfaces = [t, backward, this] {
                if( backward )
                {
                    return surfaces_around_vertices_[this->border_polygons_[t]
                                                         .v1];
                }
                else
                {
                    return surfaces_around_vertices_[this->border_polygons_[t]
                                                         .v0];
                }
            }();
            for( auto surface : surfaces )
            {
                if( !contains( cur_line_.adjacent_surfaces_, surface ) )
                {
                    return true;
                }
            }
            return false;
        }

        bool equal_to_line_adjacent_surfaces( index_t t ) const
        {
            auto input = get_adjacent_surfaces( t );
            if( input.size() != cur_line_.adjacent_surfaces_.size() )
            {
                return false;
            }
            return std::equal( input.begin(), input.end(),
                cur_line_.adjacent_surfaces_.begin() );
        }

        void add_border_polygon_vertices_to_line(
            index_t polygon_index, bool backward )
        {
            const auto& border_polygon = this->border_polygons_[polygon_index];
            add_vertices_to_line(
                border_polygon.v0, border_polygon.v1, !backward );
        }

        void add_vertices_to_line( index_t v0, index_t v1, bool at_the_end )
        {
            auto to_add = v1;
            if( at_the_end )
            {
                auto end_vertex = cur_line_.vertices_.back();
                if( v1 == end_vertex )
                {
                    to_add = v0;
                }
                else
                {
                    ringmesh_assert( v0 == end_vertex );
                }
                cur_line_.vertices_.push_back( to_add );
            }
            else
            {
                auto first_vertex = cur_line_.vertices_.front();
                if( v1 == first_vertex )
                {
                    to_add = v0;
                }
                else
                {
                    ringmesh_assert( v0 == first_vertex );
                }
                cur_line_.vertices_.insert(
                    cur_line_.vertices_.begin(), to_add );
            }
        }

        /*!
         * @brief Gets the next BorderPolygon in the same surface
         */
        index_t get_next_border_polygon( index_t from, bool backward ) const
        {
            const auto& border_polygon = this->border_polygons_[from];
            const auto& S = this->geomodel_.surface( border_polygon.surface );
            const auto& mesh = S.mesh();
            auto surface_id = S.gmme();
            const auto& geomodel_vertices = this->geomodel_.mesh.vertices;

            // Gets the next edge on border in the Surface
            auto p = border_polygon.polygon;
            auto possible_v0_id = geomodel_vertices.mesh_entity_vertex_id(
                surface_id, border_polygon.v0 );
            ringmesh_assert( !possible_v0_id.empty() );
            index_t v0_id{ NO_ID };
            for( auto id : possible_v0_id )
            {
                if( mesh.vertex_index_in_polygon( p, id ) != NO_ID )
                {
                    v0_id = id;
                }
            }
            ringmesh_assert( v0_id != NO_ID );
            auto v0_id_in_polygon = mesh.vertex_index_in_polygon( p, v0_id );
            ringmesh_assert( v0_id_in_polygon != NO_ID );

            const PolygonLocalEdge cur_polygon_local_edge{ p,
                v0_id_in_polygon };
            PolygonLocalEdge next_polygon_local_edge0_on_border =
                backward ? mesh.prev_on_border( cur_polygon_local_edge )
                         : mesh.next_on_border( cur_polygon_local_edge );
            ringmesh_assert(
                next_polygon_local_edge0_on_border.polygon_id != NO_ID );
            ringmesh_assert(
                next_polygon_local_edge0_on_border.local_edge_id != NO_ID );

            PolygonLocalEdge next_polygon_local_edge1_on_border =
                mesh.next_polygon_vertex( next_polygon_local_edge0_on_border );
            ringmesh_assert(
                next_polygon_local_edge1_on_border.polygon_id != NO_ID );
            ringmesh_assert(
                next_polygon_local_edge1_on_border.local_edge_id != NO_ID );

            // Finds the BorderPolygon that is corresponding to this
            // It must exist and there is only one
            BorderPolygon bait{ border_polygon.surface,
                next_polygon_local_edge0_on_border.polygon_id,
                geomodel_vertices.geomodel_vertex_id(
                    surface_id, next_polygon_local_edge0_on_border ),
                geomodel_vertices.geomodel_vertex_id(
                    surface_id, next_polygon_local_edge1_on_border ) };
            auto result = find_sorted( this->border_polygons_, bait );

            ringmesh_assert( result != NO_ID );
            ringmesh_assert( this->border_polygons_[result].same_edge( bait ) );
            return result;
        }

        /*!
         * @brief Marks as visited all BorderPolygons whose first edge is equal
         * to i's first edge
         */
        void visit_border_polygons_on_same_edge( index_t border_id )
        {
            visited_[border_id] = true;
            for( auto next_border_id = border_id + 1;
                 next_border_id < this->border_polygons_.size()
                 && this->have_border_polygons_same_boundary_edge(
                        border_id, next_border_id );
                 next_border_id++ )
            {
                visited_[next_border_id] = true;
            }
            for( auto prev_border_id = border_id - 1;
                 prev_border_id != NO_ID
                 && this->have_border_polygons_same_boundary_edge(
                        border_id, prev_border_id );
                 prev_border_id-- )
            {
                visited_[prev_border_id] = true;
            }
        }

        /*!
         * @brief Gets the sorted indices of the Surfaces incident to the first
         * edge of the i-th BorderPolygon
         * @note When the surface appears twice (the line is an internal border)
         * both occurrences are kept.
         */
        std::vector< index_t > get_adjacent_surfaces( index_t border_id ) const
        {
            std::vector< index_t > adjacent_surfaces;
            adjacent_surfaces.reserve( 10 );
            adjacent_surfaces.push_back(
                this->border_polygons_[border_id].surface );

            for( auto next_border_id = border_id + 1;
                 next_border_id < this->border_polygons_.size()
                 && this->have_border_polygons_same_boundary_edge(
                        border_id, next_border_id );
                 next_border_id++ )
            {
                adjacent_surfaces.push_back(
                    this->border_polygons_[next_border_id].surface );
            }

            for( auto prev_border_id = border_id - 1;
                 prev_border_id != NO_ID
                 && this->have_border_polygons_same_boundary_edge(
                        border_id, prev_border_id );
                 prev_border_id-- )
            {
                adjacent_surfaces.push_back(
                    this->border_polygons_[prev_border_id].surface );
            }
            std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() );
            return adjacent_surfaces;
        }

    private:
        /// Flag the visited border_polygons when computing the Lines
        std::vector< bool > visited_;

        /// Surfaces around vertices (only filled for boundary vertices)
        std::vector< std::vector< index_t > > surfaces_around_vertices_;

        /// Currently computed line information
        index_t cur_border_polygon_;
        LineDefinition cur_line_;
    };

    double compute_angle_at_corner(
        const Line2D& line1, const Corner2D& corner, const Line2D& line2 )
    {
        auto corner_id = corner.gmme();
        bool corner_is_line1_first_boundary =
            ( line1.boundary_gmme( 0 ) == corner_id );
        auto alongside_vertex1 = corner_is_line1_first_boundary
                                     ? line1.vertex( 1 )
                                     : line1.vertex( line1.nb_vertices() - 2 );
        Geometry::Segment2D line1_segment(
            corner.vertex( 0 ), alongside_vertex1 );

        bool corner_is_line2_first_boundary =
            ( line2.boundary_gmme( 0 ) == corner_id );
        auto alongside_vertex2 = corner_is_line2_first_boundary
                                     ? line2.vertex( 1 )
                                     : line2.vertex( line2.nb_vertices() - 2 );
        Geometry::Segment2D line2_segment(
            corner.vertex( 0 ), alongside_vertex2 );

        return Position::segment_angle( line1_segment, line2_segment );
    }

} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBuilderInfo< DIMENSION >::GeoModelBuilderInfo(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : builder_( builder ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
    }

    template < index_t DIMENSION >
    GeoModelBuilderBase< DIMENSION >::GeoModelBuilderBase(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : topology( builder, geomodel ),
          geometry( builder, geomodel ),
          geology( builder, geomodel ),
          remove( builder, geomodel ),
          info( builder, geomodel ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
    }

    template < index_t DIMENSION >
    void GeoModelBuilderBase< DIMENSION >::build_corners_from_lines()
    {
        std::vector< vecn< DIMENSION > > point_extremities;
        point_extremities.reserve( geomodel_.nb_lines() * DIMENSION );
        for( const auto& line : geomodel_.lines() )
        {
            point_extremities.push_back( line.vertex( 0 ) );
            point_extremities.push_back(
                line.vertex( line.nb_vertices() - 1 ) );
        }

        NNSearch< DIMENSION > nn_search( point_extremities );
        std::vector< index_t > index_map;
        std::vector< vecn< DIMENSION > > unique_points;
        std::tie( std::ignore, index_map, unique_points ) =
            nn_search.get_colocated_index_mapping_and_unique_points(
                geomodel_.epsilon() );

        topology.create_mesh_entities( corner_type_name_static(),
            static_cast< index_t >( unique_points.size() ) );
        for( index_t c : range( geomodel_.nb_corners() ) )
        {
            geometry.set_corner( c, unique_points[c] );
        }
        index_t index = 0;
        for( const auto& line : geomodel_.lines() )
        {
            gmme_id line_id = line.gmme();
            index_t point0 = index_map[index++];
            gmme_id corner0( corner_type_name_static(), point0 );
            index_t point1 = index_map[index++];
            gmme_id corner1( corner_type_name_static(), point1 );
            topology.add_mesh_entity_boundary_relation( line_id, corner0 );
            topology.add_mesh_entity_boundary_relation( line_id, corner1 );

            // Update line vertex extremities with corner coordinates
            geometry.set_mesh_entity_vertex(
                line_id, 0, unique_points[point0], false );
            geometry.set_mesh_entity_vertex(
                line_id, line.nb_vertices() - 1, unique_points[point1], false );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderBase<
        DIMENSION >::build_lines_and_corners_from_surfaces()
    {
        LineGeometryFromGeoModelSurfaces< DIMENSION > line_computer(
            geomodel_ );

        while( line_computer.compute_next_line_geometry() )
        {
            auto line = line_computer.current_line();

            if( line.is_closed() )
            {
                line.reoder_closed_line_vertices( geomodel_ );
            }

            auto& vertices = line.vertices_;
            auto first_corner =
                topology.find_or_create_corner( vertices.front() );
            auto second_corner =
                topology.find_or_create_corner( vertices.back() );

            const auto& adjacent_surfaces = line.adjacent_surfaces_;
            auto backup_nb_lines = geomodel_.nb_lines();
            auto line_index = topology.find_or_create_line(
                adjacent_surfaces, first_corner, second_corner );

            bool created_line{ geomodel_.nb_lines() != backup_nb_lines };
            if( created_line )
            {
                geometry.set_line( line_index.index(), vertices );

                for( auto j : adjacent_surfaces )
                {
                    gmme_id surface_index{ surface_type_name_static(), j };
                    topology.add_mesh_entity_boundary_relation(
                        surface_index, line_index );
                }
                topology.add_mesh_entity_boundary_relation(
                    line_index, first_corner );
                topology.add_mesh_entity_boundary_relation(
                    line_index, second_corner );
            }
            else
            {
                bool same_geometry{ line_equal(
                    geomodel_.line( line_index.index() ), vertices ) };
                if( !same_geometry )
                {
                    geometry.set_line( line_index.index(), vertices );
                }
            }
        }
        geometry.clear_geomodel_mesh();
    }

    template < index_t DIMENSION >
    void GeoModelBuilderBase< DIMENSION >::end_geomodel()
    {
        if( geomodel_.name().empty() )
        {
            info.set_geomodel_name( "model_default_name" );
        }

        cut_geomodel_on_internal_boundaries();
        print_geomodel( geomodel_ );
    }

    class GeoModelBuilder< 2 >::Impl
    {
    public:
        Impl( GeoModel2D& geomodel ) : geomodel_( geomodel ) {}

        struct OrientedLine
        {
            OrientedLine( index_t line_index, bool line_side )
                : index( line_index ), side( line_side )
            {
            }

            bool operator==( const OrientedLine& rhs )
            {
                return index == rhs.index && side == rhs.side;
            }

            bool operator!=( const OrientedLine& rhs )
            {
                return !operator==( rhs );
            }
            index_t index;
            bool side;
        };

        struct LineIncidentSurfacePair
        {
            LineIncidentSurfacePair(
                index_t plus_surface_id, index_t minus_surface_id )
            {
                surface_ids_[0] = plus_surface_id;
                surface_ids_[1] = minus_surface_id;
            }

            LineIncidentSurfacePair() : LineIncidentSurfacePair( NO_ID, NO_ID )
            {
            }

            bool is_undetermined() const
            {
                return surface_ids_[0] == NO_ID || surface_ids_[1] == NO_ID;
            }

            index_t plus_surface_index() const
            {
                return surface_ids_[0];
            }

            index_t minus_surface_index() const
            {
                return surface_ids_[1];
            }

            index_t side_surface_index( bool side ) const
            {
                if( side )
                {
                    return plus_surface_index();
                }
                return minus_surface_index();
            }

            void set_side_surface_index( bool side, index_t surface_id )
            {
                if( side )
                {
                    surface_ids_[0] = surface_id;
                    return;
                }
                surface_ids_[1] = surface_id;
            }

        private:
            index_t surface_ids_[2];
        };

        OrientedLine get_first_undetermined_line_side(
            const std::vector< LineIncidentSurfacePair >&
                line_indicent_surfaces ) const
        {
            for( auto line : range( line_indicent_surfaces.size() ) )
            {
                if( line_indicent_surfaces[line].side_surface_index( true )
                    == NO_ID )
                {
                    return OrientedLine( line, true );
                }
                if( line_indicent_surfaces[line].side_surface_index( false )
                    == NO_ID )
                {
                    return OrientedLine( line, false );
                }
            }
            return OrientedLine( NO_ID, false );
        }

        index_t find_line_with_minimal_angle(
            const Corner2D& pivot_corner, const index_t cur_line_id ) const
        {
            double min_angle{ max_float64() };
            const auto& cur_line = this->geomodel_.line( cur_line_id );
            index_t minimal_angle_line_id{ NO_ID };
            for( auto line_itr : range( pivot_corner.nb_incident_entities() ) )
            {
                const auto& line = pivot_corner.incident_entity( line_itr );
                if( line.index() == cur_line_id )
                {
                    continue;
                }
                auto cur_angle =
                    compute_angle_at_corner( cur_line, pivot_corner, line );
                if( cur_angle < min_angle )
                {
                    min_angle = cur_angle;
                    minimal_angle_line_id = line.index();
                }
            }
            return minimal_angle_line_id;
        }

        OrientedLine get_next_surface_boundary_line(
            const OrientedLine& cur_line_and_side ) const
        {
            const auto& cur_line =
                this->geomodel_.line( cur_line_and_side.index );
            index_t cur_corner_id{
                cur_line.boundary( cur_line_and_side.side ? 1 : 0 ).index()
            };
            const auto& cur_corner = this->geomodel_.corner( cur_corner_id );
            if( cur_corner.nb_incident_entities() == 1 )
            {
                // Case the current line is an internal border of a
                // surface
                return OrientedLine(
                    cur_line_and_side.index, !cur_line_and_side.side );
            }

            index_t next_line_id = find_line_with_minimal_angle(
                cur_corner, cur_line_and_side.index );
            ringmesh_assert( next_line_id != NO_ID );
            bool next_side =
                this->geomodel_.line( next_line_id ).boundary_gmme( 0 )
                == cur_corner.gmme();
            return OrientedLine( next_line_id, next_side );
        }

        std::vector< OrientedLine > get_surface_boundaries(
            const index_t cur_surface_id,
            std::vector< LineIncidentSurfacePair >& line_indicent_surfaces )
            const
        {
            std::vector< OrientedLine > cur_surface_boundaries;
            OrientedLine first_line_and_side =
                get_first_undetermined_line_side( line_indicent_surfaces );

            OrientedLine cur_line_and_side{ first_line_and_side };

            // From the first line, the lines are walked turning around corners
            // in the same direction (clockwise). Line after line, the
            // boundaries
            // of the currently processed surface are found until the algorithm
            // has gone back to the first line.
            // By convention, if a line is walked from its first boundary
            // towards
            // its second boundary, the surface is set as incident by the + side
            // of the line.
            do
            {
                ringmesh_assert(
                    line_indicent_surfaces[cur_line_and_side.index]
                        .side_surface_index( cur_line_and_side.side )
                    == NO_ID );
                line_indicent_surfaces[cur_line_and_side.index]
                    .set_side_surface_index(
                        cur_line_and_side.side, cur_surface_id );
                cur_surface_boundaries.emplace_back( cur_line_and_side );

                cur_line_and_side =
                    get_next_surface_boundary_line( cur_line_and_side );
            } while( cur_line_and_side != first_line_and_side );
            return cur_surface_boundaries;
        }

        void find_surfaces_boundary_lines(
            std::vector< LineIncidentSurfacePair >& line_indicent_surfaces,
            std::vector< std::vector< OrientedLine > >& surface_boundary_lines )
            const
        {
            // This vector registers for each line the index of the two incident
            // surfaces
            line_indicent_surfaces.resize( this->geomodel_.nb_lines() );
            index_t surface_counter{ 0 };
            while( std::count_if( line_indicent_surfaces.begin(),
                       line_indicent_surfaces.end(),
                       []( const LineIncidentSurfacePair& line ) {
                           return line.is_undetermined();
                       } )
                   > 0 )
            {
                auto cur_surface_boundaries = get_surface_boundaries(
                    surface_counter, line_indicent_surfaces );

                surface_boundary_lines.emplace_back(
                    std::move( cur_surface_boundaries ) );
                ++surface_counter;
            }
        }

        /*!
         * This function implement a geometric test to find the englobing
         * surface index for the current line
         * @param line_id index of the current line
         * @return the englobing surface index
         * @note you may need other inputs.
         */
        index_t find_englobing_surface( index_t line_id ) const
        {
            // implement geometrical test to detect the index of the surface
            // that englobe the line
            return 0;
        }

        /*!
         * This function iterates on lines and manage internal boundaries.
         * @param line_indicent_surfaces vector of lines and incident surfaces
         * to manage.
         * @return the number of intenal boundaries.
         */
        index_t manage_internal_boundary(
            const std::vector< LineIncidentSurfacePair >&
                line_indicent_surfaces,
            const index_t nb_found_surfaces ) const
        {
            std::vector< bool > are_surfaces_hole( nb_found_surfaces, true );
            for( auto line_id : range( geomodel_.nb_lines() ) )
            {
                if( line_indicent_surfaces[line_id].plus_surface_index()
                    != line_indicent_surfaces[line_id].minus_surface_index() )
                {
                    are_surfaces_hole[line_indicent_surfaces[line_id]
                                          .plus_surface_index()] = false;
                    are_surfaces_hole[line_indicent_surfaces[line_id]
                                          .minus_surface_index()] = false;
                }
                else
                { // if a line have the same minus and plus surface index this
                  // line is an internal boundary.
                    // TODO by Emna
                    // find the englobing surface for the current line
                    // implement the following function
                    index_t englobing_surface_id =
                        find_englobing_surface( line_id );
                    // reset the index of line incident surfaces to
                    // englobing_surface_id
                    //(line_indicent_surfaces[line_id].plus_surface_index() &
                    // line_indicent_surfaces[line_id].minus_surface_index())
                }
            }
            return { static_cast< index_t >( std::count(
                are_surfaces_hole.begin(), are_surfaces_hole.end(), true ) ) };
        }
        /*!
         * This function find all line boundary indices associated to a surface
         * @param surface_id index of the current surface
         * @return a set of line index that are either the minus or the plus
         * boundary.
         */
        std::set< index_t > get_line_boundary_indices(
            index_t surface_id ) const
        {
            // Todo implement it
            return std::set< index_t >();
        }

        void manage_intrusion_surface() const
        {
            // iterate on every surfaces
            index_t surface_id;
            std::set< index_t > line_boundary_indices =
                get_line_boundary_indices( surface_id );
            // interates on every other surfaces
            index_t upper_surface_id;
            std::set< index_t > upper_line_boundary_indices =
                get_line_boundary_indices( upper_surface_id );
            // compare the two sets: line_boundary_indices and
            // upper_line_boundary_indices  if equal do:
            // find englobing_surface_index for upper_line_boundary_indices
            // set line_indicent_surfaces corresponding to every
            // upper_line_boundary_indices to
            // englobing_surface_index (be careful it correspond to eitehr plus
            // or minus side)
            return;
        }
		
        void check_internal_intrusion_or_boundaries(
            const std::vector< LineIncidentSurfacePair >&
                line_indicent_surfaces,
            const index_t nb_found_surfaces ) const
        {
            index_t nb_floating_set_of_lines = manage_internal_boundary(
                line_indicent_surfaces, nb_found_surfaces );
            if( nb_floating_set_of_lines > 0 )
            {
                // TODO by EMNA
                // update message to say that we have taken into account the
                // internal boundary. remove exception
                throw RINGMeshException( "Surface2D",
                    "During surface from corners "
                    " and lines build, ",
                    nb_floating_set_of_lines,
                    " group(s) of lines are "
                    "floating inside a surface. This is not yet handled by the "
                    "algorithm. Aborting..." );
            }
            // TODO Manage intrusions surfaces
            manage_intrusion_surface();
        }

        std::vector< vec2 > get_surface_polygon_vertices(
            const std::vector< OrientedLine >& surface_boundaries ) const
        {
            std::vector< vec2 > polygon_vertices;
            for( auto cur_surf_boundary : surface_boundaries )
            {
                const auto& cur_line =
                    this->geomodel_.line( cur_surf_boundary.index );
                for( auto vertex : range( 1, cur_line.nb_vertices() ) )
                {
                    polygon_vertices.push_back( cur_line.vertex(
                        cur_surf_boundary.side
                            ? vertex
                            : ( cur_line.nb_vertices() - 1 ) - vertex ) );
                }
            }
            return polygon_vertices;
        }

        void build_surface_polygons( GeoModelBuilder2D& builder,
            const std::vector< std::vector< OrientedLine > >&
                surface_boundary_lines ) const
        {
            for( const auto& surface_boundaries : surface_boundary_lines )
            {
                std::vector< vec2 > polygon_vertices =
                    get_surface_polygon_vertices( surface_boundaries );
                std::vector< index_t > polygon_corners(
                    polygon_vertices.size() );
                std::iota( polygon_corners.begin(), polygon_corners.end(), 0 );
                auto surface_id = builder.topology.create_mesh_entity(
                    surface_type_name_static() );
                std::vector< index_t > polygon_vertex_ptr( 2, 0 );
                polygon_vertex_ptr[1] =
                    static_cast< index_t >( polygon_corners.size() );
                builder.geometry.set_surface_geometry( surface_id.index(),
                    polygon_vertices, polygon_corners, polygon_vertex_ptr );
            }
        }
        // TODO :: TO be removed by emna because the universe will not be
        // created anymore.
        void find_exterior_and_remove_it( GeoModelBuilder2D& builder,
            std::vector< std::vector< OrientedLine > >& surface_boundary_lines )
            const
        {
            double max_surface_area{ 0 };
            index_t exterior_surface_id{ NO_ID };
            for( const auto& surface : geomodel_.surfaces() )
            {
                double surface_area{ surface.size() };
                if( surface_area > max_surface_area )
                {
                    max_surface_area = surface_area;
                    exterior_surface_id = surface.index();
                }
            }
            std::set< gmme_id > to_remove;
            to_remove.insert(
                { surface_type_name_static(), exterior_surface_id } );
            builder.remove.remove_mesh_entities( to_remove );
            surface_boundary_lines.erase(
                surface_boundary_lines.begin() + exterior_surface_id );
        }

        void set_surface_line_boundary_relationships(
            GeoModelBuilder2D& builder,
            const std::vector< std::vector< OrientedLine > >&
                surface_boundary_lines ) const
        {
            index_t surface_id{ 0 };
            for( const auto& surface_boundaries : surface_boundary_lines )
            {
                for( const auto& cur_boundary : surface_boundaries )
                {
                    builder.topology.add_surface_line_boundary_relation(
                        surface_id, cur_boundary.index, cur_boundary.side );
                }
                ++surface_id;
            }
        }

    private:
        GeoModel2D& geomodel_;
    };

    GeoModelBuilder< 2 >::GeoModelBuilder( GeoModel2D& geomodel )
        : GeoModelBuilderBase< 2 >( *this, geomodel ), impl_( geomodel )
    {
    }

    GeoModelBuilder< 2 >::~GeoModelBuilder() {}

    template <>
    void GeoModelBuilderBase< 2 >::cut_geomodel_on_internal_boundaries()
    {
        geometry.cut_surfaces_by_internal_lines();
    }

    void GeoModelBuilder< 2 >::build_surfaces_from_corners_and_lines()
    {
        if( geomodel_.nb_surfaces() > 0 )
        {
            return;
        }
        // Each side of each Line is in one Surface(+side is first)
        std::vector< Impl::LineIncidentSurfacePair > line_incident_surfaces;
        std::vector< std::vector< Impl::OrientedLine > > surface_boundary_lines;
        impl_->find_surfaces_boundary_lines(
            line_incident_surfaces, surface_boundary_lines );

        impl_->check_internal_intrusion_or_boundaries( line_incident_surfaces,
            static_cast< index_t >( surface_boundary_lines.size() ) );

        // Generate surface polygons
        impl_->build_surface_polygons( *this, surface_boundary_lines );
        impl_->find_exterior_and_remove_it( *this, surface_boundary_lines );

        // Update topology
        impl_->set_surface_line_boundary_relationships(
            *this, surface_boundary_lines );
    }

    GeoModelBuilder< 3 >::GeoModelBuilder( GeoModel3D& geomodel )
        : GeoModelBuilderBase< 3 >( *this, geomodel )
    {
    }

    template <>
    void GeoModelBuilderBase< 3 >::cut_geomodel_on_internal_boundaries()
    {
        geometry.cut_surfaces_by_internal_lines();
        geometry.cut_regions_by_internal_surfaces();
    }

    void GeoModelBuilder< 3 >::build_regions_from_lines_and_surfaces()
    {
        std::set< gmme_id > regions_to_delete;
        for( auto r : range( geomodel_.nb_regions() ) )
        {
            const Region3D& region = geomodel_.region( r );
            regions_to_delete.insert( region.gmme() );
        }
        remove.remove_mesh_entities( regions_to_delete );

        RegionTopologyFromGeoModelSurfaces region_computer{ geomodel_ };
        region_computer.compute_region_info();

        const auto& region_info = region_computer.region_info();

        if( geomodel_.nb_surfaces() < 2 || geomodel_.nb_lines() == 0 )
        {
            throw RINGMeshException( "GeoModel",
                "You need at least 1 line and 2 surfaces to use "
                "GeoModelBuilder::build_regions_from_lines_and_surfaces" );
        }

        // Each side of each Surface is in one Region( +side is first )
        std::vector< index_t > surf_2_region(
            2 * geomodel_.nb_surfaces(), NO_ID );

        // Start with the first Surface on its + side
        std::stack< std::pair< index_t, bool > > S;
        S.emplace( 0, true );

        while( !S.empty() )
        {
            auto cur = S.top();
            S.pop();
            // This side is already assigned
            if( surf_2_region[cur.second ? 2 * cur.first : 2 * cur.first + 1]
                != NO_ID )
            {
                continue;
            }
            // Create a new region
            index_t cur_region_id{ geomodel_.nb_regions() };
            topology.create_mesh_entities( region_type_name_static(), 1 );
            // Get all oriented surfaces defining this region
            std::stack< std::pair< index_t, bool > > SR;
            SR.push( cur );
            while( !SR.empty() )
            {
                auto s = SR.top();
                SR.pop();
                index_t s_id = s.second ? 2 * s.first : 2 * s.first + 1;
                // This oriented surface has already been visited
                if( surf_2_region[s_id] != NO_ID )
                {
                    continue;
                }
                // Add the surface to the current region
                topology.add_region_surface_boundary_relation(
                    cur_region_id, s.first, s.second );
                surf_2_region[s_id] = cur_region_id;

                // Check the other side of the surface and push it in S
                index_t s_id_opp = !s.second ? 2 * s.first : 2 * s.first + 1;
                if( surf_2_region[s_id_opp] == NO_ID )
                {
                    S.emplace( s.first, !s.second );
                }
                // For each contact, push the next oriented surface that is in
                // the same region
                const auto& surface = geomodel_.surface( s.first );
                for( auto i : range( surface.nb_boundaries() ) )
                {
                    const auto& n =
                        region_info[surface.boundary_gmme( i ).index()].next(
                            s );
                    index_t n_id = n.second ? 2 * n.first : 2 * n.first + 1;

                    if( surf_2_region[n_id] == NO_ID )
                    {
                        SR.push( n );
                    }
                }
            }
        }

        // Check if all the surfaces were visited
        // If not, this means that there are additionnal regions included in
        // those built
        if( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID )
            != 0 )
        {
            Logger::err( "GeoModel",
                "Small bubble regions were skipped at geomodel building " );
            // Or, most probably, we have a problem before
            ringmesh_assert_not_reached;
            /// @todo handle the region building of small bubble regions
        }

        // We need to remove from the regions_ the one corresponding
        // to the "universe", i.e., the one with the biggest volume.
        double max_volume{ -1. };
        index_t universe_id{ NO_ID };
        for( const auto& region : geomodel_.regions() )
        {
            double cur_volume = region.size();
            if( cur_volume > max_volume )
            {
                max_volume = cur_volume;
                universe_id = region.index();
            }
        }
        const auto& cur_region = geomodel_.region( universe_id );
        std::set< gmme_id > to_erase;
        to_erase.insert( cur_region.gmme() );
        remove.remove_mesh_entities( to_erase );
    }

    template class geomodel_builder_api GeoModelBuilderBase< 2 >;
    template class geomodel_builder_api GeoModelBuilderInfo< 2 >;

    template class geomodel_builder_api GeoModelBuilderBase< 3 >;
    template class geomodel_builder_api GeoModelBuilderInfo< 3 >;

} // namespace RINGMesh
