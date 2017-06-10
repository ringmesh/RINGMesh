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
    GeoModelBuilder< DIMENSION >::GeoModelBuilder( GeoModel< DIMENSION >& geomodel )
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

    template< index_t DIMENSION >
    void GeoModelBuilder< DIMENSION >::end_geomodel()
    {
        if( geomodel_.name().empty() ) {
            info.set_geomodel_name( "model_default_name" );
        }

        geometry.cut_surfaces_by_internal_lines();
        geometry.cut_regions_by_internal_surfaces();
        topology.compute_universe();

        // Deliberate clear of the geomodel vertices used for geomodel building
        geomodel_.mesh.vertices.clear();
    }

    template< index_t DIMENSION >
    void GeoModelBuilder< DIMENSION >::build_brep_regions_from_surfaces()
    {
        if( geomodel_.nb_lines() == 0 ) {
            Logger::warn( "GeoModel", "No Line in the Geomodel ", geomodel_.name(),
                ".", "Computing Lines from Surfaces..." );
            from_surfaces.build_lines_and_corners_from_surfaces();
        }

        std::vector< GeoModelRegionFromSurfaces > region_info_( geomodel_.nb_lines() );


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
                gmme_id region_id( Region < DIMENSION > ::type_name_static(),
                    geomodel_.nb_regions() );
                builder_.topology.create_mesh_entities(
                    Region < DIMENSION > ::type_name_static(), 1 );
                gmme_id surface_id( Surface < DIMENSION > ::type_name_static(), 0 );
                builder_.topology.add_mesh_entity_boundary_relation( region_id,
                    surface_id, inside );

                // Set universe boundary
                builder_.topology.add_universe_boundary( 0, !inside );
            }
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector < index_t
                > surf_2_region( 2 * geomodel_.nb_surfaces(), NO_ID );

            // Start with the first Surface on its + side
            std::stack < std::pair< index_t, bool > > S;
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
                gmme_id cur_region_id( Region < DIMENSION > ::type_name_static(),
                    geomodel_.nb_regions() );
                builder_.topology.create_mesh_entities(
                    Region < DIMENSION > ::type_name_static(), 1 );
                // Get all oriented surfaces defining this region
                std::stack < std::pair< index_t, bool > > SR;
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
                        gmme_id( Surface < DIMENSION > ::type_name_static(),
                            s.first ), s.second );
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
            std::set < gmme_id > to_erase;
            to_erase.insert( cur_region.gmme() );
            builder_.removal.remove_mesh_entities( to_erase );
        }
        return true;
    }

//    template class RINGMESH_API GeoModelBuilder< 2 > ;
    //    template class RINGMESH_API GeoModelBuilderInfo< 2 > ;
    //    template class RINGMESH_API GeoModelBuilderCopy< 2 > ;

    template class RINGMESH_API GeoModelBuilder< 3 > ;
    template class RINGMESH_API GeoModelBuilderInfo< 3 > ;
    template class RINGMESH_API GeoModelBuilderCopy< 3 > ;

} // namespace
