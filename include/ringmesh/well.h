/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_WELL_
#define __RINGMESH_WELL_

#include <ringmesh/common.h>

namespace RINGMesh {
    class BoundaryModel ;
    class Edge ;
    class Well ;
}

namespace RINGMesh {
    class RINGMESH_API WellCorner {
public:
        WellCorner(
            const vec3& point,
            signed_index_t id = - 1 )
            : point_( point ), surface_id_( id ), resolution_( - 1 )
        {
        }

        const vec3& point() const { return point_ ;}
        void set_surface_id( signed_index_t id ) { surface_id_ = id ;}
        signed_index_t surface_id() const { return surface_id_ ;}
        double resolution() const { return resolution_ ;}
        void set_resolution( double d ) { resolution_ = d ;}

private:
        vec3 point_ ;
        signed_index_t surface_id_ ;
        double resolution_ ;
    } ;

    class RINGMESH_API WellPart {
public:
        WellPart()
            : well_( nil ), id_( - 1 )
        {
            corners_.reserve( 2 ) ;
        }

        void add_corner( signed_index_t c ) { corners_.push_back( c ) ;}
        signed_index_t corner( signed_index_t c ) const
        {
            ringmesh_debug_assert( c < corners_.size() ) ;
            return corners_[ c ] ;
        }

        void add_points( const std::vector< vec3 >& p ) { points_ = p ;}
        const vec3& point( index_t p ) const
        {
            ringmesh_debug_assert( p < points_.size() ) ;
            return points_[ p ] ;
        }

        index_t nb_points() const { return points_.size() ;}
        const std::vector< vec3 >& points() const { return points_ ;}
        double resolution( index_t p ) const
        {
            ringmesh_debug_assert( p < resolutions_.size() ) ;
            return resolutions_[ p ] ;
        }

        std::vector< double >& resolutions() { return resolutions_ ;}
        void set_well( Well* well ) { well_ = well ;}
        const Well* well() const { return well_ ;}
        void set_id( signed_index_t id ) { id_ = id ;}
        signed_index_t id() const { return id_ ;}

private:
        Well* well_ ;
        signed_index_t id_ ;
        std::vector< vec3 > points_ ;
        std::vector< signed_index_t > corners_ ;
        std::vector< double > resolutions_ ;
    } ;

    class RINGMESH_API Well {
public:
        Well() {}
        void copy_corners_and_informations( Well& well ) const ;

        void set_well_in_parts()
        {
            for( index_t p = 0; p < nb_parts(); p++ ) {
                parts_[ p ].set_well( this ) ;
            }
        }

        void add_part_points(
            index_t part,
            const std::vector< vec3 >& p )
        {
            parts_[ part ].add_points( p ) ;
        }

        void get_part_edges(
            index_t p,
            std::vector< Edge >& edges ) const ;

        void get_region_edges(
            index_t p,
            std::vector< Edge >& edges ) const ;

        void add_part(
            const WellPart& part,
            index_t r )
        {
            parts_.push_back( part ) ;
            part_region_id_.push_back( r ) ;
            parts_.back().set_id( parts_.size() - 1 ) ;
        }

        signed_index_t find_or_create_corner(
            const vec3& p,
            signed_index_t id = - 1 ) ;

        signed_index_t find_corner( const vec3& p ) const ;

        const WellCorner& corner( index_t c ) const
        {
            ringmesh_debug_assert( c < corners_.size() ) ;
            return corners_[ c ] ;
        }

        const WellPart& part( index_t part ) const
        {
            ringmesh_debug_assert( part < parts_.size() ) ;
            return parts_[ part ] ;
        }

        signed_index_t part_region_id( index_t part ) const
        {
            ringmesh_debug_assert( part < nb_parts() ) ;
            return part_region_id_[ part ] ;
        }

        double part_length( index_t part ) const ;

        void all_part_vertices(
            index_t part,
            std::vector< vec3 >& vertices ) const ;

        index_t nb_corners() const { return corners_.size() ;}
        index_t nb_parts() const { return parts_.size() ;}

        void set_name( const std::string& name ) { name_ = name ;}
        const std::string& name() const { return name_ ;}

private:
        signed_index_t add_corner(
            const vec3& p,
            signed_index_t id )
        {
            corners_.push_back( WellCorner( p, id ) ) ;
            return corners_.size() - 1 ;
        }

private:
        std::vector< WellCorner > corners_ ;
        std::vector< WellPart > parts_ ;
        std::vector< index_t > part_region_id_ ;
        std::string name_ ;
    } ;

    class RINGMESH_API WellGroup {
public:
        WellGroup() ;
        virtual ~WellGroup() ;

        void get_region_edges(
            index_t region,
            std::vector< Edge >& edges ) const ;

        BoundaryModel* model() const { return model_ ;}
        void set_model( RINGMesh::BoundaryModel* model ) { model_ = model ;}
        bool is_well_already_added( const std::string& name ) const ;

        void add_well( const Well& w ) { wells_.push_back( w ) ;}

        index_t nb_wells() const { return wells_.size() ;}
        const Well& well( index_t w ) const { return wells_[ w ] ;}
        void update_wells()
        {
            for( index_t w = 0; w < nb_wells(); w++ ) {
                wells_[ w ].set_well_in_parts() ;
            }
        }

protected:
        std::vector< Well > wells_ ;
        BoundaryModel* model_ ;
    } ;
}
#endif
