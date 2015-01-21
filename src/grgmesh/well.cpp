/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <grgmesh/well.h>
#include <grgmesh/utils.h>

#include <cmath>
#include <stack>

namespace GRGMesh {

    signed_index_t Well::find_or_create_corner( const vec3& p, signed_index_t id )
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( p == corner( c ).point() ) {
                grgmesh_debug_assert( id == corner( c ).surface_id() ) ;
                grgmesh_debug_assert( id != -1 ) ;
                return c ;
            }
        }
        return add_corner( p, id ) ;
    }

    signed_index_t Well::find_corner( const vec3& p ) const
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( p == corner( c ).point() ) {
                return c ;
            }
        }
        return -1 ;
    }

    void Well::copy_corners_and_informations( Well& well ) const
    {
        well.name_ = name_ ;
        well.corners_ = corners_ ;
        well.part_region_id_ = part_region_id_ ;
        well.parts_.resize( nb_parts() ) ;
        for( index_t p = 0; p < nb_parts(); p++ ) {
            well.parts_[p].add_corner( parts_[p].corner(0) ) ;
            well.parts_[p].add_corner( parts_[p].corner(1) ) ;
            well.parts_[p].set_id( parts_[p].id() ) ;
            well.parts_[p].set_well( &well ) ;
        }
    }

    double Well::part_length( index_t p ) const
    {
        double l = 0.0 ;
        const WellPart& well_part = part( p ) ;
        vec3 prev = corner( well_part.corner(0) ).point() ;
        for( index_t p = 0; p < well_part.nb_points(); p++ ) {
            l += length( prev - well_part.point( p ) ) ;
            prev = well_part.point( p ) ;
        }
        l += length( prev - corner( well_part.corner(1) ).point() ) ;
        return l ;
    }

    void Well::get_part_edges( index_t p, std::vector< Edge >& edges ) const
    {
        const WellPart& well_part = part( p ) ;
        vec3 prev = corner( well_part.corner(0) ).point() ;
        for( index_t p = 0; p < well_part.nb_points(); p++ ) {
            edges.push_back( Edge( prev, well_part.point( p ) ) ) ;
            prev = well_part.point( p ) ;
        }
        edges.push_back( Edge( prev, corner( well_part.corner(1) ).point() ) ) ;
    }

    void Well::get_region_edges(
        index_t region,
        std::vector< Edge >& edges ) const
    {
        for( index_t p = 0; p < nb_parts(); p++ ) {
            if( part_region_id( p ) == region ) {
                const WellPart& well_part = part( p ) ;
                vec3 prev = corner( well_part.corner( 0 ) ).point() ;
                for( index_t p = 0; p < well_part.nb_points(); p++ ) {
                    edges.push_back( Edge( prev, well_part.point( p ) ) ) ;
                    prev = well_part.point( p ) ;
                }
                edges.push_back(
                    Edge( prev, corner( well_part.corner( 1 ) ).point() ) ) ;
            }
        }
    }

    void Well::all_part_vertices(
        index_t p,
        std::vector< vec3 >& vertices ) const
    {
        const WellPart& well_part = part( p ) ;
        vertices.reserve( well_part.nb_points()+2 ) ;
        vertices.push_back( corner( well_part.corner(0) ).point() ) ;
        for( index_t v = 0; v < well_part.nb_points(); v++ ) {
            vertices.push_back( well_part.point( v ) ) ;
        }
        vertices.push_back( corner( well_part.corner(1) ).point() ) ;
    }

    WellGroup::WellGroup()
        : model_( nil )
    {
    }
    WellGroup::~WellGroup()
    {
    }

    void WellGroup::get_region_edges(
        index_t region,
        std::vector< Edge >& edges ) const
    {
        for( index_t w = 0; w < nb_wells(); w++ ) {
            const Well& cur_well = well( w ) ;
            for( index_t p = 0; p < cur_well.nb_parts(); p++ ) {
                if( cur_well.part_region_id( p ) == region ) {
                    cur_well.get_part_edges( p, edges ) ;
                }
            }
        }
    }

    bool WellGroup::is_well_already_added( const std::string& name ) const
    {
        for( index_t w = 0; w < nb_wells(); w++ ) {
            if( well( w ).name() == name ) return true ;
        }
        return false ;
    }

}

