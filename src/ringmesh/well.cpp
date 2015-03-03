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

#include <ringmesh/well.h>
#include <ringmesh/utils.h>

#include <cmath>
#include <stack>

namespace RINGMesh {
    signed_index_t Well::find_or_create_corner(
        const vec3& p,
        signed_index_t id )
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( p == corner( c ).point() ) {
                ringmesh_debug_assert( id == corner( c ).surface_id() ) ;
                ringmesh_debug_assert( id != - 1 ) ;
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

        return - 1 ;
    }


    void Well::copy_corners_and_informations( Well& well ) const
    {
        well.name_ = name_ ;
        well.corners_ = corners_ ;
        well.part_region_id_ = part_region_id_ ;
        well.parts_.resize( nb_parts() ) ;

        for( index_t p = 0; p < nb_parts(); p++ ) {
            well.parts_[ p ].add_corner( parts_[ p ].corner( 0 ) ) ;
            well.parts_[ p ].add_corner( parts_[ p ].corner( 1 ) ) ;
            well.parts_[ p ].set_id( parts_[ p ].id() ) ;
            well.parts_[ p ].set_well( &well ) ;
        }
    }


    double Well::part_length( index_t p ) const
    {
        double l = 0.0 ;
        const WellPart& well_part = part( p ) ;
        vec3 prev = corner( well_part.corner( 0 ) ).point() ;

        for( index_t p = 0; p < well_part.nb_points(); p++ ) {
            l += length( prev - well_part.point( p ) ) ;
            prev = well_part.point( p ) ;
        }

        l += length( prev - corner( well_part.corner( 1 ) ).point() ) ;
        return l ;
    }


    void Well::get_part_edges(
        index_t p,
        std::vector< Edge >& edges ) const
    {
        const WellPart& well_part = part( p ) ;
        vec3 prev = corner( well_part.corner( 0 ) ).point() ;

        for( index_t p = 0; p < well_part.nb_points(); p++ ) {
            edges.push_back( Edge( prev, well_part.point( p ) ) ) ;
            prev = well_part.point( p ) ;
        }

        edges.push_back( Edge( prev, corner( well_part.corner( 1 ) ).point() ) ) ;
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
        vertices.reserve( well_part.nb_points() + 2 ) ;
        vertices.push_back( corner( well_part.corner( 0 ) ).point() ) ;

        for( index_t v = 0; v < well_part.nb_points(); v++ ) {
            vertices.push_back( well_part.point( v ) ) ;
        }

        vertices.push_back( corner( well_part.corner( 1 ) ).point() ) ;
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
            if( well( w ).name() == name ) {return true ;}
        }

        return false ;
    }
}
