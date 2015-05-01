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
 *     Ecole Nationale Sup�rieure de G�ologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/well.h>
#include <ringmesh/utils.h>
#include <ringmesh/boundary_model.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

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
        return -1 ;
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

    void WellGroup::get_region_edges(
        index_t region,
        std::vector< std::vector< Edge > >& edges ) const
    {
        edges.clear() ;
        edges.resize( nb_wells() ) ;
        for( index_t w = 0; w < nb_wells(); w++ ) {
            const Well& cur_well = well( w ) ;
            cur_well.get_region_edges( region, edges[w] ) ;
        }
    }

    struct LineInstersection {
        LineInstersection(
            const vec3& intersection,
            index_t surface_id = NO_ID,
            index_t trgl_id = NO_ID )
            :
                intersection_( intersection ),
                surface_id_( surface_id ),
                trgl_id_( trgl_id )
        {
        }
        vec3 intersection_ ;
        index_t surface_id_ ;
        index_t trgl_id_ ;
    } ;

    class EdgeConformerAction {
    public:
        EdgeConformerAction(
            const Surface& surface,
            const vec3& v_from,
            const vec3& v_to,
            std::vector< LineInstersection >& intersections )
            :
                surface_( surface ),
                v_from_( v_from ),
                v_to_( v_to ),
                intersections_( intersections )
        {
        }

        void operator()( index_t trgl )
        {
            vec3 result ;
            if( Math::segment_triangle_intersection(
                v_from_, v_to_,
                surface_.vertex( trgl, 0 ), surface_.vertex( trgl, 1 ),
                surface_.vertex( trgl, 2 ), result ) ) {
                intersections_.push_back(
                    LineInstersection( result, surface_.id(), trgl ) ) ;
            }
        }

    private:
        const Surface& surface_ ;
        const vec3& v_from_ ;
        const vec3& v_to_ ;

        std::vector< LineInstersection >& intersections_ ;
    } ;

    void WellGroup::add_well( const GEO::Mesh& mesh, const std::string& name )
    {
        ringmesh_debug_assert( model() ) ;
        if( is_well_already_added( name ) ) return ;
        wells_.push_back( Well() ) ;
        Well& new_well = wells_.back() ;
        new_well.set_name( name ) ;

        std::vector< Box3d > boxes( model()->nb_surfaces() ) ;
        for( index_t s = 0; s < model()->nb_surfaces(); s++ ) {
            Box3d& box = boxes[s] ;
            const Surface& surface = model()->surface( s ) ;
            for( unsigned int p = 0; p < surface.nb_vertices(); p++ ) {
                box.add_point( surface.vertex( p ) ) ;
            }
        }

        bool last_sign = false ;
        LineInstersection start( mesh.vertices.point( 0 ) ) ;
        std::vector< vec3 > well_points ;
        for( index_t e = 0; e < mesh.edges.nb(); e++ ) {
            const vec3& v_from = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( e, 0 ) ) ;
            const vec3& v_to = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( e, 1 ) ) ;
            Box3d box ;
            box.add_point( v_from ) ;
            box.add_point( v_to ) ;
            std::vector< index_t > potential_surfaces ;
            for( index_t s = 0; s < model()->nb_surfaces(); s++ ) {
                if( box.bboxes_overlap( boxes[s] ) ) {
                    potential_surfaces.push_back( s ) ;
                }
            }
            std::vector< LineInstersection > intersections ;
            for( index_t s = 0; s < potential_surfaces.size(); s++ ) {
                index_t s_id = potential_surfaces[s] ;
                const Surface& surface = model_->surface( s_id ) ;
                EdgeConformerAction action( surface, v_from, v_to, intersections ) ;
                surface.tools.aabb().compute_bbox_facet_bbox_intersections( box,
                    action ) ;
            }
            if( !intersections.empty() ) {
                std::vector< index_t > indices( intersections.size() ) ;
                std::vector< double > distances( intersections.size() ) ;
                if( intersections.size() == 1 ) {
                    indices[0] = 0 ;
                } else {
                    for( index_t i = 0; i < intersections.size(); i++ ) {
                        indices[i] = i ;
                        distances[i] = length(
                            v_from - intersections[i].intersection_ ) ;
                    }
                    IndirectSort< double, index_t > sort( distances, indices ) ;
                    sort.sort() ;
                }
                for( index_t i = 0; i < intersections.size(); i++ ) {
                    vec3 v_prev =
                        well_points.empty() ?
                            start.intersection_ : vec3( well_points.back().data() ) ;
                    index_t index = indices[i] ;
                    vec3 direction = v_prev - intersections[index].intersection_ ;
                    bool sign =
                        dot( direction,
                            model_->surface( intersections[index].surface_id_ ).facet_normal(
                                intersections[index].trgl_id_ ) ) > 0 ;
                    last_sign = sign ;
                    index_t region = model_->find_region(
                        intersections[index].surface_id_, sign ) ;
                    if( region == NO_ID ) {
                        well_points.clear() ;
                    } else {
                        WellPart well_part ;
                        index_t id =
                            start.surface_id_ == NO_ID ?
                                -region - 1 : start.surface_id_ ;
                        signed_index_t corner_id = new_well.find_or_create_corner(
                            start.intersection_, id ) ;
                        well_part.add_corner( corner_id ) ;
                        corner_id = new_well.find_or_create_corner(
                            intersections[index].intersection_,
                            intersections[index].surface_id_ ) ;
                        well_part.add_corner( corner_id ) ;
                        well_part.add_points( well_points ) ;
                        new_well.add_part( well_part, region ) ;
                        well_points.clear() ;
                    }
                    start = intersections[index] ;
                }
            }
            well_points.push_back( v_to ) ;
        }
        index_t region = model_->find_region( start.surface_id_, !last_sign ) ;
        if( region != NO_ID ) {
            WellPart well_part ;
            index_t id = start.surface_id_ == -1 ? -region - 1 : start.surface_id_ ;
            signed_index_t corner_id = new_well.find_or_create_corner( start.intersection_,
                id ) ;
            well_part.add_corner( corner_id ) ;
            corner_id = new_well.find_or_create_corner( well_points.back() ) ;
            well_part.add_corner( corner_id ) ;
            well_points.resize( well_points.size() - 1 ) ;
            well_part.add_points( well_points ) ;
            new_well.add_part( well_part, region ) ;
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
