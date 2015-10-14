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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/well.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geometry.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/utils.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>

#include <cmath>
#include <stack>

namespace {
    using namespace RINGMesh ;

    /*!
    * @brief Returns the index of the region of the model neighboring the surface.
    * @param[in] model to consider
    * @param[in] surface_part_id Index of the Surface
    * @param[in] side Side of the Surface
    * @return The region index or NO_ID if none found.
    */
    index_t find_region( const GeoModel& BM, index_t surface_part_id, bool side )
    {
        ringmesh_debug_assert( surface_part_id < BM.nb_surfaces() ) ;
        GME::gme_t cur_surface( GME::SURFACE, surface_part_id ) ;
        /// @todo It would be better to directly check the region
        /// adjacent to the Surface.
        for( index_t r = 0; r < BM.nb_regions(); r++ ) {
            const Region& cur_region = BM.region( r ) ;
            for( index_t s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_id( s ) == cur_surface ) {
                    return r ;
                }
            }
        }
        return GME::NO_ID ;
    }
}

namespace RINGMesh {

    WellMesh::WellMesh( const Well* well )
        : well_( well ), mesh_( 3 )
    {
    }

    /*!
     * Gets the associated Mesh
     */
    GEO::Mesh& WellMesh::mesh() const
    {
        return const_cast< GEO::Mesh& >( mesh_ ) ;
    }

    /*!
     * Gets a point from the mesh
     * @param[in] p the point id
     * @return the corresponding point
     */
    const vec3& WellMesh::point( index_t p ) const
    {
        return mesh_.vertices.point( p ) ;
    }

    /*!
     * Gets the number of points
     */
    index_t WellMesh::nb_points() const
    {
        return mesh_.vertices.nb() ;
    }

// --------------------------------------------------------------------------

    WellCorner::WellCorner(
        const Well* well,
        const vec3& point,
        const corner_info_t& corner_info )
        : WellMesh( well )
    {
        mesh_.vertices.create_vertex( point.data() ) ;
        corner_info_.bind( mesh_.vertices.attributes(), "corner_info" ) ;
        corner_info_[0] = corner_info ;
    }

    WellCorner::~WellCorner()
    {
        if( corner_info_.is_bound() ) corner_info_.unbind() ;
    }

    /*!
     * Gets the corner info
     */
    const WellCorner::corner_info_t& WellCorner::corner_info() const {
        return corner_info_[0] ;
    }

// --------------------------------------------------------------------------

    /*!
     * Create a WellPart
     * @param[in] well the associated well
     * @param[in] id the position in the parts_ vector of the associated well
     */
    WellPart::WellPart( const Well* well, index_t id )
          : WellMesh( well ), id_( id )
    {
        corners_[0] = NO_ID ;
        corners_[1] = NO_ID ;
    }

    /*!
     * Create the associated Mesh of the part
     * @param[in] points the points of the mesh
     * @pre the points should be oriented in the order of the well path
     */
    void WellPart::set_points( const std::vector< vec3 >& points )
    {
        mesh_.vertices.create_vertices( points.size() ) ;
        for( index_t p = 0; p < points.size(); p++ ) {
            mesh_.vertices.point( p ) = points[p] ;
        }

        mesh_.edges.create_edges( points.size() - 1 ) ;
        for( index_t e = 0; e < points.size() - 1; e++ ) {
            mesh_.edges.set_vertex( e, 0, e ) ;
            mesh_.edges.set_vertex( e, 1, e + 1 ) ;
        }
    }

    /*!
     * Gets the number of edges
     */
    index_t WellPart::nb_edges() const
    {
        return mesh_.edges.nb() ;
    }

    /*!
     * Gets the length of the part
     */
    double WellPart::length() const
    {
        double l = 0.0 ;
        for( index_t e = 0; e < nb_edges(); e++ ) {
            l += ( point( e + 1 ) - point( e ) ).length() ;
        }
        return l ;
    }

// --------------------------------------------------------------------------

    Well::Well()
        : nb_edges_( NO_ID )
    {
    }

    Well::~Well()
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( corners_[c] ) delete corners_[c] ;
        }
        for( index_t part = 0; part < nb_parts(); part++ ) {
            if( parts_[part] ) delete parts_[part] ;
        }
    }

    /*!
     * Finds if a corner at a given geometric position exist
     * @param[in] p the geometric position to test
     * @return the id of the corner or NO_ID if not found any corresponding to \p p
     */
    index_t Well::find_corner( const vec3& p ) const
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( p == corner( c ).point() ) {
                return c ;
            }
        }
        return NO_ID ;
    }

    /*!
     * Copies information and resize the number of parts and corners
     * @param[in,out] well the current well information will be copied into this one
     */
    void Well::copy_corners_and_informations( Well& well ) const
    {
        well.name_ = name_ ;
        well.part_region_id_ = part_region_id_ ;

        well.corners_.reserve( nb_corners() ) ;
        for( index_t c = 0; c < nb_corners(); c++ ) {
            well.create_corner( corners_[c]->point(), corners_[c]->corner_info() ) ;
        }

        well.parts_.reserve( nb_parts() ) ;
        for( index_t p = 0; p < nb_parts(); p++ ) {
            well.create_part( part_region_id( p ) ) ;
            const WellPart& from_part = part( p ) ;
            WellPart& cur_part = well.part( p ) ;
//            cur_part.mesh().copy( from_part.mesh() ) ;
            cur_part.set_corner( 0, from_part.corner( 0 ) ) ;
            cur_part.set_corner( 1, from_part.corner( 1 ) ) ;
        }
    }

    /*!
     * Gets the number of edges of the well
     */
    index_t Well::nb_edges() const {
        if( nb_edges_ == NO_ID ) {
            index_t res = 0 ;
            for( index_t p = 0; p < nb_parts(); p++ ) {
                res += part( p ).mesh().edges.nb() ;
            }
            const_cast< Well* >( this )->nb_edges_ = res ;
        }
        return nb_edges_ ;
    }

    /*!
     * Gets the edges of a part
     * @param[in] p the part id
     * @param[out] edges the edges of the part
     */
    void Well::get_part_edges(
        index_t p,
        std::vector< Edge >& edges ) const
    {
        const WellPart& well_part = part( p ) ;
        for( index_t e = 0; e < well_part.nb_edges(); e++ ) {
            edges.push_back(
                Edge( well_part.point( e ), well_part.point( e + 1 ) ) ) ;
        }
    }

    /*!
     * Gets all the edges of a corresponding region
     * @param[in] region the region id
     * @param[out] edges the corresponding edges
     */
    void Well::get_region_edges(
        index_t region,
        std::vector< Edge >& edges ) const
    {
        for( index_t p = 0; p < nb_parts(); p++ ) {
            if( part_region_id( p ) == region ) {
                get_part_edges( p, edges ) ;
            }
        }
    }


// --------------------------------------------------------------------------


    WellGroup::WellGroup()
          : model_( nil )
    {
    }


    WellGroup::~WellGroup()
    {
    }

    /*!
     * Gets all the edges contained in a region
     * @param[in] region the region id
     * @param[out] edges the edges of the region
     */
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

    /*!
     * Gets all the edges contained in a region
     * @param[in] region the region id
     * @param[out] edges the edges of the region, one vector per well
     */
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
            if( segment_triangle_intersection(
                v_from_, v_to_,
                surface_.vertex( trgl, 0 ), surface_.vertex( trgl, 1 ),
                surface_.vertex( trgl, 2 ), result ) ) {
                intersections_.push_back(
                    LineInstersection( result, surface_.gme_id().index, trgl ) ) ;
            }
        }

    private:
        const Surface& surface_ ;
        const vec3& v_from_ ;
        const vec3& v_to_ ;

        std::vector< LineInstersection >& intersections_ ;
    } ;

    /*!
     * Creates new wells
     * @param[in] nb the number of wells to create
     */
    void WellGroup::create_wells( index_t nb )
    {
        wells_.resize( nb, nil ) ;
        for( index_t w = 0; w < nb_wells(); w++ ) {
            wells_[w] = new Well ;
        }
    }

    /*!
     * Add a well from its mesh and makes it conformal to the associated BoundarModel
     * @param[in] mesh the mesh of the well
     * @param[in] name the name of the well
     * @pre the mesh should be made of continuous edges without fork or discontinuities
     */
    void WellGroup::add_well( const GEO::Mesh& mesh, const std::string& name )
    {
        ringmesh_debug_assert( model() ) ;
        if( find_well( name ) != NO_ID ) return ;
        wells_.push_back( new Well ) ;
        Well& new_well = *wells_.back() ;
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
        well_points.push_back(  mesh.vertices.point( 0 ) ) ;
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
                    const vec3& v_prev = well_points.back() ;
                    index_t index = indices[i] ;
                    vec3 direction = v_prev - intersections[index].intersection_ ;
                    bool sign =
                        dot( direction,
                            model_->surface( intersections[index].surface_id_ ).normal(
                                intersections[index].trgl_id_ ) ) > 0 ;
                    last_sign = sign ;
                    index_t region = find_region(
                        *model_, intersections[index].surface_id_, sign ) ;
                    if( region != NO_ID ) {
                        index_t new_well_part_id = new_well.create_part( region ) ;
                        WellPart& well_part = new_well.part( new_well_part_id ) ;
                        index_t corner0 = new_well.find_corner(
                            start.intersection_ ) ;
                        if( corner0 == NO_ID ) {
                            WellCorner::corner_info_t corner_info ;
                            corner_info.is_on_surface = start.surface_id_ != NO_ID ;
                            if( corner_info.is_on_surface ) {
                                corner_info.id = start.surface_id_ ;
                            } else {
                                corner_info.id = region ;
                            }
                            corner0 = new_well.create_corner( start.intersection_,
                                corner_info ) ;
                        }
                        well_part.set_corner( 0, corner0 ) ;
                        index_t corner1 = new_well.find_corner(
                            intersections[index].intersection_ ) ;
                        if( corner1 == NO_ID ) {
                            WellCorner::corner_info_t corner_info ;
                            corner_info.is_on_surface =
                                intersections[index].surface_id_ != NO_ID ;
                            if( corner_info.is_on_surface ) {
                                corner_info.id = intersections[index].surface_id_ ;
                            } else {
                                corner_info.id = region ;
                            }
                            corner1 = new_well.create_corner( intersections[index].intersection_,
                                corner_info ) ;
                        }
                        well_part.set_corner( 1, corner1 ) ;
                        well_points.push_back( intersections[index].intersection_ ) ;
                        well_part.set_points( well_points ) ;
                    }
                    well_points.clear() ;
                    start = intersections[index] ;
                    well_points.push_back( start.intersection_ ) ;

                }
            }
            well_points.push_back( v_to ) ;
        }
        index_t region = NO_ID ;
        if( start.surface_id_ != NO_ID )
            region = find_region( *model_, start.surface_id_, !last_sign ) ;
        else {
            Box3d box ;
            for( index_t s = 0; s < model()->nb_surfaces(); s++ ) {
                const Surface& surface = model()->surface( s ) ;
                for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                    box.add_point( surface.vertex( v ) ) ;
                }
            }
            vec3 diag_box = box.max() - box.min() ;
            vec3 v_from =  well_points.back() ;
            vec3 v_to = well_points.back() + diag_box ;

            Box3d edge_box ;
            edge_box.add_point( v_from ) ;
            edge_box.add_point( v_to ) ;
            std::vector< index_t > potential_surfaces ;
            for( index_t s = 0; s < model()->nb_surfaces(); s++ ) {
                if( edge_box.bboxes_overlap( boxes[s] ) ) {
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
                index_t index = indices[0] ;
                vec3 direction = v_from - intersections[index].intersection_ ;
                bool sign = dot( direction,
                    model_->surface( intersections[index].surface_id_ ).normal(
                        intersections[index].trgl_id_ ) ) > 0 ;
                last_sign = sign ;
                region = find_region( *model_, intersections[index].surface_id_,
                    sign ) ;
            }
        }
        if( region != NO_ID ) {
//            WellPart well_part ;
//            well_part.set_well( &new_well ) ;
//            index_t id = start.surface_id_ == -1 ? -region - 1 : start.surface_id_ ;
//            signed_index_t corner_id = new_well.find_or_create_corner( start.intersection_,
//                id ) ;
//            well_part.add_corner( corner_id ) ;
//            corner_id = new_well.find_or_create_corner( well_points.back() ) ;
//            well_part.add_corner( corner_id ) ;
//            well_points.resize( well_points.size() - 1 ) ;
//            well_part.set_points( well_points ) ;
//            new_well.add_part( well_part, region ) ;


            index_t new_well_part_id = new_well.create_part( region ) ;
            WellPart& well_part = new_well.part( new_well_part_id ) ;
            index_t corner0 = new_well.find_corner(
                start.intersection_ ) ;
            if( corner0 == NO_ID ) {
                WellCorner::corner_info_t corner_info ;
                corner_info.is_on_surface = start.surface_id_ != NO_ID ;
                if( corner_info.is_on_surface ) {
                    corner_info.id = start.surface_id_ ;
                } else {
                    corner_info.id = region ;
                }
                corner0 = new_well.create_corner( start.intersection_,
                    corner_info ) ;
            }
            well_part.set_corner( 0, corner0 ) ;
//            index_t corner1 = new_well.find_corner(
//                well_points.back() ) ;
            ringmesh_debug_assert( new_well.find_corner(
                well_points.back() ) == NO_ID ) ;
//            if( corner1 == NO_ID ) {
                WellCorner::corner_info_t corner_info ;
                corner_info.is_on_surface = false ;
//                    intersections[index].surface_id_ != NO_ID ;
//                if( corner_info.is_on_surface ) {
//                    corner_info.id = intersections[index].surface_id_ ;
//                } else {
                    corner_info.id = region ;
//                }
                   index_t corner1 = new_well.create_corner( well_points.back() ,
                    corner_info ) ;
//            }
            well_part.set_corner( 1, corner1 ) ;
//            well_points.push_back( intersections[index].intersection_ ) ;
            well_part.set_points( well_points ) ;

        }

    }

    /*!
     * Finds if a well with the same name already exist
     * @param[in] name the name to test
     * @return the id of the well or NO_ID
     */
    index_t WellGroup::find_well( const std::string& name ) const
    {
        for( index_t w = 0; w < nb_wells(); w++ ) {
            if( well( w ).name() == name ) {
                return w ;
            }
        }
        return NO_ID ;
    }
}
