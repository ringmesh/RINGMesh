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

#include <ringmesh/mesh/well.h>

#include <cmath>
#include <stack>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/box3d.h>
#include <ringmesh/basic/geometry.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

#include <ringmesh/mesh/geogram_mesh.h>

/*!
 * @file Implementation of Wells
 * @author Arnaud Botella
 */

namespace {
    using namespace RINGMesh ;

    /*!
     * @brief Returns the index of the region of the geomodel neighboring the surface.
     * @param[in] geomodel to consider
     * @param[in] surface_part_id Index of the Surface
     * @param[in] side Side of the Surface
     * @return The region index or NO_ID if none found.
     */
    index_t find_region(
        const GeoModel& geomodel,
        index_t surface_part_id,
        bool side )
    {
        ringmesh_assert( surface_part_id < geomodel.nb_surfaces() ) ;
        gme_t cur_surface( Surface::type_name_static(), surface_part_id ) ;
        /// @todo It would be better to directly check the region
        /// adjacent to the Surface.
        for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
            const Region& cur_region = geomodel.region( r ) ;
            for( index_t s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_gme( s ) == cur_surface ) {
                    return r ;
                }
            }
        }
        return NO_ID ;
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

    index_t find_or_create_corner(
        Well& well,
        index_t region,
        const LineInstersection& corner )
    {
        index_t corner_id = well.find_corner( corner.intersection_ ) ;
        if( corner_id == NO_ID ) {
            bool is_on_surface = corner.surface_id_ != NO_ID ;
            index_t id = NO_ID ;
            if( is_on_surface ) {
                id = corner.surface_id_ ;
            } else {
                id = region ;
            }
            corner_id = well.create_corner( corner.intersection_, is_on_surface,
                id ) ;
        }
        return corner_id ;
    }

    void create_well_part_and_corners(
        Well& well,
        index_t region,
        const std::vector< vec3 > vertices,
        const LineInstersection& start,
        const LineInstersection& end )
    {
        index_t new_well_part_id = well.create_part( region ) ;
        WellPart& well_part = well.part( new_well_part_id ) ;

        index_t corner0 = find_or_create_corner( well, region, start ) ;
        well_part.set_corner( 0, corner0 ) ;

        index_t corner1 = find_or_create_corner( well, region, end ) ;
        well_part.set_corner( 1, corner1 ) ;
        well_part.set_points( vertices ) ;
    }

    bool get_side(
        const vec3& vertex,
        const vec3& on_surface,
        const Surface& surface,
        index_t triangle )
    {
        vec3 direction = vertex - on_surface ;
        return dot( direction, surface.facet_normal( triangle ) ) > 0 ;
    }
}

namespace RINGMesh {

    WellEntity::WellEntity( const Well* well )
        : well_( well )
    {
    }

// --------------------------------------------------------------------------

    WellCorner::WellCorner(
        const Well* well,
        const vec3& point,
        bool is_on_surface,
        index_t id )
        :
            WellEntity( well ),
            is_on_surface_( is_on_surface ),
            id_( id ),
            mesh_( Mesh0D::create_mesh( GeogramMesh0D::type_name_static() ) )
    {
        Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder( *mesh_ ) ;
        builder->create_vertex( point ) ;
    }

    WellCorner::~WellCorner()
    {
        delete mesh_ ;
    }

    const vec3& WellCorner::point() const
    {
        return mesh_->vertex( 0 ) ;
    }

// --------------------------------------------------------------------------

    /*!
     * Create a WellPart
     * @param[in] well the associated well
     * @param[in] id the position in the parts_ vector of the associated well
     */
    WellPart::WellPart( const Well* well, index_t id )
        :
            WellEntity( well ),
            id_( id ),
            mesh_( Mesh1D::create_mesh( GeogramMesh1D::type_name_static() ) )
    {
        corners_[0] = NO_ID ;
        corners_[1] = NO_ID ;
    }

    WellPart::~WellPart()
    {
        delete mesh_ ;
    }

    /*!
     * Create the associated Mesh of the part
     * @param[in] points the points of the mesh
     * @pre the points should be oriented in the order of the well path
     */
    void WellPart::set_points( const std::vector< vec3 >& points )
    {
        index_t nb_points = static_cast< index_t >( points.size() ) ;
        DEBUG( nb_points ) ;
        Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder( *mesh_ ) ;
        builder->create_vertices( nb_points ) ;
        for( index_t p = 0; p < nb_points; p++ ) {
            builder->set_vertex( p, points[p] ) ;
        }

        index_t nb_edges = nb_points - 1 ;
        builder->create_edges( nb_edges ) ;
        for( index_t e = 0; e < nb_edges; e++ ) {
            builder->set_edge_vertex( e, 0, e ) ;
            builder->set_edge_vertex( e, 1, e + 1 ) ;
        }
    }

    /*!
     * Gets the number of edges
     */
    index_t WellPart::nb_edges() const
    {
        return mesh_->nb_edges() ;
    }

    const vec3& WellPart::vertex( index_t v ) const
    {
        return mesh_->vertex( v ) ;
    }

    const vec3& WellPart::edge_vertex( index_t edge, index_t v ) const
    {
        return vertex( mesh_->edge_vertex( edge, v ) ) ;
    }
    /*!
     * Gets the length of the part
     */
    double WellPart::length() const
    {
        double l = 0.0 ;
        for( index_t e = 0; e < nb_edges(); e++ ) {
            l += ( vertex( e + 1 ) - vertex( e ) ).length() ;
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
     * @param[in] vertex the geometric position to test
     * @return the id of the corner or NO_ID if not found any corresponding to \p p
     */
    index_t Well::find_corner( const vec3& vertex ) const
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( vertex == corner( c ).point() ) {
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
            well.create_corner( corners_[c]->point(), corners_[c]->is_on_surface(),
                corners_[c]->id() ) ;
        }

        well.parts_.reserve( nb_parts() ) ;
        for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
            well.create_part( part_region_id( part_id ) ) ;
            const WellPart& from_part = part( part_id ) ;
            WellPart& cur_part = well.part( part_id ) ;
//            cur_part.mesh().copy( from_part.mesh() ) ;
            cur_part.set_corner( 0, from_part.corner( 0 ) ) ;
            cur_part.set_corner( 1, from_part.corner( 1 ) ) ;
        }
    }

    /*!
     * Gets the number of edges of the well
     */
    index_t Well::nb_edges() const
    {
        if( nb_edges_ == NO_ID ) {
            index_t nb_edges = 0 ;
            for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
                nb_edges += part( part_id ).nb_edges() ;
            }
            const_cast< Well* >( this )->nb_edges_ = nb_edges ;
        }
        return nb_edges_ ;
    }

    /*!
     * Gets the edges of a part
     * @param[in] part_id the part id
     * @param[out] edges the edges of the part
     */
    void Well::get_part_edges( index_t part_id, std::vector< Edge >& edges ) const
    {
        const WellPart& well_part = part( part_id ) ;
        for( index_t e = 0; e < well_part.nb_edges(); e++ ) {
            edges.push_back(
                Edge( well_part.vertex( e ), well_part.vertex( e + 1 ) ) ) ;
        }
    }

    /*!
     * Gets all the edges of a corresponding region
     * @param[in] region the region id
     * @param[out] edges the corresponding edges
     */
    void Well::get_region_edges( index_t region, std::vector< Edge >& edges ) const
    {
        for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
            if( part_region_id( part_id ) == region ) {
                get_part_edges( part_id, edges ) ;
            }
        }
    }

// --------------------------------------------------------------------------

    WellGroup::WellGroup()
        : geomodel_( nil )
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
            for( index_t part_id = 0; part_id < cur_well.nb_parts(); part_id++ ) {
                if( cur_well.part_region_id( part_id ) == region ) {
                    cur_well.get_part_edges( part_id, edges ) ;
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
            if( segment_triangle_intersection( v_from_, v_to_,
                surface_.mesh_element_vertex( trgl, 0 ),
                surface_.mesh_element_vertex( trgl, 1 ),
                surface_.mesh_element_vertex( trgl, 2 ), result ) ) {
                intersections_.push_back(
                    LineInstersection( result, surface_.index(), trgl ) ) ;
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

    struct OrientedEdge {
        OrientedEdge( const Mesh1D& mesh, index_t edge, index_t vertex_from )
            : edge_( edge ), vertex_from_( vertex_from )
        {
            if( mesh.edge_vertex( edge, 0 ) == vertex_from ) {
                edge_vertex_ = 0 ;
            } else {
                ringmesh_assert( mesh.edge_vertex( edge, 1 ) == vertex_from ) ;
                edge_vertex_ = 1 ;
            }
        }
        index_t edge_ ;
        index_t edge_vertex_ ;
        index_t vertex_from_ ;
    } ;

    struct OrientedEdgePart: public OrientedEdge {
        OrientedEdgePart(
            const Mesh1D& mesh,
            index_t edge,
            index_t vertex_from,
            const LineInstersection& intersection )
            : OrientedEdge( mesh, edge, vertex_from ), intersection_( intersection )
        {
        }
        OrientedEdgePart(
            const OrientedEdge& edge,
            const LineInstersection& intersection )
            : OrientedEdge( edge ), intersection_( intersection )
        {
        }
        LineInstersection intersection_ ;
    } ;

    /*!
     * Add a well from its mesh and makes it conformal to the associated GeoModel
     * @param[in] mesh the mesh of the well
     * @param[in] name the name of the well
     */
    void WellGroup::add_well( const Mesh1D& mesh, const std::string& name )
    {
        ringmesh_assert( geomodel() ) ;
        double epsilon = geomodel_->epsilon() ;
        if( find_well( name ) != NO_ID ) return ;
        wells_.push_back( new Well ) ;
        Well& new_well = *wells_.back() ;
        new_well.set_name( name ) ;

        std::vector< std::vector< index_t > > edges_around_vertices(
            mesh.nb_vertices() ) ;
        for( index_t e = 0; e < mesh.nb_edges(); e++ ) {
            for( index_t i = 0; i < 2; i++ ) {
                index_t v = mesh.edge_vertex( e, i ) ;
                edges_around_vertices[v].push_back( e ) ;
            }
        }

        std::stack< OrientedEdge > S ;
        for( index_t v = 0; v < mesh.nb_vertices(); v++ ) {
            const std::vector< index_t >& edges = edges_around_vertices[v] ;
            if( edges.size() == 1 ) {
                S.push( OrientedEdge( mesh, edges.front(), v ) ) ;
            }
        }
        DEBUG( S.size() ) ;
        DEBUG( mesh.nb_edges() ) ;
        DEBUG( mesh.nb_vertices() ) ;

        std::vector< bool > edge_visited( mesh.nb_edges(), false ) ;
        do {
            OrientedEdge cur_edge = S.top() ;
            S.pop() ;
            if( edge_visited[cur_edge.edge_] ) {
                continue ;
            }
            edge_visited[cur_edge.edge_] = true ;
            LineInstersection start(
                mesh.vertex(
                    mesh.edge_vertex( cur_edge.edge_, cur_edge.edge_vertex_ ) ) ) ;

            DEBUG( "START PART" ) ;
            std::vector< vec3 > well_part_points ;
            std::stack< OrientedEdgePart > S_part ;
            S_part.push( OrientedEdgePart( cur_edge, start ) ) ;
            do {
                OrientedEdgePart cur_edge_part = S_part.top() ;
                S_part.pop() ;
                edge_visited[cur_edge_part.edge_] = true ;
                const vec3& v_from = cur_edge_part.intersection_.intersection_ ;
                index_t v_to_id = mesh.edge_vertex( cur_edge_part.edge_,
                    ( cur_edge_part.edge_vertex_ + 1 ) % 2 ) ;
                const vec3& v_to = mesh.vertex( v_to_id ) ;
                well_part_points.push_back( v_from ) ;

                DEBUG( v_from ) ;
                DEBUG( v_to ) ;

                Box3d box ;
                box.add_point( v_from ) ;
                box.add_point( v_to ) ;
                std::vector< LineInstersection > intersections ;
                for( index_t s = 0; s < geomodel()->nb_surfaces(); s++ ) {
                    const Surface& surface = geomodel()->surface( s ) ;
                    EdgeConformerAction action( surface, v_from, v_to,
                        intersections ) ;
                    surface.facets_aabb().compute_bbox_element_bbox_intersections(
                        box, action ) ;
                }

                DEBUG( intersections.size() ) ;
                if( !intersections.empty() ) {
                    std::vector< index_t > indices( intersections.size() ) ;
                    std::vector< double > distances( intersections.size() ) ;
                    for( index_t i = 0; i < intersections.size(); i++ ) {
                        indices[i] = i ;
                        distances[i] = length(
                            v_from - intersections[i].intersection_ ) ;
                    }
                    indirect_sort( distances, indices ) ;
                    index_t i = 0 ;
                    if( distances[indices[0]] < epsilon ) {
                        start = intersections[0] ;
                        i++ ;
                    }
                    for( ; i < intersections.size(); i++ ) {
                        index_t index = indices[i] ;
                        bool sign = get_side( v_from,
                            intersections[index].intersection_,
                            geomodel()->surface( intersections[index].surface_id_ ),
                            intersections[index].trgl_id_ ) ;
                        index_t region = find_region( *geomodel_,
                            intersections[index].surface_id_, sign ) ;
                        if( region != NO_ID ) {
                            well_part_points.push_back(
                                intersections[index].intersection_ ) ;
                            create_well_part_and_corners( new_well, region,
                                well_part_points, start, intersections[index] ) ;
                        }
                        well_part_points.clear() ;
                        start = intersections[index] ;
                        well_part_points.push_back( start.intersection_ ) ;
                    }
//                    S_part.push( OrientedEdgePart( cur_edge_part, start ) ) ;
                }
//                else {
                const std::vector< index_t >& edges = edges_around_vertices[v_to_id] ;
                DEBUG( edges.size() ) ;
                if( edges.size() == 2 ) {
                    index_t count = 0 ;
                    for( index_t e = 0; e < edges.size(); e++ ) {
                        if( !edge_visited[edges[e]] ) {
                            S_part.push(
                                OrientedEdgePart( mesh, edges[e], v_to_id,
                                    LineInstersection( v_to ) ) ) ;
                            count++ ;
                        }
                    }
                    ringmesh_assert( count == 1 ) ;
                } else {
                    double best_distance = max_float64() ;
                    index_t best_surface = NO_ID ;
                    vec3 best_nearest ;
                    index_t best_triangle = NO_ID ;
                    for( index_t s = 0; s < geomodel()->nb_surfaces(); s++ ) {
                        const Surface& surface = geomodel()->surface( s ) ;
                        vec3 nearest ;
                        double distance ;
                        index_t triangle = surface.facets_aabb().closest_triangle(
                            v_to, nearest, distance ) ;
                        if( distance < best_distance ) {
                            best_distance = distance ;
                            best_nearest = nearest ;
                            best_surface = s ;
                            best_triangle = triangle ;
                        }
                    }
                    bool sign = get_side( v_to, best_nearest,
                        geomodel()->surface( best_surface ), best_triangle ) ;
                    index_t region = find_region( *geomodel_, best_surface, sign ) ;
                    well_part_points.push_back( v_to ) ;
                    create_well_part_and_corners( new_well, region, well_part_points,
                        start, LineInstersection( v_to ) ) ;
                    for( index_t e = 0; e < edges.size(); e++ ) {
                        S.push( OrientedEdge( mesh, edges[e], v_to_id ) ) ;
                    }
//                    }
                }
                DEBUG( S_part.size() ) ;
            } while( !S_part.empty() ) ;
            DEBUG( S.size() ) ;
        } while( !S.empty() ) ;

        DEBUG( nb_wells() ) ;
        DEBUG( well( 0 ).nb_corners() ) ;
        DEBUG( well( 0 ).nb_edges() ) ;
        DEBUG( well( 0 ).nb_parts() ) ;
        std::vector< vec3 > points ;
        index_t count = 0 ;
        for( index_t part_id = 0; part_id < well( 0 ).nb_parts(); part_id++ ) {
            for( index_t e = 0; e < well( 0 ).part( part_id ).nb_edges(); e++){
                const vec3& v0 = well( 0 ).part( part_id ).edge_vertex( e, 0 ) ;
                const vec3& v1 = well( 0 ).part( part_id ).edge_vertex( e, 1 ) ;
                points.push_back( 0.5*(v0+v1) ) ;
                if( length( v0 -v1 ) < epsilon ) count++ ;
            }
        }
        NNSearch nn( points, false ) ;
        std::vector< index_t > tt;
        DEBUG( nn.get_colocated_index_mapping( epsilon, tt ) ) ;
        DEBUG( count ) ;
        DEBUG( epsilon ) ;
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
