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
    using namespace RINGMesh;

    bool inexact_equal( const vec3& v1, const vec3& v2, double epsilon )
    {
        return length( v2 - v1 ) < epsilon;
    }

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
        ringmesh_assert( surface_part_id < geomodel.nb_surfaces() );
        gmme_id cur_surface( Surface::type_name_static(), surface_part_id );
        const Surface& surface = geomodel.surface( surface_part_id );
        for( index_t r = 0; r < surface.nb_incident_entities(); r++ ) {
            const Region& cur_region =
                dynamic_cast< const Region& >( surface.incident_entity( r ) );
            for( index_t s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_gmme( s ) == cur_surface ) {
                    return r;
                }
            }
        }
        return NO_ID;
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
        LineInstersection()
            : intersection_( vec3() ), surface_id_( NO_ID ), trgl_id_( NO_ID )
        {
        }
        vec3 intersection_;
        index_t surface_id_;
        index_t trgl_id_;
    };

    index_t find_or_create_corner(
        Well& well,
        index_t region,
        const LineInstersection& corner,
        double epsilon )
    {
        index_t corner_id = well.find_corner( corner.intersection_, epsilon );
        if( corner_id == NO_ID ) {
            bool is_on_surface = corner.surface_id_ != NO_ID;
            index_t id = NO_ID;
            if( is_on_surface ) {
                id = corner.surface_id_;
            } else {
                id = region;
            }
            corner_id = well.create_corner( corner.intersection_, is_on_surface,
                id );
        }
        return corner_id;
    }

    bool get_side(
        const vec3& vertex,
        const vec3& on_surface,
        const Surface& surface,
        index_t triangle )
    {
        vec3 direction = vertex - on_surface;
        return dot( direction, surface.polygon_normal( triangle ) ) > 0;
    }

    index_t find_region_from_corners(
        const GeoModel& geomodel,
        const std::vector< vec3 > vertices,
        const LineInstersection& start,
        const LineInstersection& end )
    {
        if( start.surface_id_ != NO_ID ) {
            bool sign = get_side( vertices[1], start.intersection_,
                geomodel.surface( start.surface_id_ ), start.trgl_id_ );
            return find_region( geomodel, start.surface_id_, sign );
        } else if( end.surface_id_ != NO_ID ) {
            bool sign = get_side( vertices[vertices.size() - 2], end.intersection_,
                geomodel.surface( end.surface_id_ ), end.trgl_id_ );
            return find_region( geomodel, end.surface_id_, sign );
        } else {
            double best_distance = max_float64();
            index_t best_surface = NO_ID;
            vec3 best_nearest;
            index_t best_triangle = NO_ID;
            for( index_t s = 0; s < geomodel.nb_surfaces(); s++ ) {
                const Surface& surface = geomodel.surface( s );
                vec3 nearest;
                double distance;
                index_t triangle = surface.polygon_aabb().closest_triangle(
                    start.intersection_, nearest, distance );
                if( distance < best_distance ) {
                    best_distance = distance;
                    best_nearest = nearest;
                    best_surface = s;
                    best_triangle = triangle;
                }
            }
            bool sign = get_side( start.intersection_, best_nearest,
                geomodel.surface( best_surface ), best_triangle );
            return find_region( geomodel, best_surface, sign );
        }
    }

    void create_well_part_and_corners(
        const GeoModel& geomodel,
        Well& well,
        const std::vector< vec3 > vertices,
        const LineInstersection& start,
        const LineInstersection& end )
    {
        index_t region = find_region_from_corners( geomodel, vertices, start, end );
        if( region == NO_ID ) {
            return;
        }

        index_t new_well_part_id = well.create_part( region );
        WellPart& well_part = well.part( new_well_part_id );

        index_t corner0 = find_or_create_corner( well, region, start,
            geomodel.epsilon() );
        well_part.set_corner( 0, corner0 );

        index_t corner1 = find_or_create_corner( well, region, end,
            geomodel.epsilon() );
        well_part.set_corner( 1, corner1 );
        well_part.set_points( vertices );
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
            vec3 result;
            if( segment_triangle_intersection( v_from_, v_to_,
                surface_.mesh_element_vertex( trgl, 0 ),
                surface_.mesh_element_vertex( trgl, 1 ),
                surface_.mesh_element_vertex( trgl, 2 ), result ) ) {
                intersections_.push_back(
                    LineInstersection( result, surface_.index(), trgl ) );
            }
        }

    private:
        const Surface& surface_;
        const vec3& v_from_;
        const vec3& v_to_;

        std::vector< LineInstersection >& intersections_;
    };

    struct OrientedEdge {
        OrientedEdge( const LineMesh& mesh, index_t edge, index_t vertex_from )
            : edge_( edge ), vertex_from_( vertex_from )
        {
            if( mesh.edge_vertex( edge, 0 ) == vertex_from ) {
                edge_vertex_ = 0;
            } else {
                ringmesh_assert( mesh.edge_vertex( edge, 1 ) == vertex_from );
                edge_vertex_ = 1;
            }
        }
        index_t edge_;
        index_t edge_vertex_;
        index_t vertex_from_;
    };

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
            mesh_( PointSetMesh::create_mesh( GeogramPointSetMesh::type_name_static() ) )
    {
        PointSetMeshBuilder::create_builder( *mesh_ )->create_vertex( point );
    }

    const vec3& WellCorner::point() const
    {
        return mesh_->vertex( 0 );
    }

    GEO::AttributesManager& WellCorner::vertex_attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }

// --------------------------------------------------------------------------

    WellPart::WellPart( const Well* well, index_t id )
        :
            WellEntity( well ),
            id_( id ),
            mesh_( LineMesh::create_mesh( GeogramLineMesh::type_name_static() ) )
    {
        corners_[0] = NO_ID;
        corners_[1] = NO_ID;
    }

    void WellPart::set_points( const std::vector< vec3 >& points )
    {
        index_t nb_points = static_cast< index_t >( points.size() );
        std::unique_ptr< LineMeshBuilder > builder = LineMeshBuilder::create_builder(
            *mesh_ );
        builder->create_vertices( nb_points );
        for( index_t p = 0; p < nb_points; p++ ) {
            builder->set_vertex( p, points[p] );
        }

        index_t nb_edges = nb_points - 1;
        builder->create_edges( nb_edges );
        for( index_t e = 0; e < nb_edges; e++ ) {
            builder->set_edge_vertex( e, 0, e );
            builder->set_edge_vertex( e, 1, e + 1 );
        }
    }

    index_t WellPart::nb_edges() const
    {
        return mesh_->nb_edges();
    }

    index_t WellPart::nb_vertices() const
    {
        return mesh_->nb_vertices();
    }

    const vec3& WellPart::vertex( index_t v ) const
    {
        return mesh_->vertex( v );
    }

    const vec3& WellPart::edge_vertex( index_t edge, index_t v ) const
    {
        return vertex( mesh_->edge_vertex( edge, v ) );
    }

    const NNSearch& WellPart::vertices_nn_search() const
    {
        return mesh_->vertex_nn_search();
    }

    double WellPart::length() const
    {
        double l = 0.0;
        for( index_t e = 0; e < nb_edges(); e++ ) {
            l += ( vertex( e + 1 ) - vertex( e ) ).length();
        }
        return l;
    }

    GEO::AttributesManager& WellPart::vertex_attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }
    GEO::AttributesManager& WellPart::edge_attribute_manager() const
    {
        return mesh_->edge_attribute_manager();
    }

// --------------------------------------------------------------------------

    Well::Well()
        : nb_edges_( NO_ID )
    {
    }

    index_t Well::find_corner( const vec3& vertex, double epsilon ) const
    {
        for( index_t c = 0; c < nb_corners(); c++ ) {
            if( inexact_equal( vertex, corner( c ).point(), epsilon ) ) {
                return c;
            }
        }
        return NO_ID;
    }

    void Well::copy_corners_and_informations( Well& well ) const
    {
        well.name_ = name_;
        well.part_region_id_ = part_region_id_;

        well.corners_.reserve( nb_corners() );
        for( index_t c = 0; c < nb_corners(); c++ ) {
            well.create_corner( corners_[c]->point(), corners_[c]->is_on_surface(),
                corners_[c]->id() );
        }

        well.parts_.reserve( nb_parts() );
        for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
            well.create_part( part_region_id( part_id ) );
            const WellPart& from_part = part( part_id );
            WellPart& cur_part = well.part( part_id );
//            cur_part.mesh().copy( from_part.mesh() ) ;
            cur_part.set_corner( 0, from_part.corner( 0 ) );
            cur_part.set_corner( 1, from_part.corner( 1 ) );
        }
    }

    index_t Well::nb_edges() const
    {
        if( nb_edges_ == NO_ID ) {
            index_t nb_edges = 0;
            for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
                nb_edges += part( part_id ).nb_edges();
            }
            const_cast< Well* >( this )->nb_edges_ = nb_edges;
        }
        return nb_edges_;
    }

    void Well::get_part_edges( index_t part_id, std::vector< Edge >& edges ) const
    {
        const WellPart& well_part = part( part_id );
        for( index_t e = 0; e < well_part.nb_edges(); e++ ) {
            edges.push_back(
                Edge( well_part.vertex( e ), well_part.vertex( e + 1 ) ) );
        }
    }

    void Well::get_region_edges( index_t region, std::vector< Edge >& edges ) const
    {
        for( index_t part_id = 0; part_id < nb_parts(); part_id++ ) {
            if( part_region_id( part_id ) == region ) {
                get_part_edges( part_id, edges );
            }
        }
    }

// --------------------------------------------------------------------------

    WellGroup::WellGroup()
        : geomodel_( nullptr )
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
            const Well& cur_well = well( w );
            for( index_t part_id = 0; part_id < cur_well.nb_parts(); part_id++ ) {
                if( cur_well.part_region_id( part_id ) == region ) {
                    cur_well.get_part_edges( part_id, edges );
                }
            }
        }
    }

    void WellGroup::get_region_edges(
        index_t region,
        std::vector< std::vector< Edge > >& edges ) const
    {
        edges.clear();
        edges.resize( nb_wells() );
        for( index_t w = 0; w < nb_wells(); w++ ) {
            const Well& cur_well = well( w );
            cur_well.get_region_edges( region, edges[w] );
        }
    }

    void WellGroup::create_wells( index_t nb )
    {
        wells_.resize( nb, nullptr );
        for( index_t w = 0; w < nb_wells(); w++ ) {
            wells_[w] = new Well;
        }
    }

    void WellGroup::add_well( const LineMesh& mesh, const std::string& name )
    {
        ringmesh_assert( geomodel() );
        if( find_well( name ) != NO_ID ) return;
        wells_.push_back( new Well );
        Well& new_well = *wells_.back();
        new_well.set_name( name );

        GeogramLineMesh conformal_mesh;
        compute_conformal_mesh( mesh, conformal_mesh );

        std::vector< std::vector< index_t > > edges_around_vertices(
            conformal_mesh.nb_vertices() );
        for( index_t e = 0; e < conformal_mesh.nb_edges(); e++ ) {
            for( index_t i = 0; i < 2; i++ ) {
                index_t v = conformal_mesh.edge_vertex( e, i );
                edges_around_vertices[v].push_back( e );
            }
        }

        std::stack< OrientedEdge > S;
        for( index_t v = 0; v < conformal_mesh.nb_vertices(); v++ ) {
            const std::vector< index_t >& edges = edges_around_vertices[v];
            if( edges.size() == 1 ) {
                S.emplace( conformal_mesh, edges.front(), v );
            }
        }
        if( S.empty() ) {
            throw RINGMeshException( "Well",
                "A well should have at least one starting or ending point" );
        }

        GEO::Attribute< LineInstersection > vertex_info(
            conformal_mesh.vertex_attribute_manager(), "info" );
        std::vector< bool > edge_visited( conformal_mesh.nb_edges(), false );
        do {
            OrientedEdge cur_edge = S.top();
            S.pop();
            if( edge_visited[cur_edge.edge_] ) {
                continue;
            }
            edge_visited[cur_edge.edge_] = true;

            std::vector< vec3 > well_part_points;
            std::stack< OrientedEdge > S_part;
            S_part.push( cur_edge );
            do {
                OrientedEdge cur_edge_part = S_part.top();
                S_part.pop();
                edge_visited[cur_edge_part.edge_] = true;
                const vec3& v_from = conformal_mesh.vertex(
                    cur_edge_part.vertex_from_ );
                index_t v_to_id = conformal_mesh.edge_vertex( cur_edge_part.edge_,
                    ( cur_edge_part.edge_vertex_ + 1 ) % 2 );
                const vec3& v_to = conformal_mesh.vertex( v_to_id );
                well_part_points.push_back( v_from );

                const std::vector< index_t >& edges = edges_around_vertices[v_to_id];
                if( edges.size() == 2 ) {
                    index_t count = 0;
                    for( index_t edge : edges ) {
                        if( !edge_visited[edge] ) {
                            S_part.emplace( conformal_mesh, edge, v_to_id );
                            count++;
                        }
                    }
                    ringmesh_assert( count == 1 );
                } else {
                    well_part_points.push_back( v_to );
                    create_well_part_and_corners( *geomodel(), new_well,
                        well_part_points, vertex_info[cur_edge.vertex_from_],
                        vertex_info[v_to_id] );
                    for( index_t edge : edges ) {
                        S.emplace( conformal_mesh, edge, v_to_id );
                    }
                }
            } while( !S_part.empty() );
        } while( !S.empty() );
    }

    void WellGroup::compute_conformal_mesh( const LineMesh& in, LineMesh& out )
    {
        double epsilon = geomodel_->epsilon();
        std::unique_ptr< LineMeshBuilder > builder = LineMeshBuilder::create_builder(
            out );
        builder->clear( false, false );

        GEO::Attribute< LineInstersection > vertex_info(
            out.vertex_attribute_manager(), "info" );
        builder->create_vertices( in.nb_vertices() );
        for( index_t v = 0; v < in.nb_vertices(); v++ ) {
            const vec3& vertex = in.vertex( v );
            builder->set_vertex( v, vertex );
            vertex_info[v] = LineInstersection( vertex );
        }

        for( index_t e = 0; e < in.nb_edges(); e++ ) {
            index_t from_id = in.edge_vertex( e, 0 );
            const vec3& from_vertex = in.vertex( from_id );
            index_t to_id = in.edge_vertex( e, 1 );
            const vec3& to_vertex = in.vertex( to_id );

            Box3d box;
            box.add_point( from_vertex );
            box.add_point( to_vertex );
            std::vector< LineInstersection > intersections;
            for( index_t s = 0; s < geomodel()->nb_surfaces(); s++ ) {
                const Surface& surface = geomodel()->surface( s );
                EdgeConformerAction action( surface, from_vertex, to_vertex,
                    intersections );
                surface.polygon_aabb().compute_bbox_element_bbox_intersections( box,
                    action );
            }

            std::vector< index_t > indices( intersections.size() );
            std::vector< double > distances( intersections.size() );
            for( index_t i = 0; i < intersections.size(); i++ ) {
                indices[i] = i;
                distances[i] = length(
                    from_vertex - intersections[i].intersection_ );
            }
            indirect_sort( distances, indices );
            double edge_length = length( from_vertex - to_vertex );
            index_t last_vertex = from_id;
            for( index_t i = 0; i < intersections.size(); i++ ) {
                if( distances[indices[i]] < epsilon ) {
                    vertex_info[from_id] = intersections[i];
                } else if( std::fabs( distances[indices[i]] - edge_length )
                    < epsilon ) {
                    vertex_info[to_id] = intersections[i];
                } else {
                    index_t vertex_id = builder->create_vertex(
                        intersections[i].intersection_ );
                    vertex_info[vertex_id] = intersections[i];
                    builder->create_edge( last_vertex, vertex_id );
                    last_vertex = vertex_id;
                }
            }
            builder->create_edge( last_vertex, to_id );
        }
    }

    index_t WellGroup::find_well( const std::string& name ) const
    {
        for( index_t w = 0; w < nb_wells(); w++ ) {
            if( well( w ).name() == name ) {
                return w;
            }
        }
        return NO_ID;
    }
}
