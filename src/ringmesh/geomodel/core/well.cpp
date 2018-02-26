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

#include <ringmesh/geomodel/core/well.h>

#include <cmath>
#include <numeric>
#include <stack>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/box.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>


#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>

/*!
 * @file Implementation of Wells
 * @author Arnaud Botella
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    bool inexact_equal( const vecn< DIMENSION >& v1,
        const vecn< DIMENSION >& v2,
        double epsilon )
    {
        return length( v2 - v1 ) < epsilon;
    }

    /*!
     * @brief Returns the index of the region of the geomodel neighboring the
     * surface.
     * @param[in] geomodel to consider
     * @param[in] surface_part_id Index of the Surface
     * @param[in] side Side of the Surface
     * @return The region index or NO_ID if none found.
     */
    index_t find_region(
        const GeoModel3D& geomodel, index_t surface_part_id, bool side )
    {
        ringmesh_assert( surface_part_id < geomodel.nb_surfaces() );
        gmme_id cur_surface( Surface3D::type_name_static(), surface_part_id );
        const Surface3D& surface = geomodel.surface( surface_part_id );
        for( auto r : range( surface.nb_incident_entities() ) )
        {
            const Region3D& cur_region = surface.incident_entity( r );
            for( auto s : range( cur_region.nb_boundaries() ) )
            {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_gmme( s ) == cur_surface )
                {
                    return r;
                }
            }
        }
        return NO_ID;
    }

    struct LineInstersection
    {
        explicit LineInstersection( const vec3& intersection,
            index_t surface_id = NO_ID,
            index_t trgl_id = NO_ID )
            : intersection_( intersection ),
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

    template < index_t DIMENSION >
    index_t find_or_create_corner( Well< DIMENSION >& well,
        index_t region,
        const LineInstersection& corner,
        double epsilon )
    {
        index_t corner_id = well.find_corner( corner.intersection_, epsilon );
        if( corner_id == NO_ID )
        {
            bool is_on_surface = corner.surface_id_ != NO_ID;
            index_t id = NO_ID;
            if( is_on_surface )
            {
                id = corner.surface_id_;
            }
            else
            {
                id = region;
            }
            corner_id =
                well.create_corner( corner.intersection_, is_on_surface, id );
        }
        return corner_id;
    }

    bool get_side( const vec3& vertex,
        const vec3& on_surface,
        const Surface3D& surface,
        index_t triangle )
    {
        vec3 direction = vertex - on_surface;
        return dot( direction, surface.mesh().polygon_normal( triangle ) ) > 0;
    }

    index_t find_region_from_corners( const GeoModel3D& geomodel,
        const std::vector< vec3 >& vertices,
        const LineInstersection& start,
        const LineInstersection& end )
    {
        if( start.surface_id_ != NO_ID )
        {
            bool sign = get_side( vertices[1], start.intersection_,
                geomodel.surface( start.surface_id_ ), start.trgl_id_ );
            return find_region( geomodel, start.surface_id_, sign );
        }
        if( end.surface_id_ != NO_ID )
        {
            bool sign =
                get_side( vertices[vertices.size() - 2], end.intersection_,
                    geomodel.surface( end.surface_id_ ), end.trgl_id_ );
            return find_region( geomodel, end.surface_id_, sign );
        }
        double best_distance = max_float64();
        index_t best_surface = NO_ID;
        vec3 best_nearest;
        index_t best_triangle = NO_ID;
        for( const auto& surface : geomodel.surfaces() )
        {
            index_t triangle = NO_ID;
            vec3 nearest;
            double distance = max_float64();
            std::tie( triangle, nearest, distance ) =
                surface.polygon_aabb().closest_triangle( start.intersection_ );
            if( distance < best_distance )
            {
                best_distance = distance;
                best_nearest = nearest;
                best_surface = surface.index();
                best_triangle = triangle;
            }
        }
        bool sign = get_side( start.intersection_, best_nearest,
            geomodel.surface( best_surface ), best_triangle );
        return find_region( geomodel, best_surface, sign );
    }

    template < index_t DIMENSION >
    void create_well_part_and_corners( const GeoModel3D& geomodel,
        Well< DIMENSION >& well,
        const std::vector< vec3 > vertices,
        const LineInstersection& start,
        const LineInstersection& end )
    {
        index_t region =
            find_region_from_corners( geomodel, vertices, start, end );
        if( region == NO_ID )
        {
            return;
        }

        index_t new_well_part_id = well.create_part( region );
        WellPart< DIMENSION >& well_part = well.part( new_well_part_id );

        index_t corner0 =
            find_or_create_corner( well, region, start, geomodel.epsilon() );
        well_part.set_corner( 0, corner0 );

        index_t corner1 =
            find_or_create_corner( well, region, end, geomodel.epsilon() );
        well_part.set_corner( 1, corner1 );
        well_part.set_points( vertices );
    }

    class EdgeConformerAction
    {
    public:
        EdgeConformerAction( const Surface3D& surface,
            const vec3& v_from,
            const vec3& v_to,
            std::vector< LineInstersection >& intersections )
            : surface_( surface ),
              segment_( v_from, v_to ),
              intersections_( intersections )
        {
        }

        void operator()( index_t trgl )
        {
            bool does_seg_intersect_triangle = false;
            vec3 result;
            std::tie( does_seg_intersect_triangle, result ) =
                Intersection::segment_triangle( segment_,
                    { surface_.mesh_element_vertex( { trgl, 0 } ),
                        surface_.mesh_element_vertex( { trgl, 1 } ),
                        surface_.mesh_element_vertex( { trgl, 2 } ) } );
            if( does_seg_intersect_triangle )
            {
                intersections_.emplace_back( result, surface_.index(), trgl );
            }
        }

    private:
        const Surface3D& surface_;
        Geometry::Segment3D segment_;

        std::vector< LineInstersection >& intersections_;
    };

    struct OrientedEdge
    {
        OrientedEdge(
            const LineMesh3D& mesh, index_t edge, index_t vertex_from )
            : edge_( edge ), vertex_from_( vertex_from )
        {
            if( mesh.edge_vertex( ElementLocalVertex( edge, 0 ) )
                == vertex_from )
            {
                edge_vertex_ = 0;
            }
            else
            {
                ringmesh_assert(
                    mesh.edge_vertex( ElementLocalVertex( edge, 1 ) )
                    == vertex_from );
                edge_vertex_ = 1;
            }
        }
        index_t edge_;
        index_t edge_vertex_;
        index_t vertex_from_;
    };

    std::vector< std::vector< index_t > > get_edges_around_vertices(
        const LineMesh3D& conformal_mesh )
    {
        std::vector< std::vector< index_t > > edges_around_vertices(
            conformal_mesh.nb_vertices() );
        for( auto e : range( conformal_mesh.nb_edges() ) )
        {
            for( auto i : range( 2 ) )
            {
                auto v = conformal_mesh.edge_vertex( { e, i } );
                edges_around_vertices[v].push_back( e );
            }
        }
        return edges_around_vertices;
    }

    void process_linear_edges(
        const std::vector< index_t >& edges,
        std::vector< bool >& edge_visited,
        const LineMesh3D& conformal_mesh,
        std::stack< OrientedEdge >& S_part,
        index_t v_to_id )
    {
        index_t count = 0;
        for( auto edge : edges )
        {
            if( !edge_visited[edge] )
            {
                S_part.emplace( conformal_mesh, edge, v_to_id );
                count++;
            }
        }
        ringmesh_assert( count == 1 );
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    WellEntity< DIMENSION >::WellEntity( const Well< DIMENSION >* well )
        : well_( well )
    {
    }

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    WellCorner< DIMENSION >::WellCorner( const Well< DIMENSION >* well,
        const vecn< DIMENSION >& point,
        bool is_on_surface,
        index_t id )
        : WellEntity< DIMENSION >( well ),
          is_on_surface_( is_on_surface ),
          id_( id ),
          mesh_( PointSetMesh< DIMENSION >::create_mesh() )
    {
        PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ )
            ->create_vertex( point );
    }

    template < index_t DIMENSION >
    WellCorner< DIMENSION >::~WellCorner()
    {
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& WellCorner< DIMENSION >::point() const
    {
        return mesh_->vertex( 0 );
    }

    template < index_t DIMENSION >
    GEO::AttributesManager&
        WellCorner< DIMENSION >::vertex_attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    WellPart< DIMENSION >::WellPart( const Well< DIMENSION >* well, index_t id )
        : WellEntity< DIMENSION >( well ),
          id_( id ),
          mesh_( LineMesh< DIMENSION >::create_mesh() )
    {
        corners_[0] = NO_ID;
        corners_[1] = NO_ID;
    }

    template < index_t DIMENSION >
    WellPart< DIMENSION >::~WellPart()
    {
    }

    template < index_t DIMENSION >
    void WellPart< DIMENSION >::set_points(
        const std::vector< vecn< DIMENSION > >& points )
    {
        index_t nb_points = static_cast< index_t >( points.size() );
        std::unique_ptr< LineMeshBuilder< DIMENSION > > builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        builder->create_vertices( nb_points );
        for( auto p : range( nb_points ) )
        {
            builder->set_vertex( p, points[p] );
        }

        index_t nb_edges = nb_points - 1;
        builder->create_edges( nb_edges );
        for( auto e : range( nb_edges ) )
        {
            builder->set_edge_vertex( EdgeLocalVertex( e, 0 ), e );
            builder->set_edge_vertex( EdgeLocalVertex( e, 1 ), e + 1 );
        }
    }

    template < index_t DIMENSION >
    index_t WellPart< DIMENSION >::nb_edges() const
    {
        return mesh_->nb_edges();
    }

    template < index_t DIMENSION >
    index_t WellPart< DIMENSION >::nb_vertices() const
    {
        return mesh_->nb_vertices();
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& WellPart< DIMENSION >::vertex( index_t v ) const
    {
        return mesh_->vertex( v );
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& WellPart< DIMENSION >::edge_vertex(
        const ElementLocalVertex& well_edge_local_vertex ) const
    {
        return vertex( mesh_->edge_vertex( well_edge_local_vertex ) );
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >&
        WellPart< DIMENSION >::vertices_nn_search() const
    {
        return mesh_->vertex_nn_search();
    }

    template < index_t DIMENSION >
    double WellPart< DIMENSION >::length() const
    {
        double l = 0.0;
        for( auto e : range( nb_edges() ) )
        {
            l += ( vertex( e + 1 ) - vertex( e ) ).length();
        }
        return l;
    }

    template < index_t DIMENSION >
    GEO::AttributesManager&
        WellPart< DIMENSION >::vertex_attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }
    template < index_t DIMENSION >
    GEO::AttributesManager&
        WellPart< DIMENSION >::edge_attribute_manager() const
    {
        return mesh_->edge_attribute_manager();
    }

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    index_t Well< DIMENSION >::find_corner(
        const vecn< DIMENSION >& vertex, double epsilon ) const
    {
        for( auto c : range( nb_corners() ) )
        {
            if( inexact_equal( vertex, corner( c ).point(), epsilon ) )
            {
                return c;
            }
        }
        return NO_ID;
    }

    template < index_t DIMENSION >
    void Well< DIMENSION >::copy_corners_and_informations(
        Well< DIMENSION >& well ) const
    {
        well.name_ = name_;
        well.part_region_id_ = part_region_id_;

        well.corners_.reserve( nb_corners() );
        for( auto c : range( nb_corners() ) )
        {
            well.create_corner( corners_[c]->point(),
                corners_[c]->is_on_surface(), corners_[c]->id() );
        }

        well.parts_.reserve( nb_parts() );
        for( auto part_id : range( nb_parts() ) )
        {
            well.create_part( part_region_id( part_id ) );
            const WellPart< DIMENSION >& from_part = part( part_id );
            WellPart< DIMENSION >& cur_part = well.part( part_id );
            cur_part.set_corner( 0, from_part.corner( 0 ) );
            cur_part.set_corner( 1, from_part.corner( 1 ) );
        }
    }

    template < index_t DIMENSION >
    index_t Well< DIMENSION >::nb_edges() const
    {
        if( nb_edges_ == NO_ID )
        {
            index_t nb_edges = 0;
            for( auto part_id : range( nb_parts() ) )
            {
                nb_edges += part( part_id ).nb_edges();
            }
            nb_edges_ = nb_edges;
        }
        return nb_edges_;
    }

    template < index_t DIMENSION >
    void Well< DIMENSION >::get_part_edges(
        index_t part_id, std::vector< Edge< DIMENSION > >& edges ) const
    {
        const WellPart< DIMENSION >& well_part = part( part_id );
        for( auto e : range( well_part.nb_edges() ) )
        {
            edges.emplace_back(
                well_part.vertex( e ), well_part.vertex( e + 1 ) );
        }
    }

    template < index_t DIMENSION >
    void Well< DIMENSION >::get_region_edges(
        index_t region, std::vector< Edge< DIMENSION > >& edges ) const
    {
        for( auto part_id : range( nb_parts() ) )
        {
            if( part_region_id( part_id ) == region )
            {
                get_part_edges( part_id, edges );
            }
        }
    }

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    WellGroup< DIMENSION >::WellGroup() : geomodel_( nullptr )
    {
    }

    template < index_t DIMENSION >
    void WellGroup< DIMENSION >::get_region_edges(
        index_t region, std::vector< Edge< DIMENSION > >& edges ) const
    {
        for( auto w : range( nb_wells() ) )
        {
            const Well< DIMENSION >& cur_well = well( w );
            for( auto part_id : range( cur_well.nb_parts() ) )
            {
                if( cur_well.part_region_id( part_id ) == region )
                {
                    cur_well.get_part_edges( part_id, edges );
                }
            }
        }
    }

    template < index_t DIMENSION >
    void WellGroup< DIMENSION >::get_region_edges( index_t region,
        std::vector< std::vector< Edge< DIMENSION > > >& edges ) const
    {
        edges.clear();
        edges.resize( nb_wells() );
        for( auto w : range( nb_wells() ) )
        {
            const Well< DIMENSION >& cur_well = well( w );
            cur_well.get_region_edges( region, edges[w] );
        }
    }

    template < index_t DIMENSION >
    void WellGroup< DIMENSION >::create_wells( index_t nb_wells )
    {
        wells_.resize( nb_wells, nullptr );
        for( auto w : range( this->nb_wells() ) )
        {
            wells_[w] = new Well< DIMENSION >;
        }
    }

    template <>
    void WellGroup< 2 >::compute_conformal_mesh(
        const LineMesh< 2 >& in, LineMesh< 2 >& out )
    {
        ringmesh_unused( in );
        ringmesh_unused( out );
        throw RINGMeshException(
            "Wells", "2D Wells not fully implemented yet" );
    }

    template <>
    void WellGroup< 3 >::compute_conformal_mesh(
        const LineMesh3D& in, LineMesh3D& out )
    {
        double epsilon = geomodel_->epsilon();
        std::unique_ptr< LineMeshBuilder3D > builder =
            LineMeshBuilder3D::create_builder( out );
        builder->clear( false, false );

        GEO::Attribute< LineInstersection > vertex_info(
            out.vertex_attribute_manager(), "info" );
        builder->create_vertices( in.nb_vertices() );
        for( auto v : range( in.nb_vertices() ) )
        {
            const vec3& vertex = in.vertex( v );
            builder->set_vertex( v, vertex );
            vertex_info[v] = LineInstersection( vertex );
        }

        for( auto e : range( in.nb_edges() ) )
        {
            index_t from_id = in.edge_vertex( ElementLocalVertex( e, 0 ) );
            const vec3& from_vertex = in.vertex( from_id );
            index_t to_id = in.edge_vertex( ElementLocalVertex( e, 1 ) );
            const vec3& to_vertex = in.vertex( to_id );

            Box3D box;
            box.add_point( from_vertex );
            box.add_point( to_vertex );
            std::vector< LineInstersection > intersections;
            for( const auto& surface : geomodel_->surfaces() )
            {
                EdgeConformerAction action(
                    surface, from_vertex, to_vertex, intersections );
                surface.polygon_aabb().compute_bbox_element_bbox_intersections(
                    box, action );
            }

            std::vector< index_t > indices( intersections.size() );
            std::iota( indices.begin(), indices.end(), 0 );
            std::vector< double > distances( intersections.size() );
            for( auto i : range( intersections.size() ) )
            {
                distances[i] =
                    length( from_vertex - intersections[i].intersection_ );
            }
            indirect_sort( distances, indices );
            double edge_length = length( from_vertex - to_vertex );
            index_t last_vertex = from_id;
            for( auto i : range( intersections.size() ) )
            {
                if( distances[indices[i]] < epsilon )
                {
                    vertex_info[from_id] = intersections[i];
                }
                else if( std::fabs( distances[indices[i]] - edge_length )
                         < epsilon )
                {
                    vertex_info[to_id] = intersections[i];
                }
                else
                {
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

    template <>
    void WellGroup< 2 >::add_well(
        const LineMesh< 2 >& mesh, const std::string& name )
    {
        ringmesh_unused( mesh );
        ringmesh_unused( name );
        throw RINGMeshException(
            "Wells", "2D Wells not fully implemented yet" );
    }

    template <>
    void WellGroup< 3 >::add_well(
        const LineMesh3D& mesh, const std::string& name )
    {
        ringmesh_assert( geomodel() );
        if( find_well( name ) != NO_ID )
        {
            return;
        }
        wells_.push_back( new Well3D );
        Well3D& new_well = *wells_.back();
        new_well.set_name( name );

        auto conformal_mesh = LineMesh3D::create_mesh();
        compute_conformal_mesh( mesh, *conformal_mesh );

        auto edges_around_vertices = get_edges_around_vertices( *conformal_mesh );

        std::stack< OrientedEdge > S;
        for( auto v : range( conformal_mesh->nb_vertices() ) )
        {
            const auto& edges = edges_around_vertices[v];
            if( edges.size() == 1 )
            {
                S.emplace( *conformal_mesh, edges.front(), v );
            }
        }
        if( S.empty() )
        {
            throw RINGMeshException( "Well",
                "A well should have at least one starting or ending point" );
        }

        GEO::Attribute< LineInstersection > vertex_info(
            conformal_mesh->vertex_attribute_manager(), "info" );
        std::vector< bool > edge_visited( conformal_mesh->nb_edges(), false );
        do
        {
            OrientedEdge cur_edge = S.top();
            S.pop();
            if( edge_visited[cur_edge.edge_] )
            {
                continue;
            }
            edge_visited[cur_edge.edge_] = true;

            std::vector< vec3 > well_part_points;
            std::stack< OrientedEdge > S_part;
            S_part.push( cur_edge );
            do
            {
                OrientedEdge cur_edge_part = S_part.top();
                S_part.pop();
                edge_visited[cur_edge_part.edge_] = true;
                const vec3& v_from =
                    conformal_mesh->vertex( cur_edge_part.vertex_from_ );
                index_t v_to_id = conformal_mesh->edge_vertex(
                    ElementLocalVertex( cur_edge_part.edge_,
                        ( cur_edge_part.edge_vertex_ + 1 ) % 2 ) );
                const vec3& v_to = conformal_mesh->vertex( v_to_id );
                well_part_points.push_back( v_from );

                const auto& edges = edges_around_vertices[v_to_id];
                if( edges.size() == 2 )
                {
                    process_linear_edges( edges, edge_visited, *conformal_mesh,
                        S_part, v_to_id );
                }
                else
                {
                    well_part_points.push_back( v_to );
                    create_well_part_and_corners( *geomodel(), new_well,
                        well_part_points, vertex_info[cur_edge.vertex_from_],
                        vertex_info[v_to_id] );
                    for( auto edge : edges )
                    {
                        S.emplace( *conformal_mesh, edge, v_to_id );
                    }
                }
            } while( !S_part.empty() );
        } while( !S.empty() );
    }

    template < index_t DIMENSION >
    index_t WellGroup< DIMENSION >::find_well( const std::string& name ) const
    {
        for( auto w : range( nb_wells() ) )
        {
            if( well( w ).name() == name )
            {
                return w;
            }
        }
        return NO_ID;
    }
    template class geomodel_core_api WellEntity< 2 >;
    template class geomodel_core_api WellCorner< 2 >;
    template class geomodel_core_api WellPart< 2 >;
    template class geomodel_core_api Well< 2 >;
    template class geomodel_core_api WellGroup< 2 >;

    template class geomodel_core_api WellEntity< 3 >;
    template class geomodel_core_api WellCorner< 3 >;
    template class geomodel_core_api WellPart< 3 >;
    template class geomodel_core_api Well< 3 >;
    template class geomodel_core_api WellGroup< 3 >;
} // namespace RINGMesh
