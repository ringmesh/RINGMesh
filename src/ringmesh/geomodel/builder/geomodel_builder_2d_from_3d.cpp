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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder_2d_from_3d.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh_index.h>

namespace
{
    using namespace RINGMesh;

    const MeshEntityType projectable_entity_types[3] = {
        corner_type_name_static(), line_type_name_static(),
        surface_type_name_static()
    };

    const std::map< GeologicalEntityType, GeologicalEntityType >
        geol_entity_type_2d_to_3d_map = { { Contact3D::type_name_static(),
                                              Interface2D::type_name_static() },
            { Interface3D::type_name_static(), Layer2D::type_name_static() } };

    template < typename U, typename T >
    const T& mapped_value( const std::map< U, T >& map, const U& key )
    {
        return map.find( key )->second;
    }
} // namespace

namespace RINGMesh
{
    GeoModelBuilder2DFrom3D::GeoModelBuilder2DFrom3D( GeoModel2D& geomodel2d,
        const GeoModel3D& geomodel3d_from,
        const Geometry::Plane& plane )
        : GeoModelBuilder( geomodel2d ),
          geomodel3d_from_( geomodel3d_from ),
          plane_( plane )
    {
        PlaneReferenceFrame3D plane_frame( plane );
        u_axis_ = std::move( plane_frame[0] );
        v_axis_ = std::move( plane_frame[1] );
    }

    vec2 GeoModelBuilder2DFrom3D::get_2d_coord( const vec3& coord3d )
    {
        return { dot( coord3d, u_axis_ ), dot( coord3d, v_axis_ ) };
    }

    GeoModelBuilder2DProjection::GeoModelBuilder2DProjection(
        GeoModel2D& geomodel2d,
        const GeoModel3D& geomodel3d_from,
        const Geometry::Plane& plane )
        : GeoModelBuilder2DFrom3D( geomodel2d, geomodel3d_from, plane )
    {
        info.set_geomodel_name( geomodel3d_from_.name() + "_projected" );
    }

    void GeoModelBuilder2DProjection::build_geomodel()
    {
        copy_geomodel_3d_topology();
        copy_geomodel_3d_geological_informations();
        project_geomodel_3d_mesh_entities();
        print_geomodel( geomodel_ );
    }

    void GeoModelBuilder2DProjection::copy_geomodel_3d_topology()
    {
        for( const auto& entity_type : projectable_entity_types )
        {
            topology.create_mesh_entities(
                entity_type, geomodel3d_from_.nb_mesh_entities( entity_type ) );
        }

        for( const auto& line : geomodel3d_from_.lines() )
        {
            for( auto boundary_id : range( line.nb_boundaries() ) )
            {
                topology.add_line_corner_boundary_relation(
                    line.gmme().index(), line.boundary( boundary_id ).index() );
            }
        }

        for( const auto& surface : geomodel3d_from_.surfaces() )
        {
            for( auto boundary_id : range( surface.nb_boundaries() ) )
            {
                topology.add_surface_line_boundary_relation(
                    surface.gmme().index(),
                    surface.boundary( boundary_id ).index(),
                    false ); // TODO side
            }
        }
    }

    void GeoModelBuilder2DProjection::copy_geomodel_3d_geological_informations()
    {
        for( const auto& geol_entity_id :
            range( geomodel3d_from_.nb_geological_entity_types() ) )
        {
            const auto& cur_type =
                geomodel3d_from_.geological_entity_type( geol_entity_id );
            geology.create_geological_entities(
                mapped_value( geol_entity_type_2d_to_3d_map, cur_type ),
                geomodel3d_from_.nb_geological_entities( cur_type ) );
        }
        for( const auto& geol_entity_id :
            range( geomodel3d_from_.nb_geological_entity_types() ) )
        {
            const auto& cur_type =
                geomodel3d_from_.geological_entity_type( geol_entity_id );
            for( const auto& cur_geol_entity :
                geomodel3d_from_.geol_entities( cur_type ) )
            {
                for( const auto& child_id :
                    range( cur_geol_entity.nb_children() ) )
                {
                    geology.add_parent_children_relation(
                        { mapped_value(
                              geol_entity_type_2d_to_3d_map, cur_type ),
                            cur_geol_entity.index() },
                        cur_geol_entity.child_gmme( child_id ) );
                }
            }
        }
    }

    void GeoModelBuilder2DProjection::project_geomodel_3d_mesh_entities()
    {
        for( const auto& corner : geomodel3d_from_.corners() )
        {
            auto projected_vertices = compute_projected_vertices( corner );
            ringmesh_assert( projected_vertices.size() == 1 );
            geometry.set_corner( corner.index(), projected_vertices.front() );
        }

        for( const auto& line : geomodel3d_from_.lines() )
        {
            auto projected_vertices = compute_projected_vertices( line );
            geometry.set_line( line.index(), projected_vertices );
        }

        for( const auto& surface : geomodel3d_from_.surfaces() )
        {
            auto projected_vertices = compute_projected_vertices( surface );
            std::vector< index_t > surface_polygons;
            surface_polygons.reserve( 4 * surface.nb_mesh_elements() );
            std::vector< index_t > surface_polygon_ptr( 1, 0 );
            surface_polygon_ptr.reserve( surface.nb_mesh_elements() + 1 );
            for( const auto& polygon_id : range( surface.nb_mesh_elements() ) )
            {
                for( const auto& v_id :
                    range( surface.nb_mesh_element_vertices( polygon_id ) ) )
                {
                    surface_polygons.push_back(
                        surface.mesh_element_vertex_index(
                            { polygon_id, v_id } ) );
                }
                surface_polygon_ptr.push_back(
                    surface_polygon_ptr.back()
                    + surface.nb_mesh_element_vertices( polygon_id ) );
            }
            geometry.set_surface_geometry( surface.index(), projected_vertices,
                surface_polygons, surface_polygon_ptr );
        }
    }

    std::vector< vec2 > GeoModelBuilder2DProjection::compute_projected_vertices(
        const GeoModelMeshEntity3D& entity )
    {
        std::vector< vec2 > projected_vertices;
        projected_vertices.reserve( entity.nb_vertices() );
        for( const auto v : range( entity.nb_vertices() ) )
        {
            projected_vertices.push_back( get_2d_coord( entity.vertex( v ) ) );
        }
        return projected_vertices;
    }
} // namespace RINGMesh
