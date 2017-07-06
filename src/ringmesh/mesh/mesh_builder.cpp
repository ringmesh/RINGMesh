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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh_builder.h>

#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace {
    using namespace RINGMesh;

    template< index_t DIMENSION >
    PointSetMeshBuilder< DIMENSION >* create_point_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        PointSetMeshBuilder< DIMENSION >* builder = PointMeshBuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->configure_builder(
                dynamic_cast< PointSetMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    LineMeshBuilder< DIMENSION >* create_line_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        LineMeshBuilder< DIMENSION >* builder =
            LineMeshBuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->configure_builder(
                dynamic_cast< LineMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    SurfaceMeshBuilder< DIMENSION >* create_surface_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        SurfaceMeshBuilder< DIMENSION >* builder = SurfaceMeshBuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->configure_builder(
                dynamic_cast< SurfaceMeshBase< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    VolumeMeshBuilder< DIMENSION >* create_volume_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        VolumeMeshBuilder< DIMENSION >* builder =
            VolumeMeshBuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->configure_builder(
                dynamic_cast< VolumeMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    MeshBaseBuilder< DIMENSION >* create_pointset_builder(
        MeshBase< DIMENSION >& mesh )
    {
        MeshBaseBuilder< DIMENSION >* builder = nullptr;
        PointSetMeshBuilder< DIMENSION >* builder0d = create_point_mesh_builder(
            mesh );
        if( builder0d ) {
            builder = builder0d;
        } else {
            LineMeshBuilder< DIMENSION >* builder1d = create_line_mesh_builder(
                mesh );
            if( builder1d ) {
                builder = builder1d;
            } else {
                SurfaceMeshBuilder< DIMENSION >* builder2d =
                    create_surface_mesh_builder( mesh );
                if( builder2d ) {
                    builder = builder2d;
                }
            }
        }
        return builder;
    }

}

namespace RINGMesh {

    template< index_t DIMENSION >
    void SurfaceMeshBuilder< DIMENSION >::connect_polygons()
    {
        std::vector< index_t > polygons_to_connect( surface_mesh_->nb_polygons() );
        std::iota( polygons_to_connect.begin(), polygons_to_connect.end(), 0 );
        connect_polygons( polygons_to_connect );
    }

    template< index_t DIMENSION >
    void SurfaceMeshBuilder< DIMENSION >::connect_polygons(
        const std::vector< index_t >& polygons_to_connect )
    {
        struct PolygonLocalVertex {
            PolygonLocalVertex( index_t polygon, index_t local_vertex )
                : polygon_( polygon ), local_vertex_( local_vertex )
            {
            }
            index_t polygon_ { NO_ID };
            index_t local_vertex_ { NO_ID };
        };

        index_t nb_local_vertices = 0;
        for( index_t polygon : polygons_to_connect ) {
            nb_local_vertices += surface_mesh_->nb_polygon_vertices( polygon );
        }

        std::vector< PolygonLocalVertex > polygon_vertices;
        polygon_vertices.reserve( nb_local_vertices );
        for( index_t polygon : polygons_to_connect ) {
            for( index_t v : range( surface_mesh_->nb_polygon_vertices( polygon ) ) ) {
                polygon_vertices.emplace_back( polygon, v );
            }
        }

        std::vector< index_t > next_local_vertex_around_vertex( nb_local_vertices,
            NO_ID );
        std::vector< index_t > vertex2polygon_local_vertex(
            surface_mesh_->nb_vertices(), NO_ID );
        index_t local_vertex_count = 0;
        for( index_t polygon : polygons_to_connect ) {
            for( index_t v = 0; v < surface_mesh_->nb_polygon_vertices( polygon );
                v++, local_vertex_count++ ) {
                index_t vertex = surface_mesh_->polygon_vertex( polygon, v );
                next_local_vertex_around_vertex[local_vertex_count] =
                    vertex2polygon_local_vertex[vertex];
                vertex2polygon_local_vertex[vertex] = local_vertex_count;
            }
        }

        local_vertex_count = 0;
        for( index_t polygon : polygons_to_connect ) {
            for( index_t v = 0; v < surface_mesh_->nb_polygon_vertices( polygon );
                v++, local_vertex_count++ ) {
                if( !surface_mesh_->is_edge_on_border( polygon, v ) ) {
                    continue;
                }
                index_t vertex = surface_mesh_->polygon_vertex( polygon, v );
                index_t next_vertex = surface_mesh_->polygon_vertex( polygon,
                    surface_mesh_->next_polygon_vertex( polygon, v ) );
                for( index_t local_vertex = vertex2polygon_local_vertex[next_vertex];
                    local_vertex != NO_ID; local_vertex =
                        next_local_vertex_around_vertex[local_vertex] ) {
                    if( local_vertex == local_vertex_count ) {
                        continue;
                    }
                    index_t adj_polygon = polygon_vertices[local_vertex].polygon_;
                    index_t adj_local_vertex =
                        polygon_vertices[local_vertex].local_vertex_;
                    index_t adj_next_vertex = surface_mesh_->polygon_vertex( polygon,
                        surface_mesh_->next_polygon_vertex( polygon,
                            adj_local_vertex ) );
                    if( adj_next_vertex == vertex ) {
                        set_polygon_adjacent( polygon, v, adj_polygon );
                        set_polygon_adjacent( adj_polygon, adj_local_vertex,
                            polygon );
                        break;
                    }
                }
            }
        }
    }

    template< >
    std::unique_ptr< MeshBaseBuilder< 2 > > MeshBaseBuilder< 2 >::create_builder(
        MeshBase< 2 >& mesh )
    {
        MeshBaseBuilder< 2 >* builder = create_pointset_builder( mesh );
        if( !builder ) {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh data structure: " + mesh.type_name() );
        }
        return std::unique_ptr< MeshBaseBuilder< 2 > >( builder );
    }

    template< >
    std::unique_ptr< MeshBaseBuilder< 3 > > MeshBaseBuilder< 3 >::create_builder(
        MeshBase< 3 >& mesh )
    {
        MeshBaseBuilder< 3 >* builder = create_pointset_builder( mesh );
        if( !builder ) {
            VolumeMeshBuilder< 3 >* builder3d = create_volume_mesh_builder( mesh );
            if( builder3d ) {
                builder = builder3d;
            }
        }
        if( !builder ) {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh data structure: " + mesh.type_name() );
        }
        return std::unique_ptr< MeshBaseBuilder< 3 > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< PointSetMeshBuilder< DIMENSION > > PointSetMeshBuilder<
        DIMENSION >::create_builder( PointSetMesh< DIMENSION >& mesh )
    {
        PointSetMeshBuilder< DIMENSION > *builder = PointMeshBuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_point_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh0DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "PointMeshBuilder",
                "Falling back to GeogramPointSetMeshBuilder data structure" );

            builder = new GeogramPointSetMeshBuilder< DIMENSION >;
            builder->configure_builder( mesh );
        }
        return std::unique_ptr< PointSetMeshBuilder >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< LineMeshBuilder< DIMENSION > > LineMeshBuilder< DIMENSION >::create_builder(
        LineMesh< DIMENSION >& mesh )
    {
        LineMeshBuilder< DIMENSION >* builder =
            LineMeshBuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_line_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "LineMeshBuilder",
                "Could not create mesh data structure: ", mesh.type_name() );
            Logger::warn( "LineMeshBuilder",
                "Falling back to GeogramLineMeshBuilder data structure" );

            builder = new GeogramLineMeshBuilder< DIMENSION >;
            builder->configure_builder( mesh );
        }
        return std::unique_ptr< LineMeshBuilder< DIMENSION > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > SurfaceMeshBuilder< DIMENSION >::create_builder(
        SurfaceMeshBase< DIMENSION >& mesh )
    {
        SurfaceMeshBuilder< DIMENSION >* builder = SurfaceMeshBuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_surface_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "SurfaceMeshBuilder",
                "Could not create mesh data structure: ", mesh.type_name() );
            Logger::warn( "SurfaceMeshBuilder",
                "Falling back to GeogramSurfaceMeshBuilder data structure" );

            builder = new GeogramSurfaceMeshBuilder< DIMENSION >;
            builder->configure_builder( mesh );
        }
        return std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > > VolumeMeshBuilder< DIMENSION >::create_builder(
        VolumeMesh< DIMENSION >& mesh )
    {
        VolumeMeshBuilder< DIMENSION >* builder =
            VolumeMeshBuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_volume_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "VolumeMeshBuilder",
                "Could not create mesh data structure: ", mesh.type_name() );
            Logger::warn( "VolumeMeshBuilder",
                "Falling back to GeogramVolumeMesshBuilder data structure" );

            builder = new GeogramVolumeMeshBuilder< DIMENSION >;
            builder->configure_builder( mesh );
        }
        return std::unique_ptr< VolumeMeshBuilder< DIMENSION > >( builder );
    }

    template std::unique_ptr< MeshBaseBuilder< 2 > > RINGMESH_API MeshBaseBuilder< 2 >::create_builder(
        MeshBase< 2 >& );
    template std::unique_ptr< PointSetMeshBuilder< 2 > > RINGMESH_API PointSetMeshBuilder<
        2 >::create_builder( PointSetMesh< 2 >& );
    template std::unique_ptr< LineMeshBuilder< 2 > > RINGMESH_API LineMeshBuilder< 2 >::create_builder(
        LineMesh< 2 >& );
    template std::unique_ptr< SurfaceMeshBuilder< 2 > > RINGMESH_API SurfaceMeshBuilder<
        2 >::create_builder( SurfaceMeshBase< 2 >& );

    template std::unique_ptr< MeshBaseBuilder< 3 > > RINGMESH_API MeshBaseBuilder< 3 >::create_builder(
        MeshBase< 3 >& );
    template std::unique_ptr< PointSetMeshBuilder< 3 > > RINGMESH_API PointSetMeshBuilder<
        3 >::create_builder( PointSetMesh< 3 >& );
    template std::unique_ptr< LineMeshBuilder< 3 > > RINGMESH_API LineMeshBuilder< 3 >::create_builder(
        LineMesh< 3 >& );
    template std::unique_ptr< SurfaceMeshBuilder< 3 > > RINGMESH_API SurfaceMeshBuilder<
        3 >::create_builder( SurfaceMeshBase< 3 >& );
    template std::unique_ptr< VolumeMeshBuilder< 3 > > RINGMESH_API VolumeMeshBuilder<
        3 >::create_builder( VolumeMesh< 3 >& );

    template class RINGMESH_API MeshBaseBuilder< 2 > ;
    template class RINGMESH_API PointSetMeshBuilder< 2 > ;
    template class RINGMESH_API LineMeshBuilder< 2 > ;
    template class RINGMESH_API SurfaceMeshBuilder< 2 > ;

    template class RINGMESH_API MeshBaseBuilder< 3 > ;
    template class RINGMESH_API PointSetMeshBuilder< 3 > ;
    template class RINGMESH_API LineMeshBuilder< 3 > ;
    template class RINGMESH_API SurfaceMeshBuilder< 3 > ;
    template class RINGMESH_API VolumeMeshBuilder< 3 > ;

}
