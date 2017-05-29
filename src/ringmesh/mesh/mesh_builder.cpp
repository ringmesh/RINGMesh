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
        PointSetMeshBuilder< DIMENSION >* builder =
            PointMesh2BuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< PointSetMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    LineMeshBuilder< DIMENSION >* create_line_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        LineMeshBuilder< DIMENSION >* builder =
            LineMesh2BuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< LineMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    SurfaceMeshBuilder< DIMENSION >* create_surface_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        SurfaceMeshBuilder< DIMENSION >* builder = SurfaceMesh2BuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< SurfaceMeshBase< DIMENSION >& >( mesh ) );
        }
        return builder;
    }

    template< index_t DIMENSION >
    VolumeMeshBuilder< DIMENSION >* create_volume_mesh_builder(
        MeshBase< DIMENSION >& mesh )
    {
        VolumeMeshBuilder< DIMENSION >* builder = VolumeMesh2BuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< VolumeMesh< DIMENSION >& >( mesh ) );
        }
        return builder;
    }
}

namespace RINGMesh {

    template< index_t DIMENSION >
    std::unique_ptr< MeshBaseBuilder< DIMENSION > > MeshBaseBuilder< DIMENSION >::create_builder(
        MeshBase< DIMENSION >& mesh )
    {
        MeshBaseBuilder< DIMENSION >* builder = nullptr;
        PointSetMeshBuilder< DIMENSION >* builder0d = create_point_mesh_builder(
            mesh );
        if( builder0d ) {
            builder0d->set_mesh( dynamic_cast< PointSetMesh< DIMENSION >& >( mesh ) );
            builder = builder0d;
        } else {
            LineMeshBuilder< DIMENSION >* builder1d = create_line_mesh_builder(
                mesh );
            if( builder1d ) {
                builder1d->set_mesh(
                    dynamic_cast< LineMesh< DIMENSION >& >( mesh ) );
                builder = builder1d;
            } else {
                SurfaceMeshBuilder< DIMENSION >* builder2d =
                    create_surface_mesh_builder( mesh );
                if( builder2d ) {
                    builder2d->set_mesh(
                        dynamic_cast< SurfaceMeshBase< DIMENSION >& >( mesh ) );
                    builder = builder2d;
                } else {
                    VolumeMeshBuilder< DIMENSION >* builder3d =
                        create_volume_mesh_builder( mesh );
                    if( builder3d ) {
                        builder3d->set_mesh(
                            dynamic_cast< VolumeMesh< DIMENSION >& >( mesh ) );
                        builder = builder3d;
                    }
                }
            }
        }
        if( !builder ) {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh data structure: " + mesh.type_name() );
        }
        return std::unique_ptr< MeshBaseBuilder< DIMENSION > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< PointSetMeshBuilder< DIMENSION > > PointSetMeshBuilder< DIMENSION >::create_builder(
        PointSetMesh< DIMENSION >& mesh )
    {
        PointSetMeshBuilder< DIMENSION > *builder =
            PointMesh2BuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_point_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh0DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "PointMeshBuilder",
                "Falling back to GeogramPointSetMeshBuilder data structure" );

            builder = new GeogramPointSetMeshBuilder< DIMENSION >;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< PointSetMeshBuilder >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< LineMeshBuilder< DIMENSION > > LineMeshBuilder< DIMENSION >::create_builder(
        LineMesh< DIMENSION >& mesh )
    {
        LineMeshBuilder< DIMENSION >* builder =
            LineMesh2BuilderFactory< DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_line_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "LineMeshBuilder",
                "Could not create mesh data structure: ", mesh.type_name() );
            Logger::warn( "LineMeshBuilder",
                "Falling back to GeogramLineMeshBuilder data structure" );

            builder = new GeogramLineMeshBuilder< DIMENSION >;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< LineMeshBuilder< DIMENSION > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > SurfaceMeshBuilder<
        DIMENSION >::create_builder( SurfaceMeshBase< DIMENSION >& mesh )
    {
        SurfaceMeshBuilder< DIMENSION >* builder = SurfaceMesh2BuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_surface_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "SurfaceMeshBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "SurfaceMeshBuilder",
                "Falling back to GeogramSurfaceMeshBuilder data structure" );

            builder = new GeogramSurfaceMeshBuilder< DIMENSION >;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >( builder );
    }

    template< index_t DIMENSION >
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > > VolumeMeshBuilder< DIMENSION >::create_builder(
        VolumeMesh< DIMENSION >& mesh )
    {
        VolumeMeshBuilder< DIMENSION >* builder = VolumeMesh2BuilderFactory<
            DIMENSION >::create_object( mesh.type_name() );
        if( !builder ) {
            builder = create_volume_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "VolumeMeshBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "VolumeMeshBuilder",
                "Falling back to GeogramVolumeMesshBuilder data structure" );

            builder = new GeogramVolumeMeshBuilder< DIMENSION >;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< VolumeMeshBuilder< DIMENSION > >( builder );
    }

    template class PointSetMeshBuilder< 2 >;
    template class LineMeshBuilder< 2 >;
    template class SurfaceMeshBuilder< 2 >;

    template class MeshBaseBuilder< 3 >;
    template class PointSetMeshBuilder< 3 >;
    template class LineMeshBuilder< 3 >;
    template class SurfaceMeshBuilder< 3 >;
    template class VolumeMeshBuilder< 3 >;



} // namespace
