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

    PointSetMeshBuilder* create_point_mesh_builder( MeshBase& mesh )
    {
        PointSetMeshBuilder* builder = PointSetMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< PointSetMesh& >( mesh ) );
        }
        return builder;
    }

    LineMeshBuilder* create_line_mesh_builder( MeshBase& mesh )
    {
        LineMeshBuilder* builder = LineMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< LineMesh& >( mesh ) );
        }
        return builder;
    }

    SurfaceMeshBuilder* create_surface_mesh_builder( MeshBase& mesh )
    {
        SurfaceMeshBuilder* builder = SurfaceMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< SurfaceMesh& >( mesh ) );
        }
        return builder;
    }

    VolumeMeshBuilder* create_volume_mesh_builder( MeshBase& mesh )
    {
        VolumeMeshBuilder* builder = VolumeMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< VolumeMesh& >( mesh ) );
        }
        return builder;
    }
}

namespace RINGMesh {

    std::unique_ptr< MeshBaseBuilder > MeshBaseBuilder::create_builder(
        MeshBase& mesh )
    {
        MeshBaseBuilder* builder = nullptr;
        PointSetMeshBuilder* builder0d = create_point_mesh_builder( mesh );
        if( builder0d ) {
            builder0d->set_mesh( dynamic_cast< PointSetMesh& >( mesh ) );
            builder = builder0d;
        } else {
            LineMeshBuilder* builder1d = create_line_mesh_builder( mesh );
            if( builder1d ) {
                builder1d->set_mesh( dynamic_cast< LineMesh& >( mesh ) );
                builder = builder1d;
            } else {
                SurfaceMeshBuilder* builder2d = create_surface_mesh_builder( mesh );
                if( builder2d ) {
                    builder2d->set_mesh( dynamic_cast< SurfaceMesh& >( mesh ) );
                    builder = builder2d;
                } else {
                    VolumeMeshBuilder* builder3d = create_volume_mesh_builder( mesh );
                    if( builder3d ) {
                        builder3d->set_mesh( dynamic_cast< VolumeMesh& >( mesh ) );
                        builder = builder3d;
                    }
                }
            }
        }
        if( !builder ) {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh data structure: " + mesh.type_name() );
        }
        return std::unique_ptr< MeshBaseBuilder >( builder );
    }

    std::unique_ptr< PointSetMeshBuilder > PointSetMeshBuilder::create_builder( PointSetMesh& mesh )
    {
        PointSetMeshBuilder* builder = PointSetMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_point_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh0DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "PointMeshBuilder",
                "Falling back to GeogramPointSetMeshBuilder data structure" );

            builder = new GeogramPointSetMeshBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< PointSetMeshBuilder >( builder );
    }

    std::unique_ptr< LineMeshBuilder > LineMeshBuilder::create_builder( LineMesh& mesh )
    {
        LineMeshBuilder* builder = LineMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_line_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "LineMeshBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "LineMeshBuilder",
                "Falling back to GeogramLineMeshBuilder data structure" );

            builder = new GeogramLineMeshBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< LineMeshBuilder >( builder );
    }

    std::unique_ptr< SurfaceMeshBuilder > SurfaceMeshBuilder::create_builder( SurfaceMesh& mesh )
    {
        SurfaceMeshBuilder* builder = SurfaceMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_surface_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "SurfaceMeshBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "SurfaceMeshBuilder",
                "Falling back to GeogramSurfaceMeshBuilder data structure" );

            builder = new GeogramSurfaceMeshBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< SurfaceMeshBuilder >( builder );
    }

    std::unique_ptr< VolumeMeshBuilder > VolumeMeshBuilder::create_builder( VolumeMesh& mesh )
    {
        VolumeMeshBuilder* builder = VolumeMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_volume_mesh_builder( mesh );
        }
        if( !builder ) {
            Logger::warn( "VolumeMeshBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "VolumeMeshBuilder",
                "Falling back to GeogramVolumeMesshBuilder data structure" );

            builder = new GeogramVolumeMeshBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< VolumeMeshBuilder >( builder );
    }

} // namespace
