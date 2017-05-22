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

    PointMeshBuilder* create_builder_0d( MeshBase& mesh )
    {
        PointMeshBuilder* builder = Mesh0DBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< PointMesh& >( mesh ) );
        }
        return builder;
    }

    LineMeshBuilder* create_builder_1d( MeshBase& mesh )
    {
        LineMeshBuilder* builder = LineMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< LineMesh& >( mesh ) );
        }
        return builder;
    }

    SurfaceMeshBuilder* create_builder_2d( MeshBase& mesh )
    {
        SurfaceMeshBuilder* builder = Mesh2DBuilderFactory::create_object(
            mesh.type_name() );
        if( builder ) {
            builder->set_mesh( dynamic_cast< SurfaceMesh& >( mesh ) );
        }
        return builder;
    }

    VolumeMeshBuilder* create_builder_3d( MeshBase& mesh )
    {
        VolumeMeshBuilder* builder = Mesh3DBuilderFactory::create_object(
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
        PointMeshBuilder* builder0d = create_builder_0d( mesh );
        if( builder0d ) {
            builder0d->set_mesh( dynamic_cast< PointMesh& >( mesh ) );
            builder = builder0d;
        } else {
            LineMeshBuilder* builder1d = create_builder_1d( mesh );
            if( builder1d ) {
                builder1d->set_mesh( dynamic_cast< LineMesh& >( mesh ) );
                builder = builder1d;
            } else {
                SurfaceMeshBuilder* builder2d = create_builder_2d( mesh );
                if( builder2d ) {
                    builder2d->set_mesh( dynamic_cast< SurfaceMesh& >( mesh ) );
                    builder = builder2d;
                } else {
                    VolumeMeshBuilder* builder3d = create_builder_3d( mesh );
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

    std::unique_ptr< PointMeshBuilder > PointMeshBuilder::create_builder( PointMesh& mesh )
    {
        PointMeshBuilder* builder = Mesh0DBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_builder_0d( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh0DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "Mesh0DBuilder",
                "Falling back to GeogramMesh0DBuilder data structure" );

            builder = new GeogramMesh0DBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< PointMeshBuilder >( builder );
    }

    std::unique_ptr< LineMeshBuilder > LineMeshBuilder::create_builder( LineMesh& mesh )
    {
        LineMeshBuilder* builder = LineMeshBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_builder_1d( mesh );
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
        SurfaceMeshBuilder* builder = Mesh2DBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_builder_2d( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh2DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "Mesh2DBuilder",
                "Falling back to GeogramMesh2DBuilder data structure" );

            builder = new GeogramMesh2DBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< SurfaceMeshBuilder >( builder );
    }

    std::unique_ptr< VolumeMeshBuilder > VolumeMeshBuilder::create_builder( VolumeMesh& mesh )
    {
        VolumeMeshBuilder* builder = Mesh3DBuilderFactory::create_object(
            mesh.type_name() );
        if( !builder ) {
            builder = create_builder_3d( mesh );
        }
        if( !builder ) {
            Logger::warn( "Mesh3DBuilder", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "Mesh3DBuilder",
                "Falling back to GeogramMesh3DBuilder data structure" );

            builder = new GeogramMesh3DBuilder;
        }
        builder->set_mesh( mesh );
        return std::unique_ptr< VolumeMeshBuilder >( builder );
    }

} // namespace
