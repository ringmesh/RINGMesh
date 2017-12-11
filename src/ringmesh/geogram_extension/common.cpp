/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/geogram_extension/common.h>

#include <ringmesh/basic/plugin_manager.h>

#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/geogram_extension/geogram_mesh.h>
#include <ringmesh/geogram_extension/geogram_mesh_builder.h>

namespace
{
    void register_geogram_mesh()
    {
        RINGMesh::PointSetMeshFactory2D::register_creator< RINGMesh::GeogramPointSetMesh2D >(
            RINGMesh::GeogramPointSetMesh2D::type_name_static() );
        RINGMesh::LineMeshFactory2D::register_creator< RINGMesh::GeogramLineMesh2D >(
            RINGMesh::GeogramLineMesh2D::type_name_static() );
        RINGMesh::SurfaceMeshFactory2D::register_creator< RINGMesh::GeogramSurfaceMesh2D >(
            RINGMesh::GeogramSurfaceMesh2D::type_name_static() );

        RINGMesh::PointSetMeshFactory3D::register_creator< RINGMesh::GeogramPointSetMesh3D >(
            RINGMesh::GeogramPointSetMesh3D::type_name_static() );
        RINGMesh::LineMeshFactory3D::register_creator< RINGMesh::GeogramLineMesh3D >(
            RINGMesh::GeogramLineMesh3D::type_name_static() );
        RINGMesh::SurfaceMeshFactory3D::register_creator< RINGMesh::GeogramSurfaceMesh3D >(
            RINGMesh::GeogramSurfaceMesh3D::type_name_static() );
        RINGMesh::VolumeMeshFactory3D::register_creator< RINGMesh::GeogramVolumeMesh3D >(
            RINGMesh::GeogramVolumeMesh3D::type_name_static() );

        RINGMesh::PointSetMeshBuilderFactory2D::
            register_creator< RINGMesh::GeogramPointSetMeshBuilder2D >(
                RINGMesh::GeogramPointSetMesh2D::type_name_static() );
        RINGMesh::LineMeshBuilderFactory2D::register_creator< RINGMesh::GeogramLineMeshBuilder2D >(
            RINGMesh::GeogramLineMesh2D::type_name_static() );
        RINGMesh::SurfaceMeshBuilderFactory2D::
            register_creator< RINGMesh::GeogramSurfaceMeshBuilder2D >(
                RINGMesh::GeogramSurfaceMesh2D::type_name_static() );

        RINGMesh::PointSetMeshBuilderFactory3D::
            register_creator< RINGMesh::GeogramPointSetMeshBuilder3D >(
                RINGMesh::GeogramPointSetMesh3D::type_name_static() );
        RINGMesh::LineMeshBuilderFactory3D::register_creator< RINGMesh::GeogramLineMeshBuilder3D >(
            RINGMesh::GeogramLineMesh3D::type_name_static() );
        RINGMesh::SurfaceMeshBuilderFactory3D::
            register_creator< RINGMesh::GeogramSurfaceMeshBuilder3D >(
                RINGMesh::GeogramSurfaceMesh3D::type_name_static() );
        RINGMesh::VolumeMeshBuilderFactory3D::
            register_creator< RINGMesh::GeogramVolumeMeshBuilder3D >(
                RINGMesh::GeogramVolumeMesh3D::type_name_static() );
    }

    RINGMESH_PLUGIN_INITIALIZE(
        RINGMesh_geogram_extension,
        // Plugin initialization
        RINGMesh::ringmesh_geogram_mesh_io_initialize();
        register_geogram_mesh();
    );
} // namespace
