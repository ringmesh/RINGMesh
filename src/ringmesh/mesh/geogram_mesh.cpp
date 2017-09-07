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

#include <ringmesh/mesh/geogram_mesh.h>

#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace RINGMesh {
void register_geogram_mesh() {
  PointSetMeshFactory2D::register_creator<GeogramPointSetMesh2D>(
      GeogramPointSetMesh2D::type_name_static());
  LineMeshFactory2D::register_creator<GeogramLineMesh2D>(
      GeogramLineMesh2D::type_name_static());
  SurfaceMeshFactory2D::register_creator<GeogramSurfaceMesh2D>(
      GeogramSurfaceMesh2D::type_name_static());

  PointSetMeshFactory3D::register_creator<GeogramPointSetMesh3D>(
      GeogramPointSetMesh3D::type_name_static());
  LineMeshFactory3D::register_creator<GeogramLineMesh3D>(
      GeogramLineMesh3D::type_name_static());
  SurfaceMeshFactory3D::register_creator<GeogramSurfaceMesh3D>(
      GeogramSurfaceMesh3D::type_name_static());
  VolumeMeshFactory3D::register_creator<GeogramVolumeMesh3D>(
      GeogramVolumeMesh3D::type_name_static());

  PointSetMeshBuilderFactory2D::register_creator<GeogramPointSetMeshBuilder2D>(
      GeogramPointSetMesh2D::type_name_static());
  LineMeshBuilderFactory2D::register_creator<GeogramLineMeshBuilder2D>(
      GeogramLineMesh2D::type_name_static());
  SurfaceMeshBuilderFactory2D::register_creator<GeogramSurfaceMeshBuilder2D>(
      GeogramSurfaceMesh2D::type_name_static());

  PointSetMeshBuilderFactory3D::register_creator<GeogramPointSetMeshBuilder3D>(
      GeogramPointSetMesh3D::type_name_static());
  LineMeshBuilderFactory3D::register_creator<GeogramLineMeshBuilder3D>(
      GeogramLineMesh3D::type_name_static());
  SurfaceMeshBuilderFactory3D::register_creator<GeogramSurfaceMeshBuilder3D>(
      GeogramSurfaceMesh3D::type_name_static());
  VolumeMeshBuilderFactory3D::register_creator<GeogramVolumeMeshBuilder3D>(
      GeogramVolumeMesh3D::type_name_static());
}
}  // namespace RINGMesh
