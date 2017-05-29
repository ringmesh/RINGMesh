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

#include <ringmesh/mesh/geogram_mesh.h>

#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace RINGMesh {

    void register_geogram_mesh()
    {
        ringmesh_register_point_mesh_3d( GeogramPointMesh3D );
        ringmesh_register_point_mesh_builder( GeogramPointMesh3D );
        ringmesh_register_line_mesh_3d( GeogramLineMesh3D );
        ringmesh_register_line_mesh_builder( GeogramLineMesh3D );
        ringmesh_register_surface_mesh_3d( GeogramSurfaceMesh3D );
        ringmesh_register_surface_mesh_builder( GeogramSurfaceMesh3D );
        ringmesh_register_volume_mesh_3d( GeogramVolumeMesh3D );
        ringmesh_register_volume_mesh_builder( GeogramVolumeMesh3D );
    }

    template class GeogramPointMesh< 2 >;
    template class GeogramLineMesh< 2 >;
    template class GeogramSurfaceMesh< 2 >;

    template class GeogramPointMesh< 3 >;
    template class GeogramLineMesh< 3 >;
    template class GeogramSurfaceMesh< 3 >;
    template class GeogramVolumeMesh< 3 >;
}

