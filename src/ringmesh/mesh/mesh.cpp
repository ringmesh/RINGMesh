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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh.h>

#include <numeric>
#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <stack>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/mesh_base.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    MeshSetBase< DIMENSION >::MeshSetBase()
    {
        create_point_set_mesh( "" );
        create_line_mesh( "" );
        create_well_mesh( "" );
        create_surface_mesh( "" );
    }

    template < index_t DIMENSION >
    void MeshSetBase< DIMENSION >::create_point_set_mesh( MeshType type )
    {
        point_set_mesh = PointSetMesh< DIMENSION >::create_mesh( type );
    }

    template < index_t DIMENSION >
    void MeshSetBase< DIMENSION >::create_line_mesh( MeshType type )
    {
        line_mesh = LineMesh< DIMENSION >::create_mesh( type );
    }

    template < index_t DIMENSION >
    void MeshSetBase< DIMENSION >::create_well_mesh( MeshType type )
    {
        well_mesh = LineMesh< DIMENSION >::create_mesh( type );
    }

    template < index_t DIMENSION >
    void MeshSetBase< DIMENSION >::create_surface_mesh( MeshType type )
    {
        surface_mesh = SurfaceMesh< DIMENSION >::create_mesh( type );
    }

    MeshSet< 3 >::MeshSet()
    {
        create_volume_mesh( "" );
    }

    void MeshSet< 3 >::create_volume_mesh( MeshType type )
    {
        volume_mesh = VolumeMesh3D::create_mesh( type );
    }

    template class mesh_api MeshSetBase< 2 >;
    template class mesh_api MeshSet< 2 >;

    template class mesh_api MeshSetBase< 3 >;
} // namespace RINGMesh
