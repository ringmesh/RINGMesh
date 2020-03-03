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

//#include <ringmesh/mesh/mesh.h>

#include <numeric>
#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <stack>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    std::unique_ptr< PointSetMesh< DIMENSION > >
        PointSetMesh< DIMENSION >::create_mesh( const MeshType type )
    {
        auto new_type = type;
        if( new_type.empty() )
        {
            if( !PointSetMeshFactory< DIMENSION >::has_creator(
                    "GeogramPointSetMesh" ) )
            {
                throw RINGMeshException( "PointSetMesh",
                    "Default mesh data structure not registered" );
            }
            return create_mesh( "GeogramPointSetMesh" );
        }
        auto mesh = PointSetMeshFactory< DIMENSION >::create( new_type );
        if( !mesh )
        {
            Logger::warn( "PointSetMesh",
                "Could not create mesh data structure: ", new_type );
            Logger::warn( "PointSetMesh",
                "Falling back to GeogramPointSetMesh data structure" );

            return create_mesh();
        }
        return mesh;
    }

    template < index_t DIMENSION >
    std::tuple< index_t, std::vector< index_t > >
        PointSetMesh< DIMENSION >::connected_components() const
    {
        const auto nb_compoments = this->nb_vertices();
        std::vector< index_t > components( nb_compoments );
        std::iota( components.begin(), components.end(), 0 );
        return std::make_tuple( nb_compoments, components );
    }

    template class mesh_api PointSetMesh< 2 >;
    template class mesh_api PointSetMesh< 3 >;
} // namespace RINGMesh
