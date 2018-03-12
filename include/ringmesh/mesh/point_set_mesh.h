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

#pragma once

#include <algorithm>
#include <memory>
#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/nn_search.h>
#include <ringmesh/mesh/common.h>

#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_base.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMeshBuilder );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * class base class for encapsulating Mesh structure
     * @brief encapsulate adimensional mesh functionalities in order to provide
     * an API
     * on which we base the RINGMesh algorithms
     * @note For now, we encapsulate the GEO::Mesh class.
     */
    /*!
     * class for encapsulating mesh composed of points
     */
    template < index_t DIMENSION >
    class PointSetMesh : public MeshBase< DIMENSION >
    {
        friend class PointSetMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< PointSetMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );
        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;
        bool is_mesh_valid() const override
        {
            return true;
        }

    protected:
        PointSetMesh() = default;
    };
    ALIAS_2D_AND_3D( PointSetMesh );

    template < index_t DIMENSION >
    using PointSetMeshFactory = Factory< MeshType, PointSetMesh< DIMENSION > >;
    ALIAS_2D_AND_3D( PointSetMeshFactory );
} // namespace RINGMesh
