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

#include <ringmesh/mesh/common.h>

#include <algorithm>
#include <memory>

#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/nn_search.h>

#include <ringmesh/mesh/mesh_aabb.h>

namespace GEO
{
    class AttributesManager;
} // namespace GEO

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBaseBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );

    struct EdgeLocalVertex;
    struct ElementLocalVertex;
    struct PolygonLocalEdge;
    struct CellLocalFacet;
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
// BEGINING OF MESHBASE DEFINITION
    template < index_t DIMENSION >
    class MeshBase
    {
        ringmesh_disable_copy_and_move( MeshBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class MeshBaseBuilder< DIMENSION >;

    public:
        virtual ~MeshBase() = default;

        virtual void save_mesh( const std::string& filename ) const = 0;

        virtual std::tuple< index_t, std::vector< index_t > >
            connected_components() const = 0;

        /*!
         * \name Vertex methods
         * @{
         */
        /*!
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return const reference to the point that corresponds to the vertex.
         */
        virtual const vecn< DIMENSION >& vertex( index_t v_id ) const = 0;
        /*
         * @brief Gets the number of vertices in the Mesh.
         */
        virtual index_t nb_vertices() const = 0;

        virtual GEO::AttributesManager& vertex_attribute_manager() const = 0;

        /*!
         * @brief return the NNSearch at vertices
         * @warning the NNSearch is destroyed when calling the
         * Mesh::polygons_aabb()
         * and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& vertex_nn_search() const;

        virtual MeshType type_name() const = 0;

        virtual std::string default_extension() const = 0;

        virtual bool is_mesh_valid() const = 0;

        /*!
         * @}
         */
    protected:
        MeshBase() = default;

    protected:
        mutable std::unique_ptr< NNSearch< DIMENSION > > vertex_nn_search_{};
    };
    ALIAS_2D_AND_3D( MeshBase );
// END OF MESHBASE DEFINITION
} // namespace RINGMesh
