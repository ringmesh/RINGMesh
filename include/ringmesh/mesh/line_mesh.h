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

#include <memory>
#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/nn_search.h>
#include <ringmesh/mesh/common.h>

#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_base.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMeshBuilder );

    struct ElementLocalVertex;
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
     * class for encapsulating line mesh (composed of edges)
     */
    template < index_t DIMENSION >
    class LineMesh : public MeshBase< DIMENSION >
    {
        friend class LineMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< LineMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );

        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_local_vertex index of the edge and of the
         * local index of the vertex, in {0,1}
         * @return the global index of vertex in \p edge_local_vertex.
         */
        virtual index_t edge_vertex(
            const ElementLocalVertex& edge_local_vertex ) const = 0;

        /*!
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        virtual index_t nb_edges() const = 0;

        /*!
         * @brief Gets the length of the edge \param edge_id
         */
        double edge_length( index_t edge_id ) const;

        vecn< DIMENSION > edge_barycenter( index_t edge_id ) const;

        /*!
         * @brief return the NNSearch at edges
         * @warning the NNSearch is destroyed when calling the
         * Mesh::polygons_aabb()
         * and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& edge_nn_search() const;

        /*!
         * @brief Creates an AABB tree for a Mesh edges
         */
        const LineAABBTree< DIMENSION >& edge_aabb() const;

        virtual GEO::AttributesManager& edge_attribute_manager() const = 0;

        bool is_mesh_valid() const override;

        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;

    protected:
        LineMesh() = default;

    private:
        mutable std::unique_ptr< NNSearch< DIMENSION > > edge_nn_search_{};
        mutable std::unique_ptr< LineAABBTree< DIMENSION > > edge_aabb_{};
    };
    ALIAS_2D_AND_3D( LineMesh );

    template < index_t DIMENSION >
    using LineMeshFactory = Factory< MeshType, LineMesh< DIMENSION > >;
    ALIAS_2D_AND_3D( LineMeshFactory );
} // namespace RINGMesh
