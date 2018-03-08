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
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshBuilder );

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

    /*!
     * class for encapsulating volume mesh component
     */
    template < index_t DIMENSION >
    class VolumeMesh : public MeshBase< DIMENSION >
    {
        ringmesh_template_assert_3d( DIMENSION );
        friend class VolumeMeshBuilder< DIMENSION >;

    public:
        static std::unique_ptr< VolumeMesh< DIMENSION > > create_mesh(
            const MeshType type = "" );

        /*!
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_vertex(
            const ElementLocalVertex& cell_local_vertex ) const = 0;

        /*!
         * @brief Gets a vertex index by cell and local edge and local vertex
         * index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_edge_vertex(
            index_t cell_id, index_t edge_id, index_t vertex_id ) const = 0;

        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param
         * facet_id
         * and facet_id number of facet in th cell \param cell_id
         */
        virtual index_t cell_facet_vertex(
            const CellLocalFacet& cell_local_facet,
            index_t vertex_id ) const = 0;

        /*!
         * @brief Gets a facet index by cell and local facet index.
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @return the global facet index.
         */
        virtual index_t cell_facet(
            const CellLocalFacet& cell_local_facet ) const = 0;

        /*!
         * Computes the Mesh cell edge length
         * @param[in] cell_id the facet index
         * @param[in] edge_id the edge index
         * @return the cell edge length
         */
        double cell_edge_length( index_t cell_id, index_t edge_id ) const;

        /*!
         * Computes the Mesh cell edge barycenter
         * @param[in] cell_id the facet index
         * @param[in] edge_id the edge index
         * @return the cell edge center
         */
        vecn< DIMENSION > cell_edge_barycenter(
            index_t cell_id, index_t edge_id ) const;

        /*!
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_facets( index_t cell_id ) const = 0;
        /*!
         * @brief Gets the total number of facet in a all cells
         */
        virtual index_t nb_cell_facets() const = 0;

        /*!
         * @brief Gets the number of edges in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_edges( index_t cell_id ) const = 0;

        /*!
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_local_facet index of the cell and
         * the local index of the facet in the cell
         * @return the number of vertices in the facet \param facet_id in the
         * cell \param cell_id
         */
        virtual index_t nb_cell_facet_vertices(
            const CellLocalFacet& cell_local_facet ) const = 0;

        /*!
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        virtual index_t nb_cell_vertices( index_t cell_id ) const = 0;

        /*!
         * @brief Gets the number of cells in the Mesh.
         */
        virtual index_t nb_cells() const = 0;

        virtual index_t cell_begin( index_t cell_id ) const = 0;

        virtual index_t cell_end( index_t cell_id ) const = 0;

        /*!
         * @return the index of the adjacent cell of \param cell_local_facet
         */
        virtual index_t cell_adjacent(
            const CellLocalFacet& cell_local_facet ) const = 0;

        virtual GEO::AttributesManager& cell_attribute_manager() const = 0;

        virtual GEO::AttributesManager&
            cell_facet_attribute_manager() const = 0;

        /*!
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        virtual CellType cell_type( index_t cell_id ) const = 0;

        /*!
         * @brief Tests whether all the cells are tetrahedra.
         * When all the cells are tetrahedra, storage and access is optimized.
         * @return True if all cells are tetrahedra and False otherwise.
         */
        virtual bool cells_are_simplicies() const = 0;

        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] cell_local_facet the cell index and
         * the local facet index in the cell
         * @return the cell facet center
         */
        vecn< DIMENSION > cell_facet_barycenter(
            const CellLocalFacet& cell_local_facet ) const;

        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        vecn< DIMENSION > cell_barycenter( index_t cell_id ) const;

        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_local_facet the cell index and
         * the local facet index in the cell
         * @return the cell facet normal
         */
        vecn< DIMENSION > cell_facet_normal(
            const CellLocalFacet& cell_local_facet ) const;

        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        virtual double cell_volume( index_t cell_id ) const = 0;

        std::vector< index_t > cells_around_vertex(
            index_t vertex_id, index_t cell_hint ) const;

        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const;

        bool find_cell_from_colocated_vertex_within_distance_if_any(
            const vecn< DIMENSION >& vertex_vec,
            double distance,
            index_t& cell_id,
            index_t& cell_vertex_id ) const;

        /*!
         * @brief return the NNSearch at cell facets
         * @warning the NNSearch is destroyed when calling the
         * Mesh::facets_aabb()
         *  and Mesh::cells_aabb()
         */
        const NNSearch< DIMENSION >& cell_facet_nn_search() const;

        /*!
         * @brief return the NNSearch at cells
         */
        const NNSearch< DIMENSION >& cell_nn_search() const
        {
            if( !cell_nn_search_ )
            {
                std::vector< vecn< DIMENSION > > cell_centers( nb_cells() );
                for( auto c : range( nb_cells() ) )
                {
                    cell_centers[c] = cell_barycenter( c );
                }
                cell_nn_search_.reset(
                    new NNSearch< DIMENSION >( cell_centers, true ) );
            }
            return *cell_nn_search_.get();
        }
        /*!
         * @brief Creates an AABB tree for a Mesh cells
         */
        const VolumeAABBTree< DIMENSION >& cell_aabb() const
        {
            if( !cell_aabb_ )
            {
                cell_aabb_.reset( new VolumeAABBTree< DIMENSION >( *this ) );
            }
            return *cell_aabb_.get();
        }

        bool is_mesh_valid() const override;
        std::tuple< index_t, std::vector< index_t > >
            connected_components() const final;

    protected:
        VolumeMesh() = default;

    private:
        mutable std::unique_ptr< NNSearch< DIMENSION > >
            cell_facet_nn_search_{};
        mutable std::unique_ptr< NNSearch< DIMENSION > > cell_nn_search_{};
        mutable std::unique_ptr< VolumeAABBTree< DIMENSION > > cell_aabb_{};

        void store_cells_around_vertex( GEO::index_t cell_hint,
            GEO::index_t vertex_id,
            std::vector< index_t >& result ) const;
    };

    using VolumeMesh3D = VolumeMesh< 3 >;

    template < index_t DIMENSION >
    using VolumeMeshFactory = Factory< MeshType, VolumeMesh< DIMENSION > >;
    using VolumeMeshFactory3D = VolumeMeshFactory< 3 >;
} // namespace RINGMesh
