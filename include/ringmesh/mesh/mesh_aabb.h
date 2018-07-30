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

#include <ringmesh/basic/aabb.h>
#include <ringmesh/mesh/common.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMesh );
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class mesh_api LineAABBTree : public AABBTree< DIMENSION >
    {
    public:
        explicit LineAABBTree( const LineMesh< DIMENSION >& mesh );

        /*!
         * @brief Gets the closest edge to a given point using
         * Euclidean distance (L2 norm)
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest edge index.
         * - nearest_point the nearest point on the closest edge.
         * - distance the distance between \p query and \p nearest_point.
         */
        std::tuple< index_t, vecn< DIMENSION >, double > closest_edge(
            const vecn< DIMENSION >& query ) const;

        /*!
         * @brief Gets the closest edge to a given point using the
         * given distance \p EvalDistance
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest edge index.
         * - nearest_point the nearest point on the closest edge.
         * - distance the distance between \p query and \p nearest_point.
         */
        template < typename EvalDistance >
        std::tuple< index_t, vecn< DIMENSION >, double > closest_edge(
            const vecn< DIMENSION >& query, EvalDistance action ) const
        {
            return this->closest_element_box( query, action );
        }

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box, index_t element_id ) const override;

        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and an edge
         */
        class DistanceToEdge
        {
        public:
            explicit DistanceToEdge( const LineMesh< DIMENSION >& mesh )
                : mesh_( mesh )
            {
            }

            std::tuple< double, vecn< DIMENSION > > operator()(
                const vecn< DIMENSION >& query, index_t cur_box ) const;
            double operator()( const vecn< DIMENSION >& pt1,
                const vecn< DIMENSION >& pt2 ) const;

        private:
            const LineMesh< DIMENSION >& mesh_;
        };

    private:
        const LineMesh< DIMENSION >& mesh_;
    };

    ALIAS_2D_AND_3D( LineAABBTree );

    template < index_t DIMENSION >
    class mesh_api SurfaceAABBTree : public AABBTree< DIMENSION >
    {
    public:
        explicit SurfaceAABBTree( const SurfaceMeshBase< DIMENSION >& mesh );

        /*!
         * @brief Gets the closest triangle to a given point using
         * Euclidean distance (L2 norm)
         * @pre The mesh needs to be triangulated
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest triangle index.
         * - the nearest point on the closest triangle.
         * - the distance between \p query and \p nearest_point.
         */
        std::tuple< index_t, vecn< DIMENSION >, double > closest_triangle(
            const vecn< DIMENSION >& query ) const;

        /*!
         * @brief Gets the closest triangle to a given point using the
         * given distance \p EvalDistance
         * @pre The mesh needs to be triangulated
         * @param[in] query the point to use
         * @return a tuple containing:
         * - the closest triangle index.
         * - the nearest point on the closest triangle.
         * - the distance between \p query and \p nearest_point.
         */
        template < typename EvalDistance >
        std::tuple< index_t, vecn< DIMENSION >, double > closest_triangle(
            const vecn< DIMENSION >& query, EvalDistance action ) const
        {
            return this->closest_element_box( query, action );
        }

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box, index_t element_id ) const override;
        /*!
         * This class is used as functor in closest_element_box() to compute
         * the distance between a point and a triangle
         */
        class DistanceToTriangle
        {
        public:
            explicit DistanceToTriangle(
                const SurfaceMeshBase< DIMENSION >& mesh )
                : mesh_( mesh )
            {
            }

            std::tuple< double, vecn< DIMENSION > > operator()(
                const vecn< DIMENSION >& query, index_t cur_box ) const;
            double operator()( const vecn< DIMENSION >& pt1,
                const vecn< DIMENSION >& pt2 ) const;

        private:
            const SurfaceMeshBase< DIMENSION >& mesh_;
        };

    private:
        const SurfaceMeshBase< DIMENSION >& mesh_;
    };

    ALIAS_2D_AND_3D( SurfaceAABBTree );

    template < index_t DIMENSION >
    class mesh_api VolumeAABBTree : public AABBTree< DIMENSION >
    {
        ringmesh_template_assert_3d( DIMENSION );

    public:
        explicit VolumeAABBTree( const VolumeMesh< DIMENSION >& mesh );

        /*!
         * @brief Gets the cell contining a point
         * @param[in] query the point to use
         * @return the cell index containing \p query,
         * NO_ID if no cell is corresponding
         */
        index_t containing_cell( const vecn< DIMENSION >& query ) const;

    private:
        /*!
         * @brief Gets an element point from its box
         * @details In this case, the point is the first vertex of the element
         */
        vecn< DIMENSION > get_point_hint_from_box(
            const Box< DIMENSION >& box, index_t element_id ) const override;
        index_t containing_cell_recursive( const vecn< DIMENSION >& query,
            index_t node_index,
            index_t box_begin,
            index_t box_end ) const;

    private:
        const VolumeMesh< DIMENSION >& mesh_;
    };

    ALIAS_3D( VolumeAABBTree );
} // namespace RINGMesh
