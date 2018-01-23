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

#include <ringmesh/tetrahedralize/common.h>

#ifdef RINGMESH_WITH_TETGEN

#include <memory>

#include <geogram/third_party/tetgen/tetgen.h>

/*!
 * @file Interface GEO::Mesh with Tetgen
 */
namespace GEO
{
    class Mesh;
} // namespace GEO

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshBuilder );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Tetgen wrapper
     * @author Jeanne Pellerin
     */
    class tetrahedralize_api TetgenMesher
    {
        ringmesh_disable_copy_and_move( TetgenMesher );

    public:
        TetgenMesher() = default;
        ~TetgenMesher();

        void tetrahedralize( const GEO::Mesh& input_mesh,
            VolumeMeshBuilder< 3 >& output_mesh_builder );

        void add_points_to_match_quality( double quality );

    private:
        void initialize();
        void initialize_tetgen_args();
        void tetrahedralize();

        void copy_mesh_to_tetgen_input( const GEO::Mesh& M );
        void copy_vertices_to_tetgen_input( const GEO::Mesh& M );
        void copy_edges_to_tetgen_input( const GEO::Mesh& M );
        void copy_polygons_to_tetgen_input( const GEO::Mesh& M );

        void set_regions( const std::vector< vec3 >& one_point_per_region );

        void assign_result_tetmesh_to_mesh(
            VolumeMeshBuilder< 3 >& output_mesh_builder ) const;
        std::vector< double > get_result_tetmesh_points() const;
        std::vector< index_t > get_result_tetmesh_tets() const;
        std::set< double > determine_tet_regions_to_keep() const;
        std::vector< index_t > determine_tets_to_keep() const;

    private:
        GEO_3rdParty::tetgenio tetgen_in_;
        GEO_3rdParty::tetgenio tetgen_out_;
        /*!
         * Command line options:
         * Q = Quiet (no output)
         * p = Use a piecewise linear complex
         * n = Save tetrahedron neighbors
         * Y = Do not add points on boundaries
         * AA = Save tetrahedron regions
         */
        std::string tetgen_command_line_ = std::string( "QpnYAA" );
        GEO_3rdParty::tetgenbehavior tetgen_args_;

        std::unique_ptr< GEO_3rdParty::tetgenio::polygon[] > polygons_{};
        std::unique_ptr< int[] > polygon_corners_{};
    };

    /*!
     * @brief Constrained tetrahedralize of the volumes defined by a
     * triangulated surface mesh
     * @details Does not require this mesh to be a closed manifold
     * as the equivalent in Geogram function does.
     */
    void tetrahedralize_api tetrahedralize_mesh_tetgen(
        VolumeMeshBuilder< 3 >& out_mesh_builder,
        const GEO::Mesh& in_mesh,
        bool refine,
        double quality );

} // namespace RINGMesh
#endif
