/*
* Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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


#ifndef __RINGMESH_TETGEN_MESHER__
#define __RINGMESH_TETGEN_MESHER__

#include <ringmesh/common.h>
#include <string> 
#include <vector>

#include <geogram/basic/memory.h>
#include <geogram/mesh/mesh.h>


#ifdef RINGMESH_WITH_TETGEN

#include <geogram/third_party/tetgen/tetgen.h>

/*!
 * @file Interface GEO::Mesh with Tetgen 
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
    * @brief Utility class to set Tetgen switches and check their consistency
    * @details Tetgen arguments are a mess and this class helps set the basic options
    * @todo To implement!
    *
    * Q: quiet
    * p: input data is surfacic
    * q: desired quality
    * O0: do not optimize mesh at all -> a lot of flat tets
    * V: verbose - A LOT of information
    * Y: prohibit steiner points on boundaries
    * A: generate region tags for each shell.
    *
    * Meshing with incomplete quality value "Qpq%fYA"
    */
    class TetgenCommandLine {
    public:
        const std::string command_line() const
        {
            return command_line_ ;
        }

    private:
        std::string command_line_ ;
    } ;

    /*!
    * @brief Tetgen wrapper
    * @author  Jeanne Pellerin
    */
    class TetgenMesher {
        ringmesh_disable_copy( TetgenMesher ) ;
    public:
        TetgenMesher()
            : polygons_( nil ), polygon_corners_( nil )
        {}
        ~TetgenMesher() ;

        void tetrahedralize(
            const GEO::Mesh& input_mesh,
            const std::string& command_line,
            GEO::Mesh& output_mesh ) ;

        void tetrahedralize(
            const GEO::Mesh& input_mesh,
            const std::vector< vec3 >& one_point_per_region,
            const std::string& command_line,
            GEO::Mesh& output_mesh ) ;

    private:
        void initialize() ;
        void initialize_tetgen_args() ;
        void set_command_line( const std::string& command_line ) ;
        void tetrahedralize() ;

        void copy_mesh_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_vertices_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_edges_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_facets_to_tetgen_input( const GEO::Mesh& M ) ;
        void set_regions( const std::vector< vec3 >& one_point_per_region ) ;

        void fill_region_attribute_on_mesh_cells(
            GEO::Mesh& M,
            const std::string& attribute_name ) const ;
        void assign_result_tetmesh_to_mesh( GEO::Mesh& M ) const ;
        void get_result_tetmesh_points( GEO::vector< double >& points ) const ;
        void get_result_tetmesh_tets( GEO::vector< index_t >& tets ) const ;

    private:
        GEO_3rdParty::tetgenio tetgen_in_ ;
        GEO_3rdParty::tetgenio tetgen_out_ ;
        std::string tetgen_command_line_ ;
        GEO_3rdParty::tetgenbehavior tetgen_args_ ;

        GEO_3rdParty::tetgenio::polygon* polygons_ ;
        int* polygon_corners_ ;
    } ;

    /*!
    * @brief Constrained tetrahedralize of the volumes defined by a triangulated surface mesh
    * @details Does not require this mesh to be a closed manifold
    * as the equivalent in Geogram function does.
    */
    void RINGMESH_API tetrahedralize_mesh_tetgen(
        GEO::Mesh& M,
        bool refine,
        double quality ) ;

}
#endif
#endif