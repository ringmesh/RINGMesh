/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr 
 *     Antoine.Mazuyer@univ-lorraine.fr 
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/

#ifndef __RINGMESH_TETRA_GEN__
#define __RINGMESH_TETRA_GEN__

#include <ringmesh/common.h>
#include <ringmesh/utils.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>

#include <vector>

#ifdef USE_MG_TETRA
    extern "C" {
        #include <meshgems/meshgems.h>
        #include <meshgems/tetra.h>
    }
#endif

namespace RINGMesh {
    class BoundaryModelElement ;
    class TetraGen ;
    class WellGroup ;
}

namespace RINGMesh {
    typedef GEO::SmartPointer< TetraGen > TetraGen_var ;

    static const std::vector< vec3 > vector_vec3 ;
    static const std::vector< std::vector< Edge > > vector_edge ;
    class RINGMESH_API TetraGen : public GEO::Counted {
    public:
        virtual ~TetraGen() ;
        static TetraGen_var instantiate(
            const TetraMethod& method,
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points = true,
            const std::vector< vec3 >& internal_vertices = vector_vec3,
            const WellGroup* wells = nil ) ;

        virtual bool tetrahedralize() = 0 ;

        index_t nb_points() const { return internal_vertices_ptr_ ; }
        index_t nb_internal_points() const { return nb_total_points() - internal_vertices_ptr_ ; }
        index_t nb_total_points() const { return tetmesh_.vertices.nb() ; }

    protected:
        TetraGen(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool refine,
            const std::vector< vec3 >& internal_vertices,
            const WellGroup* wells ) ;

        void initialize_storage( index_t nb_points, index_t nb_tets ) ;
        void set_point( index_t index, const double* point ) ;
        void set_tetra( index_t index, int* tet, index_t nb_lines, index_t nb_triangles ) ;

    protected:
        GEO::Mesh& tetmesh_ ;
        const BoundaryModelElement* region_ ;
        const WellGroup* wells_ ;
        index_t internal_vertices_ptr_ ;
        bool refine_ ;
    } ;


    class RINGMESH_API TetraGen_TetGen: public TetraGen {
    public:
        TetraGen_TetGen(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const WellGroup* wells ) ;
        virtual ~TetraGen_TetGen() {} ;

        virtual bool tetrahedralize() ;
    } ;

#ifdef USE_MG_TETRA
    class RINGMESH_API TetraGen_MG_Tetra: public TetraGen {
    public:
        TetraGen_MG_Tetra(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const WellGroup* wells ) ;
        virtual ~TetraGen_MG_Tetra() ;

        virtual bool tetrahedralize() ;

        static status_t my_message_cb( message_t * msg, void *user_data ) ;

    private:
        context_t* context_ ;
        mesh_t* mesh_input_ ;
        mesh_t* mesh_output_ ;
        tetra_session_t* tms_ ;
    } ;
#endif

}

#endif
