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
 *     Ecole Nationale Sup�rieure de G�ologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_MACRO_MESH__
#define __RINGMESH_MACRO_MESH__

#include <ringmesh/common.h>
#include <ringmesh/utils.h>

#include <geogram/mesh/mesh.h>

namespace GEO {
    class MeshTetsAABB ;
    class MeshFacetsAABB ;
}

namespace RINGMesh {
    class BoundaryModel ;
    class MacroMesh ;

    static std::vector< std::vector< vec3 > > empty_vertices ;

    class RINGMESH_API MacroMeshVertices {
    public:
        MacroMeshVertices()
              : initialized_( false )
        {
        }

        index_t nb_vertices( const MacroMesh& mm ) const ;

        index_t nb_vertex_indices( const MacroMesh& mm ) const ;

        index_t global_vertex_id(
            const MacroMesh& mm,
            index_t mesh,
            index_t v ) const ;

        const vec3& global_vertex(
            const MacroMesh& mm,
            index_t global_v ) const ;

    private:
        void initialize( const MacroMesh& mm ) ;

    private:
        bool initialized_ ;

        std::vector< vec3 > unique_vertices_ ;
        std::vector< index_t > global_vertex_indices_ ;
        std::vector< index_t > vertex2mesh_ ;
    } ;

    class RINGMESH_API MacroMeshFacets {
    public:
        MacroMeshFacets()
              : initialized_( false )
        {
        }

        index_t surface_facet(
            const MacroMesh& mm,
            index_t s,
            index_t f ) const ;

        index_t surface_mesh(
            const MacroMesh& mm,
            index_t s ) const ;

        index_t nb_surface_facets(
            const MacroMesh& mm,
            index_t s ) const ;

    private:
        index_t surface_begin( index_t s ) const
        {
            return surface_ptr_[ s ] ;
        }

        index_t surface_end( index_t s ) const
        {
            return surface_ptr_[ s + 1 ] ;
        }

        index_t surface_facet( index_t f ) const
        {
            return surface_facets_[ f ] ;
        }

        void initialize( const MacroMesh& mm ) ;

    private:
        bool initialized_ ;

        std::vector< index_t > surface_facets_ ;
        std::vector< index_t > surface_ptr_ ;
        std::vector< index_t > surface2mesh_ ;
    } ;

    class RINGMESH_API MacroMeshCells {
    public:
        MacroMeshCells()
              : mm_( nil ), initialized_( false ), nb_cells_( 0 )
        {
        }

        signed_index_t global_cell_adjacent(
            const MacroMesh* mm,
            index_t mesh,
            index_t c,
            index_t f ) const ;

        index_t get_local_cell_index(
            const MacroMesh* mm,
            index_t global_index ) const ;

        index_t get_mesh(
            const MacroMesh* mm,
            index_t global_index ) const ;

        index_t nb_cells( const MacroMesh* mm ) const ;

    private:
        void initialize( const MacroMesh* mm ) ;

    private:
        const MacroMesh* mm_ ;
        bool initialized_ ;

        std::vector< signed_index_t > global_cell_adjacents_ ;
        std::vector< index_t > cell2mesh_ ;
        index_t nb_cells_ ;
    } ;

    class RINGMESH_API MacroMeshTools {
    public:
        MacroMeshTools( MacroMesh& mm ) ;
        ~MacroMeshTools() ;

        const GEO::MeshFacetsAABB& facet_aabb( index_t region ) ;
        const GEO::MeshTetsAABB& tet_aabb( index_t region ) ;

    private:
        void init_facet_aabb( index_t region ) ;
        void init_tet_aabb( index_t region ) ;

    private:
        MacroMesh& mm_ ;

        std::vector< GEO::MeshFacetsAABB* > facet_aabb_ ;
        std::vector< GEO::MeshTetsAABB* > tet_aabb_ ;
    } ;

    class RINGMESH_API MacroMesh {
    public:
        MacroMesh(
            const BoundaryModel& model,
            index_t dim = 3 ) ;
        MacroMesh( const MacroMesh& mm) ;
        virtual ~MacroMesh() ;

        //    __  __     _   _            _
        //   |  \/  |___| |_| |_  ___  __| |___
        //   | |\/| / -_)  _| ' \/ _ \/ _` (_-<
        //   |_|  |_\___|\__|_||_\___/\__,_/__/
        //
        void compute_tetmesh(
            const TetraMethod& method,
            int region_id = - 1,
            bool add_steiner_points = true,
            MacroMesh* background = nil,
            std::vector< std::vector< vec3 > >& internal_vertices =
                empty_vertices ) ;

        //      _
        //     /_\  __ __ ___ _________ _ _ ___
        //    / _ \/ _/ _/ -_|_-<_-< _ \ '_(_-<
        //   /_/ \_\__\__\___/__/__|___/_| /__/
        //
        GEO::Mesh& mesh( index_t region )
        {
            return *meshes_[ region ] ;
        }

        const GEO::Mesh& mesh( index_t region ) const
        {
            return *meshes_[ region ] ;
        }

        index_t nb_meshes() const
        {
            return meshes_.size() ;
        }

        std::vector< std::vector< Edge > >& well_vertices( index_t region )
        {
            return well_vertices_[ region ] ;
        }

        const std::vector< std::vector< Edge > >& well_vertices( index_t region )
        const
        {
            return well_vertices_[ region ] ;
        }

        const BoundaryModel& model() const
        {
            return model_ ;
        }

        index_t surface_facet(
            index_t s,
            index_t f ) const
        {
            return mm_facets_.surface_facet( *this, s, f ) ;
        }

        index_t surface_mesh( index_t s ) const
        {
            return mm_facets_.surface_mesh( *this, s ) ;
        }

        index_t nb_surface_facets( index_t s ) const
        {
            return mm_facets_.nb_surface_facets( *this, s ) ;
        }

        index_t global_vertex_id(
            index_t mesh,
            index_t v ) const
        {
            return mm_vertices_.global_vertex_id( *this, mesh, v ) ;
        }

        const vec3& global_vertex( index_t global_v ) const
        {
            return mm_vertices_.global_vertex( *this, global_v ) ;
        }

        index_t nb_vertices() const
        {
            return mm_vertices_.nb_vertices( *this ) ;
        }

        index_t nb_vertex_indices() const
        {
            return mm_vertices_.nb_vertex_indices( *this ) ;
        }

        signed_index_t global_cell_adjacent(
            index_t mesh,
            index_t c,
            index_t f ) const
        {
            return mm_cells_.global_cell_adjacent( this, mesh, c, f ) ;
        }

        index_t get_local_cell_index( index_t global_index ) const
        {
            return mm_cells_.get_local_cell_index( this, global_index ) ;
        }

        index_t get_mesh( index_t global_index ) const
        {
            return mm_cells_.get_mesh( this, global_index ) ;
        }

        index_t nb_cells() const
        {
            return mm_cells_.nb_cells( this ) ;
        }


    protected:
        /// BoundaryModel representing the structural information of the mesh
        const BoundaryModel& model_ ;

        /// Vector of meshes, one by region
        std::vector< GEO::Mesh* > meshes_ ;

        /// Vector of constrained edges, one vector by region by well (well_vertices_[r][w] = edges of well w in the region r)
        std::vector< std::vector< std::vector< Edge > > > well_vertices_ ;

    private:
        MacroMeshVertices mm_vertices_ ;
        MacroMeshFacets mm_facets_ ;
        MacroMeshCells mm_cells_ ;

    public:
        MacroMeshTools tools ;

    } ;
}

#endif
