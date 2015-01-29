/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#ifndef __GRGMESH_MACRO_MESH__
#define __GRGMESH_MACRO_MESH__

#include <grgmesh/common.h>
#include <grgmesh/utils.h>

#include <geogram/mesh/mesh.h>

namespace GEO {
    class MeshTetsAABB ;
    class MeshFacetsAABB ;
}

namespace GRGMesh {

    class BoundaryModel ;
    class MacroMesh ;

    static std::vector< std::vector< vec3 > > empty_vertices ;

    class GRGMESH_API MacroMeshVertices {
    public:
        MacroMeshVertices( const MacroMesh& mm )
            : mm_( mm ), initialized_( false )
        {
        }

        index_t nb_vertices() const ;
        index_t nb_vertex_indices() const ;
        index_t global_vertex_id( index_t mesh, index_t v ) const ;
        const vec3& global_vertex( index_t global_v ) const ;

    private:
        void initialize() ;

    private:
        const MacroMesh& mm_ ;
        bool initialized_ ;

        std::vector< vec3 > unique_vertices_ ;
        std::vector< index_t > global_vertex_indices_ ;
        std::vector< index_t > vertex2mesh_ ;
    } ;

    class GRGMESH_API MacroMeshFacets {
    public:
        MacroMeshFacets( const MacroMesh& mm )
            : mm_( mm ), initialized_( false )
        {
        }
        index_t surface_facet( index_t s, index_t f ) const ;
        index_t surface_mesh( index_t s ) const ;
        index_t nb_surface_facets( index_t s ) const ;

    private:
        index_t surface_begin( index_t s ) const {
            return surface_ptr_[s] ;
        }
        index_t surface_end( index_t s ) const {
            return surface_ptr_[s + 1] ;
        }
        index_t surface_facet( index_t f ) const {
            return surface_facets_[f] ;
        }
        void initialize() ;

    private:
        const MacroMesh& mm_ ;
        bool initialized_ ;

        std::vector< index_t > surface_facets_ ;
        std::vector< index_t > surface_ptr_ ;
        std::vector< index_t > surface2mesh_ ;
    } ;


    class GRGMESH_API MacroMeshCells {
    public:
        MacroMeshCells( const MacroMesh& mm )
            : mm_( mm ), initialized_( false )
        {
        }

        signed_index_t global_cell_adjacent( index_t mesh, index_t c, index_t f ) const ;
        index_t get_local_cell_index( index_t global_index ) const ;
        index_t get_mesh( index_t global_index ) const ;

    private:
        void initialize() ;

    private:
        const MacroMesh& mm_ ;
        bool initialized_ ;

        std::vector< signed_index_t > global_cell_adjacents_ ;
        std::vector< index_t > cell2mesh_ ;
    } ;

    class GRGMESH_API MacroMesh {
    public:
        MacroMesh( const BoundaryModel& model, index_t dim = 3 ) ;
        virtual ~MacroMesh() ;

        //    __  __     _   _            _
        //   |  \/  |___| |_| |_  ___  __| |___
        //   | |\/| / -_)  _| ' \/ _ \/ _` (_-<
        //   |_|  |_\___|\__|_||_\___/\__,_/__/
        //
        void compute_tetmesh(
            const TetraMethod& method,
            int region_id = -1,
            bool add_steiner_points = true,
            MacroMesh* background = nil,
            std::vector< std::vector< vec3 > >& internal_vertices =
                empty_vertices ) ;


        const GEO::MeshFacetsAABB& facet_aabb( index_t region ) ;
        void init_facet_aabb( index_t region ) ;
        void init_all_facet_aabb() ;
        const GEO::MeshTetsAABB& tet_aabb( index_t region ) ;
        void init_tet_aabb( index_t region ) ;
        void init_all_tet_aabb() ;


        //      _
        //     /_\  __ __ ___ _________ _ _ ___
        //    / _ \/ _/ _/ -_|_-<_-< _ \ '_(_-<
        //   /_/ \_\__\__\___/__/__|___/_| /__/
        //
        GEO::Mesh& mesh( index_t region )
        {
            return *meshes_[region] ;
        }
        const GEO::Mesh& mesh( index_t region ) const
        {
            return *meshes_[region] ;
        }
        index_t nb_meshes() const
        {
            return meshes_.size() ;
        }
        std::vector< std::vector< Edge > >& well_vertices( index_t region )
        {
            return well_vertices_[region] ;
        }
        const std::vector< std::vector< Edge > >& well_vertices( index_t region ) const
        {
            return well_vertices_[region] ;
        }
        const BoundaryModel& model() const
        {
            return model_ ;
        }

        index_t surface_facet( index_t s, index_t f ) const  {
            return mm_facets_.surface_facet( s, f ) ;
        }
        index_t surface_mesh( index_t s ) const {
            return mm_facets_.surface_mesh( s ) ;
        }
        index_t nb_surface_facets( index_t s ) const {
            return mm_facets_.nb_surface_facets( s ) ;
        }

        index_t global_vertex_id( index_t mesh, index_t v ) const {
            return mm_vertices_.global_vertex_id( mesh, v ) ;
        }
        const vec3& global_vertex( index_t global_v) const {
            return mm_vertices_.global_vertex( global_v ) ;
        }
        index_t nb_vertices() const {
            return mm_vertices_.nb_vertices() ;
        }
        index_t nb_vertex_indices() const {
            return mm_vertices_.nb_vertex_indices() ;
        }

    protected:
        /// BoundaryModel representing the structural information of the mesh
        const BoundaryModel& model_ ;
        /// Vector of meshes, one by region
        std::vector< GEO::Mesh* > meshes_ ;
        /// Vector of constrained edges, one vector by region by well (well_vertices_[r][w] = edges of well w in the region r)
        std::vector< std::vector< std::vector< Edge > > > well_vertices_ ;

    private:
        std::vector< GEO::MeshFacetsAABB* > facet_aabb_ ;
        std::vector< GEO::MeshTetsAABB* > tet_aabb_ ;

        MacroMeshVertices mm_vertices_ ;
        MacroMeshFacets mm_facets_ ;

    } ;

}

#endif
