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

    class GRGMESH_API MacroMesh {
    public:
        MacroMesh( const BoundaryModel* model, uint8 dim = 3 ) ;
        virtual ~MacroMesh() ;
        void initialize_background_meshes( uint8 dim = 3 ) ;

        //    __  __     _   _            _
        //   |  \/  |___| |_| |_  ___  __| |___
        //   | |\/| / -_)  _| ' \/ _ \/ _` (_-<
        //   |_|  |_\___|\__|_||_\___/\__,_/__/
        //
        void compute_tetmesh(
            const TetraMethod& method,
            int region_id = -1,
            bool add_steiner_points = true ) ;

        void unique_points(
            std::vector< vec3 >& unique_vertices,
            std::vector< int >& indices ) const ;


        const GEO::MeshFacetsAABB& facet_aabb( uint32 region ) ;
        void init_facet_aabb( uint32 region ) ;
        void init_all_facet_aabb() ;
        const GEO::MeshTetsAABB& tet_aabb( uint32 region ) ;
        void init_tet_aabb( uint32 region ) ;
        void init_all_tet_aabb() ;

        //      _
        //     /_\  __ __ ___ _________ _ _ ___
        //    / _ \/ _/ _/ -_|_-<_-< _ \ '_(_-<
        //   /_/ \_\__\__\___/__/__|___/_| /__/
        //
        GEO::Mesh& mesh( uint32 region )
        {
            return *meshes_[region] ;
        }
        const GEO::Mesh& mesh( uint32 region ) const
        {
            return *meshes_[region] ;
        }
        GEO::Mesh* background_mesh( uint32 region )
        {
            return background_meshes_[region] ;
        }
        const GEO::Mesh* background_mesh( uint32 region ) const
        {
            return background_meshes_[region] ;
        }
        unsigned int nb_meshes() const
        {
            return meshes_.size() ;
        }
        std::vector< vec3 >& vertices( uint32 region ) {
            return vertices_[region] ;
        }
        const std::vector< vec3 >& vertices( uint32 region ) const {
            return vertices_[region] ;
        }
        std::vector< std::vector< Edge > >& well_vertices( uint32 region )
        {
            return well_vertices_[region] ;
        }
        const std::vector< std::vector< Edge > >& well_vertices( uint32 region ) const
        {
            return well_vertices_[region] ;
        }
        const BoundaryModel* model() const
        {
            return model_ ;
        }
      
        uint32 nb_vertices() ;

    protected:
        /// BoundaryModel representing the structural information of the mesh
        const BoundaryModel* model_ ;
        /// Vector of meshes, one by region
        std::vector< GEO::Mesh* > meshes_ ;
        /// Vector of background meshes, one by region
        std::vector< GEO::Mesh* > background_meshes_ ;
        /// Vector of constrained vertices, one vector by region
        std::vector< std::vector< vec3 > > vertices_ ;
        /// Vector of constrained edges, one vector by region by well (well_vertices_[r][w] = edges of well w in the region r)
        std::vector< std::vector< std::vector< Edge > > > well_vertices_ ;
    private:
        std::vector< GEO::MeshFacetsAABB* > facet_aabb_ ;
        std::vector< GEO::MeshTetsAABB* > tet_aabb_ ;
        uint32 nb_vertices_ ;
    } ;

}

#endif
