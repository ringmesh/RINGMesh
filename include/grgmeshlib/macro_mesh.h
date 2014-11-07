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

#include <grgmeshlib/common.h>
#include <grgmeshlib/mixed_mesh.h>
#include <grgmeshlib/utils.h>

namespace GRGMesh {

    class BoundaryModel ;

    template< class T >
    class GRGMESH_API MacroMesh {
    public:
        MacroMesh( const BoundaryModel* model ) ;
        virtual ~MacroMesh() {}

        T& mesh( int x )
        {
            return meshes_[x] ;
        }
        const T& mesh( int x ) const
        {
            return meshes_[x] ;
        }
        unsigned int nb_meshes() const
        {
            return meshes_.size() ;
        }

    protected:
        /// BoundaryModel representing the structural information of the mesh
        const BoundaryModel* model_ ;
        /// Vector of meshes, one by region
        std::vector< T > meshes_ ;
        /// Vector of background meshes, one by region
        std::vector< T* > background_meshes_ ;
        /// Vector of constrained vertices, one vector by region
        std::vector< std::vector< vec3 > > vertices_ ;
        /// Vector of constrained edges, one vector by region by well (well_vertices_[r][w] = edges of well w in the region r)
        std::vector< std::vector< std::vector< Edge > > > well_vertices_ ;
    } ;

    class GRGMESH_API MacroMixedMesh: public MacroMesh< MixedMesh > {
    public:
        MacroMixedMesh( const BoundaryModel* model )
            : MacroMesh< MixedMesh >( model )
        {
        }

        void compute_tetmesh(
            const Tetra_method& method,
            int region_id = -1,
            bool add_steiner_points = true ) ;

    } ;

}

#endif
