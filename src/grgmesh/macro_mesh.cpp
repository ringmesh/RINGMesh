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

#include <grgmesh/boundary_model.h>
#include <grgmesh/macro_mesh.h>
#include <grgmesh/tetra_gen.h>

namespace GRGMesh {

    template< class T >
    MacroMesh< T >::MacroMesh( const BoundaryModel* model )
        :
            model_( model ),
            meshes_( model->nb_regions() ),
            background_meshes_( model->nb_regions(), nil )
    {
    }

    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param method Mesher used
     * @param region_id Region to mesh, -1 for all
     * @param add_steiner_points if true, the mesher will add some points inside the region
     * to improve the mesh quality
     */
    void MacroMixedMesh::compute_tetmesh(
        const TetraMethod& method,
        int region_id,
        bool add_steiner_points )
    {
        if( region_id == -1 ) {
            for( unsigned int i = 0; i < nb_meshes(); i++ ) {
                TetraGen_var tetragen = TetraGen::instantiate( method, meshes_[i],
                    &model_->region( i ), add_steiner_points, vertices_[i],
                    well_vertices_[i], background_meshes_[i] ) ;
                tetragen->tetrahedralize() ;
            }
        } else {
            TetraGen_var tetragen = TetraGen::instantiate( method, meshes_[region_id],
                &model_->region( region_id ), add_steiner_points,
                vertices_[region_id], well_vertices_[region_id],
                background_meshes_[region_id] ) ;
            tetragen->tetrahedralize() ;
        }
    }

}

