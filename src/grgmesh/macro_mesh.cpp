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

#include <grgmesh/macro_mesh.h>
#include <grgmesh/boundary_model.h>
#include <grgmesh/tetra_gen.h>

namespace GRGMesh {

    MacroMesh::MacroMesh( const BoundaryModel* model, uint8 dim )
        :
            model_( model ),
            meshes_( model->nb_regions(), nil ),
            background_meshes_( model->nb_regions(), nil ),
            vertices_( model->nb_regions() ),
            well_vertices_( model->nb_regions() )
    {
        for( unsigned int r = 0; r < model_->nb_regions(); r++ ) {
            meshes_[r] = new GEO::Mesh( dim ) ;
        }
    }

    MacroMesh::~MacroMesh()
    {
        for( unsigned int r = 0; r < model_->nb_regions(); r++ ) {
            delete meshes_[r] ;
            if( background_meshes_[r] ) delete background_meshes_[r] ;
        }
    }

    uint32 MacroMesh::nb_regions() const
    {
        return model_->nb_regions() ;
    }

    void MacroMesh::initialize_background_meshes( uint8 dim )
    {
        for( unsigned int r = 0; r < model_->nb_regions(); r++ ) {
            background_meshes_[r] = new GEO::Mesh( dim ) ;
        }
    }

    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param method Mesher used
     * @param region_id Region to mesh, -1 for all
     * @param add_steiner_points if true, the mesher will add some points inside the region
     * to improve the mesh quality
     */
    void MacroMesh::compute_tetmesh(
        const TetraMethod& method,
        int region_id,
        bool add_steiner_points )
    {
        if( region_id == -1 ) {
            for( unsigned int i = 0; i < nb_meshes(); i++ ) {
                TetraGen_var tetragen = TetraGen::instantiate( method, mesh( i ),
                    &model_->region( i ), add_steiner_points, vertices( i ),
                    well_vertices( i ), background_mesh( i ) ) ;
                tetragen->tetrahedralize() ;
            }
        } else {
            TetraGen_var tetragen = TetraGen::instantiate( method, mesh( region_id ),
                &model_->region( region_id ), add_steiner_points,
                vertices( region_id ), well_vertices( region_id ),
                background_mesh( region_id ) ) ;
            tetragen->tetrahedralize() ;
        }
    }

}

