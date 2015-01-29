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

#include <geogram/basic/progress.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/algorithm.h>

namespace GRGMesh {

    MacroMesh::MacroMesh( const BoundaryModel& model, index_t dim )
        :
            model_( model ),
            meshes_( model.nb_regions(), nil ),
            well_vertices_( model.nb_regions() ),
            nb_vertices_( -1 ),
            facet_aabb_( model.nb_regions(), nil ),
            tet_aabb_( model.nb_regions(), nil )
    {
        for( unsigned int r = 0; r < model_.nb_regions(); r++ ) {
            meshes_[r] = new GEO::Mesh( dim ) ;
        }
    }

    MacroMesh::~MacroMesh()
    {
        for( unsigned int r = 0; r < model_.nb_regions(); r++ ) {
            delete meshes_[r] ;
            if( facet_aabb_[r] ) delete facet_aabb_[r] ;
            if( tet_aabb_[r] ) delete tet_aabb_[r] ;
        }
    }

    const GEO::MeshFacetsAABB& MacroMesh::facet_aabb( index_t region )
    {
        init_facet_aabb( region ) ;
        return *facet_aabb_[region] ;
    }
    void MacroMesh::init_facet_aabb( index_t region )
    {
        if( facet_aabb_[region] ) return ;
        facet_aabb_[region] = new GEO::MeshFacetsAABB( mesh( region ) ) ;
    }
    void MacroMesh::init_all_facet_aabb()
    {
        for( index_t region = 0; region < nb_meshes(); region++ ) {
            init_facet_aabb( region ) ;
        }
    }
    const GEO::MeshTetsAABB& MacroMesh::tet_aabb( index_t region )
    {
        init_tet_aabb( region ) ;
        return *tet_aabb_[region] ;
    }
    void MacroMesh::init_tet_aabb( index_t region )
    {
        if( tet_aabb_[region] ) return ;
        tet_aabb_[region] = new GEO::MeshTetsAABB( mesh( region ) ) ;
    }
    void MacroMesh::init_all_tet_aabb()
    {
        for( index_t region = 0; region < nb_meshes(); region++ ) {
            init_tet_aabb( region ) ;
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
        bool add_steiner_points,
        MacroMesh* background,
        std::vector< std::vector< vec3 > >& internal_vertices )
    {
        if( region_id == -1 ) {
            GEO::ProgressTask progress( "Compute", nb_meshes() ) ;
            for( unsigned int i = 0; i < nb_meshes(); i++ ) {
                GEO::Mesh* background_mesh =
                    background ? &background->mesh( i ) : nil ;
                const std::vector< vec3 >& vertices =
                    internal_vertices.empty() ? empty_vector : internal_vertices[i] ;
                TetraGen_var tetragen = TetraGen::instantiate( method, mesh( i ),
                    &model_.region( i ), add_steiner_points, vertices,
                    well_vertices( i ), background_mesh ) ;
                tetragen->tetrahedralize() ;
                progress.next() ;
            }
        } else {
            GEO::Mesh* background_mesh =
                background ? &background->mesh( region_id ) : nil ;
            const std::vector< vec3 >& vertices =
                internal_vertices.empty() ?
                    empty_vector : internal_vertices[region_id] ;
            TetraGen_var tetragen = TetraGen::instantiate( method, mesh( region_id ),
                &model_.region( region_id ), add_steiner_points, vertices,
                well_vertices( region_id ), background_mesh ) ;
            tetragen->tetrahedralize() ;
        }
    }

    void MacroMesh::unique_points(
        std::vector< vec3 >& unique_vertices,
        std::vector< index_t >& indices )
    {
        vertex2mesh_.resize(meshes_.size(),0) ;

        index_t nb_non_unique_vertices = 0 ;
        for( index_t i = 0; i < meshes_.size(); i++ ) {
            vertex2mesh_[i] = nb_non_unique_vertices ;
            nb_non_unique_vertices += meshes_[i]->nb_vertices() ;

        }
        std::vector< vec3 > all_vertices( nb_non_unique_vertices ) ;
        index_t index = 0 ;
        for( index_t i = 0; i < meshes_.size(); i++ ) {
            index_t nb_vertices =  meshes_[i]->nb_vertices() ;
            vertex2mesh_[i+1] = nb_vertices ;
            for( index_t j = 0; j < nb_vertices; j++ ) {
                all_vertices[index] = vec3( meshes_[i]->vertex_ptr( j )[0],
                    meshes_[i]->vertex_ptr( j )[1],
                    meshes_[i]->vertex_ptr( j )[2] ) ;
                index++ ;
            }
        }
        MakeUnique mu( all_vertices ) ;
        mu.unique() ;
        mu.unique_points( unique_vertices ) ;
        indices = mu.indices() ;
    }

    bool MacroMesh::surface_vertices_global_id(
        std::vector< index_t > surface_id,
        std::vector< index_t >& indices,
        std::vector< vec3 >& unique_vertices )
    {
        index_t nb_total_nodes ;
        for( index_t s = 0; s < surface_id.size(); s++ ) {
            nb_total_nodes += model_.surface( surface_id[s] ).nb_vertices() ;
        }
        indices.clear() ;
        indices.reserve( nb_total_nodes ) ;
        ColocaterANN ann( unique_vertices ) ;
        for( index_t s = 0; s < surface_id.size(); s++ ) {
            const Surface& surface = model_.surface( surface_id[s] ) ;
            for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                vec3 cur_v = surface.vertex( v ) ;
                std::vector< index_t > results ;
                if( !ann.get_colocated( cur_v, results, 1 ) ) {
                    GEO::Logger::err( "" )
                        << "Impossible to find colocated point mesh/model"
                        << std::endl ;
                    return false ;
                }
                indices.push_back( results[0] ) ;
            }
        }
        GEO::sort_unique( indices ) ;
        return true ;
    }
    index_t MacroMesh::nb_vertices() const
    {
        init_vertices() ;
        return unique_vertices_.size() ;
    }


    void MacroMesh::init_surfaces()
    {
        if( !surface_facets_.empty() ) return ;
        surface2mesh_.resize( model_.nb_surfaces(), Surface::NO_ID ) ;
        surface_ptr_.resize( model_.nb_surfaces() + 1, 0 ) ;

        for( index_t m = 0; m < nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mesh( m ) ;
            std::vector< signed_index_t > surface_proccessed ;
            for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                if( surface2mesh_[surface_id] != Surface::NO_ID ) continue ;
                if( !Utils::contains( surface_proccessed, surface_id ) ) {
                    surface_proccessed.push_back( surface_id ) ;
                }
                surface_ptr_[surface_id+1]++ ;
            }
            for( index_t s = 0; s < surface_proccessed.size(); s++ ) {
                surface2mesh_[surface_proccessed[s]] = m ;
            }
        }

        for( index_t s = 0; s < model_.nb_surfaces(); s++ ) {
            surface_ptr_[s+1] += surface_ptr_[s] ;
        }

        surface_facets_.resize( surface_ptr_.back() ) ;

        std::vector< index_t > surface_facet_index( model_.nb_surfaces(), 0 ) ;
        for( index_t m = 0; m < nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mesh( m ) ;
            for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                if( surface2mesh_[surface_id] != m ) continue ;
                surface_facets_[surface_ptr_[surface_id]
                    + surface_facet_index[surface_id]++ ] = f ;
            }
        }
    }

    index_t MacroMesh::surface_begin( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMesh* >( this )->init_surfaces() ;
        }
        return surface_ptr_[s] ;
    }
    index_t MacroMesh::surface_end( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMesh* >( this )->init_surfaces() ;
        }
        return surface_ptr_[s+1] ;
    }
    index_t MacroMesh::surface_mesh( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMesh* >( this )->init_surfaces() ;
        }
        return surface2mesh_[s] ;
    }

    void MacroMesh::init_vertices() const {
        if( !unique_vertices_.empty() ) return ;
        MacroMesh* not_const = const_cast< MacroMesh* >( this ) ;
        not_const->unique_points( not_const->unique_vertices_,
            not_const->global_vertex_indices_ ) ;
    }

    index_t MacroMesh::global_vertex_id(index_t mesh, index_t v) const {
        init_vertices() ;
        return global_vertex_indices_[vertex2mesh_[mesh] + v] ;
    }

    const vec3& MacroMesh::vertex( index_t global_v) const {
        return unique_vertices_[global_v] ;
    }

}

