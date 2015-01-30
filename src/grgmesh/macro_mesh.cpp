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

    void MacroMeshVertices::initialize()
    {
        vertex2mesh_.resize( mm_.nb_meshes(), 0 ) ;

        index_t nb_non_unique_vertices = 0 ;
        for( index_t i = 0; i < mm_.nb_meshes(); i++ ) {
            vertex2mesh_[i] = nb_non_unique_vertices ;
            nb_non_unique_vertices += mm_.mesh( i ).nb_vertices() ;

        }
        std::vector< vec3 > all_vertices( nb_non_unique_vertices ) ;
        index_t index = 0 ;
        for( index_t i = 0; i < mm_.nb_meshes(); i++ ) {
            index_t nb_vertices = mm_.mesh( i ).nb_vertices() ;
            vertex2mesh_[i + 1] = nb_vertices ;
            for( index_t j = 0; j < nb_vertices; j++ ) {
                all_vertices[index] = GEO::Geom::mesh_vertex( mm_.mesh( i ), j ) ;
                index++ ;
            }
        }
        MakeUnique mu( all_vertices ) ;
        mu.unique() ;
        mu.unique_points( unique_vertices_ ) ;
        global_vertex_indices_ = mu.indices() ;

        initialized_ = true ;
    }

    index_t MacroMeshVertices::nb_vertices() const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return unique_vertices_.size() ;
    }

    index_t MacroMeshVertices::nb_vertex_indices() const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return global_vertex_indices_.size() ;
    }

    index_t MacroMeshVertices::global_vertex_id( index_t mesh, index_t v ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return global_vertex_indices_[vertex2mesh_[mesh] + v] ;
    }

    const vec3& MacroMeshVertices::global_vertex( index_t global_v ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return unique_vertices_[global_v] ;
    }


    void MacroMeshFacets::initialize()
    {
        surface2mesh_.resize( mm_.model().nb_surfaces(), Surface::NO_ID ) ;
        surface_ptr_.resize( mm_.model().nb_surfaces() + 1, 0 ) ;

        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
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

        for( index_t s = 0; s < mm_.model().nb_surfaces(); s++ ) {
            surface_ptr_[s+1] += surface_ptr_[s] ;
        }

        surface_facets_.resize( surface_ptr_.back() ) ;

        std::vector< index_t > surface_facet_index( mm_.model().nb_surfaces(), 0 ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
            for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                if( surface2mesh_[surface_id] != m ) continue ;
                surface_facets_[surface_ptr_[surface_id]
                    + surface_facet_index[surface_id]++ ] = f ;
            }
        }

        initialized_ = true ;
    }

    index_t MacroMeshFacets::surface_mesh( index_t s ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface2mesh_[s] ;
    }

    index_t MacroMeshFacets::surface_facet( index_t s, index_t f ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface_facet( surface_begin( s ) + f ) ;
    }

    index_t MacroMeshFacets::nb_surface_facets( index_t s ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface_end( s ) - surface_begin( s ) ;
    }


    void MacroMeshCells::initialize()
    {
        //todo
        grgmesh_assert_not_reached ;
    }

    signed_index_t MacroMeshCells::global_cell_adjacent( index_t m, index_t c, index_t f ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize() ;
        }
        return global_cell_adjacents_[cell2mesh_[m]
            + mm_.mesh( m ).cell_adjacents_begin( c ) + f] ;
    }
    index_t MacroMeshCells::get_local_cell_index( index_t global_index ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize() ;
        }
        index_t m = get_mesh( global_index ) ;
        return global_index - cell2mesh_[m] ;
    }
    index_t MacroMeshCells::get_mesh( index_t global_index ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize() ;
        }
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            if( global_index < cell2mesh_[m+1] ) return m ;
        }
        grgmesh_assert_not_reached ;
        return dummy_index_t ;
    }



    MacroMesh::MacroMesh( const BoundaryModel& model, index_t dim )
        :
            model_( model ),
            meshes_( model.nb_regions(), nil ),
            well_vertices_( model.nb_regions() ),
            facet_aabb_( model.nb_regions(), nil ),
            tet_aabb_( model.nb_regions(), nil ),
            mm_vertices_( *this ), // The "this" pointer shouldn't be used in the initialization list
            mm_facets_( *this )    // Undefined behavior if the class is derived - Windows Warning C4355 
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

}

