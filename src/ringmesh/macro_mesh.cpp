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


#include <ringmesh/macro_mesh.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/tetra_gen.h>

#include <geogram/basic/progress.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/algorithm.h>

namespace RINGMesh {

    void MacroMeshVertices::initialize( const MacroMesh& mm )
    {
        vertex2mesh_.resize( mm.nb_meshes(), 0 ) ;

        index_t nb_non_unique_vertices = 0 ;
        for( index_t i = 0; i < mm.nb_meshes(); i++ ) {
            vertex2mesh_[i] = nb_non_unique_vertices ;
            nb_non_unique_vertices += mm.mesh( i ).vertices.nb() ;

        }
        std::vector< vec3 > all_vertices( nb_non_unique_vertices ) ;
        index_t index = 0 ;
        for( index_t i = 0; i < mm.nb_meshes(); i++ ) {
            index_t nb_vertices = mm.mesh( i ).vertices.nb() ;
            for( index_t j = 0; j < nb_vertices; j++ ) {
                all_vertices[index] = GEO::Geom::mesh_vertex( mm.mesh( i ), j ) ;
                index++ ;
            }
        }
        MakeUnique mu( all_vertices ) ;
        mu.unique() ;
        mu.unique_points( unique_vertices_ ) ;
        global_vertex_indices_ = mu.indices() ;

        initialized_ = true ;
    }

    index_t MacroMeshVertices::nb_vertices( const MacroMesh& mm ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize( mm ) ;
        }
        return unique_vertices_.size() ;
    }

    index_t MacroMeshVertices::nb_vertex_indices( const MacroMesh& mm ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize( mm ) ;
        }
        return global_vertex_indices_.size() ;
    }

    index_t MacroMeshVertices::global_vertex_id( const MacroMesh& mm, index_t mesh, index_t v ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize( mm ) ;
        }
        ringmesh_debug_assert( v < mm.mesh( mesh ).vertices.nb() ) ;
        return global_vertex_indices_[vertex2mesh_[mesh] + v] ;
    }

    const vec3& MacroMeshVertices::global_vertex( const MacroMesh& mm, index_t global_v ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshVertices* >( this )->initialize( mm ) ;
        }
        return unique_vertices_[global_v] ;
    }


    void MacroMeshFacets::initialize( const MacroMesh& mm )
    {
        surface2mesh_.resize( mm.model().nb_surfaces(), Surface::NO_ID ) ;
        surface_ptr_.resize( mm.model().nb_surfaces() + 1, 0 ) ;

        for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
            std::vector< index_t > surface_proccessed ;
            GEO::Attribute< index_t > attribute( cur_mesh.facets.attributes(), surface_att_name ) ;
            for( index_t f = 0; f < cur_mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
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

        for( index_t s = 0; s < mm.model().nb_surfaces(); s++ ) {
            surface_ptr_[s+1] += surface_ptr_[s] ;
        }

        surface_facets_.resize( surface_ptr_.back() ) ;

        std::vector< index_t > surface_facet_index( mm.model().nb_surfaces(), 0 ) ;
        for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
            GEO::Attribute< index_t > attribute( cur_mesh.facets.attributes(), surface_att_name ) ;
            for( index_t f = 0; f < cur_mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
                if( surface2mesh_[surface_id] != m ) continue ;
                surface_facets_[surface_ptr_[surface_id]
                    + surface_facet_index[surface_id]++ ] = f ;
            }
        }

        initialized_ = true ;
    }

    index_t MacroMeshFacets::surface_mesh( const MacroMesh& mm, index_t s ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize( mm ) ;
        }
        return surface2mesh_[s] ;
    }

    index_t MacroMeshFacets::surface_facet( const MacroMesh& mm, index_t s, index_t f ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize( mm ) ;
        }
        return surface_facet( surface_begin( s ) + f ) ;
    }

    index_t MacroMeshFacets::nb_surface_facets( const MacroMesh& mm, index_t s ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshFacets* >( this )->initialize( mm ) ;
        }
        return surface_end( s ) - surface_begin( s ) ;
    }

    void MacroMeshCells::initialize( const MacroMesh* mm )
    {
        mm_ = mm ;
        cell2mesh_.reserve( mm_->nb_meshes()+1 ) ;
        cell2mesh_.push_back( 0 ) ;
        index_t total_adjacents = 0 ;
        for( index_t m = 0; m < mm_->nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_->mesh( m ) ;
            nb_cells_ +=  mesh.cells.nb() ;
            cell2mesh_.push_back( nb_cells_ ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                total_adjacents += mesh.cells.nb_facets( c ) ;
            }
        }

        index_t nb_vertices = mm_->nb_vertices() ;
        std::vector< std::vector< index_t > > cells_around_vertex( nb_vertices ) ;
        global_cell_adjacents_.reserve( total_adjacents ) ;
        for( index_t m = 0; m < mm_->nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_->mesh( m ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t adj = mesh.cells.adjacent( c, f ) ;
                    if( adj != GEO::NO_CELL ) {
                        adj += cell2mesh_[m] ;
                    } else {
                        for( index_t v = 0; v < mesh.cells.facet_nb_vertices( c, f );
                            v++ ) {
                            index_t vertex_id = mm_->global_vertex_id( m,
                                mesh.cells.facet_vertex( c, f, v ) ) ;
                            cells_around_vertex[vertex_id].push_back( cell2mesh_[m] + c ) ;
                        }
                    }
                    global_cell_adjacents_.push_back( adj ) ;
                }
            }
        }

        for( index_t v = 0; v < cells_around_vertex.size(); v++ ) {
            GEO::sort_unique( cells_around_vertex[v] ) ;
        }

        for( index_t m = 0; m < mm_->nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_->mesh( m ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t adj = mesh.cells.adjacent( c, f ) ;
                    if( adj == GEO::NO_CELL ) {
                        index_t prev_vertex_id = mm_->global_vertex_id(
                            m,
                            mesh.cells.facet_vertex( c, f, 0 ) ) ;
                        std::vector< index_t >& prev_cells =
                            cells_around_vertex[prev_vertex_id] ;
                        index_t size_hint = prev_cells.size() ;
                        std::vector< index_t > intersection( size_hint ) ;
                        for( index_t v = 1; v < mesh.cells.facet_nb_vertices( c, f );
                            v++ ) {
                            index_t vertex_id = mm_->global_vertex_id( m,
                                mesh.cells.facet_vertex( c, f, v ) ) ;
                            std::vector< index_t >& cells =
                                cells_around_vertex[vertex_id] ;
                            intersection.erase( std::set_intersection( prev_cells.begin(),
                                prev_cells.end(), cells.begin(), cells.end(),
                                intersection.begin() ), intersection.end() ) ;
                            if( intersection.size() < 2 ) break ;
                            prev_cells = intersection ;
                        }

                        if( intersection.size() > 1 ) {
                            ringmesh_debug_assert( intersection.size() == 2 ) ;
                            signed_index_t new_adj =
                                intersection[0] == cell2mesh_[m] + c ?
                                    intersection[1] : intersection[0] ;
                            global_cell_adjacents_[cell2mesh_[m]
                                + mesh.cells.facets_begin( c ) + f] = new_adj ;
                        }
                    }
                }
            }
        }

        initialized_ = true ;
    }

    signed_index_t MacroMeshCells::global_cell_adjacent(
        const MacroMesh* mm,
        index_t m,
        index_t c,
        index_t f ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize( mm ) ;
        }
        return global_cell_adjacents_[cell2mesh_[m]
            + mm_->mesh( m ).cells.facets_begin( c ) + f] ;
    }
    index_t MacroMeshCells::get_local_cell_index(
        const MacroMesh* mm,
        index_t global_index ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize( mm ) ;
        }
        index_t m = get_mesh( mm, global_index ) ;
        return global_index - cell2mesh_[m] ;
    }
    index_t MacroMeshCells::get_mesh(
        const MacroMesh* mm,
        index_t global_index ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize( mm ) ;
        }
        for( index_t m = 0; m < mm_->nb_meshes(); m++ ) {
            if( global_index < cell2mesh_[m+1] ) return m ;
        }
        ringmesh_assert_not_reached ;
        return dummy_index_t ;
    }
    index_t MacroMeshCells::nb_cells( const MacroMesh* mm ) const
    {
        if( !initialized_ ) {
            const_cast< MacroMeshCells* >( this )->initialize( mm ) ;
        }
        return nb_cells_ ;
    }


    MacroMesh::MacroMesh( const BoundaryModel& model, index_t dim )
        :
            model_( model ),
            meshes_( model.nb_regions(), nil ),
            well_vertices_( model.nb_regions() ),
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

}

