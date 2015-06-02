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

#include <ringmesh/macro_mesh.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/tetra_gen.h>
#include <ringmesh/well.h>

#include <geogram/basic/progress.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/algorithm.h>

#include <stack>

namespace RINGMesh {

    /*!
     * Initializes the database of the MacroMesh vertices
     */
    void MacroMeshVertices::initialize()
    {
        vertex2mesh_.resize( mm_.nb_meshes(), 0 ) ;

        index_t nb_non_unique_vertices = 0 ;
        for( index_t i = 0; i < mm_.nb_meshes(); i++ ) {
            vertex2mesh_[i] = nb_non_unique_vertices ;
            nb_non_unique_vertices += mm_.mesh( i ).vertices.nb() ;

        }
        std::vector< vec3 > all_vertices( nb_non_unique_vertices ) ;
        index_t index = 0 ;
        for( index_t i = 0; i < mm_.nb_meshes(); i++ ) {
            index_t nb_vertices = mm_.mesh( i ).vertices.nb() ;
            for( index_t j = 0; j < nb_vertices; j++ ) {
                all_vertices[index] = GEO::Geom::mesh_vertex( mm_.mesh( i ), j ) ;
                index++ ;
            }
        }
        MakeUnique mu( all_vertices ) ;
        mu.unique() ;
        mu.unique_points( unique_vertices_ ) ;
        global_vertex_indices_ = mu.indices() ;
    }

    /*!
     * Initializes the database of the duplicated MacroMesh vertices
     */
    void MacroMeshVertices::initialize_duplication()
    {
        if( unique_vertices_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        /// 0 - Fills cell corner information
        mesh_cell_corner_ptr_.resize( mm_.nb_meshes() + 1, 0 ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ )
                mesh_cell_corner_ptr_[m + 1] += mesh.cells.nb_vertices( c ) ;
            mesh_cell_corner_ptr_[m + 1] += mesh_cell_corner_ptr_[m] ;
        }
        cell_corners_.resize( mesh_cell_corner_ptr_.back() ) ;

        std::vector< index_t > cur_cell_index_type( NB_CELL_TYPES * mm_.nb_meshes(),
            0 ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            GEO::Mesh& mesh = const_cast< GEO::Mesh& >( mm_.mesh( m ) ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    cell_corners_[mesh_cell_corner_ptr_[m]
                        + mesh.cells.corners_begin( c ) + v] = mesh.cells.vertex( c,
                        v ) ;
                }
            }
        }

        /// 1 - Gets all the corner vertices (redondant information)
        std::vector< vec3 > corner_vertices( cell_corners_.size() ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            index_t mesh_start = mesh_cell_corner_ptr_[m] ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    corner_vertices[mesh_start + mesh.cells.corners_begin( c ) + v] =
                        GEO::Geom::mesh_vertex( mesh, mesh.cells.vertex( c, v ) ) ;
                }
            }
        }

        /// 2 - Tags all corners to duplicate
        const BoundaryModel& model = mm_.model() ;
        std::vector< SurfaceAction > surface_actions( model.nb_surfaces(), SKIP ) ;
        std::vector< bool > corner_to_duplicate( corner_vertices.size(), false ) ;
        {
            ColocaterANN ann( corner_vertices ) ;
            for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                if( !is_surface_to_duplicate( s ) ) continue ;
                surface_actions[s] = TO_PROCESS ;
                const Surface& surface = model.surface( s ) ;
                for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                    std::vector< index_t > colocated_corners ;
                    ann.get_colocated( surface.vertex( v ), colocated_corners ) ;
                    for( index_t co = 0; co < colocated_corners.size(); co++ ) {
                        corner_to_duplicate[colocated_corners[co]] = true ;
                    }
                }
            }
        }
        corner_vertices.clear() ;

        /// 3 - Duplicates the corners
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
                surface_att_name ) ;
            ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    index_t co = mesh_cell_corner_ptr_[m]
                        + mesh.cells.corners_begin( c ) + v ;
                    if( !corner_to_duplicate[co] ) continue ;

                    index_t vertex_id = mesh.cells.vertex( c, v ) ;
                    std::vector< index_t > corner_used ;
                    std::vector< index_t > cell_added ;
                    std::set< surface_side > surfaces ;
                    std::stack< index_t > S ;
                    S.push( c ) ;
                    cell_added.push_back( c ) ;
                    do {
                        index_t cur_c = S.top() ;
                        S.pop() ;
                        for( index_t cur_v = 0;
                            cur_v < mesh.cells.nb_vertices( cur_c ); cur_v++ ) {
                            if( mesh.cells.vertex( cur_c, cur_v ) == vertex_id ) {
                                index_t cur_co = mesh_cell_corner_ptr_[m]
                                    + mesh.cells.corners_begin( cur_c ) + cur_v ;
                                corner_to_duplicate[cur_co] = false ;
                                corner_used.push_back( cur_co ) ;
                                break ;
                            }
                        }
                        for( index_t cur_f = 0;
                            cur_f < mesh.cells.nb_facets( cur_c ); cur_f++ ) {
                            for( index_t cur_v = 0;
                                cur_v < mesh.cells.facet_nb_vertices( cur_c, cur_f );
                                cur_v++ ) {
                                if( mesh.cells.facet_vertex( cur_c, cur_f, cur_v )
                                    == vertex_id ) {
                                    std::vector< index_t > result ;
                                    if( ann.get_colocated(
                                        Geom::mesh_cell_facet_center( mesh, cur_c,
                                            cur_f ), result ) ) {
                                        index_t surface_id = attribute[result[0]] ;
                                        vec3 facet_normal =
                                            GEO::Geom::mesh_facet_normal( mesh,
                                                result[0] ) ;
                                        vec3 cell_facet_normal =
                                            Geom::mesh_cell_facet_normal( mesh,
                                                cur_c, cur_f ) ;
                                        SurfaceAction side = SurfaceAction(
                                            dot( facet_normal, cell_facet_normal )
                                                > 0 ) ;
                                        surfaces.insert(
                                            surface_side( surface_id, side ) ) ;
                                    } else {
                                        index_t cur_adj = mesh.cells.adjacent( cur_c,
                                            cur_f ) ;
                                        if( cur_adj != GEO::NO_CELL
                                            && !Utils::contains( cell_added,
                                                cur_adj ) ) {
                                            cell_added.push_back( cur_adj ) ;
                                            S.push( cur_adj ) ;
                                        }
                                    }
                                    break ;
                                }
                            }
                        }
                    } while( !S.empty() ) ;

                    if( duplicate_corner( surfaces, surface_actions ) ) {
                        index_t duplicated_vertex_id = mm_.vertices.nb_vertices()
                            + duplicated_vertex_indices_.size() ;
                        index_t global_vertex_id = mm_.vertices.vertex_id( m,
                            vertex_id ) ;
                        duplicated_vertex_indices_.push_back( global_vertex_id ) ;
                        for( index_t cur_co = 0; cur_co < corner_used.size();
                            cur_co++ ) {
                            cell_corners_[corner_used[cur_co]] =
                                duplicated_vertex_id ;
                        }
                    }
                }
            }
        }
    }

    /*!
     * Tests is the cell corner id should be duplicated
     * @param[in] surfaces the set of surfaces and sides around the corner
     * @param[in,out] info information on the behavior to adapt according the current surface
     * @return Returns true if the corner needs to be duplicated, false otherwise
     */
    bool MacroMeshVertices::duplicate_corner(
        const std::set< surface_side >& surfaces,
        std::vector< SurfaceAction >& info )
    {
        std::vector< SurfaceAction > temp_info( info.size(), TO_PROCESS ) ;
        for( std::set< surface_side >::const_iterator it( surfaces.begin() );
            it != surfaces.end(); ++it ) {
            index_t surface_id = it->first ;
            if( info[surface_id] == SKIP || temp_info[surface_id] == SKIP )
                continue ;
            if( temp_info[surface_id] == TO_PROCESS ) {
                temp_info[surface_id] = it->second ;
            } else {
                if( temp_info[surface_id] != it->second ) {
                    // Free border
                    temp_info[surface_id] = SKIP ;
                }
            }
        }

        bool result = false ;
        for( index_t s = 0; s < info.size(); s++ ) {
            if( temp_info[s] < 0 ) continue ;
            ringmesh_debug_assert( info[s] != SKIP ) ;
            if( info[s] == TO_PROCESS ) {
                info[s] = SurfaceAction( !temp_info[s] ) ;
            } else {
                if( info[s] == temp_info[s] ) {
                    result = true ;
                }
            }
        }

        return result ;
    }

    /*!
     * Tests is the surface should be duplicate according the DuplicateMode
     * in the MacroMesh
     * @param[in] surface_id id of the surface to test
     * @return Returns true if the surface should be duplicated
     */
    bool MacroMeshVertices::is_surface_to_duplicate( index_t surface_id ) const
    {
        BoundaryModelElement::GEOL_FEATURE feature = mm_.model().surface(
            surface_id ).geological_feature() ;
        if( mm_.duplicate_mode() == ALL &&
            !mm_.model().surface( surface_id ).is_on_voi() 
           )
            return true ;
        if( mm_.duplicate_mode() == FAULT && BME::is_fault( feature ) )
            return true ;

        return false ;
    }

    /*!
     * Gets the number of vertices with different coordinates,
     * if several vertices are collocated they are only count once
     * (do no care about duplication or interface between GEO::Mesh)
     * @return the corresponding number
     */
    index_t MacroMeshVertices::nb_vertices() const
    {
        if( unique_vertices_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return unique_vertices_.size() ;
    }

    /*!
     * Gets the vertex id in the MacroMesh corresponding to a GEO::Mesh vertex id
     * @param[in] mesh id of the mesh/region
     * @param[in] v vertex id of the GEO::Mesh vertex
     * @return the id of \p v in the MacroMesh
     */
    index_t MacroMeshVertices::vertex_id( index_t mesh, index_t v ) const
    {
        if( unique_vertices_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        ringmesh_debug_assert( v < mm_.mesh( mesh ).vertices.nb() ) ;
        return global_vertex_indices_[vertex2mesh_[mesh] + v] ;
    }

    /*!
     * Gets the MacroMesh vertex coordinates
     * @param[in] global_v vertex id in the MacroMesh
     * @return the vertex coordinates
     */
    const vec3& MacroMeshVertices::vertex( index_t global_v ) const
    {
        if( unique_vertices_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize() ;
        }
        return unique_vertices_[global_v] ;
    }

    /*!
     * Gets the MacroMesh vertex coordinates
     * @param[in] mesh vertex id in the MacroMesh
     * @param[in] v vertex id of the GEO::Mesh vertex
     * @return the vertex coordinates
     */
    const vec3& MacroMeshVertices::vertex( index_t mesh, index_t v ) const
    {
        return vertex( vertex_id( mesh, v ) ) ;
    }

    const vec3& MacroMeshVertices::duplicated_vertex( index_t v ) const
    {
        if( cell_corners_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize_duplication() ;
        }
        return unique_vertices_[duplicated_vertex_indices_[v]] ;
    }

    /*!
     * Given a cell_corner id in a given mesh, return the vertex id in the MacroMesh
     * and the corresponding duplicated vertex id if the vertex is duplicated
     * @param[in] mesh id of the mesh
     * @param[in] cell_corner id of the cell corner
     * @param[out] vertex_id vertex id in the MacroMesh of the corresponding point
     * @param[out] duplicated_vertex_id duplicated vertex id if the vertex is duplicated
     * @return returns false if the vertex is duplicated (\p duplicated_vertex_id is filled),
     * true otherwise
     */
    bool MacroMeshVertices::vertex_id(
        index_t mesh,
        index_t cell_corner,
        index_t& vertex_id,
        index_t& duplicated_vertex_id ) const
    {
        if( cell_corners_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize_duplication() ;
        }
        index_t corner_value = cell_corners_[mesh_cell_corner_ptr_[mesh]
            + cell_corner] ;
        if( corner_value < mm_.vertices.nb_vertices() ) {
            vertex_id = mm_.vertices.vertex_id( mesh, corner_value ) ;
            return true ;
        } else {
            duplicated_vertex_id = corner_value - mm_.vertices.nb_vertices() ;
            vertex_id = duplicated_vertex_indices_[duplicated_vertex_id] ;
            return false ;
        }
    }

    /*!
     * Gets the number of duplicated points by the DuplicateMode algorithm
     * @return the corresponding number of duplications
     */
    index_t MacroMeshVertices::nb_duplicated_vertices() const
    {
        if( cell_corners_.empty() ) {
            const_cast< MacroMeshVertices* >( this )->initialize_duplication() ;
        }
        return duplicated_vertex_indices_.size() ;
    }

    /*!
     * Gets the total number of vertices (nb_vertices() + nb_duplicated_vertices())
     * @return the corresponding number of vertices
     */
    index_t MacroMeshVertices::nb_total_vertices() const
    {
        return nb_duplicated_vertices() + mm_.vertices.nb_vertices() ;
    }

    /*!
     * Clear the MacroMeshVertices database
     */
    void MacroMeshVertices::clear()
    {
        unique_vertices_.clear() ;
        global_vertex_indices_.clear() ;
        vertex2mesh_.clear() ;
        cell_corners_.clear() ;
        mesh_cell_corner_ptr_.clear() ;
        duplicated_vertex_indices_.clear() ;
    }

    /*!
     * Initialize the facet database of the MacroMesh
     */
    void MacroMeshFacets::initialize()
    {
        index_t facet_access[5] = { -1, -1, -1, 0, 1 } ;
        surface2mesh_.resize( mm_.model().nb_surfaces(), Surface::NO_ID ) ;
        surface_facet_ptr_.resize( NB_FACET_TYPES * mm_.model().nb_surfaces() + 1,
            0 ) ;

        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
            std::vector< index_t > surface_proccessed ;
            GEO::Attribute< index_t > attribute( cur_mesh.facets.attributes(),
                surface_att_name ) ;
            for( index_t f = 0; f < cur_mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
                if( surface2mesh_[surface_id] != Surface::NO_ID ) continue ;
                if( !Utils::contains( surface_proccessed, surface_id ) ) {
                    surface_proccessed.push_back( surface_id ) ;
                }
                surface_facet_ptr_[NB_FACET_TYPES * surface_id
                    + facet_access[cur_mesh.facets.nb_vertices( f )] + 1]++ ;
            }
            for( index_t s = 0; s < surface_proccessed.size(); s++ ) {
                surface2mesh_[surface_proccessed[s]] = m ;
            }
        }

        for( index_t s = 1; s < surface_facet_ptr_.size() - 1; s++ ) {
            surface_facet_ptr_[s + 1] += surface_facet_ptr_[s] ;
        }
        surface_facets_.resize( surface_facet_ptr_.back() ) ;

        std::vector< index_t > cur_facet_index_type(
            NB_FACET_TYPES * mm_.model().nb_surfaces(), 0 ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
            GEO::Attribute< index_t > attribute( cur_mesh.facets.attributes(),
                surface_att_name ) ;
            for( index_t f = 0; f < cur_mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
                if( surface2mesh_[surface_id] != m ) continue ;
                index_t type_access = facet_access[cur_mesh.facets.nb_vertices( f )] ;
                surface_facets_[surface_facet_ptr_[NB_FACET_TYPES * surface_id
                    + type_access]
                    + cur_facet_index_type[NB_FACET_TYPES * surface_id + type_access]++ ] =
                    f ;
            }
        }

        for( index_t s = 0; s < mm_.model().nb_surfaces(); s++ ) {
            nb_triangle_ += nb_triangle( s ) ;
            nb_quad_ += nb_quad( s ) ;
        }
        nb_facets_ = nb_triangle_ + nb_quad_ ;
    }

    /*!
     * Gets the total number of triangles in the MacroMesh
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_facets() const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return nb_facets_ ;
    }

    /*!
     * Gets the total number of triangles in the MacroMesh
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_triangle() const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return nb_triangle_ ;
    }

    /*!
     * Gets the number of triangles in the given surface
     * @param[in] s id of the surface
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_triangle( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface_facet_ptr_[NB_FACET_TYPES * s + 1]
            - surface_facet_ptr_[NB_FACET_TYPES * s] ;
    }

    /*!
     * Gets the triangle id in the MacroMesh
     * @param[in] s id of the surface
     * @param[in] t id of the triangle
     * @return the corresponding facet id
     */
    index_t MacroMeshFacets::triangle_id( index_t s, index_t t ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return facet( surface_facet_ptr_[NB_FACET_TYPES * s] + t ) ;
    }

    /*!
     * Gets the total number of quad in the MacroMesh
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_quad() const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return nb_quad_ ;
    }

    /*!
     * Gets the number of quad in the given surface
     * @param[in] s id of the surface
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_quad( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface_facet_ptr_[NB_FACET_TYPES * s + 2]
            - surface_facet_ptr_[NB_FACET_TYPES * s + 1] ;
    }

    /*!
     * Gets the quad id in the MacroMesh
     * @param[in] s id of the surface
     * @param[in] q id of the quad
     * @return the corresponding facet id
     */
    index_t MacroMeshFacets::quad_id( index_t s, index_t q ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return facet( surface_facet_ptr_[NB_FACET_TYPES * s + 1] + q ) ;
    }

    /*!
     * Gets the mesh associated to the surface id
     * @param[in] s id of the surface
     * @return id of the mesh that stores the surface \p s
     */
    index_t MacroMeshFacets::mesh( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface2mesh_[s] ;
    }

    /*!
     * Gets the GEO::Mesh facet id
     * @param[in] s id of the surface
     * @param[in] f facet id in the surface
     * @return the facet id of the corresponding facet in the GEO::Mesh
     */
    index_t MacroMeshFacets::facet( index_t s, index_t f ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return facet( surface_begin( s ) + f ) ;
    }

    /*!
     * Gets the number of facets in a surface
     * @param[in] s id of the surface
     * @return the corresponding number
     */
    index_t MacroMeshFacets::nb_facets( index_t s ) const
    {
        if( surface_facets_.empty() ) {
            const_cast< MacroMeshFacets* >( this )->initialize() ;
        }
        return surface_end( s ) - surface_begin( s ) ;
    }

    /*!
     * Clear the MacroMeshFacets database
     */
    void MacroMeshFacets::clear()
    {
        surface_facets_.clear() ;
        surface_facet_ptr_.clear() ;
        surface2mesh_.clear() ;
        nb_facets_ = 0 ;
        nb_triangle_ = 0 ;
        nb_quad_ = 0 ;
    }

    /*!
     * Initialize the cell database of the MacroMesh
     */
    void MacroMeshCells::initialize()
    {
        mesh_cell_ptr_.resize( NB_CELL_TYPES * mm_.nb_meshes() + 1, 0 ) ;
        index_t total_adjacents = 0 ;
        index_t cell_access[4] = { 0, 3, 2, 1 } ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            nb_cells_ += mesh.cells.nb() ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                mesh_cell_ptr_[NB_CELL_TYPES * m + cell_access[mesh.cells.type( c )]
                    + 1]++ ;
                total_adjacents += mesh.cells.nb_facets( c ) ;
            }
        }
        for( index_t m = 1; m < mesh_cell_ptr_.size() - 1; m++ ) {
            mesh_cell_ptr_[m + 1] += mesh_cell_ptr_[m] ;
        }
        cells_.resize( mesh_cell_ptr_.back() ) ;

        index_t nb_vertices = mm_.vertices.nb_vertices() ;
        std::vector< std::vector< index_t > > cells_around_vertex( nb_vertices ) ;
        cell_adjacents_.reserve( total_adjacents ) ;
        std::vector< index_t > cur_cell_index_type( NB_CELL_TYPES * mm_.nb_meshes(),
            0 ) ;
        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                index_t type_access = cell_access[mesh.cells.type( c )] ;
                cells_[mesh_cell_ptr_[NB_CELL_TYPES * m + type_access]
                    + cur_cell_index_type[NB_CELL_TYPES * m + type_access]++ ] = c ;

                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t adj = mesh.cells.adjacent( c, f ) ;
                    if( adj != GEO::NO_CELL ) {
                        adj += mesh_cell_ptr_[NB_CELL_TYPES * m] ;
                    } else {
                        for( index_t v = 0; v < mesh.cells.facet_nb_vertices( c, f );
                            v++ ) {
                            index_t vertex_id = mm_.vertices.vertex_id( m,
                                mesh.cells.facet_vertex( c, f, v ) ) ;
                            cells_around_vertex[vertex_id].push_back(
                                mesh_cell_ptr_[NB_CELL_TYPES * m] + c ) ;
                        }
                    }
                    cell_adjacents_.push_back( adj ) ;
                }
            }
        }

        for( index_t v = 0; v < cells_around_vertex.size(); v++ ) {
            GEO::sort_unique( cells_around_vertex[v] ) ;
        }

        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            const GEO::Mesh& mesh = mm_.mesh( m ) ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t adj = mesh.cells.adjacent( c, f ) ;
                    if( adj == GEO::NO_CELL ) {
                        index_t prev_vertex_id = mm_.vertices.vertex_id( m,
                            mesh.cells.facet_vertex( c, f, 0 ) ) ;
                        std::vector< index_t >& prev_cells =
                            cells_around_vertex[prev_vertex_id] ;
                        index_t size_hint = prev_cells.size() ;
                        std::vector< index_t > intersection( size_hint ) ;
                        for( index_t v = 1; v < mesh.cells.facet_nb_vertices( c, f );
                            v++ ) {
                            index_t vertex_id = mm_.vertices.vertex_id( m,
                                mesh.cells.facet_vertex( c, f, v ) ) ;
                            std::vector< index_t >& cells =
                                cells_around_vertex[vertex_id] ;
                            intersection.erase(
                                std::set_intersection( prev_cells.begin(),
                                    prev_cells.end(), cells.begin(), cells.end(),
                                    intersection.begin() ), intersection.end() ) ;
                            if( intersection.size() < 2 ) break ;
                            prev_cells = intersection ;
                        }

                        if( intersection.size() > 1 ) {
                            ringmesh_debug_assert( intersection.size() == 2 ) ;
                            index_t new_adj =
                                intersection[0]
                                    == mesh_cell_ptr_[NB_CELL_TYPES * m] + c ?
                                    intersection[1] : intersection[0] ;
                            cell_adjacents_[mesh_cell_ptr_[NB_CELL_TYPES * m]
                                + mesh.cells.facets_begin( c ) + f] = new_adj ;
                        }
                    }
                }
            }
        }

        for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
            nb_tet_ += nb_tet( m ) ;
            nb_pyramid_ += nb_pyramid( m ) ;
            nb_prism_ += nb_prism( m ) ;
            nb_hex_ += nb_hex( m ) ;
        }
    }

    /*!
     * Gets the cell adjacent id in the MacroMesh
     * @param[in] m id of the GEO::Mesh
     * @param[in] c id of the cell in the GEO::Mesh
     * @param[in] f id of the facet in the cell to test
     * @return the cell adjacent id or GEO::NO_ADJACENT
     */
    index_t MacroMeshCells::cell_adjacent( index_t m, index_t c, index_t f ) const
    {
        test_initialize() ;
        return cell_adjacents_[mesh_begin( m )
            + mm_.mesh( m ).cells.facets_begin( c ) + f] ;
    }

    /*!
     * Convert a MacroMesh cell id to a GEO::Mesh cell id
     * @param[in] global_index id of the MacroMesh cell
     * @param[out] mesh_id id of the GEO::Mesh containing the cell
     * @return the cell id in the GEO::Mesh whose id is \p mesh_id
     */
    index_t MacroMeshCells::cell_index_in_mesh(
        index_t global_index,
        index_t& mesh_id ) const
    {
        test_initialize() ;
        mesh_id = 0 ;
        for( ; mesh_id < mm_.nb_meshes(); mesh_id++ ) {
            if( global_index < mesh_cell_ptr_[mesh_id + 1] ) break ;
        }
        ringmesh_debug_assert( mesh_id < mm_.nb_meshes() ) ;
        return global_index - mesh_cell_ptr_[mesh_id] ;
    }

    /*!
     * Gets the total number of cells
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_cells() const
    {
        test_initialize() ;
        return nb_cells_ ;
    }

    /*!
     * Gets the number of cells in a GEO::Mesh (prefers using mesh.cells.nb())
     * @param[in] r id of the GEO::Mesh
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_cells( index_t r ) const
    {
        test_initialize() ;
        return mesh_end( r ) - mesh_begin( r ) ;
    }

    /*!
     * Gets the total number of tetrahedra
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_tet() const
    {
        test_initialize() ;
        return nb_tet_ ;
    }

    /*!
     * Gets the number of tetrahedra in a GEO::Mesh
     * @param[in] r id of the GEO::Mesh
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_tet( index_t r ) const
    {
        test_initialize() ;
        return mesh_cell_ptr_[NB_CELL_TYPES * r + 1]
            - mesh_cell_ptr_[NB_CELL_TYPES * r] ;
    }

    /*!
     * Gets the tetrahedron id in the MacroMesh
     * @param[in] r id of the GEO::Mesh
     * @param[in] t id of the terahedron
     * @return the corresponding id
     */
    index_t MacroMeshCells::tet_id( index_t r, index_t t ) const
    {
        test_initialize() ;
        return cells_[mesh_cell_ptr_[NB_CELL_TYPES * r] + t] ;
    }

    /*!
     * Gets the total number of pyramids
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_pyramid() const
    {
        test_initialize() ;
        return nb_pyramid_ ;
    }

    /*!
     * Gets the number of pyramids in a GEO::Mesh
     * @param[in] r id of the GEO::Mesh
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_pyramid( index_t r ) const
    {
        test_initialize() ;
        return mesh_cell_ptr_[NB_CELL_TYPES * r + 2]
            - mesh_cell_ptr_[NB_CELL_TYPES * r + 1] ;
    }

    /*!
     * Gets the pyramid id in the MacroMesh
     * @param[in] r id of the GEO::Mesh
     * @param[in] p id of the pyramid
     * @return the corresponding id
     */
    index_t MacroMeshCells::pyramid_id( index_t r, index_t p ) const
    {
        test_initialize() ;
        return cells_[mesh_cell_ptr_[NB_CELL_TYPES * r + 1] + p] ;
    }

    /*!
     * Gets the total number of prisms
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_prism() const
    {
        test_initialize() ;
        return nb_prism_ ;
    }

    /*!
     * Gets the number of prisms in a GEO::Mesh
     * @param[in] r id of the GEO::Mesh
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_prism( index_t r ) const
    {
        test_initialize() ;
        return mesh_cell_ptr_[NB_CELL_TYPES * r + 3]
            - mesh_cell_ptr_[NB_CELL_TYPES * r + 2] ;
    }

    /*!
     * Gets the prism id in the MacroMesh
     * @param[in] r id of the GEO::Mesh
     * @param[in] p id of the prism
     * @return the corresponding id
     */
    index_t MacroMeshCells::prism_id( index_t r, index_t p ) const
    {
        test_initialize() ;
        return cells_[mesh_cell_ptr_[NB_CELL_TYPES * r + 2] + p] ;
    }

    /*!
     * Gets the total number of hexahedra
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_hex() const
    {
        test_initialize() ;
        return nb_hex_ ;
    }

    /*!
     * Gets the number of hexahedra in a GEO::Mesh
     * @param[in] r id of the GEO::Mesh
     * @return the corresponding number
     */
    index_t MacroMeshCells::nb_hex( index_t r ) const
    {
        test_initialize() ;
        return mesh_cell_ptr_[NB_CELL_TYPES * r + 4]
            - mesh_cell_ptr_[NB_CELL_TYPES * r + 3] ;
    }

    /*!
     * Gets the hexahedron id in the MacroMesh
     * @param[in] r id of the GEO::Mesh
     * @param[in] h id of the hexahedron
     * @return the corresponding id
     */
    index_t MacroMeshCells::hex_id( index_t r, index_t h ) const
    {
        test_initialize() ;
        return cells_[mesh_cell_ptr_[NB_CELL_TYPES * r + 3] + h] ;
    }

    /*!
     * Clear the MacroMeshCells database
     */
    void MacroMeshCells::clear()
    {
        cells_.clear() ;
        cell_adjacents_.clear() ;
        mesh_cell_ptr_.clear() ;
        nb_cells_ = 0 ;
        nb_tet_ = 0 ;
        nb_pyramid_ = 0 ;
        nb_prism_ = 0 ;
        nb_hex_ = 0 ;
    }

    MacroMeshTools::MacroMeshTools( MacroMesh& mm )
        :
            mm_( mm )
    {
    }

    MacroMeshTools::~MacroMeshTools()
    {
        for( index_t r = 0; r < facet_aabb_.size(); r++ ) {
            if( facet_aabb_[r] ) delete facet_aabb_[r] ;
        }
        for( index_t r = 0; r < cell_aabb_.size(); r++ ) {
            if( cell_aabb_[r] ) delete cell_aabb_[r] ;
        }
    }

    /*!
     * Gets the GEO::MeshFacetsAABB of the given region
     * @param[in] region id of the GEO::Mesh
     * @return the const reference to the corresponding MeshFacetsAABB
     */
    const GEO::MeshFacetsAABB& MacroMeshTools::facet_aabb( index_t region ) const
    {
        init_facet_aabb( region ) ;
        return *facet_aabb_[region] ;
    }

    /*!
     * Initialize if needed the MeshFacetsAABB of the given region
     * @param[in] region id of the GEO::Mesh
     */
    void MacroMeshTools::init_facet_aabb( index_t region ) const
    {
        if( facet_aabb_.size() <= region ) {
            const_cast< MacroMeshTools* >( this )->facet_aabb_.resize( region + 1,
                nil ) ;
        }
        if( facet_aabb_[region] ) return ;
        const_cast< MacroMeshTools* >( this )->facet_aabb_[region] =
            new GEO::MeshFacetsAABB( mm_.mesh( region ) ) ;
    }

    /*!
     * Gets the GEO::MeshTetsAABB of the given region
     * @param[in] region id of the GEO::Mesh
     * @return the const reference to the corresponding MeshTetsAABB
     */
    const GEO::MeshCellsAABB& MacroMeshTools::cell_aabb( index_t region ) const
    {
        init_cell_aabb( region ) ;
        return *cell_aabb_[region] ;
    }

    /*!
     * Initialize if needed the MeshTetsAABB of the given region
     * @param[in] region id of the GEO::Mesh
     */
    void MacroMeshTools::init_cell_aabb( index_t region ) const
    {
        if( cell_aabb_.size() <= region ) {
            const_cast< MacroMeshTools* >( this )->cell_aabb_.resize( region + 1,
                nil ) ;
        }
        if( cell_aabb_[region] ) return ;
        const_cast< MacroMeshTools* >( this )->cell_aabb_[region] =
            new GEO::MeshCellsAABB( mm_.mesh( region ) ) ;
    }

    /*!
     * Clear the MacroMeshTools database
     */
    void MacroMeshTools::clear()
    {
        for( index_t r = 0; r < facet_aabb_.size(); r++ ) {
            if( facet_aabb_[r] ) {
                delete facet_aabb_[r] ;
                facet_aabb_[r] = nil ;
            }
        }
        for( index_t r = 0; r < cell_aabb_.size(); r++ ) {
            if( cell_aabb_[r] ) {
                delete cell_aabb_[r] ;
                cell_aabb_[r] = nil ;
            }
        }
    }

MacroMeshOrder::MacroMeshOrder( MacroMesh& mm )
        : mm_( mm ), nb_vertices_( 0 ), points_( 0 )
    {

    }

    MacroMeshOrder::~MacroMeshOrder()
    {

    }

    /*
     * Initialize the database by computing the new vertices of the mesh.
     * \param order -1 vertices are added per edges, the edges are divided
     * in equal parts by these vertices.
     * @param[in] order the mesh elements order
     */
    void MacroMeshOrder::initialize()
    {
        index_t offset = 0 ;
        nb_vertices_ = mm_.vertices.nb_total_vertices() ;
        index_t order = mm_.get_order() ;
        if( order != 1 ) {

            index_t nb_total_edges = 0 ;
            for( index_t r = 0; r < mm_.nb_meshes(); r++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( r ) ;
                for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                    for( index_t e = 0; e < cur_mesh.cells.nb_edges( c ); e++ ) {
                        nb_total_edges++ ;
                    }
                }
            }

            std::vector< vec3 > new_points( nb_total_edges * ( order - 1 ) ) ;
            for( index_t r = 0; r < mm_.nb_meshes(); r++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( r ) ;
                GEO::Attribute< std::vector< index_t > > order_vertices(
                    cur_mesh.cells.attributes(), "order_vertices" ) ;
                for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                    std::vector< index_t > cur_order_vertices(
                        cur_mesh.cells.nb_edges( c ) * ( order - 1 ) ) ;
                    for( index_t e = 0; e < cur_mesh.cells.nb_edges( c ); e++ ) {
                        std::vector< vec3 > new_points_in_edge ;
                        vec3 node0 = GEO::Geom::mesh_vertex( cur_mesh,
                            cur_mesh.cells.edge_vertex( c, e, 0 ) ) ;
                        vec3 node1 = GEO::Geom::mesh_vertex( cur_mesh,
                            cur_mesh.cells.edge_vertex( c, e, 1 ) ) ;
                        Geom::divide_edge_in_parts( node0, node1, order,
                            new_points_in_edge ) ;

                        for( index_t v = 0; v < new_points_in_edge.size(); v++ ) {
                            new_points[offset] = new_points_in_edge[v] ;
                            cur_order_vertices[e + v] = offset ;
                            offset++ ;
                        }
                        order_vertices[c] = cur_order_vertices ;
                    }
                }
            }

            MakeUnique uniq( new_points ) ;
            uniq.unique() ;
            std::vector< vec3 > uniq_points ;
            uniq.unique_points( uniq_points ) ;
            std::vector< index_t > map = uniq.indices() ;

            ColocaterANN ann( uniq_points ) ;
            for( index_t r = 0; r < mm_.nb_meshes(); r++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( r ) ;
                GEO::Attribute< std::vector< index_t > > order_vertices(
                    cur_mesh.cells.attributes(), "order_vertices" ) ;
                for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                    for( index_t v = 0; v < order_vertices[c].size(); v++ ) {
                        order_vertices[c][v] = map[order_vertices[c][v]]
                            + nb_vertices_ ;
                    }
                }
            }

            for( index_t s = 0; s < mm_.model().nb_surfaces(); s++ ) {
                index_t cur_mesh_id = mm_.facets.mesh( s ) ;
                const GEO::Mesh& cur_mesh = mm_.mesh( cur_mesh_id ) ;
                GEO::Attribute< std::vector< index_t > > order_vertices(
                    cur_mesh.facets.attributes(), "order_vertices" ) ;
                for( index_t f = 0; f < mm_.facets.nb_facets( s ); f++ ) {
                    index_t cur_facet = mm_.facets.facet( s, f ) ;
                    std::vector< index_t > cur_order_vertices(
                        cur_mesh.facets.nb_vertices( f ) * ( order - 1 ) ) ;
                    for( index_t e = 0; e < cur_mesh.facets.nb_vertices( cur_facet );
                        e++ ) {
                        vec3 node0 ;
                        vec3 node1 ;
                        std::vector< vec3 > new_points_in_edge ;
                        if( e == cur_mesh.facets.nb_vertices( cur_facet ) - 1 ) {
                            node0 = GEO::Geom::mesh_vertex( cur_mesh,
                                cur_mesh.facets.vertex( cur_facet, e ) ) ;
                            node1 = GEO::Geom::mesh_vertex( cur_mesh,
                                cur_mesh.facets.vertex( cur_facet, 0 ) ) ;
                        } else {
                            node0 = GEO::Geom::mesh_vertex( cur_mesh,
                                cur_mesh.facets.vertex( cur_facet, e ) ) ;
                            node1 = GEO::Geom::mesh_vertex( cur_mesh,
                                cur_mesh.facets.vertex( cur_facet, e + 1 ) ) ;
                        }
                        RINGMesh::Geom::divide_edge_in_parts( node0, node1, order,
                            new_points_in_edge ) ;
                        for( index_t v = 0; v < new_points_in_edge.size(); v++ ) {
                            std::vector< index_t > colocated_vertices ;
                            index_t real_vertex_id = ann.get_colocated(
                                new_points_in_edge[v], colocated_vertices ) ;
                            ringmesh_debug_assert( colocated_vertices.size() == 1 ) ;
                            cur_order_vertices[e + v] = colocated_vertices[0]
                                + nb_vertices_ ;
                        }
                    }
                    order_vertices[cur_facet] = cur_order_vertices ;
                }
            }
            nb_vertices_ += uniq_points.size() ;
        }

    }

    /*
     * Clear the MacroMeshOrder database
     */
    void MacroMeshOrder::clear()
    {
        nb_vertices_ = 0 ;
    }

    /*
     * Gets the mesh total number of vertices. It is the number of unique nodes
     * on the mesh plus the added nodes on the elements edges
     * @return the const number of vertices
     */
    const index_t MacroMeshOrder::nb_total_vertices() const
    {
        test_initialize() ;
        return nb_vertices_ ;
    }

//    /*
//     * Gets the id of the added node
//     * @return the const id of the node
//     */
//    const index_t MacroMeshOrder::id( const vec3& point ) const
//    {
//
//        test_initialize() ;
//        std::vector< index_t > colocated_points ;
//        ann_.get_colocated( point, colocated_points ) ;
//        ringmesh_debug_assert( colocated_points.size() == 1 ) ;
//        return mm_.vertices.nb_total_vertices() + colocated_points[0] ;
//    }

    const index_t MacroMeshOrder::nb_vertices() const
    {
        test_initialize() ;
        return nb_vertices_ - mm_.vertices.nb_total_vertices() ;
    }

    /*
     * Gets the vec3 of a added point
     * @param[in] id an id of the new created point for order > 2
     * @return the vec3 matching with the id
     */
    const vec3 MacroMeshOrder::point( const index_t id ) const
    {
        const_cast< MacroMeshOrder* >( this )->test_point_list_initialize() ;
        return points_[id] ;
    }

    /*
     * Move a added point
     * @param[in] id the id of the point
     * @param[in] u the displacement applied on this point
     */
    void MacroMeshOrder::move_point( const index_t id, const vec3& u )
    {
        const_cast< MacroMeshOrder* >( this )->test_point_list_initialize() ;
        for( index_t i = 0; i < 3; i++ ) {
            points_[id][i] += u[i] ;

        }
    }

    /*
     * Test wether the vec3 point list is initialize. If note, the point
     * list is initialize
     */
    void MacroMeshOrder::test_point_list_initialize()
    {
        index_t order = mm_.get_order() ;

        if( points_.size() == 0 && order > 1 ) {
            index_t offset = 0 ;
            index_t nb_total_edges = 0 ;
            for( index_t r = 0; r < mm_.nb_meshes(); r++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( r ) ;
                for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                    for( index_t e = 0; e < cur_mesh.cells.nb_edges( c ); e++ ) {
                        nb_total_edges++ ;
                    }
                }
            }

            std::vector< vec3 > new_points( nb_total_edges * ( order - 1 ) ) ;
            for( index_t r = 0; r < mm_.nb_meshes(); r++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( r ) ;
                for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                    std::vector< index_t > cur_order_vertices(
                        cur_mesh.cells.nb_edges( c ) * ( order - 1 ) ) ;
                    for( index_t e = 0; e < cur_mesh.cells.nb_edges( c ); e++ ) {
                        std::vector< vec3 > new_points_in_edge ;
                        vec3 node0 = GEO::Geom::mesh_vertex( cur_mesh,
                            cur_mesh.cells.edge_vertex( c, e, 0 ) ) ;
                        vec3 node1 = GEO::Geom::mesh_vertex( cur_mesh,
                            cur_mesh.cells.edge_vertex( c, e, 1 ) ) ;
                        Geom::divide_edge_in_parts( node0, node1, order,
                            new_points_in_edge ) ;

                        for( index_t v = 0; v < new_points_in_edge.size(); v++ ) {
                            new_points[offset] = new_points_in_edge[v] ;
                            offset++ ;
                        }
                    }
                }
            }

            MakeUnique uniq( new_points ) ;
            uniq.unique() ;
            uniq.unique_points( points_ ) ;
        }
    }


    MacroMesh::MacroMesh( const BoundaryModel& model )
        :
            model_( &model ),
            meshes_( model.nb_regions(), nil ),
            mode_( NONE ),
            wells_( nil ),
            order_(1),
            vertices( *this ),
            facets( *this ),
            cells( *this ),
            tools( *this ),
            order( *this )
    {
        for( index_t r = 0; r < model_->nb_regions(); r++ ) {
            meshes_[r] = new GEO::Mesh( 3 ) ;
        }
    }

    MacroMesh::MacroMesh()
        :
            model_( nil ),
            meshes_(),
            mode_( NONE ),
            wells_( nil ),
            order_(1),
            vertices( *this ),
            facets( *this ),
            cells( *this ),
            tools( *this ),
            order( *this )
    {
    }

    /*!
     * Copy a MacroMesh and its attributes
     * @param[in] rhs the MacroMesh copied
     * @param[in] keep_order To document
     * @param[in] copy_attributes tells whether or not you want to copy attributes
     * @todo Check this function - keep_order parameter is not used
     */
    void MacroMesh::copy( const MacroMesh& rhs, bool copy_attributes )
    {
        index_t dim = meshes_[0]->vertices.dimension() ;

        model_ = &rhs.model() ;
        order_ = rhs.get_order() ;
        mode_= rhs.duplicate_mode() ;
        wells_ = rhs.wells() ;
        for( index_t r = 0; r < model_->nb_regions(); r++ ) {
            meshes_[r]->copy( *rhs.meshes_[r], copy_attributes ) ;
        }

    }

    MacroMesh::~MacroMesh()
    {
        for( index_t r = 0; r < meshes_.size(); r++ ) {
#ifdef RINGMESH_DEBUG
            Utils::print_bounded_attributes( *meshes_[r] ) ;
#endif
            delete meshes_[r] ;
        }
    }

    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param[in] method Mesher used
     * @param[in] region_id Region to mesh, -1 for all
     * @param[in] add_steiner_points if true, the mesher will add some points inside the region
     * @param[out] internal_vertices TO DOCUMENT
     * to improve the mesh quality
     */
    void MacroMesh::compute_tetmesh(
        const std::string& method,
        int region_id,
        bool add_steiner_points,
        std::vector< std::vector< vec3 > >& internal_vertices )
    {
        GEO::Logger::out( "Info" ) << "Using " << method << std::endl ;
        if( region_id == -1 ) {
            GEO::ProgressTask progress( "Compute", nb_meshes() ) ;
            for( index_t i = 0; i < nb_meshes(); i++ ) {
                TetraGen_var tetragen = TetraGen::create( mesh( i ), method ) ;
                tetragen->set_boundaries( &model_->region( i ), wells() ) ;
                if( !internal_vertices.empty() ) {
                    tetragen->set_internal_points( internal_vertices[i] ) ;
                }
                GEO::Logger::instance()->set_quiet( true ) ;
                tetragen->tetrahedralize( add_steiner_points ) ;
                GEO::Logger::instance()->set_quiet( false ) ;
                progress.next() ;
            }
        } else {
            TetraGen_var tetragen = TetraGen::create( mesh( region_id ), method ) ;
            tetragen->set_boundaries( &model_->region( region_id ), wells() ) ;
            if( !internal_vertices.empty() ) {
                tetragen->set_internal_points( internal_vertices[region_id] ) ;
            }
            GEO::Logger::instance()->set_quiet( true ) ;
            tetragen->tetrahedralize( add_steiner_points ) ;
            GEO::Logger::instance()->set_quiet( false ) ;
        }
    }

    /*!
     * Associates a WellGroup to the MacroMesh
     * @param[in] wells the WellGroup
     */
    void MacroMesh::add_wells( const WellGroup* wells )
    {
        wells_ = wells ;
    }

    void MacroMesh::set_nodel( const BoundaryModel& model )
    {
        model_ = &model ;
        meshes_.resize( model_->nb_regions(), nil ) ;
        for( index_t r = 0; r < model_->nb_regions(); r++ ) {
            meshes_[r] = new GEO::Mesh( 3 ) ;
        }
    }
}

