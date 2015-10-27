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
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/

#include <ringmesh/geo_model_repair.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_editor.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/geometry.h>

#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/basic/logger.h>

#include <algorithm>

namespace {

    using namespace RINGMesh ;

    typedef GeoModelElement::gme_t gme_t ;
    typedef GeoModelMeshElement BMME ;
    typedef GeoModelMeshVertices::VertexInGME VGME ;

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
    *
    * \brief Tests whether a facet is degenerate.
    * \param[in] M the mesh that the facet belongs to
    * \param[in] f the index of the facet in \p M
    * \return true if facet \p f has duplicated vertices,
    *  false otherwise
    */
    bool facet_is_degenerate(
        const GEO::Mesh& M,
        index_t f,
        GEO::vector< index_t >& colocated_vertices )
    {
        index_t nb_vertices = M.facets.nb_vertices( f ) ;
        if( nb_vertices != 3 ) {
            index_t* vertices = (index_t*)alloca( nb_vertices * sizeof( index_t ) ) ;
            for( index_t lv = 0; lv < nb_vertices; ++lv ) {
                vertices[ lv ] = colocated_vertices[ M.facets.vertex( f, lv ) ] ;
            }
            std::sort( vertices, vertices + nb_vertices ) ;
            return std::unique( vertices, vertices + nb_vertices )
                != vertices + nb_vertices ;
        }
        index_t c1 = M.facets.corners_begin( f ) ;
        index_t c2 = c1 + 1 ;
        index_t c3 = c2 + 1 ;
        index_t v1 = colocated_vertices[ M.facet_corners.vertex( c1 ) ] ;
        index_t v2 = colocated_vertices[ M.facet_corners.vertex( c2 ) ] ;
        index_t v3 = colocated_vertices[ M.facet_corners.vertex( c3 ) ] ;
        return v1 == v2 || v2 == v3 || v3 == v1 ;
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
    */
    void mesh_detect_degenerate_facets(
        const GEO::Mesh& M,
        GEO::vector< index_t >& f_is_degenerate,
        GEO::vector< index_t >& colocated_vertices )
    {
        f_is_degenerate.resize( M.facets.nb() ) ;
        for( index_t f = 0; f < M.facets.nb(); ++f ) {
            f_is_degenerate[ f ] = facet_is_degenerate( M, f, colocated_vertices ) ;
        }
    }

    /*!
    * @brief Detect and remove degenerated facets in a Mesh
    */
    index_t detect_degenerate_facets( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_facets( M, degenerate, colocated ) ;
        return static_cast< index_t> (
            std::count( degenerate.begin(), degenerate.end(), 1 ) ) ;
    }

    bool edge_is_degenerate(
        const GEO::Mesh& M,
        index_t e,
        GEO::vector< index_t >& colocated_vertices )
    {
        index_t v1 = colocated_vertices[ M.edges.vertex( e, 0 ) ] ;
        index_t v2 = colocated_vertices[ M.edges.vertex( e, 1 ) ] ;
        return v1 == v2 ;
    }

    void mesh_detect_degenerate_edges(
        const GEO::Mesh& M,
        GEO::vector< index_t >& e_is_degenerate,
        GEO::vector< index_t >& colocated_vertices )
    {
        e_is_degenerate.resize( M.edges.nb() ) ;
        for( index_t e = 0; e < M.edges.nb(); ++e ) {
            e_is_degenerate[ e ] = edge_is_degenerate( M, e, colocated_vertices ) ;
        }
    }

    /*!
    * @brief Detect and remove degenerated edges in a Mesh
    */
    index_t repair_line_mesh( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_edges( M, degenerate, colocated ) ;
        index_t nb = static_cast< index_t >( std::count( degenerate.begin(), degenerate.end(), 1 ) ) ;
        /// We have a problem if some vertices are left isolated
        /// If we ermove them here we can kill all indices correspondances
        M.edges.delete_elements( degenerate, false ) ;
        return nb ;
    }

    /*!
    * @brief Remove degenerate facets and edges from the Surface
    *        and Line of the model.
    * @pre colocated vertices have already been removed
    */
    void remove_degenerate_facet_and_edges(
        GeoModel& BM,
        std::set< gme_t >& to_remove )
    {
        to_remove.clear() ;
        for( index_t i = 0; i < BM.nb_lines(); ++i ) {
            index_t nb = repair_line_mesh( BM.line( i ).mesh() ) ;
            if( nb > 0 ) {

                GEO::Logger::out( "GeoModel" ) << nb
                    << " degenerated edges removed in LINE " << i << std::endl ;


                // The line may be empty now - remove it from the model
                if( BM.line( i ).nb_cells() == 0 ) {
                    to_remove.insert( BM.line( i ).gme_id() ) ;
                }
            }
        }

        // The builder might be needed
        GeoModelBuilder builder( BM ) ;

        for( index_t i = 0; i < BM.nb_surfaces(); ++i ) {
            GEO::Mesh& M = BM.surface( i ).mesh() ;
            index_t nb = detect_degenerate_facets( M ) ;
            /// @todo Check if that cannot be simplified 
            if( nb > 0 ) {
                // If there are some degenerated facets 
                // We need to repair the model 
//#ifndef RINGMESH_DEBUG
                //GEO::Logger::instance()->set_quiet( true ) ;
//#endif
                // Using repair function of geogram
                // Warning - This triangulates the mesh

                if( M.vertices.nb() > 0 ) {
                    // Colocated vertices must be processed before
                    // MESH_REPAIR_DUP_F 2 ;
                    GEO::MeshRepairMode mode =
                        static_cast< GEO::MeshRepairMode >( 2 ) ;
                    GEO::mesh_repair( M, mode ) ;

                    // This might create some small components - remove them
                    /// @todo How to choose the epsilon ? and the maximum number of facets ?
                    GEO::remove_small_connected_components( M, epsilon_sq, 3 ) ;

                    // This is a bit of an overkill but I am lazy
                    if( M.vertices.nb() > 0 ) {
                        GEO::mesh_repair( M, mode ) ;
                    }
                }
//#ifndef RINGMESH_DEBUG
                //GEO::Logger::instance()->set_quiet( false ) ;
//#endif 
                if( M.vertices.nb() == 0 || M.facets.nb() == 0 ) {
                    to_remove.insert( BM.surface( i ).gme_id() ) ;
                } else {
                    // If the Surface has internal boundaries, we need to 
                    // re-cut the Surface along these lines
                    Surface& S = const_cast< Surface& >( BM.surface( i ) ) ;
                    std::set< index_t > cutting_lines ;
                    for( index_t l = 0; l < S.nb_boundaries(); ++l ) {
                        const Line& L = BM.line( S.boundary_id( l ).index ) ;
                        if( to_remove.count( L.gme_id() ) == 0 &&
                            L.is_inside_border( S )
                            ) {
                            cutting_lines.insert( L.index() ) ;
                        }
                    }
                    for( std::set< index_t >::iterator it = cutting_lines.begin();
                         it != cutting_lines.end(); ++it
                         ) {
                        // Force the recomputing of the model vertices
                        // I do not understand exactly what is happening [JP]
                        BM.mesh.vertices.clear() ;                        
                        builder.cut_surface_by_line( S, BM.line( *it ) ) ;
                    }
                }
            }
        }
    }

    /*!
    * Get the indices of the duplicated vertices that are on an inside border.
    * Only the vertex with the biggest index are added.
    * If there are more than 2 colocated vertices throws an assertion in debug mode
    */
    void vertices_on_inside_boundary(
        const GeoModelMeshElement& E,
        std::set< index_t >& vertices )
    {
        vertices.clear() ;
        if( E.type() == GME::CORNER ) {
            return ;
        }
        if( E.type() == GME::LINE ) {
            if( E.boundary( 0 ).is_inside_border( E ) ) {
                vertices.insert( E.nb_vertices() - 1 ) ;
            }
            return ;
        }
        std::vector< const GeoModelMeshElement* > inside_border ;
        for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
            if( E.boundary( i ).is_inside_border( E ) ) {
                inside_border.push_back(
                    dynamic_cast< const GeoModelMeshElement* >( &E.boundary( i ) ) ) ;
            }
        }
        if( !inside_border.empty() ) {
            // We want to get the indices of the vertices in E
            // that are colocated with those of the inside boundary
            // We assume that the model vertice are not yet computed
            ColocaterANN kdtree( E.mesh(), ColocaterANN::VERTICES ) ;

            for( index_t i = 0; i < inside_border.size(); ++i ) {
                const GEO::Mesh& m = inside_border[ i ]->mesh() ;
                for( index_t v = 0; v < m.vertices.nb(); ++v ) {
                    std::vector< index_t > colocated_indices ;
                    kdtree.get_colocated( m.vertices.point( v ),
                                          colocated_indices ) ;
                    if( colocated_indices.size() > 1 ) {
                        std::sort( colocated_indices.begin(), colocated_indices.end() ) ;
                        // Add colocated vertices except one to the duplicated vertices set
                        vertices.insert( colocated_indices.begin()+1, colocated_indices.end() ) ;
                    }
                }
            }
        }
    }

    void remove_colocated_element_vertices(
        GeoModel& BM,
        std::set< gme_t >& to_remove )
    {
        to_remove.clear() ;
        // For all Lines and Surfaces
        for( index_t t = GME::LINE; t < GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;

            for( index_t e = 0; e < BM.nb_elements( T ); ++e ) {
                const BMME& E = dynamic_cast< const BMME& >( BM.element(
                    gme_t( T, e ) ) ) ;

                GEO::Mesh& M = E.mesh() ;
                GEO::vector< index_t > colocated ;
                GEO::mesh_detect_colocated_vertices( M, colocated, epsilon ) ;

                // Get the vertices to delete
                std::set< index_t > inside_border ;
                vertices_on_inside_boundary( E, inside_border ) ;

                GEO::vector< index_t > to_delete( colocated.size(), 0 ) ;
                GEO::vector< index_t > old2new( colocated.size() ) ;
                index_t nb_todelete = 0 ;
                index_t cur = 0 ;
                for( index_t v = 0; v < colocated.size(); ++v ) {
                    if( colocated[ v ] == v
                        || inside_border.find( v ) != inside_border.end() ) {
                        // This point is kept 
                        // No colocated or on an inside boundary
                        old2new[ v ] = cur ;
                        ++cur ;
                    } else {
                        // The point is to remove
                        old2new[ v ] = old2new[ colocated[ v ] ] ;
                        to_delete[ v ] = 1 ;
                        nb_todelete++ ;
                    }
                }

                if( nb_todelete == 0 ) {
                    // Nothing to do there
                    continue ;
                } else if( nb_todelete == E.nb_vertices() ) {
                    // The complete element should be removed
                    to_remove.insert( E.gme_id() ) ;
                    continue ;
                } else {
                /*    // We need to update the VertexInGME at the model level
                    // if they exist. If not let's avoid this costly operation
                    // @todo This is bugged ?? I am not sure. JP
                    if( BM.mesh.vertices.is_initialized() ) {
                        // For all the vertices of this element
                        // we need to update the vertex ids in the gme_vertices_
                        // of the corresponding global vertex
                        for( index_t v = 0; v < to_delete.size(); ++v ) {
                            index_t model_id = E.model_vertex_id( v ) ;
                            const std::vector< VGME >& cur =
                                BM.mesh.vertices.gme_vertices( model_id ) ;
                            for( index_t i = 0; i < cur.size(); ++i ) {
                                if( cur[ i ] == VGME( E.gme_id(), v ) ) {
                                    index_t new_id = old2new[ v ] ;
                                    BM.mesh.vertices.set_gme( model_id, i,
                                                              VGME( E.gme_id(), new_id ) ) ;
                                }
                            }
                        }
                        BM.mesh.erase_invalid_vertices() ;
                    } */

                    for( index_t c = 0; c < M.facet_corners.nb(); c++ ) {
                        M.facet_corners.set_vertex( c,
                                                    colocated[ M.facet_corners.vertex( c ) ] ) ;
                    }
                    for( index_t e = 0; e < M.edges.nb(); e++ ) {
                        M.edges.set_vertex( e, 0,
                                            colocated[ M.edges.vertex( e, 0 ) ] ) ;
                        M.edges.set_vertex( e, 1,
                                            colocated[ M.edges.vertex( e, 1 ) ] ) ;
                    }
                    M.vertices.delete_elements( to_delete, false ) ;

                    GEO::Logger::out( "GeoModel" ) << nb_todelete
                        << " colocated vertices deleted in " << E.gme_id()
                        << std::endl ;

                }
            }
        }
    }
}


namespace RINGMesh {

    void geo_model_mesh_repair( GeoModel& GM ) {
        GeoModelEditor editor( GM ) ;

        // Force removal of global vertices - Bugs I do not know where [JP]
        GM.mesh.vertices.clear() ;

        // Remove colocated vertices in each element
        std::set< gme_t > empty_elements ;
        remove_colocated_element_vertices( GM, empty_elements ) ;
        if( !empty_elements.empty() ) {
            editor.get_dependent_elements( empty_elements ) ;
            editor.remove_elements( empty_elements ) ;
        }

        // Basic mesh repair for surfaces and lines
        remove_degenerate_facet_and_edges( GM, empty_elements ) ;
        if( !empty_elements.empty() ) {
            editor.get_dependent_elements( empty_elements ) ;
            editor.remove_elements( empty_elements ) ;
        }

        // This is basic requirement ! no_colocated model vertices !
        // So remove them if there are any 
        GM.mesh.remove_colocated_vertices() ;
    }

}