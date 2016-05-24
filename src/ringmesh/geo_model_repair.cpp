/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <algorithm>

#include <geogram/basic/logger.h>

#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_editor.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/geometry.h>

/*!
 * @file Implementation of repair function of the surfaces of a GeoModel
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    // using namespace RINGMesh ;

    typedef GeoModelMeshElement GMME ;

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
     *
     * \brief Tests whether a facet is degenerate.
     * \param[in] M the mesh that the facet belongs to
     * \param[in] f the index of the facet in \p M
     * \return true if facet \p f has duplicated vertices,
     *  false otherwise
     */
    bool GeoModelRepair::facet_is_degenerate(
        const Mesh& M,
        index_t f,
        GEO::vector< index_t >& colocated_vertices )
    {
        index_t nb_vertices = M.nb_facet_vertices( f ) ;
        if( nb_vertices != 3 ) {
            index_t* vertices = (index_t*) alloca( nb_vertices * sizeof(index_t) ) ;
            for( index_t lv = 0; lv < nb_vertices; ++lv ) {
                vertices[lv] = colocated_vertices[M.facet_vertex( f, lv )] ;
            }
            std::sort( vertices, vertices + nb_vertices ) ;
            return std::unique( vertices, vertices + nb_vertices )
                != vertices + nb_vertices ;
        }
        index_t v1 = colocated_vertices[M.facet_vertex( f, 0 )] ;
        index_t v2 = colocated_vertices[M.facet_vertex( f, 1 )] ;
        index_t v3 = colocated_vertices[M.facet_vertex( f, 2 )] ;
        return v1 == v2 || v2 == v3 || v3 == v1 ;
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
     */
    void GeoModelRepair::mesh_detect_degenerate_facets(
        const Mesh& M,
        GEO::vector< index_t >& f_is_degenerate,
        GEO::vector< index_t >& colocated_vertices )
    {
        f_is_degenerate.resize( M.nb_facets() ) ;
        for( index_t f = 0; f < M.nb_facets(); ++f ) {
            f_is_degenerate[f] = facet_is_degenerate( M, f, colocated_vertices ) ;
        }
    }

    /*!
     * @brief Detect and remove degenerated facets in a Mesh
     */
    index_t GeoModelRepair::detect_degenerate_facets( Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        const ColocaterANN& kdtree = M.colotater_ann( ColocaterANN::VERTICES ) ;
        kdtree.get_colocated_index_mapping( colocated ) ;
//        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        index_t nb_todelete = 0 ;
        for( index_t v = 0; v < colocated.size(); ++v ) {
            if( colocated[v] == v - nb_todelete ) {
                colocated[v] = v ;
            } else {
                nb_todelete++ ;
            }
        }

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_facets( M, degenerate, colocated ) ;
        return static_cast< index_t >( std::count( degenerate.begin(),
            degenerate.end(), 1 ) ) ;
    }

    void GeoModelRepair::mesh_detect_degenerate_edges(
        const Mesh& M,
        GEO::vector< index_t >& e_is_degenerate,
        GEO::vector< index_t >& colocated_vertices )
    {
        e_is_degenerate.resize( M.nb_edges() ) ;
        for( index_t e = 0; e < M.nb_edges(); ++e ) {
            e_is_degenerate[e] = edge_is_degenerate( M, e, colocated_vertices ) ;
        }
    }

    /*!
     * @brief Detect and remove degenerated edges in a Mesh
     */
    index_t GeoModelRepair::repair_line_mesh( Line& line )
    {
        GEO::vector< index_t > colocated ;
        const ColocaterANN& kdtree = line.vertex_colocater_ann() ;
        kdtree.get_colocated_index_mapping( colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_edges( line.mesh_, degenerate, colocated ) ;
        index_t nb = static_cast< index_t >( std::count( degenerate.begin(),
            degenerate.end(), 1 ) ) ;
        /// We have a problem if some vertices are left isolated
        /// If we remove them here we can kill all indices correspondances
        MeshBuilder builder( line.mesh_ ) ;
        builder.delete_edges( degenerate, false ) ;
        return nb ;
    }

    /*!
     * @brief Remove degenerate facets and edges from the Surface
     *        and Line of the model.
     * @param[in,out] GM Model to fix
     * @param[out] to_remove Ids of the elements of the model that are
     *  are emtpy once degenerate elements are removed
     * @pre Colocated vertices have already been removed
     */
    void GeoModelRepair::remove_degenerate_facet_and_edges(
        std::set< gme_t >& to_remove )
    {
        to_remove.clear() ;
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            Line& line = dynamic_cast<Line&>(element( gme_t( GME::LINE, i ) ));
            index_t nb = repair_line_mesh( line ) ;
            if( nb > 0 ) {
                GEO::Logger::out( "GeoModel" ) << nb
                    << " degenerated edges removed in LINE " << i << std::endl ;
                // If the Line is set it to remove
                if( model_.line( i ).nb_polytope() == 0 ) {
                    to_remove.insert( model_.line( i ).gme_id() ) ;
                }
            }
        }
        // The builder might be needed

        for( index_t i = 0; i < model_.nb_surfaces(); ++i ) {
            Surface& surface = dynamic_cast<Surface&>(element( gme_t(GME::SURFACE, i) ) );
            index_t nb = detect_degenerate_facets( surface.mesh_ ) ;
            /// @todo Check if that cannot be simplified 
            if( nb > 0 ) {
                // If there are some degenerated facets 
                // Using repair function of geogram
                // Warning - This triangulates the mesh
                if( surface.nb_vertices() > 0 ) {
                    // Colocated vertices must be processed before
                    // MESH_REPAIR_DUP_F 2 ;
                    GEO::MeshRepairMode mode =
                        static_cast< GEO::MeshRepairMode >( 2 ) ;
                    MeshBuilder builder( surface.mesh_ ) ;
                    builder.mesh_repair( mode, 0.0 ) ;

                    // This might create some small components - remove them
                    /// @todo How to choose the epsilon ? and the maximum number of facets ?
                    builder.remove_small_connected_components( epsilon_sq, 3 ) ;

                    // Alright, this is a bit of an overkill [JP]
                    if( surface.nb_vertices() > 0 ) {
                        builder.mesh_repair( mode, 0.0 ) ;
                    }
                }
                if( surface.nb_vertices() == 0 || surface.nb_polytope() == 0 ) {
                    to_remove.insert( model_.surface( i ).gme_id() ) ;
                } else {
                    // If the Surface has internal boundaries, we need to 
                    // re-cut the Surface along these lines
                    Surface& S = dynamic_cast<Surface&>(element( gme_t(GME::SURFACE, i) ) );
                    std::set< index_t > cutting_lines ;
                    for( index_t l = 0; l < S.nb_boundaries(); ++l ) {
                        const Line& L = model_.line( S.boundary_gme( l ).index ) ;
                        if( to_remove.count( L.gme_id() ) == 0
                            && L.is_inside_border( S ) ) {
                            cutting_lines.insert( L.index() ) ;
                        }
                    }
                    for( std::set< index_t >::iterator it = cutting_lines.begin();
                        it != cutting_lines.end(); ++it ) {
                        // Force the recomputing of the model vertices
                        // before performing the cut. 
                        model_.mesh.vertices.clear() ;
                        cut_surface_by_line( S, model_.line( *it ) ) ;
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
    void GeoModelRepair::vertices_on_inside_boundary(
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
            // We assume that the model vertices are not computed
            const ColocaterANN& kdtree = E.vertex_colocater_ann() ;

            for( index_t i = 0; i < inside_border.size(); ++i ) {
                for( index_t v = 0; v < inside_border[i]->nb_vertices(); ++v ) {
                    std::vector< index_t > colocated_indices ;
                    kdtree.get_colocated( inside_border[i]->vertex( v ),
                        colocated_indices ) ;
                    if( colocated_indices.size() > 1 ) {
                        std::sort( colocated_indices.begin(),
                            colocated_indices.end() ) ;
                        // Add colocated vertices except one to the duplicated vertices set
                        vertices.insert( colocated_indices.begin() + 1,
                            colocated_indices.end() ) ;
                    }
                }
            }
        }
    }

    /*!
     * @details Global GeoModel mesh is supposed to be empty
     */
    void GeoModelRepair::remove_colocated_element_vertices(
        std::set< gme_t >& to_remove )
    {
        to_remove.clear() ;
        // For all Lines and Surfaces
        for( index_t t = GME::LINE; t < GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;

            for( index_t e = 0; e < model_.nb_elements( T ); ++e ) {
                const GMME& E = model_.mesh_element( gme_t( T, e ) ) ;

                const ColocaterANN& kdtree = E.vertex_colocater_ann() ;
                GEO::vector< index_t > colocated ;
                kdtree.get_colocated_index_mapping( colocated ) ;

                // Get the vertices to delete
                std::set< index_t > inside_border ;
                vertices_on_inside_boundary( E, inside_border ) ;

                GEO::vector< index_t > to_delete( colocated.size(), 0 ) ;
                index_t nb_todelete = 0 ;
                index_t cur = 0 ;
                for( index_t v = 0; v < colocated.size(); ++v ) {
                    if( colocated[v] == v - nb_todelete
                        || inside_border.find( v ) != inside_border.end() ) {
                        // This point is kept 
                        // No colocated or on an inside boundary
                        colocated[v] = v ;
                        ++cur ;
                    } else {
                        // The point is to remove
                        to_delete[v] = 1 ;
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
                    GMME& ME = model_.modifiable_mesh_element( gme_t( T, e ) ) ;
                    MeshBuilder builder( ME.mesh_ ) ;
                    for( index_t c = 0; c < E.mesh_.nb_facet_corners(); c++ ) {
                        builder.set_facet_corner( c,
                            colocated[E.mesh_.facet_corner_vertex( c )] ) ;
                    }
                    for( index_t e = 0; e < E.mesh_.nb_edges(); e++ ) {
                        builder.set_edge_vertex( e, 0,
                            colocated[E.mesh_.edge_vertex( e, 0 )] ) ;
                        builder.set_edge_vertex( e, 1,
                            colocated[E.mesh_.edge_vertex( e, 1 )] ) ;
                    }
                    builder.delete_vertices( to_delete, false ) ;
                    GEO::Logger::out( "Repair" ) << nb_todelete
                        << " colocated vertices deleted in " << E.gme_id()
                        << std::endl ;
                }
            }
        }
    }

    void GeoModelRepair::geo_model_mesh_repair()
    {
        // Force removal of global vertices - Bugs ? I do not know where [JP]
        model_.mesh.vertices.clear() ;

        // Remove colocated vertices in each element
        std::set< gme_t > empty_elements ;
        remove_colocated_element_vertices( empty_elements ) ;
        if( !empty_elements.empty() ) {
            get_dependent_elements( empty_elements ) ;
            remove_elements( empty_elements ) ;
        }

        // Basic mesh repair for surfaces and lines
        remove_degenerate_facet_and_edges( empty_elements ) ;
        if( !empty_elements.empty() ) {
            get_dependent_elements( empty_elements ) ;
            remove_elements( empty_elements ) ;
        }

        // This is basic requirement ! no_colocated model vertices !
        // So remove them if there are any 
        model_.mesh.remove_colocated_vertices() ;
    }

}
