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

#include <ringmesh/internal_borders_builder.h>

#include <algorithm>

/*!
 * @file ringmesh/internal_borders_builder.cpp
 * @brief Class to create the internal borders of the GME. Dissociates vertices,
 * edges, facets...
 * @author Benjamin Chauvin
 */

namespace RINGMesh {

    InternalBordersBuilder::InternalBordersBuilder( GeoModel& model )
        : GeoModelBuilder( model )
    {

    }

    InternalBordersBuilder::~InternalBordersBuilder()
    {
    }

    void InternalBordersBuilder::compute_internal_borders()
    {
        //        geo_model_mesh_repair( model() ) ;
        // ========= bad copy paste from geo model repair
        for( index_t i = 0; i < model().nb_surfaces(); ++i ) {
            // If the Surface has internal boundaries, we need to
            // re-cut the Surface along these lines
            Surface& S = const_cast< Surface& >( model().surface( i ) ) ;
            std::set< index_t > cutting_lines ;
            for( index_t l = 0; l < S.nb_boundaries(); ++l ) {
                const Line& L = model().line( S.boundary_gme( l ).index ) ;
                if( /*to_remove.count( L.gme_id() ) == 0 &&*/L.is_inside_border(
                    S ) ) {
                    cutting_lines.insert( L.index() ) ;
                }
            }
            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                // Force the recomputing of the model vertices
                // before performing the cut.
                model().mesh.vertices.clear() ;
                DEBUG( "cut surface by line" ) ;
                cut_surface_by_line( S.index(), *it ) ;
            }
        }
        // ========= bad copy paste from geo model repair

        //=================== Cut the region by the surfaces
        // THAT MAY NOT WORK IF THE REGION IS ALREADY CUT. TO CHECK!!!
        for( index_t i = 0; i < model().nb_regions(); ++i ) {
            // If the Region has internal boundaries, we need to
            // re-cut the Region along these surfaces
            Region& R = const_cast< Region& >( model().region( i ) ) ;
            if( !R.is_meshed() ) {
                continue ;
            }
            // the std::set avoids that a region is cut twice by the same
            // surface (internal borner is defined by twice the region in the
            // in boundaries of the surface.
            std::set< index_t > cutting_surfaces ;
            for( index_t s = 0; s < R.nb_boundaries(); ++s ) {
                const Surface& S = model().surface( R.boundary_gme( s ).index ) ;
                if( /*to_remove.count( L.gme_id() ) == 0 &&*/S.is_inside_border(
                    R ) ) {
                    cutting_surfaces.insert( S.index() ) ;
                } /*else if( GME::is_fault( S.parent().geological_feature() ) ) {
                 cutting_surfaces.insert( S.index() ) ;
                 }*/
            }
            for( std::set< index_t >::iterator it = cutting_surfaces.begin();
                it != cutting_surfaces.end(); ++it ) {
                // Force the recomputing of the model vertices
                // before performing the cut.
                DEBUG("cut region by surface") ;
                model().mesh.vertices.clear() ;
                cut_region_by_surface( R, model().surface( *it ) ) ;
            }

            // If the region has only a line has internal border.
            // This line is shared by two boundaries (surface) of
            // the region and this line does not belong to a cutting surface.
            // The two boundaries (surface) in question seem to be horizons
            // of different interfaces. If it is really the case, the fact that
            // there are not cutting surfaces should be in assert (an horizon
            // does not cut).
            std::set< index_t > cutting_lines ;
            for( index_t surf_boun_itr = 0; surf_boun_itr < R.nb_boundaries();
                ++surf_boun_itr ) {
                const GME& cur_surf_boun = R.boundary( surf_boun_itr ) ;
                ringmesh_assert( cur_surf_boun.type() == GME::SURFACE ) ;
                if( !GME::is_stratigraphic_limit(
                    cur_surf_boun.parent().geological_feature() ) ) {
                    continue ;
                }
                if( std::find( cutting_surfaces.begin(), cutting_surfaces.end(),
                    cur_surf_boun.index() ) != cutting_surfaces.end() ) {
                    continue ;
                }
                for( index_t line_boun_itr = 0;
                    line_boun_itr < cur_surf_boun.nb_boundaries();
                    ++line_boun_itr ) {
                    const GME& cur_line_boun = cur_surf_boun.boundary(
                        line_boun_itr ) ;
                    ringmesh_assert( cur_line_boun.type() == GME::LINE ) ;
                    for( index_t surf_boun_itr2 = 0;
                        surf_boun_itr2 < R.nb_boundaries(); ++surf_boun_itr2 ) {
                        const GME& cur_surf_boun2 = R.boundary( surf_boun_itr2 ) ;
                        ringmesh_assert( cur_surf_boun2.type() == GME::SURFACE ) ;
                        if( !GME::is_stratigraphic_limit(
                            cur_surf_boun2.parent().geological_feature() ) ) {
                            continue ;
                        }
                        if( surf_boun_itr == surf_boun_itr2 ) {
                            continue ;
                        }
                        if( std::find( cutting_surfaces.begin(),
                            cutting_surfaces.end(), cur_surf_boun2.index() )
                            != cutting_surfaces.end() ) {
                            continue ;
                        }

                        // The 2 surfaces are not in the same interface? It seems... to check
                        if( cur_surf_boun.parent().index()
                            == cur_surf_boun2.parent().index() ) {
                            continue ;
                        }

                        bool found = false ;
                        for( index_t line_boun_itr2 = 0;
                            line_boun_itr2 < cur_surf_boun2.nb_boundaries();
                            ++line_boun_itr2 ) {
                            const GME& cur_line_boun2 = cur_surf_boun2.boundary(
                                line_boun_itr2 ) ;
                            ringmesh_assert( cur_line_boun2.type() == GME::LINE ) ;

                            if( cur_line_boun.index() == cur_line_boun2.index() ) {
                                DEBUG("cutting line") ;
                                DEBUG( cur_line_boun.index() ) ;
                                cutting_lines.insert( cur_line_boun.index() ) ;
                                found = true ;
                                break ;
                            }
                        }
                        if( found ) {
                            break ;
                        }
                    }
                }
            }

            for( std::set< index_t >::iterator it = cutting_lines.begin();
                it != cutting_lines.end(); ++it ) {
                DEBUG("cut region by line") ;
                model().mesh.vertices.clear() ;
                cut_region_by_line( R, model().line( *it ) ) ;
//                R.tools.delete_ann() ;
            }
            //=================== Cut the region by the lines
        }

        {
            recompute_geomodel_mesh() ;
            index_t toto = model().corner( 0 ).model_vertex_id( 0 ) ;
            const std::vector< GMEVertex >& gme_vertices =
                model().mesh.vertices.gme_vertices( toto ) ;
            index_t count = 0 ;
            for( index_t i = 0; i < gme_vertices.size(); ++i ) {
                if( gme_vertices[i].gme_id.type == GME::REGION ) {
                    ++count ;
                }
            }
            DEBUG(count) ;
        }
        //=================== Cut the region by the surfaces

        for( index_t reg_itr = 0; reg_itr < model().nb_regions(); ++reg_itr ) {
            MeshBuilder region_mesh_builder( mesh( GME::gme_t( GME::REGION, reg_itr ) ) ) ;
            region_mesh_builder.remove_isolated_vertices() ;
        }
        // Deliberate clear of the model vertices used for model building
        model().mesh.vertices.clear() ;
        {
            recompute_geomodel_mesh() ;
            index_t toto = model().corner( 0 ).model_vertex_id( 0 ) ;
            const std::vector< GMEVertex >& gme_vertices =
                model().mesh.vertices.gme_vertices( toto ) ;
            index_t count = 0 ;
            for( index_t i = 0; i < gme_vertices.size(); ++i ) {
                if( gme_vertices[i].gme_id.type == GME::REGION ) {
                    ++count ;
                }
            }
            DEBUG(count) ;
        }
    }
}
