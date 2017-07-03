/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/mesh/geogram_mesh.h>

#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace RINGMesh {

    std::vector< index_t > GeogramVolumeMesh::cells_around_vertex(
        index_t vertex_id,
        index_t cell_hint ) const
    {
        std::vector< index_t > result;

        if( cell_hint == NO_ID ) {
            cell_hint = find_first_cell_owing_vertex( vertex_id );
            if( cell_hint == NO_ID ) {
                return 0;
            }
        }

        // Flag the visited cells
        std::vector< index_t > visited;
        visited.reserve( 10 );

        // Stack of the adjacent cells
        std::stack< index_t > S;
        S.push( cell_hint );
        visited.push_back( cell_hint );

        do {
            index_t c = S.top();
            S.pop();

            bool cell_includes_vertex = false;
            for( index_t v = 0; v < nb_cell_vertices( c ); v++ ) {
                if( cell_vertex( c, v ) == vertex_id ) {
                    result.push_back( c );
                    cell_includes_vertex = true;
                    break;
                }
            }
            if( !cell_includes_vertex ) {
                continue;
            }

            for( index_t f = 0; f < nb_cell_facets( c ); f++ ) {
                for( index_t v = 0; v < nb_cell_facet_vertices( c, f ); v++ ) {
                    index_t vertex = cell_facet_vertex( c, f, v );
                    if( vertex == vertex_id ) {
                        index_t adj_P = cell_adjacent( c, f );

                        if( adj_P != NO_ID ) {
                            if( !contains( visited, adj_P ) ) {
                                S.push( adj_P );
                                visited.push_back( adj_P );
                            }
                        }
                        break;
                    }
                }
            }
        } while( !S.empty() );

        return result.size();
    }

    index_t GeogramVolumeMesh::find_first_cell_owing_vertex(
        index_t vertex_id_in_region ) const
    {
        const NNSearch& ann_cells = cells_nn_search();
        const vec3& vertex_pos = vertex( vertex_id_in_region );

        index_t nb_neighbors = std::min( index_t( 5 ), nb_cells() );
        std::vector< index_t > neighbors;
        index_t cur_neighbor = 0;
        index_t prev_neighbor = 0;
        do {
            prev_neighbor = cur_neighbor;
            cur_neighbor += nb_neighbors;
            cur_neighbor = std::min( cur_neighbor, nb_cells() );
            neighbors.resize( cur_neighbor );
            neighbors = ann_cells.get_neighbors( vertex_pos, cur_neighbor );
            // nb_neighbors can be less than cur_neighbor.
            nb_neighbors = static_cast< index_t >( neighbors.size() );
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                index_t c = neighbors[i];
                for( index_t j = 0; j < nb_cell_vertices( c ); j++ ) {
                    if( cell_vertex( c, j ) == vertex_id_in_region ) {
                        return c;
                    }
                }
            }
        } while( nb_cells() != cur_neighbor );

        return NO_ID;
    }

    void register_geogram_mesh()
    {
        ringmesh_register_point_mesh( GeogramPointSetMesh );
        ringmesh_register_point_mesh_builder( GeogramPointSetMesh );
        ringmesh_register_line_mesh( GeogramLineMesh );
        ringmesh_register_line_mesh_builder( GeogramLineMesh );
        ringmesh_register_surface_mesh( GeogramSurfaceMesh );
        ringmesh_register_surface_mesh_builder( GeogramSurfaceMesh );
        ringmesh_register_volume_mesh( GeogramVolumeMesh );
        ringmesh_register_volume_mesh_builder( GeogramVolumeMesh );
    }
}

