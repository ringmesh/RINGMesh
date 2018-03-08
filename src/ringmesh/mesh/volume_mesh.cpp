/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Francois Bonneau */

#include <numeric>
#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <stack>

#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    std::unique_ptr< VolumeMesh< DIMENSION > >
        VolumeMesh< DIMENSION >::create_mesh( const MeshType type )
    {
        auto new_type = type;
        if( new_type.empty() )
        {
            if( !PointSetMeshFactory< DIMENSION >::has_creator(
                    "GeogramPointSetMesh" ) )
            {
                throw RINGMeshException( "VolumeMesh",
                    "Default mesh data structure not registered" );
            }
            return create_mesh( "GeogramVolumeMesh" );
        }
        auto mesh = VolumeMeshFactory< DIMENSION >::create( new_type );
        if( !mesh )
        {
            Logger::warn( "VolumeMesh",
                "Could not create mesh data structure: ", new_type );
            Logger::warn( "VolumeMesh",
                "Falling back to GeogramVolumeMesh data structure" );

            return create_mesh();
        }
        return mesh;
    }

    template < index_t DIMENSION >
    double VolumeMesh< DIMENSION >::cell_edge_length(
        index_t cell_id, index_t edge_id ) const
    {
        const auto& e0 =
            this->vertex( cell_edge_vertex( cell_id, edge_id, 0 ) );
        const auto& e1 =
            this->vertex( cell_edge_vertex( cell_id, edge_id, 1 ) );
        return ( e1 - e0 ).length();
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > VolumeMesh< DIMENSION >::cell_edge_barycenter(
        index_t cell_id, index_t edge_id ) const
    {
        const auto& e0 =
            this->vertex( cell_edge_vertex( cell_id, edge_id, 0 ) );
        const auto& e1 =
            this->vertex( cell_edge_vertex( cell_id, edge_id, 1 ) );
        return ( e1 + e0 ) / 2.;
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > VolumeMesh< DIMENSION >::cell_facet_barycenter(
        const CellLocalFacet& cell_local_facet ) const
    {
        vecn< DIMENSION > result;
        index_t nb_vertices = nb_cell_facet_vertices( cell_local_facet );
        for( auto v : range( nb_vertices ) )
        {
            result += this->vertex( cell_facet_vertex( cell_local_facet, v ) );
        }
        ringmesh_assert( nb_vertices > 0 );

        return result / static_cast< double >( nb_vertices );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > VolumeMesh< DIMENSION >::cell_barycenter(
        index_t cell_id ) const
    {
        vecn< DIMENSION > result;
        ringmesh_assert( nb_cell_vertices( cell_id ) >= 1 );
        for( auto v : range( nb_cell_vertices( cell_id ) ) )
        {
            result += this->vertex( cell_vertex( { cell_id, v } ) );
        }
        return ( 1.0 / nb_cell_vertices( cell_id ) ) * result;
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > VolumeMesh< DIMENSION >::cell_facet_normal(
        const CellLocalFacet& cell_local_facet ) const
    {
        ringmesh_assert( cell_local_facet.cell_id < nb_cells() );
        ringmesh_assert( cell_local_facet.local_facet_id
                         < nb_cell_facets( cell_local_facet.cell_id ) );

        const auto& p1 =
            this->vertex( cell_facet_vertex( cell_local_facet, 0 ) );
        const auto& p2 =
            this->vertex( cell_facet_vertex( cell_local_facet, 1 ) );
        const auto& p3 =
            this->vertex( cell_facet_vertex( cell_local_facet, 2 ) );

        return cross( p2 - p1, p3 - p1 );
    }

    template < index_t DIMENSION >
    index_t VolumeMesh< DIMENSION >::find_cell_corner(
        index_t cell_id, index_t vertex_id ) const
    {
        for( auto v : range( nb_cell_vertices( cell_id ) ) )
        {
            if( cell_vertex( { cell_id, v } ) == vertex_id )
            {
                return v;
            }
        }
        return NO_ID;
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >&
        VolumeMesh< DIMENSION >::cell_facet_nn_search() const
    {
        if( !cell_facet_nn_search_ )
        {
            std::vector< vecn< DIMENSION > > cell_facet_centers(
                nb_cell_facets() );
            index_t cf = 0;
            for( auto c : range( nb_cells() ) )
            {
                for( auto f : range( nb_cell_facets( c ) ) )
                {
                    cell_facet_centers[cf] = cell_facet_barycenter( { c, f } );
                    ++cf;
                }
            }
            cell_facet_nn_search_.reset(
                new NNSearch< DIMENSION >( cell_facet_centers, true ) );
        }
        return *cell_facet_nn_search_.get();
    }

    template < index_t DIMENSION >
    std::tuple< index_t, std::vector< index_t > >
        VolumeMesh< DIMENSION >::connected_components() const
    {
        std::vector< index_t > components( nb_cells(), NO_ID );
        index_t nb_components{ 0 };
        for( auto cell : range( nb_cells() ) )
        {
            if( components[cell] == NO_ID )
            {
                std::stack< index_t > S;
                S.push( cell );
                components[cell] = nb_components;
                do
                {
                    auto cur_cell = S.top();
                    S.pop();
                    for( auto facet : range( nb_cell_facets( cur_cell ) ) )
                    {
                        auto adj_cell = cell_adjacent( { cur_cell, facet } );
                        if( adj_cell != NO_ID && components[adj_cell] == NO_ID )
                        {
                            S.push( adj_cell );
                            components[adj_cell] = nb_components;
                        }
                    }
                } while( !S.empty() );
                nb_components++;
            }
        }
        return std::make_tuple( nb_components, components );
    }

    template < index_t DIMENSION >
    bool VolumeMesh< DIMENSION >::is_mesh_valid() const
    {
        bool valid{ true };

        if( this->nb_vertices() < 4 )
        {
            Logger::warn( "VolumeMesh has less than 4 vertices " );
            valid = false;
        }
        if( nb_cells() == 0 )
        {
            Logger::warn( "VolumeMesh has no cell" );
            valid = false;
        }

        // No isolated vertices
        std::vector< index_t > nb( this->nb_vertices(), 0 );
        for( auto c : range( nb_cells() ) )
        {
            for( auto v : range( nb_cell_vertices( c ) ) )
            {
                nb[cell_vertex( { c, v } )]++;
            }
        }
        auto nb_isolated_vertices =
            static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) );
        if( nb_isolated_vertices > 0 )
        {
            Logger::warn( "VolumeMesh", "Mesh has ", nb_isolated_vertices,
                " isolated vertices " );
            valid = false;
        }

        return valid;
    }

    template < index_t DIMENSION >
    void VolumeMesh< DIMENSION >::store_cells_around_vertex( index_t cell_hint,
        index_t vertex_id,
        std::vector< index_t >& result ) const
    {
        // Flag the visited cells
        std::vector< index_t > visited;
        visited.reserve( 10 );
        // Stack of the adjacent cells
        std::stack< index_t > S;
        S.push( cell_hint );
        visited.push_back( cell_hint );
        do
        {
            auto c = S.top();
            S.pop();
            bool cell_includes_vertex{ false };
            for( auto v : range( nb_cell_vertices( c ) ) )
            {
                if( cell_vertex( { c, v } ) == vertex_id )
                {
                    result.push_back( c );
                    cell_includes_vertex = true;
                    break;
                }
            }
            if( !cell_includes_vertex )
            {
                continue;
            }
            for( auto f : range( nb_cell_facets( c ) ) )
            {
                for( auto v : range( nb_cell_facet_vertices( { c, f } ) ) )
                {
                    auto vertex = cell_facet_vertex( { c, f }, v );
                    if( vertex == vertex_id )
                    {
                        auto adj_P = cell_adjacent( { c, f } );
                        if( adj_P != NO_ID && !contains( visited, adj_P ) )
                        {
                            S.push( adj_P );
                            visited.push_back( adj_P );
                        }
                        break;
                    }
                }
            }
        } while( !S.empty() );
    }

    template < index_t DIMENSION >
    std::vector< index_t > VolumeMesh< DIMENSION >::cells_around_vertex(
        index_t vertex_id, index_t cell_hint ) const
    {
        std::vector< index_t > result;

        if( cell_hint == NO_ID )
        {
            const vecn< DIMENSION > cur_vec = this->vertex( vertex_id );
            index_t cell_vertex_not_used = NO_ID;
            bool found = find_cell_from_colocated_vertex_within_distance_if_any(
                cur_vec, global_epsilon, cell_hint, cell_vertex_not_used );
            if( !found )
            {
                return result;
            }
        }
        ringmesh_assert( cell_hint != NO_ID );

        // Flag the visited cells
        store_cells_around_vertex( cell_hint, vertex_id, result );
        return result;
    }

    template < index_t DIMENSION >
    bool VolumeMesh< DIMENSION >::
        find_cell_from_colocated_vertex_within_distance_if_any(
            const vecn< DIMENSION >& vertex_vec,
            double distance,
            index_t& cell_id,
            index_t& cell_vertex_id ) const
    {
        bool result = false;
        cell_nn_search().get_neighbors( vertex_vec,
            [this, &vertex_vec, &result, &cell_id, &cell_vertex_id, distance](
                index_t i ) {
                for( auto j : range( nb_cell_vertices( i ) ) )
                {
                    if( inexact_equal( this->vertex( cell_vertex( { i, j } ) ),
                            vertex_vec, distance ) )
                    {
                        cell_vertex_id = cell_vertex( { i, j } );
                        cell_id = i;
                        result = true;
                        break;
                    }
                }
                return result;
            } );
        return result;
    }

    template class mesh_api VolumeMesh< 3 >;
} // namespace RINGMesh
