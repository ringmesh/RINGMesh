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

//#include <ringmesh/mesh/mesh.h>

#include <numeric>
#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <stack>

#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/line_mesh.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    std::unique_ptr< LineMesh< DIMENSION > > LineMesh< DIMENSION >::create_mesh(
        const MeshType type )
    {
        MeshType new_type = type;
        if( new_type.empty() )
        {
            if( !PointSetMeshFactory< DIMENSION >::has_creator(
                "GeogramPointSetMesh" ) )
            {
                throw RINGMeshException( "LineMesh",
                    "Default mesh data structure not registered" );
            }
            return create_mesh( "GeogramLineMesh" );
        }
        auto mesh = LineMeshFactory< DIMENSION >::create( new_type );
        if( !mesh )
        {
            Logger::warn( "LineMesh", "Could not create mesh data structure: ",
                new_type );
            Logger::warn(
                "LineMesh", "Falling back to GeogramLineMesh data structure" );

            return create_mesh();
        }
        return mesh;
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& LineMesh< DIMENSION >::edge_nn_search() const
    {
        if( !edge_nn_search_ )
        {
            std::vector< vecn< DIMENSION > > edge_centers( nb_edges() );
            for( auto e : range( nb_edges() ) )
            {
                edge_centers[e] = edge_barycenter( e );
            }
            edge_nn_search_.reset(
                new NNSearch< DIMENSION >( edge_centers, true ) );
        }
        return *edge_nn_search_.get();
    }

    template < index_t DIMENSION >
    const LineAABBTree< DIMENSION >& LineMesh< DIMENSION >::edge_aabb() const
    {
        if( !edge_aabb_ )
        {
            edge_aabb_.reset( new LineAABBTree< DIMENSION >( *this ) );
        }
        return *edge_aabb_.get();
    }

    template < index_t DIMENSION >
    double LineMesh< DIMENSION >::edge_length( index_t edge_id ) const
    {
        const auto& e0 = this->vertex( edge_vertex( { edge_id, 0 } ) );
        const auto& e1 = this->vertex( edge_vertex( { edge_id, 1 } ) );
        return ( e1 - e0 ).length();
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > LineMesh< DIMENSION >::edge_barycenter(
        index_t edge_id ) const
    {
        const auto& e0 = this->vertex( edge_vertex( { edge_id, 0 } ) );
        const auto& e1 = this->vertex( edge_vertex( { edge_id, 1 } ) );
        return ( e1 + e0 ) / 2.;
    }

    template < index_t DIMENSION >
    bool LineMesh< DIMENSION >::is_mesh_valid() const
    {
        bool valid{ true };

        if( this->nb_vertices() < 2 )
        {
            Logger::err( "LineMesh", "Mesh has less than 2 vertices " );
            valid = false;
        }

        if( nb_edges() == 0 )
        {
            Logger::err( "LineMesh", "Mesh has no edge" );
            valid = false;
        }

        // No isolated vertices
        std::vector< index_t > nb( this->nb_vertices(), 0 );
        for( auto p : range( nb_edges() ) )
        {
            for( auto v : range( 2 ) )
            {
                nb[edge_vertex( { p, v } )]++;
            }
        }
        auto nb_isolated_vertices =
            static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) );
        if( nb_isolated_vertices > 0 )
        {
            Logger::warn( "LineMesh", "Mesh has ", nb_isolated_vertices,
                " isolated vertices " );
            valid = false;
        }

        return valid;
    }

    template < index_t DIMENSION >
    std::tuple< index_t, std::vector< index_t > >
        LineMesh< DIMENSION >::connected_components() const
    {
        std::vector< index_t > components( nb_edges(), NO_ID );
        std::vector< index_t > vertex_components( this->nb_vertices(), NO_ID );
        index_t nb_components{ 0 };

        for( auto edge : range( nb_edges() ) )
        {
            ringmesh_assert( components[edge] == NO_ID );
            const auto v0 = edge_vertex( { edge, 0 } );
            const auto v1 = edge_vertex( { edge, 1 } );
            if( vertex_components[v0] == NO_ID
                && vertex_components[v1] == NO_ID )
            {
                vertex_components[v0] = nb_components;
                vertex_components[v1] = nb_components;
                components[edge] = nb_components;
                ++nb_components;
            }
            else if( vertex_components[v0] != NO_ID
                     && vertex_components[v1] == NO_ID )
            {
                vertex_components[v1] = vertex_components[v0];
                components[edge] = vertex_components[v0];
            }
            else if( vertex_components[v0] == NO_ID
                     && vertex_components[v1] != NO_ID )
            {
                vertex_components[v0] = vertex_components[v1];
                components[edge] = vertex_components[v1];
            }
            else
            {
                // Case both nodes have already a connected component.
                if( vertex_components[v0] == vertex_components[v1] )
                {
                    components[edge] = vertex_components[v0];
                }
                else
                {
                    // It appears that 2 previously identified connected
                    // components
                    // correspond in fact to a unique connected component.
                    auto min_connected_components = std::min(
                        vertex_components[v0], vertex_components[v1] );
                    auto max_connected_components = std::max(
                        vertex_components[v0], vertex_components[v1] );
                    ringmesh_assert( min_connected_components != NO_ID );
                    ringmesh_assert( max_connected_components != NO_ID );
                    for( auto previous_edge : range( edge ) )
                    {
                        ringmesh_assert( components[previous_edge] != NO_ID );
                        ringmesh_assert( vertex_components[edge_vertex(
                                             { previous_edge, 0 } )]
                                         != NO_ID );
                        ringmesh_assert( vertex_components[edge_vertex(
                                             { previous_edge, 1 } )]
                                         != NO_ID );
                        if( components[previous_edge]
                            == max_connected_components )
                        {
                            components[previous_edge] =
                                min_connected_components;
                            vertex_components[edge_vertex( { previous_edge,
                                0 } )] = min_connected_components;
                            vertex_components[edge_vertex( { previous_edge,
                                1 } )] = min_connected_components;
                        }
                        else if( components[previous_edge]
                                 > max_connected_components )
                        {
                            ringmesh_assert(
                                components[previous_edge] - 1 >= 0 );
                            ringmesh_assert( vertex_components[edge_vertex(
                                                 { previous_edge, 0 } )]
                                                 - 1
                                             >= 0 );
                            ringmesh_assert( vertex_components[edge_vertex(
                                                 { previous_edge, 1 } )]
                                                 - 1
                                             >= 0 );
                            --components[previous_edge];
                            vertex_components[edge_vertex( { previous_edge,
                                0 } )] = components[previous_edge];
                            vertex_components[edge_vertex( { previous_edge,
                                1 } )] = components[previous_edge];
                            ringmesh_assert(
                                components[previous_edge] != NO_ID );
                            ringmesh_assert( vertex_components[edge_vertex(
                                                 { previous_edge, 0 } )]
                                             != NO_ID );
                            ringmesh_assert( vertex_components[edge_vertex(
                                                 { previous_edge, 1 } )]
                                             != NO_ID );
                            ringmesh_assert( components[previous_edge]
                                             == vertex_components[edge_vertex(
                                                    { previous_edge, 0 } )] );
                            ringmesh_assert( components[previous_edge]
                                             == vertex_components[edge_vertex(
                                                    { previous_edge, 1 } )] );
                        }
                    }
                    components[edge] = min_connected_components;
                    vertex_components[v0] = min_connected_components;
                    vertex_components[v1] = min_connected_components;
                    --nb_components;
                }
            }
            ringmesh_assert( components[edge] != NO_ID );
            ringmesh_assert( vertex_components[v0] != NO_ID );
            ringmesh_assert( vertex_components[v1] != NO_ID );
            ringmesh_assert( components[edge] == vertex_components[v0] );
            ringmesh_assert( components[edge] == vertex_components[v1] );
        }

        return std::make_tuple( nb_components, components );
    }

    template class mesh_api LineMesh< 2 >;

    template class mesh_api LineMesh< 3 >;
} // namespace RINGMesh
