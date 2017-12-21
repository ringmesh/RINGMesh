/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

namespace
{
    class GPRSIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        struct Pipe
        {
            Pipe( index_t v0_in, index_t v1_in ) : v0( v0_in ), v1( v1_in )
            {
            }
            index_t v0;
            index_t v1;
        };
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::string path = GEO::FileSystem::dir_name( filename );
            std::string name = GEO::FileSystem::base_name( filename );

            std::ofstream out_pipes( "pipes.in" );
            std::ofstream out_vol( "vol.in" );
            out_vol.precision( 16 );

            std::ostringstream oss_xyz;
            oss_xyz << name << ".xyz";
            std::ofstream out_xyz( oss_xyz.str().c_str() );
            out_xyz.precision( 16 );

            const GeoModelMesh3D& mesh = geomodel.mesh;
            std::deque< Pipe > pipes;
            index_t cell_offset = mesh.cells.nb();
            for( auto c : range( mesh.cells.nb() ) )
            {
                for( auto f : range( mesh.cells.nb_facets( c ) ) )
                {
                    index_t facet{ NO_ID };
                    bool not_used;
                    if( mesh.cells.is_cell_facet_on_surface(
                            c, f, facet, not_used ) )
                    {
                        pipes.emplace_back( c, facet + cell_offset );
                    }
                    else
                    {
                        index_t adj = mesh.cells.adjacent( c, f );
                        if( adj != NO_ID && adj < c )
                        {
                            pipes.emplace_back( c, adj );
                        }
                    }
                }
            }

            index_t nb_edges = 0;
            for( const auto& line : geomodel.lines() )
            {
                nb_edges += line.nb_mesh_elements();
            }
            std::vector< index_t > temp;
            temp.reserve( 3 );
            std::vector< std::vector< index_t > > edges( nb_edges, temp );
            std::vector< vec3 > edge_vertices( nb_edges );
            index_t count_edge = 0;
            for( const auto& line : geomodel.lines() )
            {
                for( index_t e : range( line.nb_mesh_elements() ) )
                {
                    edge_vertices[count_edge++] =
                        0.5 * ( line.vertex( e ) + line.vertex( e + 1 ) );
                }
            }
            NNSearch3D nn_search( edge_vertices, false );

            const GeoModelMeshPolygons3D& polygons = geomodel.mesh.polygons;
            for( index_t p : range( polygons.nb() ) )
            {
                for( index_t e : range( polygons.nb_vertices( p ) ) )
                {
                    index_t adj = polygons.adjacent( PolygonLocalEdge( p, e ) );
                    if( adj != NO_ID && adj < p )
                    {
                        pipes.emplace_back(
                            p + cell_offset, adj + cell_offset );
                    }
                    else
                    {
                        const vec3& e0 = mesh.vertices.vertex(
                            polygons.vertex( ElementLocalVertex( p, e ) ) );
                        const vec3& e1 = mesh.vertices.vertex(
                            polygons.vertex( ElementLocalVertex(
                                p, ( e + 1 ) % polygons.nb_vertices( p ) ) ) );
                        vec3 query = 0.5 * ( e0 + e1 );
                        std::vector< index_t > results =
                            nn_search.get_neighbors(
                                query, geomodel.epsilon() );
                        if( !results.empty() )
                        {
                            edges[results[0]].push_back( cell_offset + p );
                        }
                        else
                        {
                            ringmesh_assert_not_reached;
                        }
                    }
                }
            }

            auto nb_pipes = static_cast< index_t >( pipes.size() );
            for( const auto& vertices : edges )
            {
                nb_pipes +=
                    binomial_coef( static_cast< index_t >( vertices.size() ) );
            }
            out_pipes << nb_pipes << EOL;
            for( const auto& pipe : pipes )
            {
                out_pipes << pipe.v0 << SPACE << pipe.v1 << EOL;
            }
            for( const auto& vertices : edges )
            {
                for( auto v0 : range( vertices.size() - 1 ) )
                {
                    for( auto v1 : range( v0 + 1, vertices.size() ) )
                    {
                        out_pipes << vertices[v0] << SPACE << vertices[v1]
                                  << EOL;
                    }
                }
            }

            out_xyz << "Node geometry, not used by GPRS but useful to "
                       "reconstruct a pipe-network"
                    << EOL;
            for( auto c : range( mesh.cells.nb() ) )
            {
                out_xyz << mesh.cells.barycenter( c ) << EOL;
                out_vol << mesh.cells.volume( c ) << EOL;
            }
            for( auto p : range( polygons.nb() ) )
            {
                out_xyz << polygons.center( p ) << EOL;
                out_vol << polygons.area( p ) << EOL;
            }

            out_pipes << std::flush;
            out_vol << std::flush;
            out_xyz << std::flush;
        }
        index_t binomial_coef( index_t n ) const
        {
            switch( n )
            {
            case 1:
                return 0;
            case 2:
                return 1;
            case 3:
                return 3;
            case 4:
                return 6;
            case 5:
                return 10;
            case 6:
                return 15;
            case 7:
                return 21;
            case 8:
                return 28;
            case 9:
                return 36;
            case 10:
                return 45;
            default:
                ringmesh_assert_not_reached;
                return 0;
            }
        }
    };
}
