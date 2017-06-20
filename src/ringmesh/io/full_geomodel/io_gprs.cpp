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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

namespace {
    class GPRSIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        struct Pipe {
            Pipe( index_t v0_in, index_t v1_in )
                : v0( v0_in ), v1( v1_in )
            {
            }
            index_t v0;
            index_t v1;
        };
        virtual bool load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from GPRS not implemented yet" );
            return false;
        }
        virtual void save(
            const GeoModel< 3 >& geomodel,
            const std::string& filename ) final
        {
            std::string path = GEO::FileSystem::dir_name( filename );
            std::string directory = GEO::FileSystem::base_name( filename );
            if( path == "." ) {
                path = GEO::FileSystem::get_current_working_directory();
            }
            std::ostringstream oss;
            oss << path << "/" << directory;
            std::string full_path = oss.str();
            GEO::FileSystem::create_directory( full_path );

            std::ostringstream oss_pipes;
            oss_pipes << full_path << "/pipes.in";
            std::ofstream out_pipes( oss_pipes.str().c_str() );

            std::ostringstream oss_vol;
            oss_vol << full_path << "/vol.in";
            std::ofstream out_vol( oss_vol.str().c_str() );
            out_vol.precision( 16 );

            std::ostringstream oss_xyz;
            oss_xyz << full_path << "/gprs.xyz";
            std::ofstream out_xyz( oss_xyz.str().c_str() );
            out_xyz.precision( 16 );

            const GeoModelMesh< 3 >& mesh = geomodel.mesh;
            std::deque< Pipe > pipes;
            index_t cell_offset = mesh.cells.nb();
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t facet = NO_ID;
                    bool not_used;
                    if( mesh.cells.is_cell_facet_on_surface( c, f, facet,
                        not_used ) ) {
                        pipes.emplace_back( c, facet + cell_offset );
                    } else {
                        index_t adj = mesh.cells.adjacent( c, f );
                        if( adj != GEO::NO_CELL && adj < c ) {
                            pipes.emplace_back( c, adj );
                        }
                    }
                }
            }

            index_t nb_edges = 0;
            for( index_t l = 0; l < geomodel.nb_lines(); l++ ) {
                nb_edges += geomodel.line( l ).nb_mesh_elements();
            }
            std::vector< index_t > temp;
            temp.reserve( 3 );
            std::vector< std::vector< index_t > > edges( nb_edges, temp );
            std::vector< vec3 > edge_vertices( nb_edges );
            index_t count_edge = 0;
            for( index_t l = 0; l < geomodel.nb_lines(); l++ ) {
                const Line< 3 >& line = geomodel.line( l );
                for( index_t e = 0; e < line.nb_mesh_elements(); e++ ) {
                    edge_vertices[count_edge++ ] = 0.5
                        * ( line.vertex( e ) + line.vertex( e + 1 ) );
                }
            }
            NNSearch < 3 > nn_search( edge_vertices, false );

            const GeoModelMeshPolygons< 3 >& polygons = geomodel.mesh.polygons;
            for( index_t p = 0; p < polygons.nb(); p++ ) {
                for( index_t e = 0; e < polygons.nb_vertices( p ); e++ ) {
                    index_t adj = polygons.adjacent( p, e );
                    if( adj != GEO::NO_CELL && adj < p ) {
                        pipes.emplace_back( p + cell_offset, adj + cell_offset );
                    } else {
                        const vec3& e0 = mesh.vertices.vertex(
                            polygons.vertex( p, e ) );
                        const vec3& e1 = mesh.vertices.vertex(
                            polygons.vertex( p,
                                ( e + 1 ) % polygons.nb_vertices( p ) ) );
                        vec3 query = 0.5 * ( e0 + e1 );
                        std::vector< index_t > results = nn_search.get_neighbors(
                            query, geomodel.epsilon() );
                        if( !results.empty() ) {
                            edges[results[0]].push_back( cell_offset + p );
                        } else {
                            ringmesh_assert_not_reached;
                        }
                    }
                }
            }

            index_t nb_pipes = pipes.size();
            for( index_t e = 0; e < edges.size(); e++ ) {
                nb_pipes += binomial_coef( edges[e].size() );
            }
            out_pipes << nb_pipes << std::endl;
            for( index_t p = 0; p < pipes.size(); p++ ) {
                const Pipe& pipe = pipes[p];
                out_pipes << pipe.v0 << SPACE << pipe.v1 << std::endl;
            }
            for( index_t e = 0; e < edges.size(); e++ ) {
                const std::vector< index_t > vertices = edges[e];
                for( index_t v0 = 0; v0 < vertices.size() - 1; v0++ ) {
                    for( index_t v1 = v0 + 1; v1 < vertices.size(); v1++ ) {
                        out_pipes << vertices[v0] << SPACE << vertices[v1]
                            << std::endl;
                    }
                }
            }

            out_xyz
                << "Node geometry, not used by GPRS but useful to reconstruct a pipe-network"
                << std::endl;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out_xyz << mesh.cells.barycenter( c ) << std::endl;
                out_vol << mesh.cells.volume( c ) << std::endl;
            }
            for( index_t p = 0; p < polygons.nb(); p++ ) {
                out_xyz << polygons.center( p ) << std::endl;
                out_vol << polygons.area( p ) << std::endl;
            }
        }
        index_t binomial_coef( index_t n ) const
        {
            switch( n ) {
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
