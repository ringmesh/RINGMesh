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
    class TetGenIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::string directory = GEO::FileSystem::dir_name( filename );
            std::string file = GEO::FileSystem::base_name( filename );

            std::ostringstream oss_node;
            oss_node << directory << "/" << file << ".node";
            std::ofstream node( oss_node.str().c_str() );
            node.precision( 16 );

            const GeoModelMesh3D& mesh = geomodel.mesh;
            node << mesh.vertices.nb() << " 3 0 0" << EOL;
            for( auto v : range( mesh.vertices.nb() ) )
            {
                node << v << SPACE << mesh.vertices.vertex( v ) << EOL;
            }

            std::ostringstream oss_ele;
            oss_ele << directory << "/" << file << ".ele";
            std::ofstream ele( oss_ele.str().c_str() );
            std::ostringstream oss_neigh;
            oss_neigh << directory << "/" << file << ".neigh";
            std::ofstream neigh( oss_neigh.str().c_str() );

            ele << mesh.cells.nb() << " 4 1" << EOL;
            neigh << mesh.cells.nb() << " 4" << EOL;
            index_t nb_tet_exported = 0;
            for( auto m : range( geomodel.nb_regions() ) )
            {
                for( auto tet : range( mesh.cells.nb_tet( m ) ) )
                {
                    index_t cell = mesh.cells.tet( m, tet );
                    ele << nb_tet_exported << SPACE
                        << mesh.cells.vertex( ElementLocalVertex( cell, 0 ) )
                        << SPACE
                        << mesh.cells.vertex( ElementLocalVertex( cell, 1 ) )
                        << SPACE
                        << mesh.cells.vertex( ElementLocalVertex( cell, 2 ) )
                        << SPACE
                        << mesh.cells.vertex( ElementLocalVertex( cell, 3 ) )
                        << SPACE << m + 1 << EOL;
                    neigh << nb_tet_exported;
                    for( auto f : range( mesh.cells.nb_facets( tet ) ) )
                    {
                        neigh << SPACE;
                        index_t adj = mesh.cells.adjacent( cell, f );
                        if( adj == NO_ID )
                        {
                            neigh << -1;
                        }
                        else
                        {
                            neigh << adj;
                        }
                    }
                    neigh << EOL;
                    nb_tet_exported++;
                }
            }
            ele << std::flush;
            neigh << std::flush;
            node << std::flush;
        }
    };
}
