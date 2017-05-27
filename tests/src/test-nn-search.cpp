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

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <ringmesh/basic/geometry.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

void test_nn_search_ringmesh()
{
    Logger::out( "TEST", "Test NNsearch RINGMesh" );
    vec3 p1( 0, 0, 0 );
    vec3 p2( 1, 1, 1 );
    vec3 p3( 2, 2, 2 );
    vec3 p4( 3, 3, 3 );

    std::vector< vec3 > vertices( 7 );
    vertices[0] = p1;
    vertices[1] = p2;
    vertices[2] = p1;
    vertices[3] = p3;
    vertices[4] = p2;
    vertices[5] = p4;
    vertices[6] = p1;

    GEO::vector< vec3 > hardcoded_unique_vertices( 4 );
    GEO::vector< index_t > hardcoded_index_map( 7 );

    hardcoded_index_map[0] = 0;
    hardcoded_index_map[1] = 1;
    hardcoded_index_map[2] = 0;
    hardcoded_index_map[3] = 2;
    hardcoded_index_map[4] = 1;
    hardcoded_index_map[5] = 3;
    hardcoded_index_map[6] = 0;

    hardcoded_unique_vertices[0] = p1;
    hardcoded_unique_vertices[1] = p2;
    hardcoded_unique_vertices[2] = p3;
    hardcoded_unique_vertices[3] = p4;

    NNSearch nn_search( vertices );
    std::vector< vec3 > unique_vertices;
    std::vector< index_t > index_map;
    nn_search.get_colocated_index_mapping( global_epsilon, index_map,
        unique_vertices );
    for( index_t i = 0; i < index_map.size(); i++ ) {
        if( index_map[i] != hardcoded_index_map[i] ) {
            throw RINGMeshException( "TEST", "Index map found is wrong" );
        }
    }

    for( index_t v = 0; v < unique_vertices.size(); v++ ) {
        if( unique_vertices[v] != hardcoded_unique_vertices[v] ) {
            throw RINGMeshException( "TEST", "unique vertices found are wrong" );
        }
    }
}

int main()
{
    try {
        default_configure();

        test_nn_search_ringmesh();

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
