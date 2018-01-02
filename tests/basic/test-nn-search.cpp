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

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/nn_search.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

template < index_t DIMENSION >
void test_nn_search()
{
    std::vector< vecn< DIMENSION > > hardcoded_unique_vertices( 4 );
    for( index_t p : range( hardcoded_unique_vertices.size() ) )
    {
        vecn< DIMENSION >& point = hardcoded_unique_vertices[p];
        for( index_t i : range( DIMENSION ) )
        {
            point[i] = p;
        }
    }

    std::vector< vecn< DIMENSION > > vertices( 7 );
    vertices[0] = hardcoded_unique_vertices[0];
    vertices[1] = hardcoded_unique_vertices[1];
    vertices[2] = hardcoded_unique_vertices[0];
    vertices[3] = hardcoded_unique_vertices[2];
    vertices[4] = hardcoded_unique_vertices[1];
    vertices[5] = hardcoded_unique_vertices[3];
    vertices[6] = hardcoded_unique_vertices[0];

    std::vector< index_t > hardcoded_index_map( 7 );

    hardcoded_index_map[0] = 0;
    hardcoded_index_map[1] = 1;
    hardcoded_index_map[2] = 0;
    hardcoded_index_map[3] = 2;
    hardcoded_index_map[4] = 1;
    hardcoded_index_map[5] = 3;
    hardcoded_index_map[6] = 0;

    NNSearch< DIMENSION > nn_search( vertices );
    std::vector< vecn< DIMENSION > > unique_vertices;
    std::vector< index_t > index_map;
    std::tie( std::ignore, index_map, unique_vertices ) =
        nn_search.get_colocated_index_mapping_and_unique_points(
            global_epsilon );
    for( index_t i : range( index_map.size() ) )
    {
        if( index_map[i] != hardcoded_index_map[i] )
        {
            throw RINGMeshException( "TEST", "Index map found is wrong" );
        }
    }

    for( index_t v : range( unique_vertices.size() ) )
    {
        if( unique_vertices[v] != hardcoded_unique_vertices[v] )
        {
            throw RINGMeshException(
                "TEST", "Unique vertices found are wrong" );
        }
    }
}

int main()
{
    try
    {
        Logger::out( "TEST", "Test NNsearch 2D" );
        test_nn_search< 2 >();
        Logger::out( "TEST", "Test NNsearch 3D" );
        test_nn_search< 3 >();
    }
    catch( const RINGMeshException& e )
    {
        Logger::err( e.category(), e.what() );
        return 1;
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
