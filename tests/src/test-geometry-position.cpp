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

#include <geogram/basic/vecg.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

void verdict( bool condition, std::string test_name )
{
    if( !condition )
    {
        throw RINGMeshException( "TEST", test_name, ": KO" );
    }
    else
    {
        Logger::out( "TEST", test_name, ": OK" );
    }
}

template < index_t DIMENSION >
bool are_almost_equal(
    const vecn< DIMENSION >& vec0, const vecn< DIMENSION >& vec1 )
{
    return ( vec0 - vec1 ).length2() < global_epsilon_sq;
}

void test_point_plane_side()
{
    Logger::out( "TEST", "Test Point-Plane side" );
    Geometry::Plane plane{ { 1., 0., 0. }, { 1., 4., -2. } };

    // Test from each side
    Sign positive_side{ Position::point_side_to_plane(
        { 2., 2., 2. }, plane ) };
    verdict( positive_side == POSITIVE, "True point side positive" );
    Sign negative_side{ Position::point_side_to_plane(
        { -2., -2., -9. }, plane ) };
    verdict( negative_side == NEGATIVE, "True point side negative" );

    // Test on the plane
    Sign on_plane_side{ Position::point_side_to_plane(
        { 1., 6., -6. }, plane ) };
    verdict( on_plane_side == ZERO, "True point side on plane" );

    Logger::out( "TEST", " " );
}

int main()
{
    try
    {
        Logger::out( "TEST", "Test intersection algorithms" );

        test_point_plane_side();
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
