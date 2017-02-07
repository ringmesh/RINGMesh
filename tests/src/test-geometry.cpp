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

#include <ringmesh/basic/geometry.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh ;

void test_triangle_barycentric_coordinates()
{
    Logger::out( "TEST" ) << "Test triangle barycentric coordinates" << std::endl ;
    vec3 p0( 0, 0, 0 ) ;
    vec3 p1( 1, 0, 0 ) ;
    vec3 p2( 0, 1, 0 ) ;

    double lambda[3] ;
    triangle_barycentric_coordinates( vec3( 0.25, 0.25, 0 ), p0, p1, p2, lambda ) ;
    if( lambda[0] != 0.5 || lambda[1] != 0.25 || lambda[2] != 0.25 ) {
        throw RINGMeshException( "TEST",
            "Error in triangle barycentric coordinates" ) ;
    }
    triangle_barycentric_coordinates( vec3( 0.5, 0.5, 0 ), p0, p1, p2, lambda ) ;
    if( lambda[0] != 0 || lambda[1] != 0.5 || lambda[2] != 0.5 ) {
        throw RINGMeshException( "TEST",
            "Error in triangle barycentric coordinates" ) ;
    }
    triangle_barycentric_coordinates( vec3( 1, 1, 0 ), p0, p1, p2, lambda ) ;
    if( lambda[0] != -1 || lambda[1] != 1 || lambda[2] != 1 ) {
        throw RINGMeshException( "TEST",
            "Error in triangle barycentric coordinates" ) ;
    }
}

void test_point_plane_projection()
{
    Logger::out( "TEST" ) << "Test point plane projection" << std::endl ;

    vec3 test0( 1, 1, 1 ) ;
    vec3 N0( 0, 0, 2 ) ;
    vec3 O0( 0, 0, 0 ) ;
    vec3 projected0 ;
    point_plane_projection( test0, N0, O0, projected0 ) ;
    if( projected0 != vec3( 1, 1, 0 ) ) {
        throw RINGMeshException( "TEST", "Error in point plane projection" ) ;
    }

    vec3 test1( 0, 0.5, 1 ) ;
    vec3 N1( 1, 0, 0 ) ;
    vec3 O1( 1, 1, 1 ) ;
    vec3 projected1 ;
    point_plane_projection( test1, N1, O1, projected1 ) ;
    if( projected1 != vec3( 1, 0.5, 1 ) ) {
        throw RINGMeshException( "TEST", "Error in point plane projection" ) ;
    }
}

int main()
{
    using namespace RINGMesh ;

    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test geometric tools" << std::endl ;

        test_triangle_barycentric_coordinates() ;
        test_point_plane_projection() ;

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;
}
