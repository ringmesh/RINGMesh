/* Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>

/*!
 * @authors Arnaud Botella and Marine Hubert
 */

using namespace RINGMesh;

void test_triangle_barycentric_coordinates_3D()
{
    Logger::out( "TEST", "Test triangle barycentric coordinates 3D" );
    vec3 p0( 0, 0, 0 );
    vec3 p1( 1, 0, 0 );
    vec3 p2( 0, 1, 0 );
    vec3 p3( 0.9999999, 0 ,0);

    bool not_nul_area;
    std::array< double, 3 > lambdas;
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec3( 0.25, 0.25, 0 ), p0, p1, p2 );
    if( lambdas[0] != 0.5 || lambdas[1] != 0.25 || lambdas[2] != 0.25 )
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 3D" );
    }
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec3( 0.5, 0.5, 0 ), p0, p1, p2 );
    if( lambdas[0] != 0 || lambdas[1] != 0.5 || lambdas[2] != 0.5 )
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 3D" );
    }
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec3( 1, 1, 0 ), p0, p1, p2 );
    if( lambdas[0] != -1 || lambdas[1] != 1 || lambdas[2] != 1 )
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 3D" );
    }
    std::tie( std::ignore, lambdas ) =
               triangle_barycentric_coordinates( vec3( 0, 0, 0), p0, p1, p2 );
    if( lambdas[0] != 1 || lambdas[1] != 0 || lambdas[2] != 0)
    {
         throw RINGMeshException(
             "TEST", "Error in triangle barycentric coordinates 3D" );
    }
    std::tie( not_nul_area, lambdas ) =
               triangle_barycentric_coordinates( vec3( 0, 0, 0), p0, p1, p3 );
    if( not_nul_area )
    {
          throw RINGMeshException(
             "TEST", "Error in triangle barycentric coordinates 3D" );
    }
}

void test_triangle_barycentric_coordinates_2D()
{
    Logger::out( "TEST", "Test triangle barycentric coordinates 2D" );
    vec2 p0( 0, 0);
    vec2 p1( 1, 0);
    vec2 p2( 0, 1);
    vec2 p3( 0.9999999, 0 ); //triangle 0,1, 3 degenere

    bool not_nul_area;
    std::array< double, 3 > lambdas;
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec2( 0.25, 0.25), p0, p1, p2 );
    if( lambdas[0] != 0.5 || lambdas[1] != 0.25 || lambdas[2] != 0.25 )
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 2D" );
    }
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec2( 0.5, 0.5), p0, p1, p2 );
    if( lambdas[0] != 0 || lambdas[1] != 0.5 || lambdas[2] != 0.5)
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 2D" );
    }
    std::tie( std::ignore, lambdas ) =
        triangle_barycentric_coordinates( vec2( 1, 1), p0, p1, p2 );
    if( lambdas[0] != -1 || lambdas[1] != 1 || lambdas[2] != 1)
    {
        throw RINGMeshException(
            "TEST", "Error in triangle barycentric coordinates 2D" );
    }
    std::tie( std::ignore, lambdas ) =
            triangle_barycentric_coordinates( vec2( 0, 0), p0, p1, p2 );
    if( lambdas[0] != 1 || lambdas[1] != 0 || lambdas[2] != 0)
    {
         throw RINGMeshException(
             "TEST", "Error in triangle barycentric coordinates 2D" );
    }
    std::tie( not_nul_area, lambdas ) =
            triangle_barycentric_coordinates( vec2( 0, 0), p0, p1, p3 );
    if( not_nul_area )
    {
          throw RINGMeshException(
             "TEST", "Error in triangle barycentric coordinates 2D" );
    }
}

void test_point_plane_distance()
{
    Logger::out( "TEST", "Test point plane distance" );

    vec3 test0{ 1, 1, 1 };
    Geometry::Plane plane0{ { 0, 0, 2 }, { 0, 0, 0 } };
    vec3 projected0;
    std::tie( std::ignore, projected0 ) =
        Distance::point_to_plane( test0, plane0 );
    if( projected0 != vec3{ 1, 1, 0 } )
    {
        throw RINGMeshException( "TEST", "Error in point plane distance" );
    }

    vec3 test1{ 0, 0.5, 1 };
    Geometry::Plane plane1{ { 1, 0, 0 }, { 1, 1, 1 } };
    vec3 projected1;
    std::tie( std::ignore, projected1 ) =
        Distance::point_to_plane( test1, plane1 );
    if( projected1 != vec3{ 1, 0.5, 1 } )
    {
        throw RINGMeshException( "TEST", "Error in point plane distance" );
    }
}

void test_normalize_perp_scalar()
{
	Logger::out( "TEST", "Test scalar product" );
	vec2 v{1.,1.};
	vec2 vect_n { normalized_perp(v) };
	double scalar_product { v.x*vect_n.x + v.y*vect_n.y };
	if (scalar_product != 0)
		{
			throw RINGMeshException( "TEST", "Error in scalar product" );
		}
	}
void test_dot_perp()
{
	Logger::out( "TEST", "Test dot product" );
		vec2 v{1,1};
		vec2 vect_n {normalized_perp(v)};
		double scalar_product {dot(v,vect_n)};
			if (scalar_product != 0)
			{
				throw RINGMeshException( "TEST", "Error in dot product" );
	}
}

void test_operators()
{
	Logger::out( "TEST", "Test operators '==' and '!='" );
	vec3 u{1,1,1};
	vec3 v{1,0.99999999,1};
	vec3 w{1,1,1};

	if (operator==(v,u))
	{
		throw RINGMeshException( "TEST", "Error in operator ==" );
	}
	if (!operator==(u,w))
	{
		throw RINGMeshException( "TEST", "Error in operator ==" );
	}
	if (!operator!=(v,u))
	{
		throw RINGMeshException( "TEST", "Error in operator !=" );
	}
	if (operator!=(u,w))
	{
		throw RINGMeshException( "TEST", "Error in operator !=" );
	}
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test geometric tools" );

        test_triangle_barycentric_coordinates_3D();
        test_triangle_barycentric_coordinates_2D();
        test_point_plane_distance();
        test_normalize_perp_scalar();
        test_dot_perp();
        test_operators();

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
