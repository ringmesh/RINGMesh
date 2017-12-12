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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

void build_geomodel( GeoModel3D& geomodel )
{
    GeoModelBuilder3D builder( geomodel );

    vec3 v0( 0, 0, 0 );
    vec3 v1( 1, 0, 0 );
    vec3 v2( 1, 1, 0 );
    vec3 v3( 0, 1, 0 );
    vec3 v4( 0, 0, 1 );
    vec3 v5( 1, 0, 1 );
    vec3 v6( 1, 1, 1 );
    vec3 v7( 0, 1, 1 );

    std::vector< index_t > triangles;
    triangles.push_back( 0 );
    triangles.push_back( 1 );
    triangles.push_back( 2 );
    triangles.push_back( 0 );
    triangles.push_back( 2 );
    triangles.push_back( 3 );

    std::vector< index_t > surface_facet_ptr;
    surface_facet_ptr.push_back( 0 );
    surface_facet_ptr.push_back( 3 );
    surface_facet_ptr.push_back( 6 );

    std::vector< vec3 > vertices( 4 );
    gmme_id id;

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    vertices[3] = v3;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v1;
    vertices[1] = v5;
    vertices[2] = v6;
    vertices[3] = v2;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v4;
    vertices[1] = v5;
    vertices[2] = v6;
    vertices[3] = v7;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v0;
    vertices[1] = v3;
    vertices[2] = v7;
    vertices[3] = v4;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v3;
    vertices[1] = v2;
    vertices[2] = v6;
    vertices[3] = v7;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );

    id = builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v5;
    vertices[3] = v4;
    builder.geometry.set_surface_geometry(
        id.index(), vertices, triangles, surface_facet_ptr );
}

void check_vertex( const vec3& in, const vec3& result )
{
    if( in != result )
    {
        throw RINGMeshException( "Test", "Verification failed" );
    }
}

void test_translate( GeoModel3D& geomodel )
{
    Logger::out( "TEST", "Test translation" );
    vec3 translation_vector( 1., 2.5, -3.5 );
    translate( geomodel, translation_vector );

    const GeoModelMeshVertices3D& vertices = geomodel.mesh.vertices;
    check_vertex( vertices.vertex( 0 ), vec3( 1., 2.5, -3.5 ) );
    check_vertex( vertices.vertex( 1 ), vec3( 2., 2.5, -3.5 ) );
    check_vertex( vertices.vertex( 2 ), vec3( 2., 3.5, -3.5 ) );
    check_vertex( vertices.vertex( 3 ), vec3( 1., 3.5, -3.5 ) );
    check_vertex( vertices.vertex( 4 ), vec3( 2., 2.5, -2.5 ) );
    check_vertex( vertices.vertex( 5 ), vec3( 2., 3.5, -2.5 ) );
    check_vertex( vertices.vertex( 6 ), vec3( 1., 2.5, -2.5 ) );
    check_vertex( vertices.vertex( 7 ), vec3( 1., 3.5, -2.5 ) );
}

void test_rotation( GeoModel3D& geomodel )
{
    Logger::out( "TEST", "Test rotation" );
    vec3 origin( 1., 2.5, -3.5 );
    vec3 axis( 0, 0, 1 );
    rotate( geomodel, origin, axis, 90, true );

    const GeoModelMeshVertices3D& vertices = geomodel.mesh.vertices;
    check_vertex( vertices.vertex( 0 ), vec3( 1., 2.5, -3.5 ) );
    check_vertex( vertices.vertex( 1 ), vec3( 1., 3.5, -3.5 ) );
    check_vertex( vertices.vertex( 2 ), vec3( 0., 3.5, -3.5 ) );
    check_vertex( vertices.vertex( 3 ), vec3( 0., 2.5, -3.5 ) );
    check_vertex( vertices.vertex( 4 ), vec3( 1., 3.5, -2.5 ) );
    check_vertex( vertices.vertex( 5 ), vec3( 0., 3.5, -2.5 ) );
    check_vertex( vertices.vertex( 6 ), vec3( 1., 2.5, -2.5 ) );
    check_vertex( vertices.vertex( 7 ), vec3( 0., 2.5, -2.5 ) );
}

void check_matrices(
    const GEO::Matrix< 4, double >& lhs, const GEO::Matrix< 4, double >& rhs )
{
    for( index_t mat_i : range( 4 ) )
    {
        for( index_t mat_j : range( 4 ) )
        {
            double diff = lhs( mat_i, mat_j ) - rhs( mat_i, mat_j );
            if( std::fabs( diff ) > global_epsilon )
            {
                throw RINGMeshException( "Test", "Error in rotation matrix" );
            }
        }
    }
}

void test_rotation_matrix()
{
    const vec3 origin( 0, 0, 0 );
    const double pi = M_PI;
    const double step = 0.1;

    GEO::Matrix< 4, double > rot_mat_degree;
    GEO::Matrix< 4, double > rot_mat_radian;
    GEO::Matrix< 4, double > result;
    result( 0, 3 ) = 0;
    result( 1, 3 ) = 0;
    result( 2, 3 ) = 0;
    result( 3, 0 ) = 0;
    result( 3, 1 ) = 0;
    result( 3, 2 ) = 0;
    result( 3, 3 ) = 1;

    // Tests rotation along x axis
    vec3 axis( 1, 0, 0 );
    for( double angle = 0.; angle <= 360.; angle += step )
    {
        rot_mat_degree =
            rotation_matrix_about_arbitrary_axis( origin, axis, angle, true );
        double angle_rad = angle * pi / 180.;
        rot_mat_radian = rotation_matrix_about_arbitrary_axis(
            origin, axis, angle_rad, false );
        result( 0, 0 ) = 1;
        result( 0, 1 ) = 0;
        result( 0, 2 ) = 0;

        result( 1, 0 ) = 0;
        result( 1, 1 ) = std::cos( angle_rad );
        result( 1, 2 ) = -std::sin( angle_rad );

        result( 2, 0 ) = 0;
        result( 2, 1 ) = std::sin( angle_rad );
        result( 2, 2 ) = std::cos( angle_rad );

        check_matrices( rot_mat_degree, result );
        check_matrices( rot_mat_radian, result );
    }

    // Tests rotation along y axis
    axis = vec3( 0, 1, 0 );
    for( double angle = 0.; angle <= 360.; angle += step )
    {
        rot_mat_degree =
            rotation_matrix_about_arbitrary_axis( origin, axis, angle, true );
        double angle_rad = angle * pi / 180.;
        rot_mat_radian = rotation_matrix_about_arbitrary_axis(
            origin, axis, angle_rad, false );
        result( 0, 0 ) = std::cos( angle_rad );
        result( 0, 1 ) = 0;
        result( 0, 2 ) = std::sin( angle_rad );

        result( 1, 0 ) = 0;
        result( 1, 1 ) = 1;
        result( 1, 2 ) = 0;

        result( 2, 0 ) = -std::sin( angle_rad );
        result( 2, 1 ) = 0;
        result( 2, 2 ) = std::cos( angle_rad );

        check_matrices( rot_mat_degree, result );
        check_matrices( rot_mat_radian, result );
    }

    // Tests rotation along z axis
    axis = vec3( 0, 0, 1 );
    for( double angle = 0.; angle <= 360.; angle += step )
    {
        rot_mat_degree =
            rotation_matrix_about_arbitrary_axis( origin, axis, angle, true );
        double angle_rad = angle * pi / 180.;
        rot_mat_radian = rotation_matrix_about_arbitrary_axis(
            origin, axis, angle_rad, false );
        result( 0, 0 ) = std::cos( angle_rad );
        result( 0, 1 ) = -std::sin( angle_rad );
        result( 0, 2 ) = 0;

        result( 1, 0 ) = std::sin( angle_rad );
        result( 1, 1 ) = std::cos( angle_rad );
        result( 1, 2 ) = 0;

        result( 2, 0 ) = 0;
        result( 2, 1 ) = 0;
        result( 2, 2 ) = 1;

        check_matrices( rot_mat_degree, result );
        check_matrices( rot_mat_radian, result );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        GeoModel3D geomodel;
        build_geomodel( geomodel );
        test_translate( geomodel );
        test_rotation_matrix();
        test_rotation( geomodel );
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
