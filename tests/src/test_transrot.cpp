/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geomodel/geomodel_api.h>
#include <ringmesh/geomodel/geomodel_builder.h>


/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh ;

void build_geomodel( GeoModel& geomodel )
{
    GeoModelBuilder builder( geomodel ) ;

    vec3 v0( 0, 0, 0 ) ;
    vec3 v1( 1, 0, 0 ) ;
    vec3 v2( 1, 1, 0 ) ;
    vec3 v3( 0, 1, 0 ) ;
    vec3 v4( 0, 0, 1 ) ;
    vec3 v5( 1, 0, 1 ) ;
    vec3 v6( 1, 1, 1 ) ;
    vec3 v7( 0, 1, 1 ) ;

    std::vector< index_t > triangles ;
    triangles.push_back( 0 ) ;
    triangles.push_back( 1 ) ;
    triangles.push_back( 2 ) ;
    triangles.push_back( 0 ) ;
    triangles.push_back( 2 ) ;
    triangles.push_back( 3 ) ;

    std::vector< index_t > surface_facet_ptr ;
    surface_facet_ptr.push_back( 0 ) ;
    surface_facet_ptr.push_back( 3 ) ;
    surface_facet_ptr.push_back( 6 ) ;

    std::vector< vec3 > vertices( 4 ) ;
    gme_t id ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v0 ;
    vertices[1] = v1 ;
    vertices[2] = v2 ;
    vertices[3] = v3 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v1 ;
    vertices[1] = v5 ;
    vertices[2] = v6 ;
    vertices[3] = v2 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v4 ;
    vertices[1] = v5 ;
    vertices[2] = v6 ;
    vertices[3] = v7 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v0 ;
    vertices[1] = v3 ;
    vertices[2] = v7 ;
    vertices[3] = v4 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v3 ;
    vertices[1] = v2 ;
    vertices[2] = v6 ;
    vertices[3] = v7 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;

    id = builder.create_mesh_entity< Surface >() ;
    vertices[0] = v0 ;
    vertices[1] = v1 ;
    vertices[2] = v5 ;
    vertices[3] = v4 ;
    builder.set_surface_geometry( id.index, vertices, triangles,
        surface_facet_ptr ) ;
}

void check_vertex( const vec3& in, const vec3& result )
{
    if( in != result ) {
        throw RINGMeshException( "Test", "Verification failed" ) ;
    }
}

void test_translate( GeoModel& geomodel )
{
    Logger::out( "TEST" ) << "Test translation" << std::endl ;
    vec3 translation_vector( 1., 2.5, -3.5 ) ;
    translate( geomodel, translation_vector ) ;

    const GeoModelMeshVertices& vertices = geomodel.mesh.vertices ;
    check_vertex( vertices.vertex( 0 ), vec3( 1., 2.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 1 ), vec3( 2., 2.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 2 ), vec3( 2., 3.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 3 ), vec3( 1., 3.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 4 ), vec3( 2., 2.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 5 ), vec3( 2., 3.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 6 ), vec3( 1., 2.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 7 ), vec3( 1., 3.5, -2.5 ) ) ;
}


void test_rotation( GeoModel& geomodel )
{
    Logger::out( "TEST" ) << "Test rotation" << std::endl ;
    vec3 origin( 1., 2.5, -3.5 ) ;
    vec3 axis( 0, 0, 1 ) ;
    rotate( geomodel, origin, axis, 90, true ) ;

    const GeoModelMeshVertices& vertices = geomodel.mesh.vertices ;
    check_vertex( vertices.vertex( 0 ), vec3( 1., 2.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 1 ), vec3( 1., 3.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 2 ), vec3( 0., 3.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 3 ), vec3( 0., 2.5, -3.5 ) ) ;
    check_vertex( vertices.vertex( 4 ), vec3( 1., 3.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 5 ), vec3( 0., 3.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 6 ), vec3( 1., 2.5, -2.5 ) ) ;
    check_vertex( vertices.vertex( 7 ), vec3( 0., 2.5, -2.5 ) ) ;
}

int main()
{
    using namespace RINGMesh ;

    try {
        default_configure() ;

        GeoModel geomodel ;
        build_geomodel( geomodel ) ;
        test_translate( geomodel ) ;
        test_rotation( geomodel ) ;

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
