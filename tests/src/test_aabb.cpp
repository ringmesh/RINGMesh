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

#include <vector>

#include <ringmesh/mesh/aabb.h>
#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh ;

void add_vertices( Mesh2DBuilder* builder, index_t size )
{
    builder->create_vertices( size * size ) ;
    index_t id = 0 ;
    for( index_t i = 0; i < size; i++ ) {
        for( index_t j = 0; j < size; j++ ) {
            builder->set_vertex( id++, vec3( i, j, 0 ) ) ;
        }
    }
}

void add_triangles( Mesh2DBuilder* builder, index_t size )
{
    builder->create_facet_triangles( ( size - 1 ) * ( size - 1 ) * 2 ) ;
    index_t id = 0 ;
    for( index_t i = 0; i < size - 1; i++ ) {
        for( index_t j = 0; j < size - 1; j++ ) {
            builder->set_facet_vertex( id, 0, i * size + j ) ;
            builder->set_facet_vertex( id, 1, i * size + j + 1 ) ;
            builder->set_facet_vertex( id, 2, ( i + 1 ) * size + j ) ;
            id++ ;
            builder->set_facet_vertex( id, 0, i * size + j + 1 ) ;
            builder->set_facet_vertex( id, 1, ( i + 1 ) * size + j + 1 ) ;
            builder->set_facet_vertex( id, 2, ( i + 1 ) * size + j ) ;
            id++ ;
        }
    }
}

void check_tree( const AABBTree2D& tree, index_t size )
{
    double offset = 0.2 ;
    index_t id = 0 ;
    for( index_t i = 0; i < size - 1; i++ ) {
        for( index_t j = 0; j < size - 1; j++ ) {
            vec3 query1( i + offset, j + offset, 0 ) ;
            vec3 nearest_point1 ;
            double distance1 ;
            index_t triangle1 = tree.closest_triangle( query1, nearest_point1,
                distance1 ) ;
            if( triangle1 != id++ ) {
                throw RINGMeshException( "TEST", "Not the correct triangle found" ) ;
            }
            if( nearest_point1 != vec3( i + offset, j + offset, 0 ) ) {
                throw RINGMeshException( "TEST",
                    "Not the correct nearest point found" ) ;
            }

            vec3 query2( i + 1 - offset, j + 1 - offset, offset ) ;
            vec3 nearest_point2 ;
            double distance2 ;
            index_t triangle2 = tree.closest_triangle( query2, nearest_point2,
                distance2 ) ;
            if( triangle2 != id++ ) {
                throw RINGMeshException( "TEST", "Not the correct triangle found" ) ;
            }
            if( nearest_point2 != vec3( i + 1 - offset, j + 1 - offset, 0 ) ) {
                throw RINGMeshException( "TEST",
                    "Not the correct nearest point found" ) ;
            }
        }
    }

    vec3 query( 0, 0, 0 ) ;
    vec3 nearest_point ;
    double distance ;
    index_t triangle = tree.closest_triangle( query, nearest_point, distance ) ;
    if( triangle != 0 ) {
        throw RINGMeshException( "TEST", "Not the correct triangle found" ) ;
    }
    if( nearest_point != vec3( 0, 0, 0 ) ) {
        throw RINGMeshException( "TEST", "Not the correct nearest point found" ) ;
    }
}

int main()
{
    using namespace RINGMesh ;

    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test AABB" << std::endl ;
        GeogramMesh2D geogram_mesh ;
        Mesh2DBuilder* builder = geogram_mesh.get_mesh2d_builder() ;

        index_t size = 10 ;
        add_vertices( builder, size ) ;
        add_triangles( builder, size ) ;

        AABBTree2D tree( geogram_mesh ) ;
        tree.save_tree( "tree" ) ;
        check_tree( tree, size ) ;

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
