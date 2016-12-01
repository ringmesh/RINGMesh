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

void add_vertices( Mesh3DBuilder* builder, index_t size )
{
    builder->create_vertices( size * size * size ) ;
    index_t id = 0 ;
    for( index_t i = 0; i < size; i++ ) {
        for( index_t j = 0; j < size; j++ ) {
            for( index_t k = 0; k < size; k++ ) {
                builder->set_vertex( id++, vec3( i, j, k ) ) ;
            }
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

void add_hexs( Mesh3DBuilder* builder, index_t size )
{
    builder->create_cells( ( size - 1 ) * ( size - 1 ) * ( size - 1 ),
        GEO::MESH_HEX ) ;
    index_t id = 0 ;
    for( index_t i = 0; i < ( size - 1 ); i++ ) {
        for( index_t j = 0; j < ( size - 1 ); j++ ) {
            for( index_t k = 0; k < ( size - 1 ); k++ ) {
                index_t corner = i + j * size + k * size * size ;
                builder->set_cell_vertex( id, 0, corner ) ;
                builder->set_cell_vertex( id, 4, corner + size * size ) ;
                builder->set_cell_vertex( id, 6, corner + size * size + 1 ) ;
                builder->set_cell_vertex( id, 2, corner + 1 ) ;
                builder->set_cell_vertex( id, 1, corner + size ) ;
                builder->set_cell_vertex( id, 5, corner + size * size + size ) ;
                builder->set_cell_vertex( id, 7, corner + size * size + size + 1 ) ;
                builder->set_cell_vertex( id, 3, corner + size + 1 ) ;
                id++ ;
            }

        }
    }
    builder->connect_cells() ;
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

void add_hex_C5( Mesh3DBuilder& builder, const GeogramMesh3D& mesh_hex, index_t hex )
{
    std::vector< index_t > vertices_in_hex( 8 ) ;
    for( index_t v = 0; v < 8; v++ ) {
        vertices_in_hex[v] = mesh_hex.cell_vertex( hex, v ) ;
    }
    builder.set_cell_vertex( 5 * hex, 0, vertices_in_hex[0] ) ;
    builder.set_cell_vertex( 5 * hex, 1, vertices_in_hex[4] ) ;
    builder.set_cell_vertex( 5 * hex, 2, vertices_in_hex[5] ) ;
    builder.set_cell_vertex( 5 * hex, 3, vertices_in_hex[6] ) ;

    builder.set_cell_vertex( 5 * hex + 1, 0, vertices_in_hex[0] ) ;
    builder.set_cell_vertex( 5 * hex + 1, 1, vertices_in_hex[2] ) ;
    builder.set_cell_vertex( 5 * hex + 1, 2, vertices_in_hex[3] ) ;
    builder.set_cell_vertex( 5 * hex + 1, 3, vertices_in_hex[6] ) ;

    builder.set_cell_vertex( 5 * hex + 2, 0, vertices_in_hex[7] ) ;
    builder.set_cell_vertex( 5 * hex + 2, 1, vertices_in_hex[6] ) ;
    builder.set_cell_vertex( 5 * hex + 2, 2, vertices_in_hex[3] ) ;
    builder.set_cell_vertex( 5 * hex + 2, 3, vertices_in_hex[5] ) ;

    builder.set_cell_vertex( 5 * hex + 3, 0, vertices_in_hex[1] ) ;
    builder.set_cell_vertex( 5 * hex + 3, 1, vertices_in_hex[0] ) ;
    builder.set_cell_vertex( 5 * hex + 3, 2, vertices_in_hex[5] ) ;
    builder.set_cell_vertex( 5 * hex + 3, 3, vertices_in_hex[3] ) ;

    builder.set_cell_vertex( 5 * hex + 4, 0, vertices_in_hex[0] ) ;
    builder.set_cell_vertex( 5 * hex + 4, 1, vertices_in_hex[5] ) ;
    builder.set_cell_vertex( 5 * hex + 4, 2, vertices_in_hex[6] ) ;
    builder.set_cell_vertex( 5 * hex + 4, 3, vertices_in_hex[3] ) ;
}
void decompose_in_tet(
    const GeogramMesh3D& hex_mesh,
    GeogramMesh3D& tet_mesh,
    index_t size )
{
    Mesh3DBuilder* builder = tet_mesh.get_mesh3d_builder() ;
    builder->create_cells( hex_mesh.nb_cells() * 5, GEO::MESH_TET ) ;
    add_vertices( builder, size ) ;
    for( index_t hex = 0; hex < hex_mesh.nb_cells(); hex++ ) {
        add_hex_C5( *builder, hex_mesh, hex ) ;
    }

}
void test_AABB2D()
{
    Logger::out( "TEST" ) << "Test AABB 2D" << std::endl ;
    GeogramMesh2D geogram_mesh ;
    Mesh2DBuilder* builder = geogram_mesh.get_mesh2d_builder() ;

    index_t size = 10 ;
    add_vertices( builder, size ) ;
    add_triangles( builder, size ) ;

    AABBTree2D tree( geogram_mesh ) ;
    tree.save_tree( "tree" ) ;
    check_tree( tree, size ) ;

}

void test_AABB3D()
{
    Logger::out( "TEST" ) << "Test AABB 3D" << std::endl ;
    GeogramMesh3D geogram_mesh_hex ;
    Mesh3DBuilder* builder = geogram_mesh_hex.get_mesh3d_builder() ;

    index_t size = 3 ;
    add_vertices( builder, size ) ;
    add_hexs( builder, size ) ;
    GeogramMesh3D geogram_mesh_tet ;

    decompose_in_tet( geogram_mesh_hex, geogram_mesh_tet, size ) ;
    for( index_t c = 0; c < geogram_mesh_tet.nb_cells(); c++ ) {
        vec3 barycenter = geogram_mesh_tet.cell_barycenter( c ) ;
        const AABBTree3D& aabb3D = geogram_mesh_tet.cells_aabb() ;
        index_t containing_cell = aabb3D.containing_cell( barycenter ) ;
        if( containing_cell != c ) {
            throw RINGMeshException( "TEST", "Not the correct cell found" ) ;
        }
    }
}

int main()
{
    using namespace RINGMesh ;

    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test AABB" << std::endl ;
        test_AABB2D() ;
        test_AABB3D() ;

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
