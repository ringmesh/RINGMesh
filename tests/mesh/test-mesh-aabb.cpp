/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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
#include <ringmesh/basic/matrix.h>

#include <ringmesh/geogram_extension/geogram_mesh.h>

#include <ringmesh/mesh/mesh_aabb.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/mesh_set.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

template < index_t DIMENSION >
vecn< DIMENSION > create_vertex( double i, double j );

template <>
vecn< 2 > create_vertex( double i, double j )
{
    return vec2( i, j );
}

template <>
vecn< 3 > create_vertex( double i, double j )
{
    return vec3( i, j, 0 );
}

template < index_t DIMENSION >
void add_vertices( LineMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_vertices( size );
    for( index_t i : range( size ) )
    {
        builder->set_vertex( i, create_vertex< DIMENSION >( i, i + 1 ) );
    }
}

template < index_t DIMENSION >
void add_vertices( SurfaceMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_vertices( size * size );
    index_t id = 0;
    for( index_t i : range( size ) )
    {
        for( index_t j : range( size ) )
        {
            builder->set_vertex( id++, create_vertex< DIMENSION >( i, j ) );
        }
    }
}

template < index_t DIMENSION >
void add_vertices( VolumeMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_vertices( size * size * size );
    index_t id = 0;
    for( index_t i : range( size ) )
    {
        for( index_t j : range( size ) )
        {
            for( index_t k : range( size ) )
            {
                builder->set_vertex( id++, vecn< DIMENSION >( i, j, k ) );
            }
        }
    }
}

template < index_t DIMENSION >
void add_edges( LineMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_edges( size - 1 );
    for( index_t i : range( size - 1 ) )
    {
        builder->set_edge_vertex( EdgeLocalVertex( i, 0 ), i );
        builder->set_edge_vertex( EdgeLocalVertex( i, 1 ), i + 1 );
    }
}

template < index_t DIMENSION >
void add_triangles( SurfaceMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_triangles( ( size - 1 ) * ( size - 1 ) * 2 );
    index_t id = 0;
    for( index_t i : range( size - 1 ) )
    {
        for( index_t j : range( size - 1 ) )
        {
            builder->set_polygon_vertex( { id, 0 }, i * size + j );
            builder->set_polygon_vertex( { id, 1 }, i * size + j + 1 );
            builder->set_polygon_vertex( { id, 2 }, ( i + 1 ) * size + j );
            id++;
            builder->set_polygon_vertex( { id, 0 }, i * size + j + 1 );
            builder->set_polygon_vertex( { id, 1 }, ( i + 1 ) * size + j + 1 );
            builder->set_polygon_vertex( { id, 2 }, ( i + 1 ) * size + j );
            id++;
        }
    }
}

template < index_t DIMENSION >
void add_hexs( VolumeMeshBuilder< DIMENSION >* builder, index_t size )
{
    builder->create_cells(
        ( size - 1 ) * ( size - 1 ) * ( size - 1 ), CellType::HEXAHEDRON );
    index_t id = 0;
    for( index_t i : range( size - 1 ) )
    {
        for( index_t j : range( size - 1 ) )
        {
            for( index_t k : range( size - 1 ) )
            {
                index_t corner = i + j * size + k * size * size;
                builder->set_cell_vertex( ElementLocalVertex( id, 0 ), corner );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 4 ), corner + size * size );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 6 ), corner + size * size + 1 );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 2 ), corner + 1 );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 1 ), corner + size );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 5 ), corner + size * size + size );
                builder->set_cell_vertex( ElementLocalVertex( id, 7 ),
                    corner + size * size + size + 1 );
                builder->set_cell_vertex(
                    ElementLocalVertex( id, 3 ), corner + size + 1 );
                id++;
            }
        }
    }
    builder->connect_cells();
}

template < index_t DIMENSION >
void check_tree( const SurfaceAABBTree< DIMENSION >& tree, index_t size )
{
    double offset = 0.2;
    index_t id = 0;
    for( index_t i : range( size - 1 ) )
    {
        for( index_t j : range( size - 1 ) )
        {
            vecn< DIMENSION > query1 =
                create_vertex< DIMENSION >( i + offset, j + offset );
            index_t triangle1 = NO_ID;
            vecn< DIMENSION > nearest_point1;
            std::tie( triangle1, nearest_point1, std::ignore ) =
                tree.closest_triangle( query1 );
            if( triangle1 != id++ )
            {
                throw RINGMeshException(
                    "TEST", "Not the correct triangle found" );
            }
            if( nearest_point1 != query1 )
            {
                throw RINGMeshException(
                    "TEST", "Not the correct nearest point found" );
            }

            vecn< DIMENSION > query2 =
                create_vertex< DIMENSION >( i + 1 - offset, j + 1 - offset );
            index_t triangle2 = NO_ID;
            vecn< DIMENSION > nearest_point2;
            std::tie( triangle2, nearest_point2, std::ignore ) =
                tree.closest_triangle( query2 );
            if( triangle2 != id++ )
            {
                throw RINGMeshException(
                    "TEST", "Not the correct triangle found" );
            }
            if( nearest_point2 != query2 )
            {
                throw RINGMeshException(
                    "TEST", "Not the correct nearest point found" );
            }
        }
    }

    vecn< DIMENSION > query;
    index_t triangle = NO_ID;
    vecn< DIMENSION > nearest_point;
    std::tie( triangle, nearest_point, std::ignore ) =
        tree.closest_triangle( query );
    if( triangle != 0 )
    {
        throw RINGMeshException( "TEST", "Not the correct triangle found" );
    }
    if( nearest_point != vecn< DIMENSION >() )
    {
        throw RINGMeshException(
            "TEST", "Not the correct nearest point found" );
    }
}

template < index_t DIMENSION >
void create_5_tets_from_hex( VolumeMeshBuilder< DIMENSION >& builder,
    const VolumeMesh< DIMENSION >& mesh_hex,
    index_t hex )
{
    std::vector< index_t > vertices_in_hex( 8 );
    for( index_t v : range( 8 ) )
    {
        vertices_in_hex[v] =
            mesh_hex.cell_vertex( ElementLocalVertex( hex, v ) );
    }
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex, 0 ), vertices_in_hex[0] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex, 1 ), vertices_in_hex[4] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex, 2 ), vertices_in_hex[5] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex, 3 ), vertices_in_hex[6] );

    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 1, 0 ), vertices_in_hex[0] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 1, 1 ), vertices_in_hex[2] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 1, 2 ), vertices_in_hex[3] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 1, 3 ), vertices_in_hex[6] );

    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 2, 0 ), vertices_in_hex[7] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 2, 1 ), vertices_in_hex[6] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 2, 2 ), vertices_in_hex[3] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 2, 3 ), vertices_in_hex[5] );

    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 3, 0 ), vertices_in_hex[1] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 3, 1 ), vertices_in_hex[0] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 3, 2 ), vertices_in_hex[5] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 3, 3 ), vertices_in_hex[3] );

    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 4, 0 ), vertices_in_hex[0] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 4, 1 ), vertices_in_hex[5] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 4, 2 ), vertices_in_hex[6] );
    builder.set_cell_vertex(
        ElementLocalVertex( 5 * hex + 4, 3 ), vertices_in_hex[3] );
}

template < index_t DIMENSION >
void decompose_in_tet( const VolumeMesh< DIMENSION >& hex_mesh,
    VolumeMesh< DIMENSION >& tet_mesh,
    index_t size )
{
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > > builder =
        VolumeMeshBuilder< DIMENSION >::create_builder( tet_mesh );
    builder->create_cells( hex_mesh.nb_cells() * 5, CellType::TETRAHEDRON );
    add_vertices( builder.get(), size );
    for( index_t hex : range( hex_mesh.nb_cells() ) )
    {
        create_5_tets_from_hex( *builder, hex_mesh, hex );
    }
}

template < index_t DIMENSION >
void test_SurfaceAABB()
{
    Logger::out( "TEST", "Test Surface AABB ", DIMENSION, "D" );
    auto mesh = SurfaceMesh< DIMENSION >::create_mesh();
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > builder =
        SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh );

    index_t size = 10;
    add_vertices( builder.get(), size );
    add_triangles( builder.get(), size );

    SurfaceAABBTree< DIMENSION > tree( *mesh );
    check_tree( tree, size );
}

template < index_t DIMENSION >
void test_locate_cell_on_3D_mesh( const VolumeMesh< DIMENSION >& mesh )
{
    for( index_t c : range( mesh.nb_cells() ) )
    {
        vecn< DIMENSION > barycenter = mesh.cell_barycenter( c );
        const VolumeAABBTree< DIMENSION >& aabb3D = mesh.cell_aabb();
        index_t containing_cell = aabb3D.containing_cell( barycenter );
        if( containing_cell != c )
        {
            throw RINGMeshException( "TEST", "Not the correct cell found" );
        }
    }
}

template < index_t DIMENSION >
void test_VolumeAABB()
{
    Logger::out( "TEST", "Test Volume AABB ", DIMENSION, "D" );
    auto mesh_hex = VolumeMesh< DIMENSION >::create_mesh();
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > > builder =
        VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_hex );

    index_t size = 10;
    add_vertices( builder.get(), size );
    add_hexs( builder.get(), size );

    auto mesh_tet = VolumeMesh< DIMENSION >::create_mesh();
    decompose_in_tet( *mesh_hex, *mesh_tet, size );
    test_locate_cell_on_3D_mesh( *mesh_tet );
}

template < index_t DIMENSION >
void test_locate_edge_on_1D_mesh( const LineMesh< DIMENSION >& mesh )
{
    for( index_t e : range( mesh.nb_edges() ) )
    {
        vecn< DIMENSION > barycenter = mesh.edge_barycenter( e );
        const LineAABBTree< DIMENSION >& aabb1D = mesh.edge_aabb();
        index_t closest_edge = NO_ID;
        std::tie( closest_edge, std::ignore, std::ignore ) =
            aabb1D.closest_edge( barycenter );
        if( closest_edge != e )
        {
            throw RINGMeshException( "TEST", "Not the correct edge found" );
        }
    }
}

template < index_t DIMENSION >
void test_LineAABB()
{
    Logger::out( "TEST", "Test Line AABB ", DIMENSION, "D" );
    auto mesh = LineMesh< DIMENSION >::create_mesh();
    std::unique_ptr< LineMeshBuilder< DIMENSION > > builder =
        LineMeshBuilder< DIMENSION >::create_builder( *mesh );

    index_t size = 10;
    add_vertices( builder.get(), size );
    add_edges( builder.get(), size );
    test_locate_edge_on_1D_mesh( *mesh );
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test AABB" );

        test_LineAABB< 2 >();
        test_LineAABB< 3 >();
        test_SurfaceAABB< 2 >();
        test_SurfaceAABB< 3 >();
        test_VolumeAABB< 3 >();
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
