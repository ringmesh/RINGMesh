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

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>

/*!
 * @author Benjamin Chauvin
 */

using namespace RINGMesh;

namespace
{
    void point_set_mesh_connected_component_test()
    {
        auto point_set = PointSetMesh3D::create_mesh();
        auto point_set_builder =
            PointSetMeshBuilder3D::create_builder( *point_set );

        index_t nb_connected_components{ NO_ID };
        std::vector< index_t > connected_components;
        std::tie( nb_connected_components, connected_components ) =
            point_set->connected_components();
        if( nb_connected_components != 0 || connected_components.size() != 0 )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Point set should have 0 connected component." );
        }
        const index_t max_itr{ 5 };
        for( auto i : range( 1, max_itr + 1 ) )
        {
            point_set_builder->create_vertex( { static_cast< double >( i ),
                static_cast< double >( i ), static_cast< double >( i ) } );
            nb_connected_components = NO_ID;
            connected_components.clear();
            std::tie( nb_connected_components, connected_components ) =
                point_set->connected_components();
            if( nb_connected_components != i
                || connected_components.size() != i )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Point set should have ", i, " connected components." );
            }

            std::vector< index_t > solution( i );
            std::iota( solution.begin(), solution.end(), 0 );
            if( connected_components != solution )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Point set connected components are not correct." );
            }
        }
    }

    void line_mesh_connected_component_test()
    {
        auto line_mesh = LineMesh3D::create_mesh();
        auto line_mesh_builder =
            LineMeshBuilder3D::create_builder( *line_mesh );

        index_t nb_connected_components{ NO_ID };
        std::vector< index_t > connected_components;
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        if( nb_connected_components != 0 || connected_components.size() != 0 )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 0 connected component." );
        }

        nb_connected_components = NO_ID;
        auto v0 = line_mesh_builder->create_vertex( { 0., 0., 0. } );
        auto v1 = line_mesh_builder->create_vertex( { 0., 0., 1. } );
        line_mesh_builder->create_edge( v0, v1 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        std::vector< index_t > solution( 1, 0 );
        if( nb_connected_components != 1 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 1 connected component with index 0." );
        }

        nb_connected_components = NO_ID;
        auto v2 = line_mesh_builder->create_vertex( { 0., 0., 2. } );
        line_mesh_builder->create_edge( v1, v2 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution.push_back( 0 );
        if( nb_connected_components != 1 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 1 connected component with index 0." );
        }

        nb_connected_components = NO_ID;
        auto v3 = line_mesh_builder->create_vertex( { 0., 0., 3. } );
        auto v4 = line_mesh_builder->create_vertex( { 0., 0., 4. } );
        line_mesh_builder->create_edge( v3, v4 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution.push_back( 1 );
        if( nb_connected_components != 2 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 2 connected components: 0 0 1." );
        }

        nb_connected_components = NO_ID;
        auto v5 = line_mesh_builder->create_vertex( { 1., 0., 0. } );
        line_mesh_builder->create_edge( v5, v0 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution.push_back( 0 );
        if( nb_connected_components != 2 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 2 connected components: 0 0 1 0." );
        }

        nb_connected_components = NO_ID;
        auto v6 = line_mesh_builder->create_vertex( { 0., 0., 5. } );
        auto v7 = line_mesh_builder->create_vertex( { 0., 0., 6. } );
        line_mesh_builder->create_edge( v6, v7 );
        line_mesh_builder->create_edge( v4, v6 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution.push_back( 1 );
        solution.push_back( 1 );
        if( nb_connected_components != 2 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 2 connected components: 0 0 1 0 1 1." );
        }

        nb_connected_components = NO_ID;
        line_mesh_builder->create_edge( v2, v3 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution[2] = 0;
        solution[4] = 0;
        solution[5] = 0;
        solution.push_back( 0 );
        if( nb_connected_components != 1 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Line mesh should have 1 connected component: 0 0 0 0 0 0 0." );
        }

        nb_connected_components = NO_ID;
        auto v8 = line_mesh_builder->create_vertex( { 0., 0., 7. } );
        auto v9 = line_mesh_builder->create_vertex( { 0., 0., 8. } );
        line_mesh_builder->create_edge( v8, v9 );
        auto v10 = line_mesh_builder->create_vertex( { 0., 0., 9. } );
        auto v11 = line_mesh_builder->create_vertex( { 0., 0., 10. } );
        line_mesh_builder->create_edge( v10, v11 );
        line_mesh_builder->create_edge( v7, v8 );
        std::tie( nb_connected_components, connected_components ) =
            line_mesh->connected_components();
        solution.push_back( 0 );
        solution.push_back( 1 );
        solution.push_back( 0 );
        if( nb_connected_components != 2 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test", "Line mesh should have 2 "
                                                      "connected components: 0 "
                                                      "0 0 0 0 0 0 0 1 0." );
        }
    }

    void surface_mesh_connected_component_test()
    {
        auto surface_mesh = SurfaceMesh3D::create_mesh();
        auto surface_mesh_builder =
            SurfaceMeshBuilder3D::create_builder( *surface_mesh );

        index_t nb_connected_components{ NO_ID };
        std::vector< index_t > connected_components;
        std::tie( nb_connected_components, connected_components ) =
            surface_mesh->connected_components();
        if( nb_connected_components != 0 || connected_components.size() != 0 )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Surface mesh should have 0 connected component." );
        }

        nb_connected_components = NO_ID;
        auto v0 = surface_mesh_builder->create_vertex( { 0., 0., 0. } );
        auto v1 = surface_mesh_builder->create_vertex( { 1., 0., 0. } );
        auto v2 = surface_mesh_builder->create_vertex( { 0., 1., 0. } );
        surface_mesh_builder->create_polygon( { v0, v1, v2 } );
        surface_mesh_builder->connect_polygons();
        std::tie( nb_connected_components, connected_components ) =
            surface_mesh->connected_components();
        std::vector< index_t > solution( 1, 0 );
        if( nb_connected_components != 1 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Surface mesh should have 1 connected "
                "component with an index at 0." );
        }

        nb_connected_components = NO_ID;
        auto v3 = surface_mesh_builder->create_vertex( { 1., 1., 1. } );
        surface_mesh_builder->create_polygon( { v1, v3, v2 } );
        auto v4 = surface_mesh_builder->create_vertex( { 5., 5., 0. } );
        auto v5 = surface_mesh_builder->create_vertex( { 5., 6., 0. } );
        auto v6 = surface_mesh_builder->create_vertex( { 6., 5., 0. } );
        surface_mesh_builder->create_polygon( { v4, v5, v6 } );
        surface_mesh_builder->connect_polygons();
        std::tie( nb_connected_components, connected_components ) =
            surface_mesh->connected_components();
        solution.push_back( 0 );
        solution.push_back( 1 );
        if( nb_connected_components != 2 || connected_components != solution )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Surface mesh should have 2 connected "
                "component with an index at 0 0 1." );
        }
    }

    void volume_mesh_connected_component_test()
    {
        auto volume_mesh = VolumeMesh3D::create_mesh();
        auto volume_mesh_builder =
            VolumeMeshBuilder3D::create_builder( *volume_mesh );

        index_t nb_connected_components{ NO_ID };
        std::vector< index_t > connected_components;
        std::tie( nb_connected_components, connected_components ) =
            volume_mesh->connected_components();
        if( nb_connected_components != 0 || connected_components.size() != 0 )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Volume mesh should have 0 connected component." );
        }
    }

    void run_tests()
    {
        point_set_mesh_connected_component_test();
        line_mesh_connected_component_test();
        surface_mesh_connected_component_test();
        volume_mesh_connected_component_test();
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test connected components" );
        run_tests();
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
