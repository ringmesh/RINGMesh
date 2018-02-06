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

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/plugin_manager.h>

#include <ringmesh/geogram_extension/geogram_mesh.h>

#include <ringmesh/mesh/mesh_builder.h>

/*!
 * @file Test if geogram factories are loaded
 * @author Pierre Anquez
 */

using namespace RINGMesh;

void throw_error_key(
    const std::string& factory_name, const std::string& key_name )
{
    throw RINGMeshException( "RINGMesh Test", "Factory ", factory_name,
        " has no register for the key: ", key_name );
}

#define test_key_in_factory( Factory, Key )                                    \
    if( !Factory::has_creator( Key ) )                                         \
    {                                                                          \
        throw_error_key( #Factory, Key );                                      \
    }

void test_geogram_factory_2D()
{
    test_key_in_factory(
        PointSetMeshFactory2D, GeogramPointSetMesh2D::type_name_static() );
    test_key_in_factory(
        LineMeshFactory2D, GeogramLineMesh2D::type_name_static() );
    test_key_in_factory(
        SurfaceMeshFactory2D, GeogramSurfaceMesh2D::type_name_static() );

    test_key_in_factory( PointSetMeshBuilderFactory2D,
        GeogramPointSetMesh2D::type_name_static() );
    test_key_in_factory(
        LineMeshBuilderFactory2D, GeogramLineMesh2D::type_name_static() );
    test_key_in_factory(
        SurfaceMeshBuilderFactory2D, GeogramSurfaceMesh2D::type_name_static() );
}

void test_geogram_factory_3D()
{
    test_key_in_factory(
        PointSetMeshFactory3D, GeogramPointSetMesh3D::type_name_static() );
    test_key_in_factory(
        LineMeshFactory3D, GeogramLineMesh3D::type_name_static() );
    test_key_in_factory(
        SurfaceMeshFactory3D, GeogramSurfaceMesh3D::type_name_static() );
    test_key_in_factory(
        VolumeMeshFactory3D, GeogramVolumeMesh3D::type_name_static() );

    test_key_in_factory( PointSetMeshBuilderFactory3D,
        GeogramPointSetMesh3D::type_name_static() );
    test_key_in_factory(
        LineMeshBuilderFactory3D, GeogramLineMesh3D::type_name_static() );
    test_key_in_factory(
        SurfaceMeshBuilderFactory3D, GeogramSurfaceMesh3D::type_name_static() );
    test_key_in_factory(
        VolumeMeshBuilderFactory3D, GeogramVolumeMesh3D::type_name_static() );
}

void test_geograph_mesh_io_initialize()
{
    test_key_in_factory( GEO::MeshIOHandlerFactory, "ts" );
    test_key_in_factory( GEO::MeshIOHandlerFactory, "lin" );
}

int main()
{
    try
    {
        std::string geo_ext = "RINGMesh_geogram_extension";
#ifdef RINGMESH_DEBUG
        geo_ext += "d";
#endif
        PluginManager::load_plugin( geo_ext );

        Logger::out( "TEST", "Is geogram plugin well loaded?" );
        // Test geogram mesh register
        test_geogram_factory_2D();
        test_geogram_factory_3D();
        // Test geogram mesh IO initialize (.ts and .lin)
        test_geograph_mesh_io_initialize();
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
