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

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/logger.h>

#include <ringmesh/geogram_extension/geogram_mesh.h>

#include <ringmesh/geomodel/builder/geomodel_builder_gocad.h>
#include <ringmesh/geomodel/core/entity_type.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>

#include <ringmesh/io/io.h>

#include <ringmesh/mesh/mesh_builder.h>

#include <ringmesh/tetrahedralize/tetra_gen.h>

/*!
 * @file Test default configure
 * @author Pierre Anquez
 */

using namespace RINGMesh;

void throw_error_empty( const std::string& factory_name )
{
    throw RINGMeshException(
        "RINGMesh Test", "Factory ", factory_name, " has no register." );
}

void throw_error_key(
    const std::string& factory_name, const std::string& key_name )
{
    throw RINGMeshException( "RINGMesh Test", "Factory ", factory_name,
        " has no register for the key: ", key_name );
}

void test_mesh_initialize()
{
    if( GeoModelInputHandlerFactory2D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelInputHandler2D" );
    }
    if( GeoModelInputHandlerFactory3D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelInputHandler3D" );
    }
    if( GeoModelOutputHandlerFactory2D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelOutputHandler2D" );
    }
    if( GeoModelOutputHandlerFactory3D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelOutputHandler3D" );
    }
    if( WellGroupIOHandlerFactory::list_creators().empty() )
    {
        throw_error_empty( "WellGroupIOHandler" );
    }
}

void test_tetragen_initialize()
{
#ifdef RINGMESH_WITH_TETGEN
    if( !TetraGenFactory::has_creator( "TetGen" ) )
    {
        throw_error_key( "TetraGenFactory", "TetGen" );
    }
#endif

#ifdef USE_MG_TETRA
    if( !TetraGenFactory::has_creator( "MG_Tetra" ) )
    {
        throw_error_key( "TetraGenFactory", "MG_Tetra" );
    }
#endif
}

void test_geomodel_geological_entity_factories()
{
    if( GeoModelGeologicalEntityFactory2D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelGeologicalEntityFactory2D" );
    }
    if( GeoModelGeologicalEntityFactory3D::list_creators().empty() )
    {
        throw_error_empty( "GeoModelGeologicalEntityFactory3D" );
    }
}

void test_ringmesh_mesh_io_initialize()
{
    if( !GEO::MeshIOHandlerFactory::has_creator( "ts" ) )
    {
        throw_error_key( "GEO::MeshIOHandlerFactory", "ts" );
    }
    if( !GEO::MeshIOHandlerFactory::has_creator( "lin" ) )
    {
        throw_error_key( "GEO::MeshIOHandlerFactory", "lin" );
    }
}

void test_gocad_import_factories()
{
    if( GocadLineFactory::list_creators().empty() )
    {
        throw_error_empty( "GocadLineFactory" );
    }
    if( TSolidLineFactory::list_creators().empty() )
    {
        throw_error_empty( "TSolidLineFactory" );
    }
    if( MLLineFactory::list_creators().empty() )
    {
        throw_error_empty( "MLLineFactory" );
    }
}

void test_geogram_factory_2D()
{
    if( !PointSetMeshFactory2D::has_creator(
            GeogramPointSetMesh2D::type_name_static() ) )
    {
        throw_error_key( "PointSetMeshFactory2D",
            GeogramPointSetMesh2D::type_name_static() );
    }
    if( !LineMeshFactory2D::has_creator(
            GeogramLineMesh2D::type_name_static() ) )
    {
        throw_error_key(
            "LineMeshFactory2D", GeogramLineMesh2D::type_name_static() );
    }
    if( !SurfaceMeshFactory2D::has_creator(
            GeogramSurfaceMesh2D::type_name_static() ) )
    {
        throw_error_key(
            "SurfaceMeshFactory2D", GeogramSurfaceMesh2D::type_name_static() );
    }
    if( !PointSetMeshBuilderFactory2D::has_creator(
            GeogramPointSetMesh2D::type_name_static() ) )
    {
        throw_error_key( "PointSetMeshBuilderFactory2D",
            GeogramPointSetMesh2D::type_name_static() );
    }
    if( !LineMeshBuilderFactory2D::has_creator(
            GeogramLineMesh2D::type_name_static() ) )
    {
        throw_error_key(
            "LineMeshBuilderFactory2D", GeogramLineMesh2D::type_name_static() );
    }
    if( !SurfaceMeshBuilderFactory2D::has_creator(
            GeogramSurfaceMesh2D::type_name_static() ) )
    {
        throw_error_key( "SurfaceMeshBuilderFactory2D",
            GeogramSurfaceMesh2D::type_name_static() );
    }
}

void test_geogram_factory_3D()
{
    if( !PointSetMeshFactory3D::has_creator(
            GeogramPointSetMesh3D::type_name_static() ) )
    {
        throw_error_key( "PointSetMeshFactory3D",
            GeogramPointSetMesh3D::type_name_static() );
    }
    if( !LineMeshFactory3D::has_creator(
            GeogramLineMesh3D::type_name_static() ) )
    {
        throw_error_key(
            "LineMeshFactory3D", GeogramLineMesh3D::type_name_static() );
    }
    if( !SurfaceMeshFactory3D::has_creator(
            GeogramSurfaceMesh3D::type_name_static() ) )
    {
        throw_error_key(
            "SurfaceMeshFactory3D", GeogramSurfaceMesh3D::type_name_static() );
    }
    if( !VolumeMeshFactory3D::has_creator(
            GeogramVolumeMesh3D::type_name_static() ) )
    {
        throw_error_key(
            "VolumeMeshFactory3D", GeogramVolumeMesh3D::type_name_static() );
    }
    if( !PointSetMeshBuilderFactory3D::has_creator(
            GeogramPointSetMesh3D::type_name_static() ) )
    {
        throw_error_key( "PointSetMeshBuilderFactory3D",
            GeogramPointSetMesh3D::type_name_static() );
    }
    if( !LineMeshBuilderFactory3D::has_creator(
            GeogramLineMesh3D::type_name_static() ) )
    {
        throw_error_key(
            "LineMeshBuilderFactory3D", GeogramLineMesh3D::type_name_static() );
    }
    if( !SurfaceMeshBuilderFactory3D::has_creator(
            GeogramSurfaceMesh3D::type_name_static() ) )
    {
        throw_error_key( "SurfaceMeshBuilderFactory3D",
            GeogramSurfaceMesh3D::type_name_static() );
    }
    if( !VolumeMeshBuilderFactory3D::has_creator(
            GeogramVolumeMesh3D::type_name_static() ) )
    {
        throw_error_key( "VolumeMeshBuilderFactory3D",
            GeogramVolumeMesh3D::type_name_static() );
    }
}

int main()
{
    try
    {
        default_configure();

        Logger::out( "TEST", "Is RINGMesh correctly configured?" );

        // Test mesh initialize
        test_mesh_initialize();

        // Test TetraGen initialize
        test_tetragen_initialize();

        // Test Geological entity initialize
        test_geomodel_geological_entity_factories();

        // Test RINGMesh mesh IO initialize
        test_ringmesh_mesh_io_initialize();

        // Test Gocad import factory initialize
        test_gocad_import_factories();

        // Test geogram mesh register
        test_geogram_factory_2D();
        test_geogram_factory_3D();
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
