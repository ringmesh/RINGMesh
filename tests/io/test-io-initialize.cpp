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
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/plugin_manager.h>

#include <ringmesh/io/io.h>

#include <ringmesh/io/geomodel_builder_gocad.h>

/*!
 * @file Test if geogram factories are loaded
 * @author Pierre Anquez
 */

using namespace RINGMesh;

void throw_error_empty( const std::string& factory_name )
{
    throw RINGMeshException(
        "RINGMesh Test", "Factory ", factory_name, " has no register." );
}

#define test_factory_not_empty( Factory )                                      \
    if( Factory::list_creators().empty() )                                     \
    {                                                                          \
        throw_error_empty( #Factory );                                         \
    }

void test_mesh_initialize()
{
    test_factory_not_empty( GeoModelInputHandlerFactory2D );
    test_factory_not_empty( GeoModelInputHandlerFactory3D );
    test_factory_not_empty( GeoModelOutputHandlerFactory2D );
    test_factory_not_empty( GeoModelOutputHandlerFactory3D );
    test_factory_not_empty( WellGroupIOHandlerFactory );
}

void test_gocad_import_factories()
{
    test_factory_not_empty( GocadLineFactory );
    test_factory_not_empty( TSolidLineFactory );
    test_factory_not_empty( MLLineFactory );
}

int main()
{
    try
    {
        std::string io = "RINGMesh_io";
#ifdef RINGMESH_DEBUG
        io += "d";
#endif
        PluginManager::load_plugin( io );

        Logger::out( "TEST", "Is io plugin well loaded?" );
        // Test mesh initialize
        test_mesh_initialize();
        // Test Gocad import factory initialize
        test_gocad_import_factories();
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
