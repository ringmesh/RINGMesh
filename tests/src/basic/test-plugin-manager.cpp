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

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>

#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/plugin_manager.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

void load_plugin_from_command_line()
{
    GEO::CmdLine::set_arg( "sys:plugins", "RINGMesh_geomodel_tools;Bar" );
    if( PluginManager::load_plugins() )
    {
        throw RINGMeshException( "TEST", "Not supposed to be able load Bar" );
    }
    GEO::CmdLine::set_arg( "sys:plugins", "" );
}

void load_plugin_from_file()
{
    std::ofstream config( PluginManager::configuration_file );
    config << "RINGMesh_io";
    config.close();
    if( !PluginManager::load_plugins() )
    {
        throw RINGMeshException( "TEST", "Failed to load RINGMesh_io" );
    }
    GEO::FileSystem::delete_file( PluginManager::configuration_file );
}

void load_plugin_from_code()
{
    auto status = PluginManager::load_plugin( "RINGMesh_geomodel_builder" );
    if( !status )
    {
        throw RINGMeshException(
            "TEST", "Failed to load RINGMesh_geomodel_builder" );
    }

    status = PluginManager::load_plugin( "RINGMesh_geomodel_core" );
    if( !status )
    {
        throw RINGMeshException(
            "TEST", "Failed to load RINGMesh_geomodel_core" );
    }

    status = PluginManager::load_plugin( "RINGMesh_geomodel_core" );
    if( status )
    {
        throw RINGMeshException(
            "TEST", "Not supposed to load RINGMesh_geomodel_core twice" );
    }

    status = PluginManager::load_plugin( "Foo" );
    if( status )
    {
        throw RINGMeshException( "TEST", "Not supposed to be able load Foo" );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test Plugin Manager" );

        load_plugin_from_code();
        load_plugin_from_command_line();
        load_plugin_from_file();
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
