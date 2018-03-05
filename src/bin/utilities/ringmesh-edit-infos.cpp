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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/io/io.h>

/*!
 * @author Pierre Anquez
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    void edit_geomodel_name( const std::string& geomodel_in_path,
        const std::string& geomodel_new_name,
        const std::string& geomodel_out_path )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, geomodel_in_path );
        GeoModelBuilder< DIMENSION > builder( geomodel );
        builder.info.set_geomodel_name( geomodel_new_name );
        if( geomodel_out_path == geomodel_in_path )
        {
            GEO::FileSystem::delete_file( geomodel_in_path );
        }
        geomodel_save( geomodel, geomodel_out_path );
    }

    template < index_t DIMENSION >
    void edit_surface_name( const std::string& geomodel_in_path,
        const std::string& surface_old_name,
        const std::string& surface_new_name,
        const std::string& geomodel_out_path )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, geomodel_in_path );
        GeoModelBuilder< DIMENSION > builder( geomodel );
        for( const auto& surface : geomodel.surfaces() )
        {
            if( surface.name() == surface_old_name )
            {
                builder.info.set_mesh_entity_name(
                    surface.gmme(), surface_new_name );
            }
        }
        if( geomodel_out_path == geomodel_in_path )
        {
            GEO::FileSystem::delete_file( geomodel_in_path );
        }
        geomodel_save( geomodel, geomodel_out_path );
    }

    void show_usage_example()
    {
        Logger::div( "Example" );
        Logger::out( "",
            "ringmesh-edit-infos in:geomodel=path/to/input/geomodel.ext ",
            "edit:name=new_geomodel_name ",
            "[out:geomodel=path/to/input/new_geomodel_file.ext]" );
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        print_header_information();
        Logger::div( "RINGMesh" );
        Logger::out( "", "Welcome to RINGMesh-edit-infos !" );

        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        GEO::CmdLine::declare_arg_group( "edit", "Edit GeoModel infos" );
        GEO::CmdLine::declare_arg(
            "edit:geomodel_name", "", "New GeoModel name" );
        GEO::CmdLine::declare_arg(
            "edit:surface_old_name", "", "Old Surface name" );
        GEO::CmdLine::declare_arg(
            "edit:surface_new_name", "", "New Surface name" );
        if( argc == 1 )
        {
            GEO::CmdLine::show_usage();
            show_usage_example();
            return 0;
        }

        std::vector< std::string > filenames;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) )
        {
            show_usage_example();
            return 1;
        }

        std::string geomodel_in_file = GEO::CmdLine::get_arg( "in:geomodel" );
        if( geomodel_in_file.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }

        std::string geomodel_new_name =
            GEO::CmdLine::get_arg( "edit:geomodel_name" );
        std::string surface_old_name =
            GEO::CmdLine::get_arg( "edit:surface_old_name" );
        std::string surface_new_name =
            GEO::CmdLine::get_arg( "edit:surface_new_name" );
        if( geomodel_new_name.empty() && surface_old_name.empty() )
        {
            throw RINGMeshException( "I/O",
                "Give at least a new GeoModel name in edit:geomodel_name or a "
                "Old Surface Name in edit:surface_old_name" );
        }

        std::string geomodel_out_file = GEO::CmdLine::get_arg( "out:geomodel" );
        if( geomodel_out_file.empty() )
        {
            geomodel_out_file = geomodel_in_file;
        }

        index_t dimension = find_geomodel_dimension( geomodel_in_file );
        if( dimension == 2 )
        {
            if( !geomodel_new_name.empty() )
            {
                edit_geomodel_name< 2 >(
                    geomodel_in_file, geomodel_new_name, geomodel_out_file );
            }
            if( !surface_old_name.empty() )
            {
                edit_surface_name< 2 >( geomodel_in_file, surface_old_name,
                    surface_new_name, geomodel_out_file );
            }
        }
        else if( dimension == 3 )
        {
            if( !geomodel_new_name.empty() )
            {
                edit_geomodel_name< 3 >(
                    geomodel_in_file, geomodel_new_name, geomodel_out_file );
            }
            if( !surface_old_name.empty() )
            {
                edit_surface_name< 3 >( geomodel_in_file, surface_old_name,
                    surface_new_name, geomodel_out_file );
            }
        }
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
    return 0;
}
