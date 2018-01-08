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
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/core/entity_type_manager.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/io/io.h>

/*!
 * @author Arnaud Botella
 */

namespace
{
    using namespace RINGMesh;

    void import_arg_group_stats()
    {
        GEO::CmdLine::declare_arg_group( "stats", "Statistics options" );
        GEO::CmdLine::declare_arg(
            "stats:volume", false, "Print statistics on the volume" );
        GEO::CmdLine::declare_arg(
            "stats:nb", true, "Print statistics on the number of entities" );
    }

    void print_geomodel2d_stats( const std::string& model_name )
    {
        GeoModel< 2 > geomodel;
        geomodel_load( geomodel, model_name );
        if( GEO::CmdLine::get_arg_bool( "stats:nb" ) )
        {
            print_geomodel_mesh_stats( geomodel );
        }
    }

    void print_geomodel3d_stats( const std::string& model_name )
    {
        GeoModel< 3 > geomodel;
        geomodel_load( geomodel, model_name );
        if( GEO::CmdLine::get_arg_bool( "stats:nb" ) )
        {
            print_geomodel_mesh_stats( geomodel );
        }
        if( GEO::CmdLine::get_arg_bool( "stats:volume" ) )
        {
            print_geomodel_mesh_cell_volumes( geomodel );
        }
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        print_header_information();
        Logger::div( "RINGMesh-Stats" );
        Logger::out( "", "Welcome to RINGMesh-Stats !" );

        CmdLine::import_arg_group( "in" );
        import_arg_group_stats();

        if( argc == 1 )
        {
            GEO::CmdLine::show_usage();
            return 0;
        }

        std::vector< std::string > filenames;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) )
        {
            return 1;
        }

        GEO::Stopwatch total( "Total time" );

        std::string model_name = GEO::CmdLine::get_arg( "in:geomodel" );
        if( model_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }

        index_t dimension = find_geomodel_dimension( model_name );
        if( dimension == 2 )
        {
            print_geomodel2d_stats( model_name );
        }
        else if( dimension == 3 )
        {
            print_geomodel3d_stats( model_name );
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
