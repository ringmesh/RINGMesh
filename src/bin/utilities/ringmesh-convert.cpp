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
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @author Arnaud Botella
 */

namespace
{
    using namespace RINGMesh;

    void convert_mesh( const std::string& mesh_in_name )
    {
        GEO::Mesh mesh;
        GEO::mesh_load( mesh_in_name, mesh );
        std::string mesh_out_name = GEO::CmdLine::get_arg( "out:mesh" );
        if( mesh_out_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give the parameter out:mesh to save the mesh" );
        }
        GEO::mesh_save( mesh, mesh_out_name );
    }

    template < index_t DIMENSION >
    void convert_geomodel( const std::string& geomodel_in_name )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, geomodel_in_name );
        std::string geomodel_out_name = GEO::CmdLine::get_arg( "out:geomodel" );
        if( geomodel_out_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give the parameter out:geomodel to save the geomodel" );
        }
        geomodel_save( geomodel, geomodel_out_name );
    }

    void show_usage_example()
    {
        Logger::div( "Example" );
        Logger::out( "",
            "ringmesh-convert in:geomodel=path/to/input/geomodel.ext ",
            "out:geomodel=path/to/output/geomodel.ext" );
    }
}

namespace RINGMesh
{
    namespace CmdLine
    {
        void import_more_in_out()
        {
            GEO::CmdLine::declare_arg(
                "in:mesh", "", "Filename of the input mesh" );
            GEO::CmdLine::declare_arg( "out:mesh", "", "Saves the mesh" );
        }
    }
}
int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        print_header_information();
        Logger::div( "RINGMesh-Convert" );
        Logger::out( "", "Welcome to RINGMesh-Convert !" );

        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        CmdLine::import_more_in_out();
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

        GEO::Stopwatch total( "Total time" );

        std::string geomodel_in_name = GEO::CmdLine::get_arg( "in:geomodel" );
        std::string mesh_in_name = GEO::CmdLine::get_arg( "in:mesh" );
        if( geomodel_in_name.empty() && mesh_in_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel or in:mesh" );
        }

        if( geomodel_in_name.empty() )
        {
            convert_mesh( mesh_in_name );
        }
        else
        {
            index_t dimension = find_geomodel_dimension( geomodel_in_name );
            if( dimension == 2 )
            {
                convert_geomodel< 2 >( geomodel_in_name );
            }
            else if( dimension == 3 )
            {
                convert_geomodel< 3 >( geomodel_in_name );
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
