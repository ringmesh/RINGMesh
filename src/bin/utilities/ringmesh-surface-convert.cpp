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

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_io.h>
#include <ringmesh/basic/common.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/io/io.h>

/*!
 * @file ringmesh_surface_convert/main.cpp
 * @brief Executable to High level functions on GeoModel
 * @author Gautier Laurent
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        // welcome
        print_header_information();
        Logger::div( "RINGMesh-Surface-Convert" );
        Logger::out( "", "Welcome to RINGMesh-Surface-Convert !" );

        // help
        if( argc == 1 )
        {
            Logger::div( "Help" );
            Logger::out( "", "usage: ", GEO::FileSystem::base_name( argv[0] ),
                " [out_format]" );
            Logger::out( "",
                "out_format: a non empty list of output format amongst:",
                "obj mesh meshb ply off stl " );
            Logger::out( "",
                "Should be launched in the directory that contains the .ts ",
                "files to convert" );
            Logger::out( "",
                "This will create a directory for each selected output format ",
                "in the directory one level above the current one,",
                " and create a new file with the corresponding format for ",
                "each .ts in the current directory" );
            return 0;
        }

        // parsing the output formats
        std::vector< std::string > output_formats;
        GEO::CmdLine::parse( argc, argv, output_formats );
        if( output_formats.empty() )
        {
            throw RINGMeshException( "I/O", "Give at least one output format "
                                            "amongst: obj mesh meshb ply off "
                                            "stl" );
        }

        GEO::Stopwatch total( "Total time" );

        // get files in current directory
        std::string starting_directory =
            GEO::FileSystem::get_current_working_directory();
        std::vector< std::string > input_file_names;
        GEO::FileSystem::get_files( starting_directory, input_file_names );

        // filter only the files with .ts extension
        std::vector< std::string > input_ts_names;
        for( const std::string& file : input_file_names )
        {
            if( GEO::FileSystem::extension( file ) == "ts" )
            {
                input_ts_names.push_back( file );
            }
        }

        if( input_ts_names.empty() )
        {
            throw RINGMeshException( "I/O",
                "Run this command in a folder with at least one .ts file." );
        }

        // create the output format folders
        if( !GEO::FileSystem::set_current_working_directory( ".." ) )
        {
            throw RINGMeshException( "I/O", "Can't access parent directory." );
        }

        std::string saving_directory =
            GEO::FileSystem::get_current_working_directory();
        for( const std::string& format : output_formats )
        {
            if( !GEO::FileSystem::create_directory(
                    saving_directory + "/" + format ) )
            {
                throw RINGMeshException(
                    "I/O", "Can't create ", format, " directory." );
            }
        }

        // for each .ts file save it in each appropriate format
        for( const std::string& ts_name : input_ts_names )
        {
            Logger::out( "", "Processing: ", ts_name );
            GEO::FileSystem::set_current_working_directory(
                starting_directory );

            // load the tsurf
            GEO::Mesh mesh_surface_in;
            if( !GEO::mesh_load( ts_name, mesh_surface_in ) )
            {
                throw RINGMeshException( "I/O", "Can't load: ", ts_name );
            }

            // get the basename
            std::string surface_in_basename =
                GEO::FileSystem::base_name( ts_name );

            // for each format save it
            for( const std::string& output_format : output_formats )
            {
                GEO::FileSystem::set_current_working_directory( ".." );
                GEO::FileSystem::set_current_working_directory( output_format );

                std::string surface_output_name =
                    surface_in_basename + "." + output_format;
                if( !GEO::mesh_save( mesh_surface_in, surface_output_name ) )
                {
                    throw RINGMeshException(
                        "I/O", "Can't save to: ", surface_output_name );
                }
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
