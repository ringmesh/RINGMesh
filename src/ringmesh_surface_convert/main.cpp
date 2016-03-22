/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/common.h>

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/io.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_io.h>

/*!
 * @file ringmesh_surface_convert/main.cpp
 * @brief Executable to High level functions on GeoModel
 * @author Gautier Laurent
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;


    try {
        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        // welcome
        GEO::Logger::div( "RINGMeshConvert" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshSurfaceConvert !"
            << std::endl ;
        GEO::Logger::out( "" ) << "People working on this project in RING"
            << std::endl ;
        GEO::Logger::out( "" ) << "Gautier Laurent<g.laurent.research@gmail.com> "
            << std::endl ;

        // help
        if( argc == 1 ) {
            GEO::Logger::div( "Help" ) ;
            GEO::Logger::out( "" ) << "usage: "
                << GEO::FileSystem::base_name( argv[0] ) << " [out_format]"
                << std::endl
                << "out_format: a non empty list of output format amongst:"
                << "obj mesh meshb ply off stl " << std::endl
                << "Should be launched in the directory that contains the .ts " 
                << "files to convert" << std::endl
                << "This will create a directory for each selected output format "
                << "in the directory one level above the current one,"
                << " and create a new file with the corresponding format for "
                << "each .ts in the current directory" << std::endl ;
            return 0 ;
        }

        // parsing the output formats
        std::vector< std::string > output_formats ;
        GEO::CmdLine::parse( argc, argv, output_formats ) ;
        if( output_formats.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least one output format amongst: obj mesh meshb ply off stl" ) ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        // get files in current directory
        std::string starting_directory =
            GEO::FileSystem::get_current_working_directory() ;
        std::vector< std::string > input_file_names ;
        GEO::FileSystem::get_files( starting_directory, input_file_names ) ;

        // filter only the files with .ts extension
        std::vector< std::string > input_ts_names ;
        for( index_t file_itr = 0; file_itr < input_file_names.size(); ++file_itr ) {
            if( GEO::FileSystem::extension( input_file_names[file_itr] ) == "ts" ) {
                input_ts_names.push_back( input_file_names[file_itr] ) ;
            }
        }

        if( input_ts_names.empty() ) {
            throw RINGMeshException( "I/O",
                "Run this command in a folder with at least one .ts file." ) ;
        }

        // create the output format folders
        if( !GEO::FileSystem::set_current_working_directory( ".." ) ) {
            throw RINGMeshException( "I/O", "Can't access parent directory." ) ;
        }
        
        std::string saving_directory =
            GEO::FileSystem::get_current_working_directory() ;
        for( index_t format_itr = 0; format_itr < output_formats.size();
            ++format_itr ) {
            if( !GEO::FileSystem::create_directory(
                saving_directory + "/" + output_formats[format_itr] ) ) {
                throw RINGMeshException( "I/O",
                    "Can't create " + output_formats[format_itr] + " directory." ) ;
            }
        }

        // for each .ts file save it in each appropriate format
        for( index_t ts_itr = 0; ts_itr < input_ts_names.size(); ++ts_itr ) {
            GEO::Logger::out( "" ) << "Processing: " << input_ts_names[ts_itr]
                << std::endl ;
            GEO::FileSystem::set_current_working_directory( starting_directory ) ;

            // load the tsurf
            GEO::Mesh mesh_surface_in ;
            if( !GEO::mesh_load( input_ts_names[ts_itr], mesh_surface_in ) ) {
                throw RINGMeshException( "I/O", "Can't load: " + input_ts_names[ts_itr] ) ;
            }

            // get the basename
            std::string surface_in_basename = GEO::FileSystem::base_name(
                input_ts_names[ts_itr] ) ;

            // for each format save it
            for( index_t format_itr = 0; format_itr < output_formats.size();
                ++format_itr ) {
                GEO::FileSystem::set_current_working_directory( ".." ) ;
                GEO::FileSystem::set_current_working_directory(
                    output_formats[format_itr] ) ;

                std::string surface_output_name = surface_in_basename + "."
                    + output_formats[format_itr] ;
                if( !GEO::mesh_save( mesh_surface_in, surface_output_name ) ) {
                    throw RINGMeshException( "I/O",
                        "Can't save to: " + surface_output_name ) ;
                }
            }

        }
    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
