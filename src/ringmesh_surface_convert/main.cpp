/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *     g.laurent.research@gmail.com
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


int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::div( "RINGMeshConvert" ) ;
    GEO::Logger::out( "" ) << "Welcome to RINGMeshSurfaceConvert !" << std::endl ;
    GEO::Logger::out( "" ) << "People working on this project in RING" << std::endl ;
    GEO::Logger::out( "" ) << "Gautier Laurent<g.laurent.research@gmail.com> "
        << std::endl ;
	
	// help
    if( argc == 1 ) {
		GEO::Logger::div( "Help" ) ;
        GEO::Logger::out( "" ) << "usage: " << argv[0] << " [out_format]" << std::endl ;
        GEO::Logger::out( "" ) << "out_format: a non empty list of output format amongst: obj mesh meshb ply off stl " << std::endl ;
        GEO::Logger::out( "" ) << "This will create a directory for each selected output format in the directory one level above the current one,"
			<< " and create a new file with the corresponding format for each .ts in the current directory." << std::endl ;
        return 0 ;
    }

	// parsing the output formats
    std::vector< std::string > output_formats ;
	GEO::CmdLine::parse( argc, argv, output_formats );
    if( output_formats.empty() ) {
        GEO::Logger::err( "I/O" ) << "Give at least one output format amongst: obj mesh meshb ply off stl"
            << std::endl ;
        return 1 ;
    }

    GEO::Stopwatch total( "Total time" ) ;

	// get current directory .ts files
	std::string starting_directory = GEO::FileSystem::get_current_working_directory();
	std::vector<std::string> input_file_names;
	GEO::FileSystem::get_files( starting_directory , input_file_names );
	
	std::vector<std::string> input_ts_names;
	for( std::vector<std::string>::iterator file_itr = input_file_names.begin(); file_itr < input_file_names.end(); ++file_itr ){
		if( GEO::FileSystem::extension( *file_itr ) == "ts" ){
			input_ts_names.push_back( *file_itr );
		}
	}

    if( input_ts_names.empty() ) {
        GEO::Logger::err( "I/O" ) << "Run this command in a folder with at least on .ts file."
            << std::endl ;
        return 1 ;
    }
	
	// create the output format folders
	if( !GEO::FileSystem::set_current_working_directory("..") ){
        GEO::Logger::err( "I/O" ) << "Can't access parent directory."
            << std::endl ;
        return 1 ;
	}
	for( std::vector<std::string>::iterator format_itr = output_formats.begin(); format_itr < output_formats.end(); ++format_itr ){
		if( !GEO::FileSystem::create_directory( *format_itr ) ){
			GEO::Logger::err( "I/O" ) << "Can't create " << *format_itr << " directory."
				<< std::endl ;
			return 1 ;
		}
	}

	// for each .ts file save it in each appropriate format
	for( std::vector<std::string>::iterator ts_itr = input_ts_names.begin() ; ts_itr < input_ts_names.end(); ++ts_itr ){
							
        GEO::Logger::out( "" ) << "Processing: " << (*ts_itr) << std::endl ;
		GEO::FileSystem::set_current_working_directory(starting_directory) ;

		// load the tsurf
		GEO::Mesh mesh_surface_in ;
		if( !load_ts_file( mesh_surface_in, *ts_itr ) ){
			GEO::Logger::err( "I/O" ) << "Can't load: " << *ts_itr << std::endl ;
			continue ;
		}
		// get the basename
		std::string surface_in_basename = GEO::FileSystem::base_name( *ts_itr );
        GEO::Logger::out( "" ) << " ... basename: " << surface_in_basename << std::endl ;

		// for each format save it
		for( std::vector<std::string>::iterator format_itr = output_formats.begin(); format_itr < output_formats.end(); ++format_itr ){
			
            GEO::Logger::out( "" ) << " ... going to: " << "..\\" << *format_itr << std::endl ;
		    GEO::FileSystem::set_current_working_directory("..") ;
		    GEO::FileSystem::set_current_working_directory(*format_itr) ;

            GEO::Logger::out( "" ) << " ... saving to: " << surface_in_basename + "." + (*format_itr) << std::endl ;
			if( !GEO::mesh_save( mesh_surface_in, surface_in_basename + "." + (*format_itr) ) ){
				GEO::Logger::err( "I/O" ) << "Can't save to: " << *ts_itr << std::endl ;
				continue ;
			}
		}

	}

    return 0 ;
}
