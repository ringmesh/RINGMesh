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

    CmdLine::import_arg_group( "in" ) ;
    CmdLine::import_arg_group( "out" ) ;

    if( argc == 1 ) {
        GEO::CmdLine::show_usage() ;
        return 0 ;
    }

    std::vector< std::string > filenames ;
    if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
        return 1 ;
    }

    GEO::Stopwatch total( "Total time" ) ;

    std::string surface_in_name = GEO::CmdLine::get_arg( "in:model" ) ;
    if( surface_in_name  == "" ) {
        GEO::Logger::err( "I/O" ) << "Give at least a filename in in:model"
            << std::endl ;
        return 1 ;
    }
    GEO::Mesh mesh_surface_in ;
    if( !load_ts_file( mesh_surface_in, surface_in_name ) ){
        return 1 ;
	}
	
    std::string surface_out_name = GEO::CmdLine::get_arg( "out:model" ) ;
    if( surface_out_name != "" ) {
        if( !GEO::mesh_save( mesh_surface_in, surface_out_name ) )
            return 1 ;
    }
	
    return 0 ;
}
