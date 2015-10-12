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
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/common.h>

#include <ringmesh/command_line.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/io.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>


int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::div( "RINGMeshConvert" ) ;
    GEO::Logger::out( "" ) << "Welcome to RINGMeshConvert !" << std::endl ;
    GEO::Logger::out( "" ) << "People working on the project in RING" << std::endl ;
    GEO::Logger::out( "" ) << "Arnaud Botella <arnaud.botella@univ-lorraine.fr> "
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

    std::string model_in_name = GEO::CmdLine::get_arg( "in:model" ) ;
    if( model_in_name == "" ) {
        GEO::Logger::err( "I/O" ) << "Give at least a filename in in:model"
            << std::endl ;
        return 1 ;
    }
    BoundaryModel model_in ;
    if( !model_load( model_in_name, model_in ) )
        return 1 ;

    std::string mesh_in_name = GEO::CmdLine::get_arg( "in:mesh" ) ;
    if( mesh_in_name != "" ) {
        if( !mesh_load( mesh_in_name, model_in ) ) return 1 ;
    }

    std::string model_out_name = GEO::CmdLine::get_arg( "out:model" ) ;
    if( model_out_name != "" ) {
        if( !model_save( model_in, model_out_name ) )
            return 1 ;
    }

    std::string mesh_out_name = GEO::CmdLine::get_arg( "out:mesh" ) ;
    if( mesh_out_name != "" ) {
        if( !mesh_save( model_in, mesh_out_name ) )
            return 1 ;
    }

    return 0 ;
}
