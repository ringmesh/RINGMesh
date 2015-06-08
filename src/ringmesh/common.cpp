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

#ifdef WIN32

  #pragma section(".CRT$XCU",read)
  #define INITIALIZER( f ) \
    static void __cdecl f( void ) ; \
    __declspec( allocate( ".CRT$XCU" ) ) void( __cdecl * f ## _ ) (void) = f ; \
    static void __cdecl f( void )

#elif defined( __GNUC__ )

  #define INITIALIZER( f ) \
    static void f( void ) __attribute__( ( constructor ) ) ; \
    static void f( void )

#endif

#include <ringmesh/attribute.h>
#include <ringmesh/io.h>
#include <ringmesh/tetra_gen.h>

#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#ifdef RINGMESH_WITH_GRAPHICS
#   include <geogram_gfx/basic/common.h>
#endif

static bool inited = false ;
INITIALIZER( initialize ) {
    if( !inited ) {
        inited = true ;
        GEO::initialize() ;
        GEO::CmdLine::import_arg_group( "sys" ) ;
        GEO::CmdLine::set_arg( "sys:assert", "abort" ) ;
        GEO::CmdLine::set_arg( "sys:FPE", false ) ;
        GEO::CmdLine::import_arg_group( "algo" ) ;
        GEO::CmdLine::set_arg( "algo:predicates", "exact" ) ;
        GEO::CmdLine::import_arg_group( "log" ) ;
        GEO::CmdLine::set_arg( "sys:use_doubles", true ) ;

#ifdef RINGMESH_WITH_GRAPHICS
        GEO::Graphics::initialize();
        GEO::CmdLine::import_arg_group( "gfx" ) ;
#endif

        RINGMesh::RINGMeshIO::initialize() ;
        RINGMesh::TetraGen::initialize() ;

        // Initialization for BoundaryModel attribute serialization
//        RINGMesh::AttributeSerializer::initialize() ;
//        RINGMesh::ringmesh_register_attribute_type< int           > ( "int" ) ;
//        RINGMesh::ringmesh_register_attribute_type< unsigned int  > ( "index" ) ;
//        RINGMesh::ringmesh_register_attribute_type< double        > ( "double" ) ;
//        RINGMesh::ringmesh_register_attribute_type< float         > ( "float" ) ;
//        RINGMesh::ringmesh_register_attribute_type< bool          > ( "bool" ) ;
//        atexit( RINGMesh::AttributeSerializer::terminate ) ;
    }
}
