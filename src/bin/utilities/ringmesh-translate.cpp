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
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/io/io.h>

/*!
 * @author Benjamin Chauvin
 */

namespace
{
    using namespace RINGMesh;

    void hello()
    {
        print_header_information();
        Logger::div( "RINGMesh-Translate" );
        Logger::out( "", "Welcome to RINGMesh-Translate !" );
    }

    void import_arg_group_translation()
    {
        GEO::CmdLine::declare_arg_group(
            "translation", "Options to translate a GeoModel" );
        GEO::CmdLine::declare_arg( "translation:vector", "0 0 0",
            "Translation vector to be written between quotation marks" );
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        import_arg_group_translation();
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > extract_coords_from_string(
        const std::string& coords_in_string )
    {
        std::vector< std::string > split_coords;
        split_coords.reserve( DIMENSION );
        GEO::String::split_string( coords_in_string, ' ', split_coords, true );
        if( split_coords.size() != DIMENSION )
        {
            throw RINGMeshException( "I/O", "Vector (", coords_in_string,
                ") has not exactly ", DIMENSION, " components" );
        }
        vecn< DIMENSION > coords_vec;
        for( auto split_coords_itr : range( DIMENSION ) )
        {
            coords_vec[split_coords_itr] =
                GEO::String::to_double( split_coords[split_coords_itr] );
        }
        return coords_vec;
    }

    template < index_t DIMENSION >
    void translate_geomodel( const std::string& input_geomodel_name )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, input_geomodel_name );

        std::string translation_vector_string =
            GEO::CmdLine::get_arg( "translation:vector" );
        vecn< DIMENSION > translation_vector =
            extract_coords_from_string< DIMENSION >(
                translation_vector_string );

        translate( geomodel, translation_vector );

        std::string output_geomodel_name =
            GEO::CmdLine::get_arg( "out:geomodel" );
        if( output_geomodel_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in out:geomodel" );
        }
        geomodel_save( geomodel, output_geomodel_name );
    }

    void run()
    {
        GEO::Stopwatch total( "Total time" );

        std::string input_geomodel_name =
            GEO::CmdLine::get_arg( "in:geomodel" );
        if( input_geomodel_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }

        index_t dimension = find_geomodel_dimension( input_geomodel_name );
        if( dimension == 2 )
        {
            translate_geomodel< 2 >( input_geomodel_name );
        }
        else if( dimension == 3 )
        {
            translate_geomodel< 3 >( input_geomodel_name );
        }
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        hello();
        import_arg_groups();

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
        run();
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
