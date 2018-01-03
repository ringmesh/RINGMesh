/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/ringmesh_tests_config.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/core/geomodel.h>

#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel saving
 * @author Arnaud Botella
 */

using namespace RINGMesh;

const std::string ringmesh_test_save_path = ringmesh_test_path + "io/data/save/";

const index_t MAX_DIFF_DIGIT{ 1 };

void get_line( GEO::LineInput& in )
{
    in.get_line();
    in.get_fields();
}

bool is_double( const std::string& field )
{
    try
    {
        std::stod( field );
    }
    catch( const std::invalid_argument& )
    {
        return false;
    }
    return true;
}

void resize( std::string& word1, std::string& word2 )
{
    if( word1.size() < word2.size() )
    {
        word1.resize( word2.size(), '0' );
    }
    else if( word2.size() < word1.size() )
    {
        word2.resize( word1.size(), '0' );
    }
}

index_t remove_point( std::string& word )
{
    auto position = find( word, '.' );
    if( position != NO_ID )
    {
        word.erase( word.begin() + position );
    }
    return position;
}

bool compare_double( std::string word1, std::string word2 )
{
    if( word1.size() != word2.size() )
    {
        resize( word1, word2 );
    }
    auto position1 = remove_point( word1 );
    auto position2 = remove_point( word2 );
    if( position1 != position2 )
    {
        return false;
    }
    try
    {
        auto value1 = std::stoll( word1 );
        auto value2 = std::stoll( word2 );
        return std::abs( value1 - value2 ) <= MAX_DIFF_DIGIT;
    }
    catch( const std::invalid_argument& )
    {
        return false;
    }
}

bool compare_words( const std::string& word1, const std::string& word2 )
{
    if( is_double( word1 ) )
    {
        if( !compare_double( word1, word2 ) )
        {
            return false;
        }
    }
    else
    {
        if( std::strcmp( word1.c_str(), word2.c_str() ) != 0 )
        {
            return false;
        }
    }
    return true;
}

bool compare_output_files( const std::string& file1, const std::string& file2 )
{
    GEO::LineInput in1{ file1 };
    GEO::LineInput in2{ file2 };

    while( !in1.eof() && !in2.eof() && in1.get_line() && in2.get_line() )
    {
        in1.get_fields();
        in2.get_fields();
        if( in1.nb_fields() != in2.nb_fields() )
        {
            DEBUG( in1.nb_fields() );
            DEBUG( in2.nb_fields() );
            return false;
        }
        for( auto f : range( in1.nb_fields() ) )
        {
            if( !compare_words( in1.field( f ), in2.field( f ) ) )
            {
                return false;
            }
        }
    }

    return true;
}

void check_files( const std::string& file1, const std::string& file2 )
{
    if( !compare_output_files( file1, file2 ) )
    {
        throw RINGMeshException( "TEST",
            "Output file " + file1 + " does not match template file " + file2 );
    }
}

void check_output( GEO::LineInput& in )
{
    while( !in.eof() && in.get_line() )
    {
        in.get_fields();
        std::string template_output{ ringmesh_test_save_path + in.field( 0 ) };
        std::string new_output{ ringmesh_test_output_path + in.field( 0 ) };
        check_files( new_output, template_output );
    }
}

template < index_t DIMENSION >
void io_geomodel(
    const std::string& geomodel_file, const std::string& extension )
{
    GeoModel< DIMENSION > geomodel;
    geomodel_load( geomodel, geomodel_file );
    geomodel_save( geomodel, ringmesh_test_output_path + "geomodel"
                                 + std::to_string( DIMENSION ) + "d."
                                 + extension );
}

template < index_t DIMENSION >
void process_extension( const std::string& extension )
{
    if( extension == "gm" )
    {
        return;
    }
#if defined( WIN32 ) || defined( __APPLE__ )
    if( extension == "adeli" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // ADELI file format. Going to be fixed...
        return;
    }
    if( extension == "csmp" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // csmp file format. Going to be fixed...
        return;
    }
    if( extension == "fem" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // fem file format. Going to be fixed...
        return;
    }
    if( extension == "inp" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // inp file format. Going to be fixed...
        return;
    }
    if( extension == "mail" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // mail file format. Going to be fixed...
        return;
    }
    if( extension == "mfem" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // mfem file format. Going to be fixed...
        return;
    }
    if( extension == "so" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // so file format. Going to be fixed...
        return;
    }
    if( extension == "stl" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // stl file format. Going to be fixed...
        return;
    }
    if( extension == "tetgen" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // tetgen file format. Going to be fixed...
        return;
    }
    if( extension == "vtk" )
    {
        // @todo Temporary switch off the test for saving GeoModel into
        // vtk file format. Going to be fixed...
        return;
    }
#endif
    std::string info{ ringmesh_test_save_path + extension
                      + std::to_string( DIMENSION ) + "d.txt" };
    GEO::LineInput in{ info };
    if( !in.OK() )
    {
        throw RINGMeshException( "TEST", "Failed to load file: ", info );
    }
    get_line( in );
    io_geomodel< DIMENSION >(
        ringmesh_test_data_path + in.field( 0 ), extension );
    check_output( in );
    Logger::out( "TEST", "Format ", extension, " OK" );
}

template < index_t DIMENSION >
void test_output_geomodel()
{
    Logger::out( "TEST", "Save GeoModel", DIMENSION, "D files" );
    auto extensions =
        GeoModelOutputHandlerFactory< DIMENSION >::list_creators();
    for( const auto& extension : extensions )
    {
        process_extension< DIMENSION >( extension );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        GEO::CmdLine::set_arg( "validity:do_not_check", "A" );
        test_output_geomodel< 2 >();
        test_output_geomodel< 3 >();
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
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
