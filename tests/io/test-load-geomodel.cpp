/* * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <geogram/basic/line_stream.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel loading
 * @author Arnaud Botella
 */

using namespace RINGMesh;

const std::string ringmesh_test_load_path = ringmesh_test_path + "io/data/load/";

template < index_t DIMENSION >
void load_input_geomodel(
    GeoModel< DIMENSION >& geomodel, const std::string& file )
{
    auto loaded_model_is_valid =
        geomodel_load( geomodel, ringmesh_test_data_path + file );
    if( !loaded_model_is_valid )
    {
        throw RINGMeshException(
            "RINGMesh Test", "Failed when loading model " + geomodel.name()
                                 + ": the loaded model is not valid." );
    }
}

std::map< std::string, index_t > load_reference_info( const std::string& file )
{
    std::map< std::string, index_t > reference;
    std::string info{ ringmesh_test_load_path + file };
    GEO::LineInput in( info );
    if( !in.OK() )
    {
        throw RINGMeshException( "TEST", "Failed to load file: ", info );
    }
    while( !in.eof() && in.get_line() )
    {
        in.get_fields();
        reference[in.field( 0 )] = in.field_as_uint( 1 );
    }
    return reference;
}

template < index_t DIMENSION >
void throw_error(
    const GeoModel< DIMENSION >& geomodel, const std::string& entity )
{
    throw RINGMeshException( "RINGMesh Test", "Failed when loading model ",
        geomodel.name(), ": the loaded model as not the correct number of ",
        entity );
}

template < index_t DIMENSION >
void check_geomodel_base( const GeoModel< DIMENSION >& geomodel,
    const std::map< std::string, index_t >& reference )
{
    const auto& manager = geomodel.entity_type_manager();
    for( const auto& type : manager.mesh_entity_manager.mesh_entity_types() )
    {
        if( reference.at( type.string() ) != geomodel.nb_mesh_entities( type ) )
        {
            throw_error( geomodel, type.string() );
        }
    }

    for( const auto& type :
        manager.geological_entity_manager.geological_entity_types() )
    {
        if( reference.at( type.string() )
            != geomodel.nb_geological_entities( type ) )
        {
            throw_error( geomodel, type.string() );
        }
    }

    if( geomodel.mesh.vertices.nb() != reference.at( "Vertex" ) )
    {
        throw_error( geomodel, "vertices" );
    }
    if( geomodel.mesh.polygons.nb() != reference.at( "Polygon" ) )
    {
        throw_error( geomodel, "polygons" );
    }
}

template < index_t DIMENSION >
void check_geomodel( const GeoModel< DIMENSION >& geomodel,
    const std::map< std::string, index_t >& reference );

template <>
void check_geomodel( const GeoModel2D& geomodel,
    const std::map< std::string, index_t >& reference )
{
    check_geomodel_base( geomodel, reference );
}

template <>
void check_geomodel( const GeoModel3D& geomodel,
    const std::map< std::string, index_t >& reference )
{
    check_geomodel_base( geomodel, reference );
    if( geomodel.mesh.cells.nb() != reference.at( "Cell" ) )
    {
        throw_error( geomodel, "cells" );
    }
}

template < index_t DIMENSION >
void process_extension( const std::string& extension )
{
    std::string info{ ringmesh_test_load_path + extension
                      + std::to_string( DIMENSION ) + "d.txt" };
    GEO::LineInput in{ info };
    if( !in.OK() )
    {
        throw RINGMeshException( "TEST", "Failed to load file: ", info );
    }
    while( !in.eof() && in.get_line() )
    {
        in.get_fields();
        std::string file{ in.field( 0 ) };
        GeoModel< DIMENSION > geomodel;
        load_input_geomodel( geomodel, file );
        check_geomodel( geomodel, load_reference_info( file + ".txt" ) );
        Logger::out( "TEST", "Import GeoModel from ", file, " OK" );
    }
}

template < index_t DIMENSION >
void test_input_geomodels()
{
    Logger::out( "TEST", "Load GeoModel", DIMENSION, "D files" );
    auto extensions = GeoModelInputHandlerFactory< DIMENSION >::list_creators();
    for( const auto& extension : extensions )
    {
        process_extension< DIMENSION >( extension );
    }
}

int main()
{
    try
    {
        Logger::out( "TEST", "Import GeoModel files" );
        test_input_geomodels< 2 >();
        test_input_geomodels< 3 >();
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
