/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/ringmesh_tests_config.h>

#include <geogram/basic/line_stream.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel building
 * @author Arnaud Botella
 */

using namespace RINGMesh;

void load_geomodel( GeoModel& geomodel, const std::string& file )
{
    bool loaded_model_is_valid = geomodel_load( geomodel, file );
    if( !loaded_model_is_valid ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when loading model " + geomodel.name()
                + ": the loaded model is not valid." );
    }
}

void get_line( GEO::LineInput& in )
{
    in.get_line();
    in.get_fields();
}

void throw_error( const GeoModel& geomodel, const std::string& entity )
{
    throw RINGMeshException( "RINGMesh Test",
        "Failed when loading model " + geomodel.name()
            + ": the loaded model as the correct number of " + entity );
}

void check_geomodel( const GeoModel& geomodel, const std::string& result )
{
    GEO::LineInput in( result );
    get_line( in );
    index_t nb_corners = in.field_as_uint( 1 );
    if( geomodel.nb_corners() != nb_corners ) {
        throw_error( geomodel, "corners" );
    }
    get_line( in );
    index_t nb_lines = in.field_as_uint( 1 );
    if( geomodel.nb_lines() != nb_lines ) {
        throw_error( geomodel, "lines" );
    }
    get_line( in );
    index_t nb_surfaces = in.field_as_uint( 1 );
    if( geomodel.nb_surfaces() != nb_surfaces ) {
        throw_error( geomodel, "surfaces" );
    }
    get_line( in );
    index_t nb_regions = in.field_as_uint( 1 );
    if( geomodel.nb_regions() != nb_regions ) {
        throw_error( geomodel, "regions" );
    }
    get_line( in );
    index_t nb_vertices = in.field_as_uint( 1 );
    if( geomodel.mesh.vertices.nb() != nb_vertices ) {
        throw_error( geomodel, "vertices" );
    }
    get_line( in );
    index_t nb_polygons = in.field_as_uint( 1 );
    if( geomodel.mesh.polygons.nb() != nb_polygons ) {
        throw_error( geomodel, "polygons" );
    }
    get_line( in );
    index_t nb_cells = in.field_as_uint( 1 );
    if( geomodel.mesh.cells.nb() != nb_cells ) {
        throw_error( geomodel, "cells" );
    }
    while( !in.eof() && in.get_line() ) {
        in.get_fields();
        std::string geol_type = in.field( 0 );
        index_t nb_geol_entities = in.field_as_uint( 1 );
        if( geomodel.nb_geological_entities( geol_type ) != nb_geol_entities ) {
            throw_error( geomodel, geol_type );
        }
    }
}

void process_file( const std::string& file )
{
    Logger::out( "TEST", "Import GeoModel from ", file );
    GeoModel geomodel;
    load_geomodel( geomodel, file );

    std::string result = file + ".txt";
    check_geomodel( geomodel, result );
    Logger::out( "TEST", "Import GeoModel from ", file, " OK" );
}

int main()
{
    using namespace RINGMesh;

    try {
        default_configure();

        Logger::out( "TEST", "Import GeoModel files" );
        std::vector< std::string > files;
        std::string data_path = ringmesh_test_data_path + "/load";
        GEO::FileSystem::get_files( data_path, files );

        for( const std::string& file : files ) {
            if( GEO::FileSystem::extension( file ) == "txt" ) {
                continue;
            }
            process_file( file );
        }

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
