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

#include <ringmesh/ringmesh_tests_config.h>

#include <future>

#include <geogram/basic/command_line.h>

#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*! Tests the GeoModel invalidity tracking.
 * Loads a .ml file, check its validity. Then, we alter the geomodel in order to
 * break its validity and we check that we detect the invalidities.
 * @returns 0 if success or an error code if not.
 * @author Pierre Anquez
 */

using namespace RINGMesh;

void make_geomodel_copy(
    const GeoModel3D& from, const std::string& name, GeoModel3D& to )
{
    GeoModelBuilder3D geomodel_breaker2( to );
    copy_geomodel( from, to );
    geomodel_breaker2.info.set_geomodel_name( name );
}

void verdict( const GeoModel3D& invalid_model,
    const std::string& feature,
    const ValidityCheckMode& validity_check_mode )
{
    if( is_geomodel_valid( invalid_model, validity_check_mode ) )
    {
        throw RINGMeshException( "RINGMesh Test", "Fail to ", feature );
    }
    else
    {
        Logger::out( "TEST", "Succeed to ", feature );
    }
}

int main()
{
    try
    {
        // Load a model without region : 6 surfaces defining a cube with holes
        // between surfaces (all surface borders are free borders)
        GeoModel3D not_sealed_cube_geomodel;
        std::string input_cube_model_file_name( ringmesh_test_data_path );
        input_cube_model_file_name += "not_sealed_cube.ml";
        bool is_cube_valid = geomodel_load(
            not_sealed_cube_geomodel, input_cube_model_file_name );

        if( is_cube_valid )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Failed to detect invalidities in loaded model ",
                not_sealed_cube_geomodel.name() );
        }

        // No validity checks at loading
        GEO::CmdLine::set_arg( "validity:do_not_check", "A" );

        std::vector< std::future< void > > futures;

        futures.emplace_back( std::async( std::launch::async, [] {
            std::string input_model_file_name{ ringmesh_test_data_path
                                               + "modelA6.ml" };
            GeoModel3D geomodel;
            geomodel_load( geomodel, input_model_file_name );

            Logger::out( "TEST", "Break geomodel:" );

            GeoModel3D invalid_model;
            make_geomodel_copy( geomodel, "broken model 1", invalid_model );
            GeoModelBuilder3D geomodel_breaker( invalid_model );
            geomodel_breaker.topology.create_mesh_entity(
                RINGMesh::Surface3D::type_name_static() );
            if( is_geomodel_valid(
                    invalid_model, ValidityCheckMode::TOPOLOGY ) )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Fail to detect addition of an isolated "
                    "GeoModelMeshEntity" );
            }
        } ) );

        futures.emplace_back( std::async( std::launch::async, [] {
            GeoModel3D cloudspin;
            geomodel_load(
                cloudspin, ringmesh_test_data_path + "CloudSpin.ml" );
            if( is_geomodel_valid( cloudspin,
                    ValidityCheckMode::SURFACE_LINE_MESH_CONFORMITY ) )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Fail to SURFACE_LINE_MESH_CONFORMITY on CloudSpin" );
            }
        } ) );

        for( auto& future : futures )
        {
            future.wait();
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
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
