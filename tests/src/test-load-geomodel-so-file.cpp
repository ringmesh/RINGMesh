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

#include <future>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel building from a mesh loaded from a .so file
 * @author Pierre Anquez
 */

void test_mesh( const std::string& file_name )
{
    using namespace RINGMesh;

    GeoModel3D model;
    bool loaded_model_is_valid { geomodel_load( model, file_name ) };

    if( !loaded_model_is_valid ) {
        throw RINGMeshException( "RINGMesh Test", "Failed when loading model ",
            model.name(), ": the loaded model is not valid." );
    }

    // Check number of entities of the imported GeoModel (from TSolid or LightTSolid file)
    if( model.nb_corners() != 52 || model.nb_lines() != 98
        || model.nb_surfaces() != 55 || model.nb_regions() != 8
        || model.nb_geological_entities( Interface3D::type_name_static() ) != 11
        || model.nb_geological_entities( Contact3D::type_name_static() ) != 38
        || model.mesh.vertices.nb() != 6691 || model.mesh.polygons.nb() != 10049
        || model.mesh.cells.nb() != 34540 ) {
        throw RINGMeshException( "RINGMesh Test", "Failed when loading model ",
            model.name(), ": wrong number of entities." );
    }
}

int main()
{
    using namespace RINGMesh;

    try {
        default_configure();

        Logger::out( "TEST",
            "Import two meshed GeoModels (1 TSolid & 1 LightTSolid) from .so" );

        std::vector< std::future< void > > futures;

        futures.emplace_back(
            std::async( std::launch::async, &test_mesh,
                ringmesh_test_data_path + "modelA4.so" ) );

        futures.emplace_back(
            std::async( std::launch::async, &test_mesh,
                ringmesh_test_data_path + "modelA4_lts.so" ) );

        for( auto& future : futures ) {
            future.get();
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
