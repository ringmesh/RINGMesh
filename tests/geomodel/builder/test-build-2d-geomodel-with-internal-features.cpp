/* * Copyright (c) 2012-2019, Association Scientifique pour la Geologie et ses
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
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel loading
 * @author Pierre Anquez
 */

using namespace RINGMesh;

const std::string ringmesh_test_load_path =
    ringmesh_test_path + "io/data/load/";

void throw_error( const std::string& feature )
{
    throw RINGMeshException( "RINGMesh Test", "Failed when loading model: ",
        "the loaded model as not the correct number of ", feature );
}

void check_geomodel( const GeoModel2D& geomodel )
{
    if( geomodel.nb_corners() != 38 )
    {
        throw_error( "corners" );
    }

    if( geomodel.nb_lines() != 36 )
    {
        throw_error( "lines" );
    }

    if( geomodel.nb_surfaces() != 6 )
    {
        throw_error( "surfaces" );
    }

    index_t max_boundaries{ 0 };
    for( const auto& surface : geomodel.surfaces() )
    {
        max_boundaries = std::max( max_boundaries, surface.nb_boundaries() );
    }

    if( max_boundaries != 34 )
    {
        throw_error( "number of boundaries" );
    }
}

void test_geomodel_2D_with_internal_features()
{
    Logger::out( "TEST", "Load GeoModel2D with internal features" );
    GeoModel2D geomodel;
    geomodel_load(
        geomodel, ringmesh_test_data_path + "model_with_internal_features.gm" );
    check_geomodel( geomodel );
}

int main()
{
    try
    {
        Logger::out( "TEST", "Import GeoModel files" );
        test_geomodel_2D_with_internal_features();
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
