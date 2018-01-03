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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/builder/geomodel_builder_2d_from_3d.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*! Tests the creation of a GeoModel2D from the projection of a
 * GeoModel3D in a plane.
 * @author Pierre Anquez
 */
int main()
{
    using namespace RINGMesh;

    try
    {
        std::string input_geomodel3d_file_name =
            ringmesh_test_data_path + "seg_overthrust_afault.gm";

        Logger::out( "TEST", "Loading GeoModel3D input file ",
            input_geomodel3d_file_name );

        GeoModel3D geomodel3d;
        geomodel_load( geomodel3d, input_geomodel3d_file_name );

        vec3 plane_normal{ 987., 0., 2150. };
        vec3 plane_origin{ 6300., 10500., -3200. };
        Geometry::Plane projection_plane( plane_normal, plane_origin );
        PlaneReferenceFrame3D plane_frame( projection_plane );
        GeoModel2D projection_geomodel2d( plane_frame );
        GeoModelBuilder2DProjection geomodel2d_builder(
            projection_geomodel2d, geomodel3d, projection_plane );
        geomodel2d_builder.build_geomodel();

        std::vector< std::future< void > > checks;
        checks.emplace_back(
            std::async( std::launch::async, [&projection_geomodel2d] {
                if( !is_geomodel_valid( projection_geomodel2d ) )
                {
                    throw RINGMeshException(
                        "TEST", "FAILED : built GeoModel2D is not valid" );
                }
            } ) );

        std::string output_model_file_name( ringmesh_test_output_path );
        output_model_file_name +=
            projection_geomodel2d.name() + "_saved_out.gm";
        geomodel_save( projection_geomodel2d, output_model_file_name );

        GeoModel2D reloaded_geomodel2d;
        auto is_reloaded_model_valid =
            geomodel_load( reloaded_geomodel2d, output_model_file_name );
        if( !is_reloaded_model_valid )
        {
            std::string output_model_file_name_bis( ringmesh_test_output_path );
            output_model_file_name_bis +=
                reloaded_geomodel2d.name() + "_saved_out_bis.gm";
            geomodel_save( reloaded_geomodel2d, output_model_file_name_bis );
            throw RINGMeshException(
                "TEST", "FAILED : reloaded GeoModel2D is not valid" );
        }

        for( auto& check : checks )
        {
            check.wait();
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
