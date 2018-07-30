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

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/io/io.h>

/*! Tests the GeoModel copy API function.
 * Loads a .gm file, copy it, and save it.
 * @author Pierre Anquez
 */

using namespace RINGMesh;

int main()
{
    try {
        // 3D
        GeoModel3D input;
        std::string input_filename( ringmesh_test_data_path );
        input_filename += "modelA1_version2.gm";
        bool is_valid = geomodel_load( input, input_filename );

        if( !is_valid ) {
            throw RINGMeshException( "RINGMesh Test",
                "Input 3D GeoModel must be valid." );
        }

        //Copy it
        GeoModel3D copy;
        copy_geomodel( input, copy );

        if( !is_geomodel_valid( copy ) ) {
            throw RINGMeshException( "RINGMesh Test",
                "Copied 3D GeoModel is not valid." );
        }

        std::string output_filename( ringmesh_test_output_path );
        output_filename += "copied_geomodel3D.gm";
        geomodel_save( copy, output_filename );

        // 2D
        GeoModel2D input2d;
        std::string input2d_filename( ringmesh_test_data_path );
        input2d_filename += "model_2d_version2.gm";
        bool is_valid2d = geomodel_load( input2d, input2d_filename );

        if( !is_valid2d ) {
            throw RINGMeshException( "RINGMesh Test",
                "Input 2D GeoModel must be valid." );
        }

        //Copy it
        GeoModel2D copy2d;
        copy_geomodel( input2d, copy2d );

        if( !is_geomodel_valid( copy2d ) ) {
            throw RINGMeshException( "RINGMesh Test",
                "Copied 2D GeoModel is not valid." );
        }

        std::string output2d_filename( ringmesh_test_output_path );
        output2d_filename += "copied_geomodel2D.gm";
        geomodel_save( copy2d, output2d_filename );
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
