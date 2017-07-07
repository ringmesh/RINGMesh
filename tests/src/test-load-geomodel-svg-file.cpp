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

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*! Tests the loading a svg file.
 * @author Arnaud Botella
 */
int main()
{
    using namespace RINGMesh;

    try {
        default_configure();
        std::string input_model_file_name = ringmesh_test_data_path
            + "sketch_conform.svg";

        Logger::out( "TEST", "Geomodel input test. Loading file ",
            input_model_file_name );

        GeoModel2D in;
        /*bool loaded_model_is_valid =*/geomodel_load( in, input_model_file_name );

//        if( !loaded_model_is_valid ) {
//            throw RINGMeshException( "RINGMesh Test",
//                "Failed when loading model " + in.name()
//                    + ": the loaded model is not valid." );
//        }

        std::string output_model_file_name( ringmesh_test_output_path );
        output_model_file_name += in.name() + "_saved_out.gm";
        geomodel_save( in, output_model_file_name );

        GeoModel2D in2;
        /*bool reloaded_model_is_valid =*/geomodel_load( in2,
            output_model_file_name );

//        if( !reloaded_model_is_valid ) {
//            throw RINGMeshException( "RINGMesh Test",
//                "Failed when reloading model " + in2.name()
//                    + ": the reloaded model is not valid." );
//        }

        std::string output_model_file_name_bis( ringmesh_test_output_path );
        output_model_file_name_bis += in.name() + "_saved_out_bis.gm";
        geomodel_save( in2, output_model_file_name_bis );

        if( !compare_files( output_model_file_name, output_model_file_name_bis ) ) {
            throw RINGMeshException( "TEST", "FAILED" );
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
