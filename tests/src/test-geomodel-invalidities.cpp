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
#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*! Tests the GeoModel invalidity tracking.
 * Loads a .ml file, check its validity. Then, we alter the geomodel in order to
 * break its validity and we check that we detect the invalidities.
 * @returns 0 if success or an error code if not.
 * @author Pierre Anquez
 */
int main()
{
    using namespace RINGMesh ;

    try {
        default_configure() ;

        // Set an output log file
        std::string log_file( ringmesh_test_output_path ) ;
        log_file += "log.txt" ;
        GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
        Logger::instance()->register_client( file_logger ) ;

        std::string input_model_file_name = ringmesh_test_data_path + "modelA6.ml" ;

        Logger::out( "TEST", "Geomodel input. Loading file ",
            input_model_file_name ) ;

        GeoModel in ;
        bool loaded_model_is_valid = geomodel_load( in, input_model_file_name ) ;

        if( !loaded_model_is_valid ) {
            throw RINGMeshException( "RINGMesh Test",
                "Failed when loading model " + in.name()
                    + ": the loaded model is not valid." ) ;
        }


        GeoModel invalid_model ;
        GeoModelBuilder geomodel_copier( invalid_model ) ;
        geomodel_copier.copy.copy_geomodel( in ) ;
        geomodel_copier.info.set_geomodel_name( "broken model 1" ) ;

        GeoModelBuilder geomodel_breaker( invalid_model ) ;
        geomodel_breaker.topology.add_mesh_entity_boundary(
            invalid_model.surface( 0 ).gmme(), 4 ) ;
        bool broken_model_is_valid = is_geomodel_valid(invalid_model) ;
        DEBUG(broken_model_is_valid) ;

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() ) ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() ) ;
        return 1 ;
    }
    Logger::out( "TEST", "SUCCESS" ) ;
    return 0 ;

}
