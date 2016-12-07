/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geomodel/geo_model_repair.h>
#include <ringmesh/geomodel/geo_model_validity.h>
#include <ringmesh/io/io.h>

/*! 
 * Load and fix a given structural model file.
 * @author Jeanne Pellerin
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

        GeoModel M ;
        std::string file_name( ringmesh_test_data_path ) ;
        file_name += "annot.ml" ;
        
        Logger::out( "RINGMesh Test" ) << "Loading and fixing structural model:"
            << file_name << std::endl ;

        // Set the debug directory for the validity checks
        set_validity_errors_directory( ringmesh_test_output_path ) ;

        // Load the model
        bool init_model_is_valid = geomodel_load( M, file_name ) ;
        if( init_model_is_valid ) {
            throw RINGMeshException( "RINGMesh Test",
                "Input test model " + M.name()
                    + " must be invalid to check the repair functionalities." ) ;
        }

        Logger::out( "RINGMesh Test" ) << "Repairing "
            << std::endl << std::endl << std::endl ;
        // Repair the model

        GeoModelRepair model_repair( M ) ;
        model_repair.repair( GeoModelRepair::ALL ) ;

        // Test the validity again
        if( is_geomodel_valid( M ) ) {
            std::string fixed_file_name( ringmesh_test_output_path ) ;
            fixed_file_name += M.name() + "_repaired.ml" ;
            geomodel_save( M, fixed_file_name ) ;
            Logger::out( "RINGMesh Test" ) << "Invalid geological model "
                << M.name()
                << " has been successfully fixed and is saved under: "
                << fixed_file_name << std::endl ;
            Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
            return 0 ;
        } else {
            throw RINGMeshException( "RINGMesh Test",
                "Fixing the invalid geological model " + M.name()
                    + " failed." ) ;
        }

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
}
