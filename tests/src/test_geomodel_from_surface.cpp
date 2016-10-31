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

#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/logger.h>

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_api.h>
#include <ringmesh/geomodel/geo_model_validity.h>
#include <ringmesh/geomodel/geo_model_builder_from_mesh.h>
#include <ringmesh/io/io.h>

#include <chrono>

/*! 
 * Test the creation of a GeoModel from a conformal surface mesh 
 * @todo Test on other datasets: nested spheres.
 * @author Jeanne Pellerin
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        std::string file_name = ringmesh_test_data_path ;
        file_name += "modelA6.mesh" ;

        // Set an output log file
        std::string log_file( ringmesh_test_output_path ) ;
        log_file += "log.txt" ;
        GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
        Logger::instance()->register_client( file_logger ) ;

        Logger::out( "TEST" ) << "Test GeoModel building from Surface"
            << std::endl ;

		auto t00 = std::chrono::steady_clock::now();
        GEO::Mesh in ;
        GEO::mesh_load( file_name, in ) ;


		auto t0 = std::chrono::steady_clock::now();

        GeoModel model ;
        GeoModelBuilderSurfaceMesh BB( model, in ) ;
        BB.build_polygonal_surfaces_from_connected_components() ;
        BB.build_model_from_surfaces() ;

		
		auto t1 = std::chrono::steady_clock::now();
		auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(t0 - t00);
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
		Logger::out("TIMING") << "Mesh load: " << duration0.count() << " milliseconds" << std::endl;
		Logger::out("TIMING") << "Model construction: " << duration.count() << " milliseconds" << std::endl;

        print_geomodel( model ) ;
        //GEO::CmdLine::set_arg( "in:intersection_check", false ) ;
        is_geomodel_valid( model, true ) ;

		std::string output_file_name(ringmesh_test_output_path);
		output_file_name += model.name() + "_reconstructed.gm";
		geomodel_save(model, output_file_name);	
    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;

}
