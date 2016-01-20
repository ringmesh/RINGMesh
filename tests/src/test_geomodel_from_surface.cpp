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
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
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

#include <ringmesh/io.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geo_model_builder.h>


/*! 
 * Test the creation of a GeoModel from a conformal surface mesh 
 * @todo Test on other datasets: nested spheres.
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    std::string file_name = ringmesh_test_data_path ;
    file_name += "modelA6.mesh" ;

    // Set an output log file
    std::string log_file( ringmesh_test_output_path ) ;
    log_file += "log.txt" ;
    GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
    GEO::Logger::instance()->register_client( file_logger ) ;

    GEO::Logger::out( "TEST" ) << "Test GeoModel building from Surface" << std::endl ;

    GEO::Mesh in ;
    GEO::mesh_load( file_name, in ) ;
    GeoModel model ;
	
    GeoModelBuilderSurfaceMesh BB( model, in ) ;
    BB.build_polygonal_surfaces_from_connected_components() ;
    if( !BB.build_model_from_surfaces() ) {
		GEO::Logger::out("TEST") << "FAILED" << std::endl ;	
		return 1 ;
	}    
    print_geomodel( model ) ;
    is_geomodel_valid( model, false ) ;
	GEO::Logger::out("TEST") << "SUCCESS" << std::endl ;
	return 0 ;
   
 }
