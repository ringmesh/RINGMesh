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

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_api.h>
#include <ringmesh/geomodel/geo_model_builder_from_mesh.h>
#include <ringmesh/geomodel/geo_model_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel building from a Mesh
 * @author Jeanne Pellerin
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {
        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        // Set an output log file
        std::string log_file( ringmesh_test_output_path + "log.txt" ) ;
        GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
        Logger::instance()->register_client( file_logger ) ;

        {
            std::string file_name( ringmesh_test_data_path ) ;
            file_name += "split_cube.eobj" ;
            std::string result_file_name( ringmesh_test_output_path ) ;
            result_file_name += "split_cube_rebuilt.gm" ;

            // GeoModel from Surface with Attribute
            GEO::Mesh mesh ;
            GEO::MeshIOFlags mesh_io_flags ;
            mesh_io_flags.set_attribute( GEO::MESH_FACET_REGION ) ;
            // Warning: In an .eobj file Geogram loads only the facet integer attribute named "chart"
            // Used to fill the attribute called "region" on the mesh.
            GEO::mesh_load( file_name, mesh, mesh_io_flags ) ;

            GeoModel geomodel ;
            GeoModelBuilderMesh builder( geomodel, mesh, "region", "" ) ;
            RINGMesh::GeoModelBuildingFlags building_flags ;
            building_flags.compute_regions_brep = true ;
            builder.set_options( building_flags ) ;

            builder.create_and_build_surfaces() ;
            builder.build_model_from_surfaces() ;

            print_geomodel( geomodel ) ;

            // Checking if building has been successfully done
            if( !is_geomodel_valid( geomodel ) ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Failed when building model " + geomodel.name()
                        + ": the loaded model is not valid." ) ;
            }
            if( geomodel.mesh.vertices.nb() != mesh.vertices.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of vertices "
                        "than input mesh." ) ;
            }
            if( geomodel.mesh.facets.nb() != mesh.facets.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of facets "
                        "than input mesh." ) ;
            }
            if( geomodel.mesh.cells.nb() != mesh.cells.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of cells "
                        "than input mesh." ) ;
            }
            geomodel_save( geomodel, result_file_name ) ;
        }

        {
            // GeoModel from Surface connected components
            std::string file_name( ringmesh_test_data_path + "split_cube.obj" ) ;
            GEO::Mesh mesh ;
            GEO::mesh_load( file_name, mesh ) ;

            GeoModel geomodel ;
            GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components(
                mesh, "region" ) ;
            GeoModelBuilderMesh builder( geomodel, mesh, "region", "" ) ;
            RINGMesh::GeoModelBuildingFlags building_flags ;
            building_flags.compute_regions_brep = true ;
            builder.set_options( building_flags ) ;

            builder.create_and_build_surfaces() ;
            builder.build_model_from_surfaces() ;

            print_geomodel( geomodel ) ;

            // Checking if building has been successfully done
            if( !is_geomodel_valid( geomodel ) ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Failed when building model " + geomodel.name()
                        + ": the loaded model is not valid." ) ;
            }
            if( geomodel.mesh.vertices.nb() != mesh.vertices.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of vertices "
                        "than input mesh." ) ;
            }
            if( geomodel.mesh.facets.nb() != mesh.facets.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of facets "
                        "than input mesh." ) ;
            }
            if( geomodel.mesh.cells.nb() != mesh.cells.nb() ) {
                throw RINGMeshException( "RINGMesh Test",
                    "Error when building model: not same number of cells "
                        "than input mesh." ) ;
            }
        }

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
