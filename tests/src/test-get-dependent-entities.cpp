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

#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*! 
 * Test of the method GeoModelBuilderTopology::get_dependent_entities.
 * @author Benjamin Chauvin
 */

using namespace RINGMesh;

void run_tests( GeoModel& geomodel )
{
    GeoModelBuilder model_builder( geomodel );

    std::set< gmme_id > in_mesh_entities;
    in_mesh_entities.insert( gmme_id( Region::type_name_static(), 4 ) );
    std::set< gmge_id > in_geological_entities;
    model_builder.topology.get_dependent_entities( in_mesh_entities,
        in_geological_entities );

    // Solution:
    // Corners: 31, 33, 54, 55, 56, 57, 58, 93, 118, 128, 129.
    // Lines: 41, 43, 68, 69, 70, 71, 72, 73, 131, 135, 144, 177, 182, 203, 205, 207, 210, 233, 234, 238.
    // Surfaces: 11, 37, 40, 60, 85, 91, 99, 110, 114.
    // Region: 4.
    std::vector< gmme_id > solution_gmme_id = { gmme_id( Corner::type_name_static(), 31 ),
                                                gmme_id( Corner::type_name_static(), 33 ),
                                                gmme_id( Corner::type_name_static(), 54 ),
                                                gmme_id( Corner::type_name_static(), 55 ),
                                                gmme_id( Corner::type_name_static(), 56 ),
                                                gmme_id( Corner::type_name_static(), 57 ),
                                                gmme_id( Corner::type_name_static(), 58 ),
                                                gmme_id( Corner::type_name_static(), 93 ),
                                                gmme_id( Corner::type_name_static(), 118 ),
                                                gmme_id( Corner::type_name_static(), 128 ),
                                                gmme_id( Corner::type_name_static(), 129 ),
                                                gmme_id( Line::type_name_static(), 41 ),
                                                gmme_id( Line::type_name_static(), 43 ),
                                                gmme_id( Line::type_name_static(), 68 ),
                                                gmme_id( Line::type_name_static(), 69 ),
                                                gmme_id( Line::type_name_static(), 70 ),
                                                gmme_id( Line::type_name_static(), 71 ),
                                                gmme_id( Line::type_name_static(), 72 ),
                                                gmme_id( Line::type_name_static(), 73 ),
                                                gmme_id( Line::type_name_static(), 131 ),
                                                gmme_id( Line::type_name_static(), 135 ),
                                                gmme_id( Line::type_name_static(), 144 ),
                                                gmme_id( Line::type_name_static(), 177 ),
                                                gmme_id( Line::type_name_static(), 182 ),
                                                gmme_id( Line::type_name_static(), 203 ),
                                                gmme_id( Line::type_name_static(), 205 ),
                                                gmme_id( Line::type_name_static(), 207 ),
                                                gmme_id( Line::type_name_static(), 210 ),
                                                gmme_id( Line::type_name_static(), 233 ),
                                                gmme_id( Line::type_name_static(), 234 ),
                                                gmme_id( Line::type_name_static(), 238 ),
                                                gmme_id( Surface::type_name_static(), 11 ),
                                                gmme_id( Surface::type_name_static(), 37 ),
                                                gmme_id( Surface::type_name_static(), 40 ),
                                                gmme_id( Surface::type_name_static(), 60 ),
                                                gmme_id( Surface::type_name_static(), 85 ),
                                                gmme_id( Surface::type_name_static(), 91 ),
                                                gmme_id( Surface::type_name_static(), 99 ),
                                                gmme_id( Surface::type_name_static(), 110 ),
                                                gmme_id( Surface::type_name_static(), 114 ),
                                                gmme_id( Region::type_name_static(), 4 )
    };
    // Contacts: 26, 27, 28, 78, 79, 83.
    // Interface: 21.
    // Layer: 0.
    std::vector< gmge_id > solution_gmge_id = { gmge_id( Contact::type_name_static(), 26 ),
                                                gmge_id( Contact::type_name_static(), 27 ),
                                                gmge_id( Contact::type_name_static(), 28 ),
                                                gmge_id( Contact::type_name_static(), 78 ),
                                                gmge_id( Contact::type_name_static(), 79 ),
                                                gmge_id( Contact::type_name_static(), 83 ),
                                                gmge_id( Interface::type_name_static(), 21 ),
                                                gmge_id( Layer::type_name_static(), 0 ),
    };

    for( const gmme_id& cur_gmme_id : in_mesh_entities ) {
        //std::cout<<std::string(cur_gmme_id.type()) + " " + std::to_string( cur_gmme_id.index() )<<std::endl;
        if( std::find( solution_gmme_id.begin(), solution_gmme_id.end(),
            cur_gmme_id ) == solution_gmme_id.end() ) {
            throw RINGMeshException( "RINGMesh Test",
                std::string(cur_gmme_id.type()) + " " + std::to_string( cur_gmme_id.index() )
                    + " is not in the solution." );
        }
    }

    for( const gmge_id& cur_gmge_id : in_geological_entities ) {
        if( std::find( solution_gmge_id.begin(), solution_gmge_id.end(),
            cur_gmge_id ) == solution_gmge_id.end() ) {
            throw RINGMeshException( "RINGMesh Test",
                std::string(cur_gmge_id.type()) + " " + std::to_string( cur_gmge_id.index() )
                    + " is not in the solution." );
        }
    }

    for( const gmme_id& cur_gmme_id : solution_gmme_id ) {
        if( std::find( in_mesh_entities.begin(), in_mesh_entities.end(),
            cur_gmme_id ) == in_mesh_entities.end() ) {
            throw RINGMeshException( "RINGMesh Test",
                std::string(cur_gmme_id.type()) + " " + std::to_string( cur_gmme_id.index() )
                    + " is not in the output." );
        }
    }

    for( const gmge_id& cur_gmge_id : solution_gmge_id ) {
        if( std::find( in_geological_entities.begin(), in_geological_entities.end(),
            cur_gmge_id ) == in_geological_entities.end() ) {
            throw RINGMeshException( "RINGMesh Test",
                std::string(cur_gmge_id.type()) + " " + std::to_string( cur_gmge_id.index() )
                    + " is not in the output." );
        }
    }
}

void load_geomodel( GeoModel& geomodel )
{
    std::string file_name( ringmesh_test_data_path );
    file_name += "CloudSpin.ml";

    // Load the model
    bool init_model_is_valid = geomodel_load( geomodel, file_name );
    if( !init_model_is_valid ) {
        throw RINGMeshException( "RINGMesh Test",
            "Input test model " + geomodel.name() + " must be valid." );
    }
}

int main()
{
    try {
        default_configure();
        Logger::out( "TEST",
            "Test GeoModelBuilderTopology::get_dependent_entities" );

        GeoModel geomodel;
        load_geomodel( geomodel );
        run_tests( geomodel );

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
