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

#include <geogram/basic/command_line.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * Test of the method GeoModelBuilderTopology::get_dependent_entities.
 * @author Benjamin Chauvin
 */

using namespace RINGMesh;

template< typename GME >
void check_element_of_a_set_are_in_another_set(
    const std::set< GME >& to_compare,
    const std::set< GME >& with,
    const std::string& set_name )
{
    if( to_compare != with ) {
        // To debug it is nice to know which entity fails...
        for( const GME& cur_gme_id : to_compare ) {
            if( find( with.begin(), with.end(), cur_gme_id ) == with.end() ) {
                throw RINGMeshException( "RINGMesh Test", cur_gme_id.type(), " ",
                    cur_gme_id.index(), " is not in the ", set_name, "." );
            }
        }
        ringmesh_assert_not_reached;
    }
}

void test_template(
    GeoModel3D& geomodel,
    const std::set< gmme_id >& solution_gmme_id,
    const std::set< gmge_id >& solution_gmge_id,
    const std::string& to_insert_type,
    index_t to_insert_id )
{
    const GeoModelBuilder3D model_builder( geomodel );
    std::set< gmme_id > in_mesh_entities;
    std::set< gmge_id > in_geological_entities;

    const MeshEntityType mesh_type( to_insert_type );
    if( geomodel.entity_type_manager().mesh_entity_manager.is_valid_type(
        mesh_type ) ) {
        in_mesh_entities.insert( { mesh_type, to_insert_id } );
    } else {
        const GeologicalEntityType geological_type( to_insert_type );
        ringmesh_assert(
            geomodel.entity_type_manager().geological_entity_manager.is_valid_type(
                geological_type ) );
        in_geological_entities.insert( { geological_type, to_insert_id } );
    }

    const GeoModelBuilder3D geomodel_builder( geomodel );
    geomodel_builder.topology.get_dependent_entities( in_mesh_entities,
        in_geological_entities );
    check_element_of_a_set_are_in_another_set< gmme_id >( in_mesh_entities,
        solution_gmme_id, "solution" );
    check_element_of_a_set_are_in_another_set< gmge_id >( in_geological_entities,
        solution_gmge_id, "solution" );
    check_element_of_a_set_are_in_another_set< gmme_id >( solution_gmme_id,
        in_mesh_entities, "output" );
    check_element_of_a_set_are_in_another_set< gmge_id >( solution_gmge_id,
        in_geological_entities, "output" );
}

template< class ENTITY, typename T >
void add_entities_in_set( std::set< T >& set )
{
    // Do nothing
    ringmesh_unused( set );
}

template< class ENTITY, typename T >
void add_entity_in_set( std::set< T >& set, int entity_id )
{
    set.emplace( ENTITY::type_name_static(), entity_id );
}

template< class ENTITY, typename T, typename ... Args >
void add_entities_in_set(
    std::set< T >& set,
    const int& first_id,
    const Args&... other_ids )
{
    add_entity_in_set< ENTITY >( set, first_id );
    add_entities_in_set< ENTITY >( set, other_ids... );
}

void test_on_top_region( GeoModel3D& geomodel )
{
    std::set< gmme_id > solution_gmme_id;
    add_entities_in_set< Corner3D >( solution_gmme_id, 31, 33, 54, 55, 56, 57, 58,
        99, 125, 138, 139 );
    add_entities_in_set< Line3D >( solution_gmme_id, 41, 43, 68, 69, 70, 71, 72, 73,
        137, 141, 151, 184, 189, 213, 215, 217, 220, 243, 244, 248 );
    add_entities_in_set< Surface3D >( solution_gmme_id, 11, 37, 40, 60, 85, 91, 99,
        110, 114 );
    add_entities_in_set< Region3D >( solution_gmme_id, 4 );

    std::set< gmge_id > solution_gmge_id;
    add_entities_in_set< Contact3D >( solution_gmge_id, 26, 27, 28, 78, 79, 83 );
    add_entities_in_set< Interface3D >( solution_gmge_id, 21 );
    add_entities_in_set< Layer3D >( solution_gmge_id, 0 );

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Region3D::type_name_static().string(), 4 );
}

void test_on_surface_within_bottom_region_partially_connected_to_voi(
    GeoModel3D& geomodel )
{
    std::set< gmme_id > solution_gmme_id = { { Line3D::type_name_static(), 103 }, {
        Surface3D::type_name_static(), 24 } };

    std::set< gmge_id > solution_gmge_id = { { Contact3D::type_name_static(), 36 }, {
        Interface3D::type_name_static(), 3 } };

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Surface3D::type_name_static().string(), 24 );
}

void test_on_fault_not_connected_to_any_surface( GeoModel3D& geomodel )
{
    std::set< gmme_id > solution_gmme_id;
    add_entities_in_set< Line3D >( solution_gmme_id, 170, 172, 173, 175 );
    add_entities_in_set< Surface3D >( solution_gmme_id, 52, 53, 54, 55, 56, 57 );

    std::set< gmge_id > solution_gmge_id = { { Contact3D::type_name_static(), 55 }, {
        Interface3D::type_name_static(), 8 } };

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Interface3D::type_name_static().string(), 8 );
}

void test_on_corner_on_botom_corner_voi( GeoModel3D& geomodel )
{
    std::set< gmme_id > solution_gmme_id = { { Corner3D::type_name_static(), 135 } };

    // No entity
    std::set< gmge_id > solution_gmge_id;

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Corner3D::type_name_static().string(), 135 );
}

void test_on_top_layer( GeoModel3D& geomodel )
{

    std::set< gmme_id > solution_gmme_id;
    add_entities_in_set< Corner3D >( solution_gmme_id, 31, 33, 54, 55, 56, 57, 58,
        99, 125, 138, 139 );
    add_entities_in_set< Line3D >( solution_gmme_id, 41, 43, 68, 69, 70, 71, 72, 73,
        137, 141, 151, 184, 189, 213, 215, 217, 220, 243, 244, 248 );
    add_entities_in_set< Surface3D >( solution_gmme_id, 11, 37, 40, 60, 85, 91, 99,
        110, 114 );
    add_entities_in_set< Region3D >( solution_gmme_id, 4 );

    std::set< gmge_id > solution_gmge_id;
    add_entities_in_set< Contact3D >( solution_gmge_id, 26, 27, 28, 78, 79, 83 );
    add_entities_in_set< Interface3D >( solution_gmge_id, 21 );
    add_entities_in_set< Layer3D >( solution_gmge_id, 0 );

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Layer3D::type_name_static().string(), 0 );
}

void run_tests( GeoModel3D& geomodel )
{
    test_on_top_region( geomodel );
    test_on_surface_within_bottom_region_partially_connected_to_voi( geomodel );
    test_on_fault_not_connected_to_any_surface( geomodel );
    test_on_corner_on_botom_corner_voi( geomodel );
    test_on_top_layer( geomodel ); // Should be the same as test_on_top_region
}

void load_geomodel( GeoModel3D& geomodel )
{
    std::string file_name { ringmesh_test_data_path };
    file_name += "CloudSpin_fixed.ml";

    // No validity checks at loading
    GEO::CmdLine::set_arg( "validity:do_not_check", "A" );

    // Load the model
    geomodel_load( geomodel, file_name );
}

int main()
{
    try {
        Logger::out( "TEST",
            "Test GeoModelBuilderTopology::get_dependent_entities" );

        GeoModel3D geomodel;
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
