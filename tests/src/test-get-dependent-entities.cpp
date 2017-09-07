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

#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * Test of the method GeoModelBuilderTopology::get_dependent_entities.
 * @author Benjamin Chauvin
 */

using namespace RINGMesh;

template < typename GME >
void check_element_of_a_set_are_in_another_set(
    const std::set< GME >& to_compare,
    const std::set< GME >& with,
    const std::string& set_name )
{
    if( to_compare != with )
    {
        // To debug it is nice to know which entity fails...
        for( const GME& cur_gme_id : to_compare )
        {
            if( find( with.begin(), with.end(), cur_gme_id ) == with.end() )
            {
                throw RINGMeshException( "RINGMesh Test", cur_gme_id.type(),
                    " ", cur_gme_id.index(), " is not in the ", set_name, "." );
            }
        }
        ringmesh_assert_not_reached;
    }
}

void test_template( GeoModel3D& geomodel,
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
            mesh_type ) )
    {
        in_mesh_entities.insert( { mesh_type, to_insert_id } );
    }
    else
    {
        const GeologicalEntityType geological_type( to_insert_type );
        ringmesh_assert(
            geomodel.entity_type_manager()
                .geological_entity_manager.is_valid_type( geological_type ) );
        in_geological_entities.insert( { geological_type, to_insert_id } );
    }

    const GeoModelBuilder3D geomodel_builder( geomodel );
    geomodel_builder.topology.get_dependent_entities(
        in_mesh_entities, in_geological_entities );
    check_element_of_a_set_are_in_another_set< gmme_id >(
        in_mesh_entities, solution_gmme_id, "solution" );
    check_element_of_a_set_are_in_another_set< gmge_id >(
        in_geological_entities, solution_gmge_id, "solution" );
    check_element_of_a_set_are_in_another_set< gmme_id >(
        solution_gmme_id, in_mesh_entities, "output" );
    check_element_of_a_set_are_in_another_set< gmge_id >(
        solution_gmge_id, in_geological_entities, "output" );
}

void test_on_top_region( GeoModel3D& geomodel )
{
    // Solution:
    // Corners: 31, 33, 54, 55, 56, 57, 58, 93, 118, 128, 129.
    // Lines: 41, 43, 68, 69, 70, 71, 72, 73, 131, 135, 144, 177, 182, 203, 205,
    // 207, 210, 233, 234, 238.
    // Surfaces: 11, 37, 40, 60, 85, 91, 99, 110, 114.
    // Region: 4.
    std::set< gmme_id > solution_gmme_id = { { Corner3D::type_name_static(),
                                                 31 },
        { Corner3D::type_name_static(), 33 },
        { Corner3D::type_name_static(), 54 },
        { Corner3D::type_name_static(), 55 },
        { Corner3D::type_name_static(), 56 },
        { Corner3D::type_name_static(), 57 },
        { Corner3D::type_name_static(), 58 },
        { Corner3D::type_name_static(), 93 },
        { Corner3D::type_name_static(), 118 },
        { Corner3D::type_name_static(), 128 },
        { Corner3D::type_name_static(), 129 },
        { Line3D::type_name_static(), 41 }, { Line3D::type_name_static(), 43 },
        { Line3D::type_name_static(), 68 }, { Line3D::type_name_static(), 69 },
        { Line3D::type_name_static(), 70 }, { Line3D::type_name_static(), 71 },
        { Line3D::type_name_static(), 72 }, { Line3D::type_name_static(), 73 },
        { Line3D::type_name_static(), 131 },
        { Line3D::type_name_static(), 135 },
        { Line3D::type_name_static(), 144 },
        { Line3D::type_name_static(), 177 },
        { Line3D::type_name_static(), 182 },
        { Line3D::type_name_static(), 203 },
        { Line3D::type_name_static(), 205 },
        { Line3D::type_name_static(), 207 },
        { Line3D::type_name_static(), 210 },
        { Line3D::type_name_static(), 233 },
        { Line3D::type_name_static(), 234 },
        { Line3D::type_name_static(), 238 },
        { Surface3D::type_name_static(), 11 },
        { Surface3D::type_name_static(), 37 },
        { Surface3D::type_name_static(), 40 },
        { Surface3D::type_name_static(), 60 },
        { Surface3D::type_name_static(), 85 },
        { Surface3D::type_name_static(), 91 },
        { Surface3D::type_name_static(), 99 },
        { Surface3D::type_name_static(), 110 },
        { Surface3D::type_name_static(), 114 },
        { Region3D::type_name_static(), 4 } };

    // Solution:
    // Contacts: 26, 27, 28, 78, 79, 83.
    // Interface: 21.
    // Layer: 0.
    std::set< gmge_id > solution_gmge_id = {
        { Contact3D::type_name_static(), 26 },
        { Contact3D::type_name_static(), 27 },
        { Contact3D::type_name_static(), 28 },
        { Contact3D::type_name_static(), 78 },
        { Contact3D::type_name_static(), 79 },
        { Contact3D::type_name_static(), 83 },
        { Interface3D::type_name_static(), 21 },
        { Layer3D::type_name_static(), 0 },
    };

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Region3D::type_name_static().string(), 4 );
}

void test_on_surface_within_bottom_region_partially_connected_to_voi(
    GeoModel3D& geomodel )
{
    // Solution:
    // Corner: none.
    // Line: 98.
    // Surface: 24.
    // Region: none.
    std::set< gmme_id > solution_gmme_id = { { Line3D::type_name_static(), 98 },
        { Surface3D::type_name_static(), 24 } };

    // Solution:
    // Contact: 36.
    // Interface: 3.
    // Layer: none.
    std::set< gmge_id > solution_gmge_id = { { Contact3D::type_name_static(),
                                                 36 },
        { Interface3D::type_name_static(), 3 } };

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Surface3D::type_name_static().string(), 24 );
}

void test_on_fault_not_connected_to_any_surface( GeoModel3D& geomodel )
{
    // Solution:
    // Corner: none.
    // Lines: 163, 165, 166, 168.
    // Surfaces: 52, 53, 54, 55, 56, 57.
    // Region: none.
    std::set< gmme_id > solution_gmme_id = { { Line3D::type_name_static(),
                                                 163 },
        { Line3D::type_name_static(), 165 },
        { Line3D::type_name_static(), 166 },
        { Line3D::type_name_static(), 168 },
        { Surface3D::type_name_static(), 52 },
        { Surface3D::type_name_static(), 53 },
        { Surface3D::type_name_static(), 54 },
        { Surface3D::type_name_static(), 55 },
        { Surface3D::type_name_static(), 56 },
        { Surface3D::type_name_static(), 57 } };

    // Solution:
    // Contact: 55.
    // Interface: 8.
    // Layer: none.
    std::set< gmge_id > solution_gmge_id = { { Contact3D::type_name_static(),
                                                 55 },
        { Interface3D::type_name_static(), 8 } };

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Interface3D::type_name_static().string(), 8 );
}

void test_on_corner_on_botom_corner_voi( GeoModel3D& geomodel )
{
    // Solution:
    // Corner: 135.
    // Line: none.
    // Surface: none.
    // Region: none.
    std::set< gmme_id > solution_gmme_id = { { Corner3D::type_name_static(),
        135 } };

    // Solution:
    // Contact: none.
    // Interface: none.
    // Layer: none.
    std::set< gmge_id > solution_gmge_id;

    test_template( geomodel, solution_gmme_id, solution_gmge_id,
        Corner3D::type_name_static().string(), 135 );
}

void test_on_top_layer( GeoModel3D& geomodel )
{
    // Solution:
    // Corners: 31, 33, 54, 55, 56, 57, 58, 93, 118, 128, 129.
    // Lines: 41, 43, 68, 69, 70, 71, 72, 73, 131, 135, 144, 177, 182, 203, 205,
    // 207, 210, 233, 234, 238.
    // Surfaces: 11, 37, 40, 60, 85, 91, 99, 110, 114.
    // Region: 4.
    std::set< gmme_id > solution_gmme_id = { { Corner3D::type_name_static(),
                                                 31 },
        { Corner3D::type_name_static(), 33 },
        { Corner3D::type_name_static(), 54 },
        { Corner3D::type_name_static(), 55 },
        { Corner3D::type_name_static(), 56 },
        { Corner3D::type_name_static(), 57 },
        { Corner3D::type_name_static(), 58 },
        { Corner3D::type_name_static(), 93 },
        { Corner3D::type_name_static(), 118 },
        { Corner3D::type_name_static(), 128 },
        { Corner3D::type_name_static(), 129 },
        { Line3D::type_name_static(), 41 }, { Line3D::type_name_static(), 43 },
        { Line3D::type_name_static(), 68 }, { Line3D::type_name_static(), 69 },
        { Line3D::type_name_static(), 70 }, { Line3D::type_name_static(), 71 },
        { Line3D::type_name_static(), 72 }, { Line3D::type_name_static(), 73 },
        { Line3D::type_name_static(), 131 },
        { Line3D::type_name_static(), 135 },
        { Line3D::type_name_static(), 144 },
        { Line3D::type_name_static(), 177 },
        { Line3D::type_name_static(), 182 },
        { Line3D::type_name_static(), 203 },
        { Line3D::type_name_static(), 205 },
        { Line3D::type_name_static(), 207 },
        { Line3D::type_name_static(), 210 },
        { Line3D::type_name_static(), 233 },
        { Line3D::type_name_static(), 234 },
        { Line3D::type_name_static(), 238 },
        { Surface3D::type_name_static(), 11 },
        { Surface3D::type_name_static(), 37 },
        { Surface3D::type_name_static(), 40 },
        { Surface3D::type_name_static(), 60 },
        { Surface3D::type_name_static(), 85 },
        { Surface3D::type_name_static(), 91 },
        { Surface3D::type_name_static(), 99 },
        { Surface3D::type_name_static(), 110 },
        { Surface3D::type_name_static(), 114 },
        { Region3D::type_name_static(), 4 } };

    // Solution:
    // Contacts: 26, 27, 28, 78, 79, 83.
    // Interface: 21.
    // Layer: 0.
    std::set< gmge_id > solution_gmge_id = {
        { Contact3D::type_name_static(), 26 },
        { Contact3D::type_name_static(), 27 },
        { Contact3D::type_name_static(), 28 },
        { Contact3D::type_name_static(), 78 },
        { Contact3D::type_name_static(), 79 },
        { Contact3D::type_name_static(), 83 },
        { Interface3D::type_name_static(), 21 },
        { Layer3D::type_name_static(), 0 },
    };

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
    std::string file_name{ ringmesh_test_data_path };
    file_name += "CloudSpin_fixed.ml";

    // Load the model
    auto init_model_is_valid = geomodel_load( geomodel, file_name );
    if( !init_model_is_valid )
    {
        throw RINGMeshException( "RINGMesh Test", "Input test model ",
            geomodel.name(), " must be valid." );
    }
}

int main()
{
    try
    {
        default_configure();
        Logger::out(
            "TEST", "Test GeoModelBuilderTopology::get_dependent_entities" );

        GeoModel3D geomodel;
        load_geomodel( geomodel );
        run_tests( geomodel );
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
