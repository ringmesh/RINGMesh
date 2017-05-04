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
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/basic/common.h>
#include <ringmesh/io/io.h>
#include <iostream>
#include <ringmesh/geomodel/stratigraphic_column.h>
#include <ringmesh/geomodel/stratigraphic_column_builder.h>

using namespace RINGMesh;

void test_rock_feature()
{
    RINGMesh::Logger::out( "RockFeature" ) << "Test RockFeature building and editing"
        << std::endl;

    RockFeature rock( "rock", NONE );
    if( rock.get_rock_type() != NONE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing RockFeature::get_rock_type()" );
    }
    rock.set_rock_type( MULTIPLE );
    if( rock.get_rock_type() != MULTIPLE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when changing RockFeature type" );
    }

    RockFeature rock_two( "rock2" );
    if( rock_two.get_name() != "rock2" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing RockFeature::get_name() " );
    }
    rock_two.set_name( "rock2_renamed" );
    if( rock_two.get_name() != "rock2_renamed" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when renaming the RockFeature" );
    }
}

void test_stratigraphic_unit( const GeoModel& in )
{
    RINGMesh::Logger::out( "StratigraphicUnit" ) << "Test StratigraphicUnit building"
        << std::endl;

    RockFeature rock( "rock", NONE );
    StratigraphicUnit test_strat_unit( "strat unit",
        in.geological_entity( Interface::type_name_static(), 1 ),
        in.geological_entity( Interface::type_name_static(), 0 ),
        in.geological_entity( Layer::type_name_static(), 0 ), CONFORMABLE,
        CONFORMABLE, rock, 0, 10 );
    if( test_strat_unit.get_name() != "strat unit" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::get_name()" );
    }
    if( !test_strat_unit.is_conformable_base() ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::is_conformable_base()" );
    }
    if( !test_strat_unit.is_conformable_top() ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::is_conformable_top()" );
    }

    if( test_strat_unit.get_relation_base() != CONFORMABLE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::get_relation_base()" );
    }
    if( test_strat_unit.get_relation_top() != CONFORMABLE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::get_relation_top()" );
    }
    if( test_strat_unit.get_min_thick() != 0 ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::get_min_thick()" );
    }
    if( test_strat_unit.get_max_thick() != 10 ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicUnit::get_max_thick()" );
    }
}

void test_stratigraphic_column_building( const GeoModel& in )
{
    RINGMesh::Logger::out( "StratigraphicColumn" )
        << "Test StratigraphicColumn building" << std::endl;
    RINGMesh::Logger::out( "StratigraphicColumn" )
        << "Init RockFeature and StratigraphicUnit" << std::endl;

    RockFeature rock( "rock", NONE );

    std::string one_name = "one";
    std::string two_name = "two";
    std::string three_name = "three";
    std::string four_name = "four";

    StratigraphicUnit one( one_name,
        in.geological_entity( Interface::type_name_static(), 1 ),
        in.geological_entity( Interface::type_name_static(), 0 ),
        in.geological_entity( Layer::type_name_static(), 0 ), CONFORMABLE,
        CONFORMABLE, rock, 0, 10 );
    StratigraphicUnit two( two_name,
        in.geological_entity( Interface::type_name_static(), 2 ),
        in.geological_entity( Interface::type_name_static(), 1 ),
        in.geological_entity( Layer::type_name_static(), 1 ), CONFORMABLE,
        CONFORMABLE, rock, 0, 10 );
    StratigraphicUnit three( three_name,
        in.geological_entity( Interface::type_name_static(), 3 ),
        in.geological_entity( Interface::type_name_static(), 2 ),
        in.geological_entity( Layer::type_name_static(), 2 ), CONFORMABLE,
        CONFORMABLE, rock, 0, 10 );
    StratigraphicUnit four( four_name,
        in.geological_entity( Interface::type_name_static(), 11 ),
        in.geological_entity( Interface::type_name_static(), 3 ),
        in.geological_entity( Layer::type_name_static(), 3 ), CONFORMABLE,
        CONFORMABLE, rock, 0, 10 );

    RINGMesh::Logger::out( "StratigraphicColumn" )
        << "First building with a vector of StratigraphicUnit" << std::endl;

    std::vector< const StratigraphicUnit* > units;
    units.push_back( &one );
    units.push_back( &two );
    units.push_back( &three );
    units.push_back( &four );

    StratigraphicColumn test1( "test 1", units, CHRONOSTRATIGRAPHIC );
    if( test1.get_unit_above( two )->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_unit_above()" );
    }

    if( test1.get_unit_below( three )->get_name() != "four" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_unit_below()" );
    }

    test1.remove_unit( two );
    if( test1.get_unit( 2 )->get_name() != "three" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::remove_unit()" );
    }

    test1.insert_unit_below( one, two );
    if( test1.get_unit( 2 )->get_name() != "two" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::insert_unit_below()" );
    }

    if( test1.get_top_unit()->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_top_unit()" );
    };

    if( test1.get_base_unit()->get_name() != "four" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_base_unit()" );
    }

    RINGMesh::Logger::out( "StratigraphicColumn" )
        << "Second building vith a vector of StratigraphicUnit" << std::endl;

    StratigraphicColumn test2( "test 2" );
    test2.set_paradigm( CHRONOSTRATIGRAPHIC );
    if( test2.get_paradigm() != CHRONOSTRATIGRAPHIC ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::set_paradigm()" );
    }

    test2.insert_top_unit( one );
    if( test2.get_top_unit()->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::insert_top_unit()" );
    }

    test2.insert_base_unit( four );
    if( test2.get_base_unit()->get_name() != "four" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::insert_base_unit()" );
    }

    if( test2.find_unit( "one" )->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::find_unit()" );
    }

    RINGMesh::Logger::out( "StratigraphicColumn" )
        << "Third building a mixed StratigraphicColumn" << std::endl;

    std::vector< const StratigraphicUnit* > mixed;
    mixed.push_back( &one );
    mixed.push_back( &two );
    mixed.push_back( &test1 );
    std::cout << "apres vec de StratUnit mixed" << std::endl;

    StratigraphicColumn mix( "mix", mixed, CHRONOSTRATIGRAPHIC );
    std::cout << "apres creation de StratiColumn mixed" << std::endl;
    if( !mix.is_conformable_base() ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::is_conformable_base()" );
    }

    if( !mix.is_conformable_top() ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::is_conformable_top()" );
    }

    if( mix.get_relation_base() != CONFORMABLE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_relation_base()" );
    }

    if( mix.get_relation_top() != CONFORMABLE ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_relation_top()" );
    }

    if( mix.get_min_thick() != 0 ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_min_thick()" );
    }

    if( mix.get_max_thick() != 60 ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_max_thick()" );
    }
    std::cout << "get_max_thick" << std::endl;

    if( mix.get_unit_above( two )->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_unit_above()" );
    }
    std::cout << "get_unit_above" << std::endl;

    if( mix.get_unit_below( two )->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_unit_below()" );
    }
    std::cout << "get_unit_below" << std::endl;

    mix.remove_unit( two );
    if( mix.get_sub_column( 2 )->get_name() != "test 1" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::remove_unit()" );
    }
    std::cout << "remove_unit" << std::endl;

    mix.insert_unit_below( one, two );
    if( mix.get_unit( 2 )->get_name() != "two" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::insert_unit_below()" );
    }
    std::cout << "insert_unit_below" << std::endl;

    if( mix.get_top_unit()->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_top_unit()" );
    }
    std::cout << "get_top_unit" << std::endl;

    if( mix.get_base_unit()->get_name() != "four" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::get_base_unit()" );
    }
    std::cout << "get_base_unit" << std::endl;

    if( mix.find_unit( "one" )->get_name() != "one" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::find_unit()" );
    }
    std::cout << "find_unit" << std::endl;

    if( mix.find_sub_column( "test 1" )->get_name() != "test 1" ) {
        throw RINGMeshException( "RINGMesh Test",
            "Failed when testing StratigraphicColumn::find_sub_unit()" );
    }
    std::cout << "find_sub_unit" << std::endl;
}

void test_load_from_gocad_xml_file()
{

    std::string input_geomodel_file_name( ringmesh_test_data_path );
    input_geomodel_file_name += "CloudSpin.ml";
    GeoModel in;
    geomodel_load( in, input_geomodel_file_name );

    std::string input_column_file_name( RINGMesh::ringmesh_test_data_path );
    input_column_file_name += "CloudSpin.xml";

    StratigraphicColumn column( "test" );

    StratigraphicColumnBuilderXML sc_builder( column, in, input_column_file_name );
    sc_builder.load_file();

}

int main()
{
    using namespace RINGMesh;

    try {
        default_configure();

        // Set an output log file
        std::string log_file( ringmesh_test_output_path );
        log_file += "log.txt";
        GEO::FileLogger* file_logger = new GEO::FileLogger( log_file );
        Logger::instance()->register_client( file_logger );

        //build geomodel
        std::string input_model_file_name( ringmesh_test_data_path );
        input_model_file_name += "corbi_layers.ml";
        GeoModel in;
        geomodel_load( in, input_model_file_name );

        index_t nb_interface = in.nb_geological_entities(
            Interface::type_name_static() );
        std::cout << "nb_interface: " << nb_interface << std::endl;

        //test RockFeature
        test_rock_feature();

        //test Stratigraphic Unit building
        test_stratigraphic_unit( in );

        //test StratigraphicColumn building
        test_stratigraphic_column_building( in );

        //load StratigraphicColumn from gocad XML file
        test_load_from_gocad_xml_file();
        system( "Pause" );
    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl;
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl;
        return 1;
    }
    Logger::out( "TEST" ) << "SUCCESS" << std::endl;
    return 0;
}

