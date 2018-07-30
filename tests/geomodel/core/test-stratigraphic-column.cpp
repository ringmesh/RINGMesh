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

#include <iostream>

#include <geogram/basic/command_line.h>

#include <ringmesh/basic/common.h>
#include <ringmesh/geomodel/builder/stratigraphic_column_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/stratigraphic_column.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

using namespace RINGMesh;

void test_rock_feature()
{
    RINGMesh::Logger::out(
        "RockFeature", "Test RockFeature building and editing" );

    RockFeature rock( "rock", ROCKTYPE::NONE );
    if( rock.get_rock_type() != ROCKTYPE::NONE )
    {
        throw RINGMeshException(
            "Test", "Failed when testing RockFeature::get_rock_type()" );
    }
    rock.set_rock_type( ROCKTYPE::MULTIPLE );
    if( rock.get_rock_type() != ROCKTYPE::MULTIPLE )
    {
        throw RINGMeshException(
            "Test", "Failed when changing RockFeature type" );
    }

    RockFeature rock_two( "rock2" );
    if( rock_two.get_name() != "rock2" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing RockFeature::get_name() " );
    }
    rock_two.set_name( "rock2_renamed" );
    if( rock_two.get_name() != "rock2_renamed" )
    {
        throw RINGMeshException(
            "Test", "Failed when renaming the RockFeature" );
    }
}

void test_stratigraphic_unit( const GeoModel3D& in )
{
    RINGMesh::Logger::out(
        "StratigraphicUnit", "Test StratigraphicUnit building" );

    RockFeature rock( "rock", ROCKTYPE::NONE );
    UnsubdividedStratigraphicUnit* test_strat_unit =
        new UnsubdividedStratigraphicUnit( "strat unit",
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 1 ) ),
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 0 ) ),
            dynamic_cast< const Layer3D& >(
                in.geological_entity( Layer3D::type_name_static(), 0 ) ),
            RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0, 10 );
    if( test_strat_unit->get_name() != "strat unit" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicUnit::get_name()" );
    }
    if( !test_strat_unit->is_conformable_base() )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicUnit::is_conformable_base()" );
    }
    if( !test_strat_unit->is_conformable_top() )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicUnit::is_conformable_top()" );
    }

    if( test_strat_unit->get_relation_base() != RELATION::CONFORMABLE )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicUnit::get_relation_base()" );
    }
    if( test_strat_unit->get_relation_top() != RELATION::CONFORMABLE )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicUnit::get_relation_top()" );
    }
    if( test_strat_unit->get_min_thick() != 0 )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicUnit::get_min_thick()" );
    }
    if( test_strat_unit->get_max_thick() != 10 )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicUnit::get_max_thick()" );
    }
}

void test_stratigraphic_column_building( const GeoModel3D& in )
{
    RINGMesh::Logger::out(
        "StratigraphicColumn", "Test StratigraphicColumn building" );
    RINGMesh::Logger::out(
        "StratigraphicColumn", "Init RockFeature and StratigraphicUnit" );

    RockFeature rock( "rock", ROCKTYPE::NONE );

    std::string one_name = "one";
    std::string two_name = "two";
    std::string three_name = "three";
    std::string four_name = "four";

    std::shared_ptr< const StratigraphicUnit > one(
        new UnsubdividedStratigraphicUnit( one_name,
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 1 ) ),
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 0 ) ),
            dynamic_cast< const Layer3D& >(
                in.geological_entity( Layer3D::type_name_static(), 0 ) ),
            RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0, 10 ) );
    std::shared_ptr< const StratigraphicUnit > two(
        new UnsubdividedStratigraphicUnit( two_name,
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 2 ) ),
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 1 ) ),
            dynamic_cast< const Layer3D& >(
                in.geological_entity( Layer3D::type_name_static(), 1 ) ),
            RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0, 20 ) );
    std::shared_ptr< const StratigraphicUnit > three(
        new UnsubdividedStratigraphicUnit( three_name,
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 3 ) ),
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 2 ) ),
            dynamic_cast< const Layer3D& >(
                in.geological_entity( Layer3D::type_name_static(), 2 ) ),
            RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0, 30 ) );
    std::shared_ptr< const StratigraphicUnit > four(
        new UnsubdividedStratigraphicUnit( four_name,
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 11 ) ),
            dynamic_cast< const Interface3D* >(
                &in.geological_entity( Interface3D::type_name_static(), 3 ) ),
            dynamic_cast< const Layer3D& >(
                in.geological_entity( Layer3D::type_name_static(), 3 ) ),
            RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0, 40 ) );

    RINGMesh::Logger::out( "StratigraphicColumn",
        "First building with a vector of StratigraphicUnit" );

    std::vector< std::shared_ptr< const StratigraphicUnit > > units;
    units.push_back( one );
    units.push_back( two );
    units.push_back( three );
    units.push_back( four );

    StratigraphicColumn test1(
        "test 1", units, STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC );
    if( test1.get_unit_above( *two )->get_name() != "one" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_unit_above()" );
    }

    if( test1.get_unit_below( *three )->get_name() != "four" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_unit_below()" );
    }

    test1.remove_unit( *two );
    if( test1.get_unit( 1 )->get_name() != "three" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicColumn::remove_unit()" );
    }

    test1.insert_unit_below( *one, two );
    if( test1.get_unit( 1 )->get_name() != "two" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::insert_unit_below()" );
    }

    if( test1.get_top_unit()->get_name() != "one" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicColumn::get_top_unit()" );
    };

    if( test1.get_base_unit()->get_name() != "four" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_base_unit()" );
    }

    RINGMesh::Logger::out( "StratigraphicColumn",
        "Second building with a vector of StratigraphicUnit" );

    StratigraphicColumn test2( "test 2",
        NestedStratigraphicUnit::StratigraphicUnits(),
        STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC );

    test2.insert_top_unit( one );
    if( test2.get_top_unit()->get_name() != "one" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::insert_top_unit()" );
    }

    test2.insert_base_unit( four );
    if( test2.get_base_unit()->get_name() != "four" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::insert_base_unit()" );
    }

    if( test2.get_unit( "one" )->get_name() != "one" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicColumn::get_unit()" );
    }

    RINGMesh::Logger::out(
        "StratigraphicColumn", "Third building a mixed StratigraphicColumn" );

    std::vector< std::shared_ptr< const StratigraphicUnit > > mixed;
    mixed.push_back( one );
    mixed.push_back( two );

    std::vector< std::shared_ptr< const StratigraphicUnit > > sub_units;
    RockFeature rocks( "several", ROCKTYPE::MULTIPLE );
    sub_units.push_back( three );
    sub_units.push_back( four );

    std::shared_ptr< const StratigraphicUnit > subdivided_unit(
        new NestedStratigraphicUnit( "subdivided", rocks, sub_units ) );
    mixed.push_back( subdivided_unit );

    StratigraphicColumn mix(
        "mix", mixed, STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC );
    if( !mix.is_conformable_base() )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::is_conformable_base()" );
    }

    if( !mix.is_conformable_top() )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::is_conformable_top()" );
    }

    if( mix.get_relation_base() != RELATION::CONFORMABLE )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_relation_base()" );
    }

    if( mix.get_relation_top() != RELATION::CONFORMABLE )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_relation_top()" );
    }

    if( mix.get_min_thick() != 0 )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_column_min_thick()" );
    }

    if( mix.get_max_thick() != 100 )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_column_max_thick()" );
    }

    if( mix.get_unit_above( *two )->get_name() != "one" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_unit_above()" );
    }

    if( mix.get_unit_below( *two )->get_name() != "subdivided" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_unit_below()" );
    }

    mix.remove_unit( *two );
    if( mix.get_unit( 1 )->get_name() != "subdivided" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicColumn::remove_unit()" );
    }

    mix.insert_unit_below( *one, two );
    if( mix.get_unit( 1 )->get_name() != "two" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::insert_unit_below()" );
    }

    if( mix.get_top_unit()->get_name() != "one" )
    {
        throw RINGMeshException(
            "Test", "Failed when testing StratigraphicColumn::get_top_unit()" );
    }

    if( mix.get_base_unit()->get_name() != "subdivided" )
    {
        throw RINGMeshException( "Test",
            "Failed when testing StratigraphicColumn::get_base_unit()" );
    }
}

void test_load_from_gocad_xml_file()
{
    std::string input_geomodel_file_name( ringmesh_test_data_path );
    input_geomodel_file_name += "CloudSpin_fixed.ml";
    GeoModel3D in;
    geomodel_load( in, input_geomodel_file_name );

    std::string input_column_file_name( RINGMesh::ringmesh_test_data_path );
    input_column_file_name += "CloudSpin.xml";

    StratigraphicColumn column( "test",
        NestedStratigraphicUnit::StratigraphicUnits(),
        STRATIGRAPHIC_PARADIGM::UNSPECIFIED );

    StratigraphicColumnBuilderXML sc_builder(
        column, in, input_column_file_name );
    sc_builder.load_file();
}

int main()
{
    using namespace RINGMesh;

    try
    {
        // Set an output log file
        std::string log_file( ringmesh_test_output_path );
        log_file += "log.txt";
        GEO::FileLogger* file_logger = new GEO::FileLogger( log_file );
        Logger::instance()->register_client( file_logger );

        // No validity checks at loading
        GEO::CmdLine::set_arg( "validity:do_not_check", "A" );

        // Load geomodel
        std::string input_model_file_name( ringmesh_test_data_path );
        input_model_file_name += "CloudSpin_fixed.ml";
        GeoModel3D in;
        geomodel_load( in, input_model_file_name );

        // test RockFeature
        test_rock_feature();

        // test Stratigraphic Unit building
        test_stratigraphic_unit( in );

        // test StratigraphicColumn building
        test_stratigraphic_column_building( in );

        // load StratigraphicColumn from gocad XML file
        test_load_from_gocad_xml_file();
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
