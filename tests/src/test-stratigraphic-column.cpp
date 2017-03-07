/*
 * test-stratigraphic-column.cpp
 *
 *  Created on: Feb 14, 2017
 *      Author: sirvent1u
 */

#include <ringmesh/ringmesh_tests_config.h>
#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/basic/common.h>
#include <ringmesh/io/io.h>
#include <iostream>
#include <ringmesh/geomodel/stratigraphic_column.h>

int main() {
	using namespace RINGMesh;

	try {
		default_configure();
		std::string input_model_file_name(ringmesh_test_data_path);
		input_model_file_name += "corbi_layers.ml";

		GeoModel in;
		bool loaded_model_is_valid = geomodel_load(in, input_model_file_name);

		index_t nb_interface = in.nb_geological_entities(
				Interface::type_name_static());
		std::cout << "nb_interface: " << nb_interface << std::endl;
		RockFeature rock("rock");
		RockFeature rock_two("rock", NONE);
		if (rock.get_name() != "rock") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing RockFeature::get_name() ");
		}
		if (rock_two.get_rock_type() != NONE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing RockFeature::get_rock_type()");
		}

		rock.set_rock_type(NONE);
		if (rock.get_rock_type() != NONE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing RockFeature::get_rock_type()");
		}

		std::string one_name = "one";
		std::string two_name = "two";
		std::string three_name = "three";
		std::string four_name = "four";

		StratigraphicUnit one(one_name,
				in.geological_entity(Interface::type_name_static(), 1),
				in.geological_entity(Interface::type_name_static(), 0),
				in.geological_entity(Layer::type_name_static(), 0), CONFORMABLE,
				CONFORMABLE, rock, 0, 10);
		StratigraphicUnit two(two_name,
				in.geological_entity(Interface::type_name_static(), 2),
				in.geological_entity(Interface::type_name_static(), 1),
				in.geological_entity(Layer::type_name_static(), 1), CONFORMABLE,
				CONFORMABLE, rock, 0, 10);
		StratigraphicUnit three(three_name,
				in.geological_entity(Interface::type_name_static(), 3),
				in.geological_entity(Interface::type_name_static(), 2),
				in.geological_entity(Layer::type_name_static(), 2), CONFORMABLE,
				CONFORMABLE, rock, 0, 10);
		StratigraphicUnit four(four_name,
				in.geological_entity(Interface::type_name_static(), 11),
				in.geological_entity(Interface::type_name_static(), 3),
				in.geological_entity(Layer::type_name_static(), 3), CONFORMABLE,
				CONFORMABLE, rock, 0, 10);
		if (one.get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::get_name()");
		}
		if (!one.is_conformable_base()) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::is_conformable_base()");
		}
		if (!one.is_conformable_top()) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::is_conformable_top()");
		}

		if (one.get_relation_base() != CONFORMABLE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::get_relation_base()");
		}
		if (one.get_relation_top() != CONFORMABLE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::get_relation_top()");
		}
		if (one.get_min_thick() != 0) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::get_min_thick()");
		}
		if (one.get_max_thick() != 10) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicUnit::get_max_thick()");
		}
		std::vector<const StratigraphicUnit*> units;
		units.push_back(&one);
		units.push_back(&two);
		units.push_back(&three);
		units.push_back(&four);

		StratigraphicColumn test_two("test_two", units, CHRONOSTRATIGRAPHIC);
		StratigraphicColumn test("test");
		test.set_paradigm(CHRONOSTRATIGRAPHIC);
		if (test.get_paradigm() != CHRONOSTRATIGRAPHIC) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::set_paradigm()");
		}

		//StratigraphicColumn test_copy(test_two);

		if (test_two.get_unit_above(two)->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_unit_above()");
		}

		if (test_two.get_unit_below(three)->get_name() != "four") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_unit_below()");
		}

		test_two.remove_unit(two);
		if (test_two.get_unit(2)->get_name() != "three") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::remove_unit()");
		}

		test_two.insert_unit_below(one, two);
		if (test_two.get_unit(2)->get_name() != "two") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::insert_unit_below()");
		}

		if (test_two.get_top_unit()->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_top_unit()");
		};

		if (test_two.get_base_unit()->get_name() != "four") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_base_unit()");
		}

		test.insert_top_unit(one);
		if (test.get_top_unit()->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::insert_top_unit()");
		}

		test.insert_base_unit(four);
		if (test.get_base_unit()->get_name() != "four") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::insert_base_unit()");
		}

		if (test.find_unit("one")->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::find_unit()");
		}

		test.get_units();
		//    test.find_unit_from_rock_feature(rock);
		//  test.find_unit_from_rock_feature_name("rock");

		//test_two.get_units_between(&one, &three);
		std::cout << "avant vec de StratUnit mixed" << std::endl;
		std::vector<const StratigraphicUnit*> mixed;
		mixed.push_back(&one);
		mixed.push_back(&two);
		mixed.push_back(&test_two);
		std::cout << "apres vec de StratUnit mixed" << std::endl;

		StratigraphicColumn mix("mix", mixed, CHRONOSTRATIGRAPHIC);
		std::cout << "apres creation de StratiColumn mixed" << std::endl;
		if (!mix.is_conformable_base()) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::is_conformable_base()");
		}

		if (!mix.is_conformable_top()) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::is_conformable_top()");
		}

		if (mix.get_relation_base() != CONFORMABLE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_relation_base()");
		}

		if (mix.get_relation_top() != CONFORMABLE) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_relation_top()");
		}

		if (mix.get_min_thick() != 0) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_min_thick()");
		}

		if (mix.get_max_thick() != 60) {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_max_thick()");
		}
		std::cout << "get_max_thick" << std::endl;

		if (mix.get_unit_above(two)->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_unit_above()");
		}
		std::cout << "get_unit_above" << std::endl;

		if (mix.get_unit_below(two)->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_unit_below()");
		}
		std::cout << "get_unit_below" << std::endl;

		mix.remove_unit(two);
		if (mix.get_sub_column(2)->get_name() != "test_two") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::remove_unit()");
		}
		std::cout << "remove_unit" << std::endl;

		mix.insert_unit_below(one, two);
		if (mix.get_unit(2)->get_name() != "two") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::insert_unit_below()");
		}
		std::cout << "insert_unit_below" << std::endl;

		if (mix.get_top_unit()->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_top_unit()");
		}
		std::cout << "get_top_unit" << std::endl;

		if (mix.get_base_unit()->get_name() != "four") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::get_base_unit()");
		}
		std::cout << "get_base_unit" << std::endl;

		if (mix.find_unit("one")->get_name() != "one") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::find_unit()");
		}
		std::cout << "find_unit" << std::endl;

		if (mix.find_sub_column("test_two")->get_name() != "test_two") {
			throw RINGMeshException("RINGMesh Test",
					"Failed when testing StratigraphicColumn::find_sub_unit()");
		}
		std::cout << "find_sub_unit" << std::endl;

		mix.get_units();

	} catch (const RINGMeshException& e) {
		Logger::err(e.category()) << e.what() << std::endl;
		return 1;
	} catch (const std::exception& e) {
		Logger::err("Exception") << e.what() << std::endl;
		return 1;
	}
	Logger::out("TEST") << "SUCCESS" << std::endl;
	return 0;
}

