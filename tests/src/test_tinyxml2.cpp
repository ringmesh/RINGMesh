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

#include <third_party/tinyxml2/tinyxml2.h>
#include <ringmesh/basic/common.h>
#include <ringmesh/ringmesh_tests_config.h>

using namespace RINGMesh;

int main() {

	tinyxml2::XMLDocument column;
	std::string input_model_file_name(ringmesh_test_data_path);
	input_model_file_name += "column_gocad.xml";
	tinyxml2::XMLError Result = column.LoadFile(input_model_file_name.c_str());
	if (Result != tinyxml2::XML_SUCCESS) {
		std::cout << "Error while loading file" << std::endl;
	}
	tinyxml2::XMLNode *root = column.FirstChild();
	if (root == nil) {
		std::cout << "Error when getting root" << std::endl;
	}

	tinyxml2::XMLElement* local = root->FirstChildElement("LocalStratigraphicColumn");
	tinyxml2::XMLElement* name_column = local->FirstChildElement("name");
	const char *name_of_column = name_column->GetText();
	std::cout << "name of column:" << name_of_column << std::endl;

	tinyxml2::XMLElement* paradigm = local->FirstChildElement("classification_type");
	std::cout << "paradigm:" << paradigm->GetText() << std::endl;

	tinyxml2::XMLElement* units = local->FirstChildElement("units");
	tinyxml2::XMLElement* unit = units->FirstChildElement("unit");

	std::vector<const char*> unitList;
	while (unit != nil) {
		tinyxml2::XMLElement* name = unit->FirstChildElement("name");
		//const char* NameValue;
		//const char * TopValue;
		const char* BaseValue;
		std::cout << "name of unit:" << name->GetText() << std::endl;
		unitList.push_back(name->GetText());
		tinyxml2::XMLElement* top = unit->FirstChildElement("top");
		if (top != nil) {
			tinyxml2::XMLElement* name_top = top->FirstChildElement("name");
			unitList.push_back(name_top->GetText());
			std::cout << "top of unit:" << name_top->GetText() << std::endl;
		} else {
			unitList.push_back("none");
			std::cout << "top = none" << std::endl;
		}
		tinyxml2::XMLElement* base = unit->FirstChildElement("base");
		if (base != nil) {
			tinyxml2::XMLElement* name_base = base->FirstChildElement("name");
			unitList.push_back(name_base->GetText());
			std::cout << "base of unit:" << name_base->GetText() << std::endl;
		} else {
			unitList.push_back("none");
			std::cout << "base = none" << std::endl;
		}
		unit = unit->NextSiblingElement("unit");

	}

	return 0;
}
