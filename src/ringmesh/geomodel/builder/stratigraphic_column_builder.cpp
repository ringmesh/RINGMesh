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

#include <ringmesh/geomodel/builder/stratigraphic_column_builder.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <tinyxml2.h>

namespace RINGMesh
{
    StratigraphicColumnBuilderFile::StratigraphicColumnBuilderFile(
        StratigraphicColumn& column, GeoModel3D& model, std::string filename )
        : StratigraphicColumnBuilder( column, model ),
          filename_( std::move( filename ) )
    {
    }

    StratigraphicColumnBuilder::StratigraphicColumnBuilder(
        StratigraphicColumn& column, GeoModel3D& model )
        : column_( column ), model_( model )
    {
        if( model_.nb_geological_entities( GeologicalEntityType( "Layer" ) )
            == 0 )
        {
            throw RINGMeshException(
                "I/O", "The GeoModel have to be defined with layers." );
        }
    }

    void StratigraphicColumnBuilderXML::load_file()
    {
        tinyxml2::XMLDocument column;
        tinyxml2::XMLError Result = column.LoadFile( filename_.c_str() );
        if( Result != tinyxml2::XML_SUCCESS )
        {
            throw RINGMeshException(
                "I/O", "Error while loading Stratigraphic Column XML file." );
        }
        tinyxml2::XMLNode* root = column.FirstChild();
        if( root == nullptr )
        {
            throw RINGMeshException( "I/O",
                "Error while getting root of Stratigraphic Column XML file." );
        }
        tinyxml2::XMLElement* local =
            root->FirstChildElement( "LocalStratigraphicColumn" );
        tinyxml2::XMLElement* name_column = local->FirstChildElement( "name" );
        const std::string name_of_column = name_column->GetText();

        tinyxml2::XMLElement* paradigm =
            local->FirstChildElement( "classification_type" );
        std::string paradigm_str = paradigm->GetText();

        tinyxml2::XMLElement* units = local->FirstChildElement( "units" );
        tinyxml2::XMLElement* unit = units->FirstChildElement( "unit" );

        std::vector< std::string > unitList;
        while( unit != nullptr )
        {
            tinyxml2::XMLElement* name = unit->FirstChildElement( "name" );
            unitList.emplace_back( name->GetText() );
            tinyxml2::XMLElement* top = unit->FirstChildElement( "top" );
            if( top != nullptr )
            {
                tinyxml2::XMLElement* name_top =
                    top->FirstChildElement( "name" );
                unitList.emplace_back( name_top->GetText() );
            }
            else
            {
                unitList.emplace_back( "none" );
            }
            tinyxml2::XMLElement* base = unit->FirstChildElement( "base" );
            if( base != nullptr )
            {
                tinyxml2::XMLElement* name_base =
                    base->FirstChildElement( "name" );
                unitList.emplace_back( name_base->GetText() );
            }
            else
            {
                unitList.emplace_back( "none" );
            }
            unit = unit->NextSiblingElement( "unit" );
        }

        // Creation of StratigraphicUnit

        std::vector< std::shared_ptr< const StratigraphicUnit > >
            units_vec_construction;
        for( index_t i = 0; i < unitList.size(); i += 3 )
        {
            const std::string& name_of_unit = unitList[i];
            if( name_of_unit != "none" )
            {
                index_t layer_id = find_geological_entity_id_from_name(
                    model_, GeologicalEntityType( "Layer" ), name_of_unit );
                const Layer< 3 >* layer = dynamic_cast< const Layer< 3 >* >(
                    &( model_.geological_entity(
                        GeologicalEntityType( "Layer" ), layer_id ) ) );
                ringmesh_assert( layer != nullptr );
                const Interface< 3 >* top_interface = nil;
                const Interface< 3 >* base_interface = nil;
                RockFeature rock( name_of_unit );
                if( unitList[i + 1] != "none" )
                {
                    std::string name_of_interface_top = unitList[i + 1];
                    index_t top_interface_id =
                        find_geological_entity_id_from_name( model_,
                            GeologicalEntityType( "Interface" ),
                            name_of_interface_top );
                    top_interface = dynamic_cast< const Interface< 3 >* >(
                        &( model_.geological_entity(
                            GeologicalEntityType( "Interface" ),
                            top_interface_id ) ) );
                    ringmesh_assert( layer != nullptr );
                }
                if( unitList[i + 2] != "none" )
                {
                    std::string name_of_interface_base = unitList[i + 2];
                    index_t base_interface_id =
                        find_geological_entity_id_from_name( model_,
                            GeologicalEntityType( "Interface" ),
                            name_of_interface_base );
                    base_interface = dynamic_cast< const Interface< 3 >* >(
                        &( model_.geological_entity(
                            GeologicalEntityType( "Interface" ),
                            base_interface_id ) ) );
                    ringmesh_assert( layer != nullptr );
                }
                std::shared_ptr< const StratigraphicUnit > unit(
                    new UnsubdividedStratigraphicUnit( name_of_unit,
                        top_interface, base_interface, *layer,
                        RELATION::CONFORMABLE, RELATION::CONFORMABLE, rock, 0,
                        std::numeric_limits< double >::max() ) );
                units_vec_construction.push_back( unit );
            }
        }
        const std::vector< std::shared_ptr< const StratigraphicUnit > >
            units_vec = units_vec_construction;
        STRATIGRAPHIC_PARADIGM paradigm_upper;
        if( paradigm_str == "chronostratigraphy" )
        {
            paradigm_upper = STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC;
        }
        else if( paradigm_str == "lithostratigraphy" )
        {
            paradigm_upper = STRATIGRAPHIC_PARADIGM::LITHOSTRATIGRAPHIC;
        }
        else
        {
            paradigm_upper = STRATIGRAPHIC_PARADIGM::BIOSTRATIGRAPHIC;
        }
        column_ =
            StratigraphicColumn( name_of_column, units_vec, paradigm_upper );
    }

} // namespace RINGMesh
