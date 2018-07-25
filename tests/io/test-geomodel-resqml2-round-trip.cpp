/* * Copyright (c)2018, Association Scientifique pour la Geologie et ses
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

#include <geogram/basic/line_stream.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/core/stratigraphic_column.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel RESQML2 round trip
 * @author Wan-Chiu Li
 */

using namespace RINGMesh;
using GMGE = GeoModelGeologicalEntity< 3 >;

const std::string ringmesh_test_save_path =
    ringmesh_test_path + "io/data/save/";

void load_geomodel( GeoModel3D& geomodel, const std::string& file )
{
    bool loaded_model_is_valid = geomodel_load( geomodel, file );
    if( !loaded_model_is_valid )
    {
        throw RINGMeshException(
            "RINGMesh Test", "Failed when loading model from file " + file
                                 + ": the loaded model is not valid." );
    }

    if( geomodel.stratigraphic_column() == nullptr )
    {
        throw RINGMeshException(
            "TEST", "GeoModel does not contain a stratigraphic column." );
    }
}

void compare_geomodels(
    const GeoModel3D& geomodel1, const GeoModel3D& geomodel2 )
{
    const auto& all_mesh_entity_types =
        geomodel1.entity_type_manager().mesh_entity_manager.mesh_entity_types();
    for( const auto& cur_entity_type : all_mesh_entity_types )
    {
        if( geomodel1.nb_mesh_entities( cur_entity_type )
            != geomodel2.nb_mesh_entities( cur_entity_type ) )
        {
            throw RINGMeshException( "TEST", "Number of entities of type ",
                cur_entity_type.string(), " of two GeoModels not equal: ",
                geomodel1.nb_mesh_entities( cur_entity_type ), " ",
                geomodel2.nb_mesh_entities( cur_entity_type ) );
        }
    }

    const auto& all_geological_entity_types =
        geomodel1.entity_type_manager()
            .geological_entity_manager.geological_entity_types();

    for( const auto& cur_entity_type : all_geological_entity_types )
    {
        if( geomodel1.nb_geological_entities( cur_entity_type )
            != geomodel2.nb_geological_entities( cur_entity_type ) )
        {
            throw RINGMeshException( "TEST", "Number of entities of type ",
                cur_entity_type.string(), " of two GeoModels not equal: ",
                geomodel1.nb_geological_entities( cur_entity_type ), " ",
                geomodel2.nb_geological_entities( cur_entity_type ) );
        }
    }

    const StratigraphicColumn* column1 = geomodel1.stratigraphic_column();
    const StratigraphicColumn* column2 = geomodel2.stratigraphic_column();

    std::vector< std::shared_ptr< const StratigraphicUnit > > all_units1 =
        column1->get_all_units();
    std::vector< std::shared_ptr< const StratigraphicUnit > > all_units2 =
        column2->get_all_units();

    if( column1->get_name() != column2->get_name() )
    {
        throw RINGMeshException( "TEST",
            "Names of the stratigraphic column of the two GeoModels "
            "not equal: ",
            column1->get_name(), " vs ", column2->get_name() );
    }

    if( all_units1.size() != all_units2.size() )
    {
        throw RINGMeshException( "TEST",
            "Number of units in the stratigraphic column of the two GeoModels "
            "not equal: ",
            all_units1.size(), " vs ", all_units2.size() );
    }

    bool identical = true;
    for( auto unit_index : range( all_units1.size() ) )
    {
        if( all_units1[unit_index]->is_conformable_base()
                != all_units2[unit_index]->is_conformable_base()
            || all_units1[unit_index]->is_conformable_top()
                   != all_units2[unit_index]->is_conformable_top()
            || all_units1[unit_index]->get_relation_base()
                   != all_units2[unit_index]->get_relation_base()
            || all_units1[unit_index]->get_relation_top()
                   != all_units2[unit_index]->get_relation_top() )
        {
            identical = false;
            break;
        }

        std::string feat_name_base1 =
            ( all_units1[unit_index]->get_interface_base() != nullptr )
                ? all_units1[unit_index]->get_interface_base()->name()
                : "";

        std::string feat_name_base2 =
            ( all_units2[unit_index]->get_interface_base() != nullptr )
                ? all_units2[unit_index]->get_interface_base()->name()
                : "";

        std::string feat_name_top1 =
            ( all_units1[unit_index]->get_interface_top() != nullptr )
                ? all_units1[unit_index]->get_interface_top()->name()
                : "";

        std::string feat_name_top2 =
            ( all_units2[unit_index]->get_interface_top() != nullptr )
                ? all_units2[unit_index]->get_interface_top()->name()
                : "";

        if( feat_name_base1 != feat_name_base2
            || feat_name_top1 != feat_name_top2 )
        {
            identical = false;
            break;
        }

        GMGE::GEOL_FEATURE geo_feat_base1 =
            ( all_units1[unit_index]->get_interface_base() != nullptr )
                ? all_units1[unit_index]
                      ->get_interface_base()
                      ->geological_feature()
                : GMGE::GEOL_FEATURE::NO_GEOL;

        GMGE::GEOL_FEATURE geo_feat_base2 =
            ( all_units2[unit_index]->get_interface_base() != nullptr )
                ? all_units2[unit_index]
                      ->get_interface_base()
                      ->geological_feature()
                : GMGE::GEOL_FEATURE::NO_GEOL;

        GMGE::GEOL_FEATURE geo_feat_top1 =
            ( all_units1[unit_index]->get_interface_top() != nullptr )
                ? all_units1[unit_index]
                      ->get_interface_top()
                      ->geological_feature()
                : GMGE::GEOL_FEATURE::NO_GEOL;

        GMGE::GEOL_FEATURE geo_feat_top2 =
            ( all_units2[unit_index]->get_interface_top() != nullptr )
                ? all_units2[unit_index]
                      ->get_interface_top()
                      ->geological_feature()
                : GMGE::GEOL_FEATURE::NO_GEOL;

        if( geo_feat_base1 != geo_feat_base2 || geo_feat_top1 != geo_feat_top2 )
        {
            identical = false;
            break;
        }

        if( all_units1[unit_index]->get_layer() == nullptr
            || all_units2[unit_index]->get_layer() == nullptr
            || all_units1[unit_index]->get_layer()->geological_feature()
                   != all_units2[unit_index]->get_layer()->geological_feature()
            || all_units1[unit_index]->get_layer()->name()
                   != all_units2[unit_index]->get_layer()->name() )
        {
            identical = false;
            break;
        }
    }

    if( !identical )
    {
        throw RINGMeshException( "TEST",
            "Units in the stratigraphic column of the two GeoModels "
            "not equal" );
    }
}

void do_test()
{
    GeoModel3D origin_geomodel;
    load_geomodel( origin_geomodel,
        ringmesh_test_data_path + "geomodel3d_w_mesh_w_prop_w_column.epc" );

    const std::string saved_filename(
        ringmesh_test_output_path + "geomodel_resqml2_round_trip.epc" );
    geomodel_save( origin_geomodel, saved_filename );

    GeoModel3D reloaded_geomodel;
    load_geomodel( reloaded_geomodel, saved_filename );

    compare_geomodels( origin_geomodel, reloaded_geomodel );
}

int main()
{
    try
    {
        Logger::out( "TEST", "GeoModel3D RESQML2 Round-trip" );
        do_test();
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
