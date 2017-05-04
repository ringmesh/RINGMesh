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

#include <ringmesh/geomodel/stratigraphic_column.h>

namespace RINGMesh {

    RockFeature::RockFeature( const std::string& name )
        : name_( name )
    {
        ROCKTYPE n = NONE;
        type_ = n;
    }

    RockFeature::RockFeature( const std::string& name, ROCKTYPE type )
        : name_( name ), type_( type )
    {
    }

    RockFeature::~RockFeature()
    {
    }

    const ROCKTYPE& RockFeature::get_rock_type() const
    {
        return type_;
    }

    void RockFeature::set_rock_type( ROCKTYPE type )
    {
        type_ = type;
    }

    StratigraphicUnit::StratigraphicUnit()
        :
            name_( "none" ),
            interface_top_( nullptr ),
            interface_base_( nullptr ),
            layer_( nullptr ),
            relation_top_( CONFORMABLE ),
            relation_base_( CONFORMABLE ),
            rock_( RockFeature( "none" ) ),
            min_thick_( 0 ),
            max_thick_( max_float64() )

    {
    }
    StratigraphicUnit::StratigraphicUnit(
        const std::string name,
        const GeoModelGeologicalEntity& interface_base,
        const GeoModelGeologicalEntity& interface_top,
        const GeoModelGeologicalEntity& layer,
        RELATION relation_top,
        RELATION relation_base,
        const RockFeature& rock,
        double min_thick,
        double max_thick )
        :
            name_( name ),
            interface_top_( &interface_top ),
            interface_base_( &interface_base ),
            layer_( &layer ),
            relation_top_( relation_top ),
            relation_base_( relation_base ),
            rock_( rock ),
            min_thick_( min_thick ),
            max_thick_( max_thick )

    {
    }

    void StratigraphicColumn::insert_unit_below(
        const StratigraphicUnit& above,
        const StratigraphicUnit& to_add )
    {
        index_t index = get_index( above.get_name() );
        ringmesh_assert( index != NO_ID );
        const StratigraphicUnit* ptr_add = &to_add;
        units_.insert( units_.begin() + index + 1, ptr_add );
    }

    void StratigraphicColumn::insert_top_unit( const StratigraphicUnit& to_add )
    {
        const StratigraphicUnit* ptr_add = &to_add;
        units_.insert( units_.begin(), ptr_add );
    }

    void StratigraphicColumn::insert_base_unit( const StratigraphicUnit& to_add )
    {
        const StratigraphicUnit* ptr_add = &to_add;
        units_.push_back( ptr_add );
    }

    void StratigraphicColumn::remove_unit( const StratigraphicUnit& unit )
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index != NO_ID );
        units_.erase( units_.begin() + index );
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit_above(
        const StratigraphicUnit& unit ) const
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index > 0 && index != NO_ID );
        return units_[index - 1];
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit_below(
        const StratigraphicUnit& unit ) const
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index < units_.size() - 1 && index != NO_ID );
        return units_[index + 1];
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit( const std::string& name ) const
    {
        index_t index = get_index( name );
        ringmesh_assert( index != NO_ID );
        return units_[index];
    }

    double StratigraphicColumn::get_column_min_thick() const
    {
        double sum = 0;
        for( index_t i = 0; i < units_.size(); ++i ) {
            sum += units_[i]->get_min_thick();
        }
        return sum;
    }

    double StratigraphicColumn::get_column_max_thick() const
    {
        double sum = 0;
        for( index_t i = 0; i < units_.size(); ++i ) {
            sum += units_[i]->get_max_thick();
        }
        return sum;
    }

    index_t StratigraphicColumn::get_index( const std::string& name ) const
    {
        for( index_t i = 0; i < units_.size(); ++i ) {
            if( units_[i]->get_name() == name ) {
                return i;
            }
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }
}
