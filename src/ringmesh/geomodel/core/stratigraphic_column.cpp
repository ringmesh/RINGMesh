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

#include <ringmesh/geomodel/core/stratigraphic_column.h>

namespace RINGMesh
{
    StratigraphicUnit::StratigraphicUnit()
        : name_( "none" ), rock_( RockFeature( "none" ) ), layer_( nullptr )
    {
    }

    StratigraphicUnit::StratigraphicUnit(
        std::string name, RockFeature rock, const Layer3D* layer )
        : name_( std::move( name ) ),
          rock_( std::move( rock ) ),
          layer_( layer )
    {
    }

    UnsubdividedStratigraphicUnit::UnsubdividedStratigraphicUnit(
        std::string name,
        const Interface3D* interface_base,
        const Interface3D* interface_top,
        const Layer3D& layer,
        RELATION relation_top,
        RELATION relation_base,
        RockFeature rock,
        double min_thick,
        double max_thick )
        : StratigraphicUnit( std::move( name ), std::move( rock ), &layer ),
          interface_top_( interface_top ),
          interface_base_( interface_base ),
          relation_top_( relation_top ),
          relation_base_( relation_base ),
          min_thick_( min_thick ),
          max_thick_( max_thick )
    {
    }

    void NestedStratigraphicUnit::insert_unit_below(
        const StratigraphicUnit& above,
        std::shared_ptr< const StratigraphicUnit > unit_to_add )
    {
        index_t index = get_index( above.get_name() );
        ringmesh_assert( index != NO_ID );
        units_.insert( units_.begin() + index + 1, unit_to_add );
    }

    void NestedStratigraphicUnit::insert_top_unit(
        std::shared_ptr< const StratigraphicUnit > to_add )
    {
        units_.insert( units_.begin(), to_add );
    }

    void NestedStratigraphicUnit::insert_base_unit(
        std::shared_ptr< const StratigraphicUnit > to_add )
    {
        units_.push_back( to_add );
    }

    void NestedStratigraphicUnit::remove_unit( const StratigraphicUnit& unit )
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index != NO_ID );
        units_.erase( units_.begin() + index );
    }

    const StratigraphicUnit* NestedStratigraphicUnit::get_unit_above(
        const StratigraphicUnit& unit ) const
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index > 0 && index != NO_ID );
        return units_[index - 1].get();
    }

    const StratigraphicUnit* NestedStratigraphicUnit::get_unit_below(
        const StratigraphicUnit& unit ) const
    {
        index_t index = get_index( unit.get_name() );
        ringmesh_assert( index < units_.size() - 1 && index != NO_ID );
        return units_[index + 1].get();
    }

    const StratigraphicUnit* NestedStratigraphicUnit::get_unit(
        const std::string& name ) const
    {
        index_t index = get_index( name );
        ringmesh_assert( index != NO_ID );
        return units_[index].get();
    }

    index_t NestedStratigraphicUnit::get_index(
        const std::string& unit_name ) const
    {
        for( auto i : range( units_.size() ) )
        {
            if( units_[i]->get_name() == unit_name )
            {
                return i;
            }
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    unsigned int NestedStratigraphicUnit::get_max_rank() const
    {
        unsigned int max_rank = 0;
        for( const auto& unit : units_ )
        {
            const NestedStratigraphicUnit* nested_unit =
                dynamic_cast< const NestedStratigraphicUnit* >( unit.get() );
            if( nested_unit == nullptr )
            {
                continue;
            }
            ringmesh_assert( nested_unit->get_units().size() > 1 );

            std::max( max_rank, nested_unit->get_max_rank() );
        }
        return ( max_rank > 0 ) ? ( 1 + max_rank ) : 0;
    }

    /*!
     * @return a vector of the ranked units.
     */
    NestedStratigraphicUnit::RankedUnits
        NestedStratigraphicUnit::get_units_with_rank( unsigned int rank ) const
    {
        RankedUnits units;

        const unsigned int max_rank = get_max_rank();
        if( rank > max_rank )
        {
            rank = max_rank;
        }

        for( const auto& unit : units_ )
        {
            const NestedStratigraphicUnit* nested_unit =
                dynamic_cast< const NestedStratigraphicUnit* >( unit.get() );
            if( nested_unit == nullptr || rank == 0 )
            {
                units.push_back( unit.get() );
            }
            else
            {
                RankedUnits nested_units =
                    nested_unit->get_units_with_rank( rank - 1 );
                units.insert(
                    units.end(), nested_units.begin(), nested_units.end() );
            }
        }

        return units;
    }
} // namespace RINGMesh
