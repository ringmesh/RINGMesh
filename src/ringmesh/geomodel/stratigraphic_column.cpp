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
        ROCKTYPE n = NONE ;
        type_ = n ;
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
        return type_ ;
    }

    void RockFeature::set_rock_type( ROCKTYPE type )
    {
        type_ = type ;
    }


    StratigraphicUnit::StratigraphicUnit()
        :
            name_( "none" ),
            interface_top_( nil ),
            interface_base_( nil ),
            layer_( nil ),
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

    StratigraphicUnit::~StratigraphicUnit()
    {

    }

    const std::string& StratigraphicUnit::get_name() const
    {
        const StratigraphicColumn* cast_strat_col =
            dynamic_cast< const StratigraphicColumn* >( this ) ;
        if( cast_strat_col != nil ) {
            return cast_strat_col->get_name() ;
        } //StratigraphicColumn
        else {
            return name_ ;
        } //StratigraphicUnit

    }

    const GeoModelGeologicalEntity& StratigraphicUnit::get_interface_base() const
    {
        return *interface_base_ ;
    }

    const GeoModelGeologicalEntity& StratigraphicUnit::get_interface_top() const
    {
        return *interface_top_ ;
    }

    const RockFeature& StratigraphicUnit::get_rock_feature( RockFeature& out ) const
    {
        const StratigraphicColumn* cast_strat_col =
            dynamic_cast< const StratigraphicColumn* >( this ) ;
        if( cast_strat_col != nil ) {
            out.set_rock_type( MULTIPLE ) ;
            out.set_name( name_ ) ;
            return out ;
            //RockFeature rock_column("rock_column", multiple);
        } //StratigraphicColumn
        else {
            out.set_rock_type( rock_.get_rock_type() ) ;
            out.set_name( rock_.get_name() ) ;
            return out ;
        } //StratigraphicUnit

    }

    bool StratigraphicUnit::is_conformable_base() const
    {
        return ( relation_base_ == CONFORMABLE ) ;
    }

    bool StratigraphicUnit::is_conformable_top() const
    {
        return ( relation_top_ == CONFORMABLE ) ;
    }

    RELATION StratigraphicUnit::get_relation_top() const
    {
        return relation_top_ ;
    }

    RELATION StratigraphicUnit::get_relation_base() const
    {
        return relation_base_ ;
    }

    double StratigraphicUnit::get_min_thick() const
    {
        const StratigraphicColumn* cast_strat_col =
            dynamic_cast< const StratigraphicColumn* >( this ) ;
        if( cast_strat_col != nil ) {
            return cast_strat_col->sum_min_thick() ;
        } //StratigraphicColumn
        else {
            return min_thick_ ;
        } //StratigraphicUnit
    }

    double StratigraphicUnit::get_max_thick() const
    {
        const StratigraphicColumn* cast_strat_col =
            dynamic_cast< const StratigraphicColumn* >( this ) ;
        if( cast_strat_col != nil ) {
            return cast_strat_col->sum_max_thick() ;
        } //StratigraphicColumn
        else {
            return max_thick_ ;
        } //StratigraphicUnit
    }

    StratigraphicColumn::StratigraphicColumn( const std::string& name )
        : name_( name ),
          layers_(),
          type_(CHRONOSTRATIGRAPHIC)
    {

    }

    StratigraphicColumn::StratigraphicColumn(
        const std::string& name,
        const std::vector< const StratigraphicUnit* >& layers,
        STRATIGRAPHIC_PARADIGM type )
        : name_( name ), layers_( layers ), type_( type )
    {

    }

    StratigraphicColumn::~StratigraphicColumn()
    {
    }

    void StratigraphicColumn::set_paradigm( STRATIGRAPHIC_PARADIGM type )
    {
        type_ = type ;
    }

    index_t StratigraphicColumn::get_index( const StratigraphicUnit& unit ) const
    {
        return get_index( unit.get_name() ) ;
    }

    index_t StratigraphicColumn::get_index( const std::string& name ) const
    {
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_name() == name ) {
                return i;
            }
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit_above(
        const StratigraphicUnit& unit )
    {
        index_t index = get_index( unit ) ;
        ringmesh_assert( index > 0 && index < layers_.size() + 1 ) ;
        const StratigraphicColumn* cast_above =
            dynamic_cast< const StratigraphicColumn* >( layers_[index - 1] ) ;
        if( cast_above != nil ) {
            return cast_above->get_base_unit() ;
        } else {
            return layers_[index - 1] ;
        }
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit_below(
        const StratigraphicUnit& unit )
    {
        index_t index = get_index( unit ) ;
        ringmesh_assert( index < layers_.size() ) ;
        const StratigraphicColumn* cast_below =
            dynamic_cast< const StratigraphicColumn* >( layers_[index + 1] ) ;
        if( cast_below != nil ) {
            return cast_below->get_top_unit() ;
        } else {
            return layers_[index + 1] ;
        }

    }

    void StratigraphicColumn::remove_unit( const StratigraphicUnit& unit )
    {
        index_t index = get_index( unit ) ;
        ringmesh_assert( index < layers_.size() ) ;
        layers_.erase( layers_.begin() + index ) ;
    }

    void StratigraphicColumn::insert_unit_below(
        const StratigraphicUnit& above,
        const StratigraphicUnit& to_add )
    {
        index_t index = get_index( above ) ;
        ringmesh_assert( index < layers_.size() ) ;
        const StratigraphicUnit* ptadd = &to_add ;
        layers_.insert( layers_.begin() + index + 1, ptadd ) ;
    }

    void StratigraphicColumn::insert_top_unit( const StratigraphicUnit& to_add )
    {
        const StratigraphicUnit* ptadd = &to_add ;
        layers_.insert( layers_.begin(), ptadd ) ;
    }

    void StratigraphicColumn::insert_base_unit( const StratigraphicUnit& to_add )
    {
        const StratigraphicUnit* ptadd = &to_add ;
        layers_.push_back( ptadd ) ;
    }

    bool StratigraphicColumn::check_if_column( const std::string& name )
    {
        index_t index = get_index( name ) ;
        ringmesh_assert( index < layers_.size() ) ;
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index] ) ;
        return ( cast_sub_column != nil ) ;
    }

    const StratigraphicUnit* StratigraphicColumn::find_unit(
        const std::string& name )
    {
        index_t index = get_index( name ) ;
        ringmesh_assert( index < layers_.size() ) ;
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index] ) ;
        ringmesh_assert( cast_sub_column == nil ) ;
        return layers_[index] ;
    }

    const StratigraphicColumn* StratigraphicColumn::find_sub_column(
        const std::string& name )
    {
        index_t index = get_index( name ) ;
        ringmesh_assert( index < layers_.size() ) ;
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index] ) ;
        ringmesh_assert( cast_sub_column != nil ) ;
        return cast_sub_column ;
    }

    const StratigraphicUnit* StratigraphicColumn::get_top_unit() const
    {
        const StratigraphicColumn* cast_top =
            dynamic_cast< const StratigraphicColumn* >( layers_[0] ) ;
        if( cast_top != nil ) {
            return cast_top->get_top_unit() ;
        } else {
            return layers_[0] ;
        }

    }

    const StratigraphicUnit* StratigraphicColumn::get_base_unit() const
    {
        const StratigraphicColumn* cast_base =
            dynamic_cast< const StratigraphicColumn* >( layers_[layers_.size() - 1] ) ;
        if( cast_base != nil ) {
            return cast_base->get_base_unit() ;
        } else {
            return layers_[layers_.size() - 1] ;
        }
    }

    bool StratigraphicColumn::check_if_column( index_t index )
    {
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index - 1] ) ;
        if( cast_sub_column == nil ) {
            return false ;
        } else {
            return true ;
        }
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit( index_t index )
    {
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index - 1] ) ;
        ringmesh_assert( cast_sub_column == nil ) ;
        return layers_[index - 1] ;
    }

    const StratigraphicColumn* StratigraphicColumn::get_sub_column( index_t index )
    {
        const StratigraphicColumn* cast_sub_column =
            dynamic_cast< const StratigraphicColumn* >( layers_[index - 1] ) ;
        ringmesh_assert( cast_sub_column != nil ) ;
        return cast_sub_column ;
    }
    const std::vector< const StratigraphicUnit* >& StratigraphicColumn::get_units() const
    {
        return layers_ ;
    }

    const StratigraphicUnit* StratigraphicColumn::find_unit_from_rock_feature(
        const RockFeature& feature )
    {

        RockFeature out( "temporary_name" ) ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_rock_feature( out ).get_name()
                == feature.get_name() ) {
                return layers_[i] ;
            }
        }
        ringmesh_assert_not_reached ;
        return nil;
    }

    const StratigraphicUnit* StratigraphicColumn::find_unit_from_rock_feature_name(
        const std::string& name )
    {
        RockFeature out( "temporary_name" ) ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_rock_feature( out ).get_name() == name ) {
                return layers_[i] ;
            }
        }
        ringmesh_assert_not_reached ;
        return nil ;
    }

    void StratigraphicColumn::get_units_between(
        const StratigraphicUnit& top,
        const StratigraphicUnit& base,
        std::vector< const StratigraphicUnit* >& units )
    {
        ringmesh_assert( units.empty() ) ;
        index_t index_top = get_index( top ) ;
        index_t index_base = get_index( base ) ;
        ringmesh_assert( index_top < layers_.size() ) ;
        ringmesh_assert( index_base < layers_.size() ) ;
        if( index_top > index_base ) {
            index_t tmp = index_base ;
            index_base = index_top ;
            index_top = tmp ;
        }
        for( index_t i = index_top; i < ( index_base - index_top ); ++i ) {
            units.push_back( layers_[i] ) ;
        }

    }

    STRATIGRAPHIC_PARADIGM StratigraphicColumn::get_paradigm() const
    {
        return type_ ;
    }

    bool StratigraphicColumn::is_conformable_base() const
    {
        return ( layers_.back()->is_conformable_base() ) ;
    }

    bool StratigraphicColumn::is_conformable_top() const
    {
        return ( layers_[0]->is_conformable_top() ) ;
    }

    RELATION StratigraphicColumn::get_relation_base()
    {
        return ( layers_[layers_.size() - 1]->get_relation_base() ) ;
    }

    RELATION StratigraphicColumn::get_relation_top()
    {
        return ( layers_[0]->get_relation_top() ) ;
    }

    const GeoModelGeologicalEntity& StratigraphicColumn::get_interface_base() const
    {
        return ( layers_[layers_.size() - 1]->get_interface_base() ) ;
    }

    const GeoModelGeologicalEntity& StratigraphicColumn::get_interface_top() const
    {
        return ( layers_[0]->get_interface_top() ) ;
    }

    double StratigraphicColumn::sum_min_thick() const
    {
        double sum = 0 ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            sum += layers_[i]->get_min_thick() ;
        }
        return sum ;
    }

    double StratigraphicColumn::sum_max_thick() const
    {
        double sum = 0 ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            sum += layers_[i]->get_max_thick() ;
        }
        return sum ;
    }
}
