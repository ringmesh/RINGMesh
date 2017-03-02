/*
 * stratigraphic_column.cpp
 *
 *  Created on: Feb 14, 2017
 *      Author: sirvent1u
 */

#include <ringmesh/geomodel/stratigraphic_column.h>

namespace RINGMesh {
////// class RockFeature/////////

    RockFeature::RockFeature( const std::string& name )
        : name_( name )
    {
        ROCKTYPE n = none ;
        type_ = n ;
    }

    RockFeature::RockFeature( const std::string& name, const ROCKTYPE& type )
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

    void RockFeature::set_rock_type( const ROCKTYPE& type )
    {
        type_ = type ;
    }

///////class StratigraphicUnit//////////////////

    StratigraphicUnit::StratigraphicUnit()
        :
            name_( "none" ),
            interface_top_( nil ),
            interface_base_( nil ),
            layer_( nil ),
            relation_top_( conformable ),
            relation_base_( conformable ),
            rock_( RockFeature( "none" ) ),
            min_thick_( 0 ),
            max_thick_( std::numeric_limits< double >::max() )

    {
    }
    StratigraphicUnit::StratigraphicUnit(
        const std::string name,
        const RINGMesh::GeoModelGeologicalEntity& interface_base,
        const RINGMesh::GeoModelGeologicalEntity& interface_top,
        const RINGMesh::GeoModelGeologicalEntity& layer,
        const RELATION& relation_top,
        const RELATION& relation_base,
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

    const RINGMesh::GeoModelGeologicalEntity& StratigraphicUnit::get_interface_base() const
    {
        return *interface_base_ ;
    }

    const RINGMesh::GeoModelGeologicalEntity& StratigraphicUnit::get_interface_top() const
    {
        return *interface_top_ ;
    }

    const RockFeature& StratigraphicUnit::get_rock_feature( RockFeature& out ) const
    {
        const StratigraphicColumn* cast_strat_col =
            dynamic_cast< const StratigraphicColumn* >( this ) ;
        if( cast_strat_col != nil ) {
            out.set_rock_type( multiple ) ;
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
        return ( relation_base_ == conformable ) ;
    }

    bool StratigraphicUnit::is_conformable_top() const
    {
        return ( relation_top_ == conformable ) ;
    }

    const RELATION& StratigraphicUnit::get_relation_top() const
    {
        return relation_top_ ;
    }

    const RELATION& StratigraphicUnit::get_relation_base() const
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
/////////////////class StratigraphicColumn////////////////

    StratigraphicColumn::StratigraphicColumn( const std::string& name )
        : name_( name )
    {
        type_ = chronostratigraphic ;
        std::vector< const StratigraphicUnit* > layers ;
        layers_ = layers ;
    }

    StratigraphicColumn::StratigraphicColumn(
        const std::string& name,
        const std::vector< const StratigraphicUnit* >& layers,
        const STRATIGRAPHIC_PARADIGM& type )
        : name_( name ), layers_( layers ), type_( type )
    {

    }

    StratigraphicColumn::~StratigraphicColumn()
    {
    }

    void StratigraphicColumn::set_paradigm( const STRATIGRAPHIC_PARADIGM& type )
    {
        type_ = type ;
    }

    index_t StratigraphicColumn::get_index( StratigraphicUnit* unit )
    {
        return get_index( unit->get_name() ) ;
    }

    index_t StratigraphicColumn::get_index( const std::string& name )
    {
        index_t index = 0 ;
        bool found = false ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_name() == name ) {
                index = i ;
                found = true ;
                return index ;
            }
        }
        ringmesh_assert( found ) ;
        return index ;
    }

    const StratigraphicUnit* StratigraphicColumn::get_unit_above(
        StratigraphicUnit* unit )
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
        StratigraphicUnit* unit )
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

    void StratigraphicColumn::remove_unit( StratigraphicUnit* unit )
    {
        index_t index = get_index( unit ) ;
        ringmesh_assert( index < layers_.size() ) ;
        layers_.erase( layers_.begin() + index ) ;
    }

    void StratigraphicColumn::insert_unit_below(
        StratigraphicUnit* above,
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
        index_t index = 0 ;
        bool found = false ;
        RockFeature out( "temporary_name" ) ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_rock_feature( out ).get_name()
                == feature.get_name() ) {
                index = i ;
                found = true ;
                return layers_[index] ;
            }
        }
        ringmesh_assert( found ) ;
        return layers_[index] ;
    }

    const StratigraphicUnit* StratigraphicColumn::find_unit_from_rock_feature_name(
        const std::string& name )
    {
        index_t index = 0 ;
        bool found = false ;
        RockFeature out( "temporary_name" ) ;
        for( index_t i = 0; i < layers_.size(); ++i ) {
            if( layers_[i]->get_rock_feature( out ).get_name() == name ) {
                index = i ;
                found = true ;
                return layers_[index] ;
            }
        }
        ringmesh_assert( found ) ;
        return layers_[index] ;
    }

    void StratigraphicColumn::get_units_between(
        StratigraphicUnit* top,
        StratigraphicUnit* base,
        std::vector< const StratigraphicUnit* > units )
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

    const STRATIGRAPHIC_PARADIGM& StratigraphicColumn::get_paradigm() const
    {
        return type_ ;
    }

    bool StratigraphicColumn::is_conformable_base()
    {
        return ( layers_.back()->is_conformable_base() ) ;
    }

    bool StratigraphicColumn::is_conformable_top()
    {
        return ( layers_[0]->is_conformable_top() ) ;
    }

    const RELATION& StratigraphicColumn::get_relation_base()
    {
        return ( layers_[layers_.size() - 1]->get_relation_base() ) ;
    }

    const RELATION& StratigraphicColumn::get_relation_top()
    {
        return ( layers_[0]->get_relation_top() ) ;
    }

    const RINGMesh::GeoModelGeologicalEntity& StratigraphicColumn::get_interface_base() const
    {
        return ( layers_[layers_.size() - 1]->get_interface_base() ) ;
    }

    const RINGMesh::GeoModelGeologicalEntity& StratigraphicColumn::get_interface_top() const
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
