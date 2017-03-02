/*
 * stratigraphic_column.cpp
 *
 *  Created on: Feb 14, 2017
 *      Author: sirvent1u
 */

#include <ringmesh/geomodel/stratigraphic_column.h>
#include <assert.h>

////// class RockFeature/////////

RockFeature::RockFeature( const std::string& name )
    : name_( name )
{
    RockType n = none ;
    type_ = n ;
}

RockFeature::RockFeature(
    const std::string& name,
    const RockType& type )
    : name_( name ), type_( type )
{
}

RockFeature::~RockFeature()
{
}

const RockType& RockFeature::get_rock_type() const
{
    return type_ ;
}

void RockFeature::set_rock_type( const RockType& type )
{
    type_ = type ;
}



///////class StratigraphicUnit//////////////////

StratigraphicUnit::StratigraphicUnit()
    :
        name_( "none" ),
        interface_top_( NULL ),
        interface_base_( NULL ),
        layer_( NULL ),
        relation_top_ ( conformable ),
        relation_base_( conformable ),
        rock_( RockFeature("none") ),
        min_thick_( 0 ),
        max_thick_( 10^6 )

{
}
StratigraphicUnit::StratigraphicUnit(
    const std::string name,
    const RINGMesh::GeoModelGeologicalEntity& interface_base,
    const RINGMesh::GeoModelGeologicalEntity& interface_top,
    const RINGMesh::GeoModelGeologicalEntity& layer,
    const relation& relation_top,
    const relation& relation_base,
    const RockFeature& rock ,
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
        min_thick_ ( min_thick ),
        max_thick_ ( max_thick )

{
}

StratigraphicUnit::~StratigraphicUnit()
{

}

const std::string& StratigraphicUnit::get_name() const
{
    const StratigraphicColumn* test = dynamic_cast<const StratigraphicColumn*>(this);
    if (test != NULL)
    {
        return test->get_name();
    } //StratigraphicColumn
    else
    {
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

const RockFeature& StratigraphicUnit::get_rock_feature( RockFeature& out) const
{
    const StratigraphicColumn* test = dynamic_cast<const StratigraphicColumn*>(this);
    if (test != NULL)
    {
        out.set_rock_type(multiple);
        out.set_name(name_);
        return out;
        //RockFeature rock_column("rock_column", multiple);
    } //StratigraphicColumn
    else
    {
        out.set_rock_type(rock_.get_rock_type());
        out.set_name(rock_.get_name());
        return out;
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

const relation& StratigraphicUnit::get_relation_top() const
{
    return relation_top_ ;
}

const relation& StratigraphicUnit::get_relation_base() const
{
    return relation_base_ ;
}

double StratigraphicUnit::get_min_thick() const
{
    const StratigraphicColumn* test = dynamic_cast<const StratigraphicColumn*>(this);
    if (test != NULL)
    {
        return test->sum_min_thick();
    } //StratigraphicColumn
    else
    {
        return min_thick_ ;
    } //StratigraphicUnit
}

double StratigraphicUnit::get_max_thick() const
{
    const StratigraphicColumn* test = dynamic_cast<const StratigraphicColumn*>(this);
    if (test != NULL)
    {
        return test->sum_max_thick();
    } //StratigraphicColumn
    else
    {
        return max_thick_ ;
    } //StratigraphicUnit
}
/////////////////class StratigraphicColumn////////////////
typedef unsigned int uint ;

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
    const StratigraphicParadigm& type )
    : name_( name ), layers_(layers), type_( type )
{

}


StratigraphicColumn::~StratigraphicColumn()
{
}

void StratigraphicColumn::set_paradigm( const StratigraphicParadigm& type )
{
    type_ = type ;
}

uint StratigraphicColumn::get_index( StratigraphicUnit* unit )
{
    return get_index( unit->get_name() ) ;
}

uint StratigraphicColumn::get_index( const std::string& name )
{
    uint index = 0 ;
    bool found = false ;
    for( uint i = 0; i < layers_.size(); ++i ) {
        if( layers_[i]->get_name() == name ) {
            index = i ;
            found = true ;
        }
    }
    assert( found == true ) ;
    return index ;
}

const StratigraphicUnit* StratigraphicColumn::get_unit_above(
    StratigraphicUnit* unit )
{
    uint index = get_index( unit ) ;
    assert( index > 0 && index < layers_.size() + 1 ) ;
    const StratigraphicColumn* cast_above = dynamic_cast<const StratigraphicColumn*>(layers_[index - 1]);
    if(cast_above != NULL)
    {
        return cast_above->get_base_unit();
    }
    else
    {
        return layers_[index - 1] ;
    }

}

const StratigraphicUnit* StratigraphicColumn::get_unit_below(
    StratigraphicUnit* unit )
{
    uint index = get_index( unit ) ;
    assert( index < layers_.size() ) ;
    const StratigraphicColumn* cast_below = dynamic_cast<const StratigraphicColumn*>(layers_[index + 1]);
    if(cast_below != NULL)
    {
        return cast_below->get_top_unit();
    }
    else
    {
        return layers_[index + 1] ;
    }

}

void StratigraphicColumn::remove_unit( StratigraphicUnit* unit )
{
    uint index = get_index( unit ) ;
    assert( index < layers_.size() ) ;
    layers_.erase( layers_.begin() + index ) ;
}

void StratigraphicColumn::insert_unit_below(
    StratigraphicUnit* above,
    const StratigraphicUnit& to_add )
{
    uint index = get_index( above ) ;
    assert( index < layers_.size() ) ;
    const StratigraphicUnit* ptadd = &to_add;
    layers_.insert( layers_.begin() + index + 1, ptadd ) ;
}

void StratigraphicColumn::insert_top_unit( const StratigraphicUnit& to_add )
{
    const StratigraphicUnit* ptadd = &to_add;
    layers_.insert( layers_.begin(), ptadd ) ;
}

void StratigraphicColumn::insert_base_unit( const StratigraphicUnit& to_add )
{
    const StratigraphicUnit* ptadd = &to_add;
    layers_.push_back( ptadd ) ;
}

bool StratigraphicColumn::check_if_column(const std::string& name)
{
    uint index = get_index( name ) ;
    assert( index < layers_.size() ) ;
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index]) ;
    if(cast_sub_column == NULL)
    {
        return false;
    }
    else
    {
        return true;
    }
}

const StratigraphicUnit* StratigraphicColumn::find_unit( const std::string& name )
{
    uint index = get_index( name ) ;
    assert( index < layers_.size() ) ;
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index]) ;
    assert(cast_sub_column == NULL);
    return layers_[index] ;
}

const StratigraphicColumn* StratigraphicColumn::find_sub_column(const std::string& name)
{
    uint index = get_index( name ) ;
    assert( index < layers_.size() ) ;
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index]) ;
    assert(cast_sub_column != NULL);
    return cast_sub_column;
}

const StratigraphicUnit* StratigraphicColumn::get_top_unit() const
{
    const StratigraphicColumn* cast_top = dynamic_cast<const StratigraphicColumn*>(layers_[0]) ;
    if (cast_top != NULL)
    {
        return cast_top->get_top_unit();
    }
    else
    {
        return layers_[0] ;
    }

}

const StratigraphicUnit* StratigraphicColumn::get_base_unit() const
{
    const StratigraphicColumn* cast_base = dynamic_cast<const StratigraphicColumn*>(layers_[layers_.size() - 1] ) ;
    if(cast_base != NULL)
    {
        return cast_base->get_base_unit();
    }
    else
    {
        return layers_[layers_.size() - 1] ;
    }
}

bool StratigraphicColumn::check_if_column(uint index)
{
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index - 1] ) ;
    if(cast_sub_column == NULL)
    {
        return false;
    }
    else
    {
        return true;
    }
}

const StratigraphicUnit* StratigraphicColumn::get_unit(uint index)
{
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index - 1] ) ;
    assert(cast_sub_column == NULL);
    return layers_[index - 1];
}

const StratigraphicColumn* StratigraphicColumn::get_sub_column(uint index)
{
    const StratigraphicColumn* cast_sub_column = dynamic_cast<const StratigraphicColumn*>(layers_[index - 1] ) ;
    assert(cast_sub_column != NULL);
    return cast_sub_column;
}
const std::vector< const StratigraphicUnit* >& StratigraphicColumn::get_units() const
{
    return layers_ ;
}

const StratigraphicUnit* StratigraphicColumn::find_unit_from_rock_feature(
    const RockFeature& feature )
{
    uint index = 0 ;
    bool found = false ;
    RockFeature out("temporary_name");
    for( uint i = 0; i < layers_.size(); ++i ) {
        if( layers_[i]->get_rock_feature(out).get_name() == feature.get_name() ) {
            index = i ;
            found = true ;
        }
    }
    assert( found == true ) ;
    return layers_[index] ;
}

const StratigraphicUnit* StratigraphicColumn::find_unit_from_rock_feature_name(
    const std::string& name )
{
    uint index = 0 ;
    bool found = false ;
    RockFeature out("temporary_name");
    for( uint i = 0; i < layers_.size(); ++i ) {
        if( layers_[i]->get_rock_feature(out).get_name() == name ) {
            index = i ;
            found = true ;
        }
    }
    assert( found == true ) ;
    return layers_[index] ;
}

std::vector< const StratigraphicUnit* > StratigraphicColumn::get_units_between(
    StratigraphicUnit* top,
    StratigraphicUnit* base )
{
    std::vector< const StratigraphicUnit* > result ;
    uint index_top = get_index( top ) ;
    uint index_base = get_index( base ) ;
    assert( index_top < layers_.size() ) ;
    assert( index_base < layers_.size() ) ;
    if( index_top > index_base ) {
        uint tmp = index_base ;
        index_base = index_top ;
        index_top = tmp ;
    }
    for( uint i = index_top; i < ( index_base - index_top ); ++i ) {
        result.push_back( layers_[i] ) ;
    }
    return result ;
}

const StratigraphicParadigm& StratigraphicColumn::get_paradigm() const
{
    return type_ ;
}

bool StratigraphicColumn::is_conformable_base()
{
    return (layers_.back()->is_conformable_base() );
}

bool StratigraphicColumn::is_conformable_top()
{
    return (layers_[0]->is_conformable_top() );
}

const relation& StratigraphicColumn::get_relation_base()
{
    return ( layers_[layers_.size() - 1]->get_relation_base() );
}

const relation& StratigraphicColumn::get_relation_top()
{
    return ( layers_[0]->get_relation_top() );
}

const RINGMesh::GeoModelGeologicalEntity& StratigraphicColumn::get_interface_base() const
{
    return ( layers_[layers_.size() - 1]->get_interface_base() );
}

const RINGMesh::GeoModelGeologicalEntity& StratigraphicColumn::get_interface_top() const
{
    return ( layers_[0]->get_interface_top() );
}

double StratigraphicColumn::sum_min_thick() const
{
    double sum = 0;
    for(uint i = 0; i < layers_.size(); ++i)
    {
        sum += layers_[i]->get_min_thick();
    }
    return sum;
}

double StratigraphicColumn::sum_max_thick() const
{
    double sum = 0;
    for(uint i = 0; i < layers_.size(); ++i)
    {
        sum += layers_[i]->get_max_thick();
    }
    return sum;
}
