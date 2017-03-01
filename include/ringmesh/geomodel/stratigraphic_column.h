/*
 * stratigraphic_column.hpp
 *
 *  Created on: Feb 14, 2017
 *      Author: sirvent1u
 */

#ifndef INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_H_
#define INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_H_

#include <ringmesh/basic/common.h>

#include <string>
#include <vector>

/*!
 * @file ringmesh/stratigraphic_column.h
 * @brief Class representing a stratigraphic column
 */

namespace RINGMesh {
    class GeoModelGeologicalEntity ;
}

enum RockType {
    none,
    multiple
} ;

/*!
 * @brief Manages the RockFeature, which contains a RockType and more informations
 */
class RockFeature {
public:
    /*!
     * @brief Complete constructor of a RockFeature
     * @param[in] name Name of the feature
     * @param[in] type Rocktype
     */
    RockFeature(
        const std::string& name,
        const RockType& type) ;
    /*!
     * @brief Simple constructor of RockFeature
     * @param[in] name Name of the feature
     */
    RockFeature( const std::string& name ) ;
    /*!
     *@return name of the feature
     */
    const std::string get_name() const
    {
        return name_ ;
    }
    void set_name(const std::string& name) {name_ = name;} ;
    ~RockFeature() ;
    const RockType& get_rock_type() const ;
    void set_rock_type( const RockType& type ) ;



private:

    //RockFeature(const RockFeature&);

    std::string name_ ;
    RockType type_ ;

} ;

enum relation {
    conformable = 0,
    eroded = 10,
    truncation = 11,
    toplap = 12,
    baselap = 20,
    onlap = 21,
    downlap = 22
} ;

/*!
 * @brief Representing Stratigraphic Units
 * Each Unit has a name, two delimiting interfaces with two corresponding relations, a layer, a RockFeature,
 * a minimum thickness and a maximum thickness. A StratigraphicColumn can be a StratigraphicUnit.
 */
class StratigraphicUnit {
public:

    /*!
     * @brief Complete Constructor of StratigraphicUnit
     * @param[in] name Name of the unit
     * @param[in] interface_base Interface delimiting the bottom of the unit
     * @param[in] interface_top Interface delimiting the top of the unit
     * @param[in] layer Layer delimiting the volume of the unit
     * @param[in] relation_top Relation characterizing the upper contact of the unit
     * @param[in] relation_base Relation characterizing the lower contact of the unit
     * @param[in] rock RockFeature characterizing among others the rock type associated with the unit
     * @param[in] min_thick Minimum thickness of the layer
     * @param[in] max_thick Maximum thickness of the layer
     */
    StratigraphicUnit(
        const std::string name,
        const RINGMesh::GeoModelGeologicalEntity& interface_base, // TODO interface
        const RINGMesh::GeoModelGeologicalEntity& interface_top,
        const RINGMesh::GeoModelGeologicalEntity& layer,
        const relation& relation_top,
        const relation& relation_base,
        const RockFeature& rock,
        double min_thick,
        double max_thick ) ;
    virtual ~StratigraphicUnit() ;
    /*!
     * @return name of the unit or name of the column if the unit is a StratigraphicColumn
     */
    virtual const std::string& get_name() const ;
    /*!
     * @param[out] out RockFeature returned after modification in the function
     * @return the RockFeature of the unit, if StratigraphicColumn, return a RockFeature with "multiple" rocktype
     */
    const RockFeature& get_rock_feature( RockFeature& out ) const ;

    /*!
     * @brief Check if the relation_top is conformable or not, if StratigraphicColumn, relation_top of first unit
     */
    virtual bool is_conformable_base() const ;
    /*!
     * @brief Check if the relation_base is conformable or not, if StratigraphicColumn relation_base of last unit
     */
    virtual bool is_conformable_top() const ;

    /*!
     * @return relation_base of unit, if StratigraphicColumn, return relation_base of last unit
     */
    virtual const relation& get_relation_base() const ;
    /*!
     * @return relation_top of unit, if StratigraphicColumn, return relation_top of first unit
     */
    virtual const relation& get_relation_top() const ;

    /*!
     * @return interface_base of unit, if StratigraphicColumn, return interface_base of last unit
     */
    virtual const RINGMesh::GeoModelGeologicalEntity& get_interface_base() const ;
    /*!
     * @return interface_top of unit, if StratigraphicColumn, return interface_top of first unit
     */
    virtual const RINGMesh::GeoModelGeologicalEntity& get_interface_top() const ;

    /*!
     * @return min_thick_ of unit, if StratigraphicColumn, return sum of min_thick_ on all units
     */
    double get_min_thick() const ;
    /*!
     * @return max_thick_ of unit, if StratigraphicColumn, return sum of max_thick_ on all units
     */
    double get_max_thick() const ;


protected:
    /*!
     * @brief Constructor of StratigraphicUnit
     */
    StratigraphicUnit() ;


private:

    std::string name_ ;
    const RINGMesh::GeoModelGeologicalEntity* interface_top_ ;
    const RINGMesh::GeoModelGeologicalEntity* interface_base_ ;
    const RINGMesh::GeoModelGeologicalEntity* layer_ ;
    relation relation_top_ ;
    relation relation_base_ ;
    RockFeature rock_ ;
    double min_thick_ ;
    double max_thick_ ;
} ;

enum StratigraphicParadigm {
    chronostratigraphic, biostratigraphic, lithostratigraphic, unspecified
} ;

typedef unsigned int uint ;

/*!
 * @brief Manages the Stratigraphic Column
 */
class StratigraphicColumn : public StratigraphicUnit{
public:
    /*!
     * @brief Complete constructor of StratigraphicColumn
     * @param[in] name Name of the stratigraphic column
     * @param[in] layers Vector of the StratigraphicUnit constituting the StratigraphicColumn
     * @param[in] type Chronostratigraphic, Lithostratigraphic or Biostratigraphic
     */
    StratigraphicColumn(
        const std::string& name,
        const std::vector< const StratigraphicUnit* >& layers,
        const StratigraphicParadigm& type ) ;
    /*!
     * @brief Simple Constructor of StratigraphicColumn
     * @param[in] name Name of the unit
     */
    StratigraphicColumn( const std::string& name ) ;
    ~StratigraphicColumn() ;
    /*!
     * @param[in] type Chronostratigraphic, Lithostratigraphic or Biostratigraphic
     */
    void set_paradigm( const StratigraphicParadigm& type ) ;

    /*!
     * @param[in] unit Reference unit
     * @return the StratigraphicUnit which position in the StratigraphicColumn is just above the reference unit
     */
    const StratigraphicUnit* get_unit_above( StratigraphicUnit* unit ) ;
    /*!
     * @param[in] unit Reference unit
     * @return the StratigraphicUnit which position in the StratigraphicColumn is just below the reference unit
     */
    const StratigraphicUnit* get_unit_below( StratigraphicUnit* unit ) ;
    void remove_unit( StratigraphicUnit* unit ) ;
    /*!
     * @param[in] above Reference unit, the new unit will be added below it
     * @param[in] unit_to_add Unit you want to add to the column
     */
    void insert_unit_below(
        StratigraphicUnit* above,
        const StratigraphicUnit& unit_to_add ) ;
    /*!
     * @param[in] to_add Unit to add at the top of the column
     */
    void insert_top_unit( const StratigraphicUnit& to_add ) ;
    /*!
     * @param[in] to_add Unit to add at the bottom of the column
     */
    void insert_base_unit( const StratigraphicUnit& to_add ) ;

    /*!
     * @return true if the unit tested is actually a StratigraphicColumn
     */
    bool check_if_column(const std::string& name) ;
    /*!
     * @brief Cannot be used on a StratigraphicColumn, using check_if_column(name) first is needed
     * @param[in] name Name of the unit
     */
    const StratigraphicUnit* find_unit( const std::string& name ) ;
    /*!
     * @brief Only available for a StratigraphicColumn, using check_if_column(name) first is needed
     * @param[in] name Name of the unit
     */
    const StratigraphicColumn* find_sub_column( const std::string& name ) ;
    /*!
     * @return the top unit of the column
     */
    const StratigraphicUnit* get_top_unit() const;
    /*!
     * @return the bottom unit of the column
     */
    const StratigraphicUnit* get_base_unit() const ;

    /*!
     * @return true if the unit is a sub-column
     */
    bool check_if_column(uint index);
    /*!
     * @brief Cannot be used if the StratigraphicUnit is a sub-column, using check_if_column(index) first is needed
     */
    const StratigraphicUnit* get_unit(uint index) ;
    /*!
     * @brief Only available if the unit is a sub-column, using check_if_column(index) first is needed
     */
    const StratigraphicColumn* get_sub_column(uint index) ;
    /*!
     * @return a vector of all the units of the column
     */
    const std::vector< const StratigraphicUnit* >& get_units() const ;
    const StratigraphicParadigm& get_paradigm() const ;

    /*!
     * @brief is_conformable_base for the Stratigraphic Column
     * @return true if the base of the last unit of the Stratigraphic Column is conformable
     */
    bool is_conformable_base();
    /*!
     * @brief is_conformable_top for the Stratigraphic Column
     * @return true if the top of the first unit of the Stratigraphic Column is conformable
     */
    bool is_conformable_top();

    /*!
     * @brief get_relation_base for the Stratigraphic Column
     * @return the relation of the base of the first unit of the StratigraphicColumn
     */
    const relation& get_relation_base();
    /*!
     * @brief get_relation_top for the Stratigraphic Column
     * @return the relation of the top of the first unit of the StratigraphicColumn
     */
    const relation& get_relation_top();

    /*!
     * @brief get_interface_base for the Stratigraphic Column
     * @return the base interface of the last unit in the Stratigraphic Column
     */
    const RINGMesh::GeoModelGeologicalEntity& get_interface_base() const ;
    /*!
     * @brief get_interface_top for the Stratigraphic Column
     * @return the top interface of the first unit in the Stratigraphic Column
     */
    const RINGMesh::GeoModelGeologicalEntity& get_interface_top() const ;

    /*!
     * @brief get_min_thick for the Stratigraphic Column
     * @return the minimum thickness for the whole Column (i.e sum on layers)
     */
    double sum_min_thick() const ;
    /*!
     * @brief get_max_thick for the Stratigraphic Column
     * @return the maximum thickness for the whole Column (i.e sum on layers)
     */
     double sum_max_thick() const ;

     const std::string& get_name() const {return name_;};

private:

    StratigraphicColumn( const StratigraphicColumn& to_copy ) ;
    /*!
     * @return the position of a unit in the stratigraphic column
     */
    uint get_index( StratigraphicUnit* unit ) ;
    /*!
     *@param[in] name of the unit to find
     * @return the position of a unit in the stratigraphic column
     */
    uint get_index( const std::string& name ) ;

    /*!
     * @param[in] feature RockFeature used to find the corresponding unit in the column
     */
    const StratigraphicUnit* find_unit_from_rock_feature(
        const RockFeature& feature ) ;
    /*!
     * @param[in] name Name of the RockFeature used to find the corresponding unit in the column
     */
    const StratigraphicUnit* find_unit_from_rock_feature_name(
        const std::string& name ) ;
    /*!
     * @return a vector of the StratigraphicUnit found between the two given units in the stratigraphic column
     */
    std::vector< const StratigraphicUnit* > get_units_between(
        StratigraphicUnit* top,
        StratigraphicUnit* base ) ;

    std::string name_ ;
    std::vector< const StratigraphicUnit* > layers_ ;
    StratigraphicParadigm type_ ;

} ;

#endif /* INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_H_ */
