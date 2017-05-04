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

#pragma once

#include <ringmesh/basic/common.h>

#include <string>
#include <vector>

/*!
 * @file ringmesh/stratigraphic_column.h
 * @brief Class representing a stratigraphic column
 */

namespace RINGMesh {
    class GeoModelGeologicalEntity;

    enum ROCKTYPE {
        NONE, MULTIPLE
    };

    /*!
     * @brief Manages the RockFeature, which contains a RockType and more informations
     */
    class RINGMESH_API RockFeature {
    public:
        /*!
         * @brief Complete constructor of a RockFeature
         * @param[in] name Name of the feature
         * @param[in] type Rocktype
         */
        RockFeature( const std::string& name, ROCKTYPE type );
        /*!
         * @brief Simple constructor of RockFeature
         * @param[in] name Name of the feature
         */
        RockFeature( const std::string& name );
        /*!
         *@return name of the feature
         */
        const std::string& get_name() const
        {
            return name_;
        }
        void set_name( const std::string& name )
        {
            name_ = name;
        }
        ~RockFeature();
        const ROCKTYPE& get_rock_type() const;
        void set_rock_type( ROCKTYPE type );
    private:
        std::string name_;
        ROCKTYPE type_;

    };

    enum RELATION {
        CONFORMABLE = 0,
        ERODED = 10,
        TRUNCATION = 11,
        TOPLAP = 12,
        BASELAP = 20,
        ONLAP = 21,
        DOWNLAP = 22
    };

    /*!
     * @brief Representing Stratigraphic Units
     * Each Unit has a name, two delimiting interfaces with two corresponding relations, a layer, a RockFeature,
     * a minimum thickness and a maximum thickness. A StratigraphicColumn can be a StratigraphicUnit.
     */
    class RINGMESH_API StratigraphicUnit {
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
            const GeoModelGeologicalEntity& interface_base, // TODO interface
            const GeoModelGeologicalEntity& interface_top,
            const GeoModelGeologicalEntity& layer,
            RELATION relation_top,
            RELATION relation_base,
            const RockFeature& rock,
            double min_thick,
            double max_thick );
        virtual ~StratigraphicUnit();
        /*!
         * @return name of the unit or name of the column if the unit is a StratigraphicColumn
         */
        virtual const std::string& get_name() const;
        /*!
         * @param[out] out RockFeature returned after modification in the function
         * @return the RockFeature of the unit,
         * @return if StratigraphicColumn, return a RockFeature with "multiple" rocktype
         */
        const RockFeature& get_rock_feature( RockFeature& out ) const;

        /*!
         * @brief Check if the relation_top is conformable or not, if StratigraphicColumn, relation_top of first unit
         */
        virtual bool is_conformable_base() const;
        /*!
         * @brief Check if the relation_base is conformable or not, if StratigraphicColumn relation_base of last unit
         */
        virtual bool is_conformable_top() const;

        /*!
         * @return relation_base of unit, if StratigraphicColumn, return relation_base of last unit
         */
        virtual RELATION get_relation_base() const;
        /*!
         * @return relation_top of unit, if StratigraphicColumn, return relation_top of first unit
         */
        virtual RELATION get_relation_top() const;

        /*!
         * @return interface_base of unit, if StratigraphicColumn, return interface_base of last unit
         */
        virtual const GeoModelGeologicalEntity& get_interface_base() const;
        /*!
         * @return interface_top of unit, if StratigraphicColumn, return interface_top of first unit
         */
        virtual const GeoModelGeologicalEntity& get_interface_top() const;

        /*!
         * @return min_thick_ of unit, if StratigraphicColumn, return sum of min_thick_ on all units
         */
        double get_min_thick() const;
        /*!
         * @return max_thick_ of unit, if StratigraphicColumn, return sum of max_thick_ on all units
         */
        double get_max_thick() const;

    protected:
        /*!
         * @brief Constructor of StratigraphicUnit
         */
        StratigraphicUnit();

    private:

        std::string name_;
        const GeoModelGeologicalEntity* interface_top_;
        const GeoModelGeologicalEntity* interface_base_;
        const GeoModelGeologicalEntity* layer_;
        RELATION relation_top_;
        RELATION relation_base_;
        RockFeature rock_;
        double min_thick_;
        double max_thick_;
    };

    enum STRATIGRAPHIC_PARADIGM {
        CHRONOSTRATIGRAPHIC, BIOSTRATIGRAPHIC, LITHOSTRATIGRAPHIC, UNSPECIFIED
    };

    /*!
     * @brief Manages the Stratigraphic Column
     */
    class RINGMESH_API StratigraphicColumn: public StratigraphicUnit {
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
            STRATIGRAPHIC_PARADIGM type );
        /*!
         * @brief Simple Constructor of StratigraphicColumn
         * @param[in] name Name of the unit
         */
        StratigraphicColumn( const std::string& name );
        virtual ~StratigraphicColumn();
        /*!
         * @param[in] type Chronostratigraphic, Lithostratigraphic or Biostratigraphic
         */
        void set_paradigm( STRATIGRAPHIC_PARADIGM type );

        /*!
         * @param[in] unit Reference unit
         * @return the StratigraphicUnit which position in the StratigraphicColumn is just above the reference unit
         */
        const StratigraphicUnit* get_unit_above( const StratigraphicUnit& unit );
        /*!
         * @param[in] unit Reference unit
         * @return the StratigraphicUnit which position in the StratigraphicColumn is just below the reference unit
         */
        const StratigraphicUnit* get_unit_below( const StratigraphicUnit& unit );
        void remove_unit( const StratigraphicUnit& unit );
        /*!
         * @param[in] above Reference unit, the new unit will be added below it
         * @param[in] unit_to_add Unit you want to add to the column
         */
        void insert_unit_below(
            const StratigraphicUnit& above,
            const StratigraphicUnit& unit_to_add );
        /*!
         * @param[in] to_add Unit to add at the top of the column
         */
        void insert_top_unit( const StratigraphicUnit& to_add );
        /*!
         * @param[in] to_add Unit to add at the bottom of the column
         */
        void insert_base_unit( const StratigraphicUnit& to_add );

        /*!
         * @return true if the unit tested is actually a StratigraphicColumn
         */
        bool check_if_column( const std::string& name );
        /*!
         * @brief Cannot be used on a StratigraphicColumn, using check_if_column(name) first is needed
         * @param[in] name Name of the unit
         */
        const StratigraphicUnit* find_unit( const std::string& name );
        /*!
         * @brief Only available for a StratigraphicColumn, using check_if_column(name) first is needed
         * @param[in] name Name of the unit
         */
        const StratigraphicColumn* find_sub_column( const std::string& name );
        /*!
         * @return the top unit of the column
         */
        const StratigraphicUnit* get_top_unit() const;
        /*!
         * @return the bottom unit of the column
         */
        const StratigraphicUnit* get_base_unit() const;

        /*!
         * @return true if the unit is a sub-column
         */
        bool check_if_column( index_t index );
        /*!
         * @brief Cannot be used if the StratigraphicUnit is a sub-column, using check_if_column(index) first is needed
         */
        const StratigraphicUnit* get_unit( index_t index );
        /*!
         * @brief Only available if the unit is a sub-column, using check_if_column(index) first is needed
         */
        const StratigraphicColumn* get_sub_column( index_t index );
        /*!
         * @return a vector of all the units of the column
         */
        const std::vector< const StratigraphicUnit* >& get_units() const;
        STRATIGRAPHIC_PARADIGM get_paradigm() const;

        /*!
         * @brief is_conformable_base for the Stratigraphic Column
         * @return true if the base of the last unit of the Stratigraphic Column is conformable
         */
        virtual bool is_conformable_base() const;
        /*!
         * @brief is_conformable_top for the Stratigraphic Column
         * @return true if the top of the first unit of the Stratigraphic Column is conformable
         */
        virtual bool is_conformable_top() const;

        /*!
         * @brief get_relation_base for the Stratigraphic Column
         * @return the relation of the base of the first unit of the StratigraphicColumn
         */
        virtual RELATION get_relation_base();
        /*!
         * @brief get_relation_top for the Stratigraphic Column
         * @return the relation of the top of the first unit of the StratigraphicColumn
         */
        virtual RELATION get_relation_top();

        /*!
         * @brief get_interface_base for the Stratigraphic Column
         * @return the base interface of the last unit in the Stratigraphic Column
         */
        virtual const GeoModelGeologicalEntity& get_interface_base() const;
        /*!
         * @brief get_interface_top for the Stratigraphic Column
         * @return the top interface of the first unit in the Stratigraphic Column
         */
        virtual const GeoModelGeologicalEntity& get_interface_top() const;

        /*!
         * @brief get_min_thick for the Stratigraphic Column
         * @return the minimum thickness for the whole Column (i.e sum on layers)
         */
        double sum_min_thick() const;
        /*!
         * @brief get_max_thick for the Stratigraphic Column
         * @return the maximum thickness for the whole Column (i.e sum on layers)
         */
        double sum_max_thick() const;

        virtual const std::string& get_name() const
        {
            return name_;
        }
        ;

    private:

        /*!
         * @return the position of a unit in the stratigraphic column
         */
        index_t get_index( const StratigraphicUnit& unit ) const;
        /*!
         *@param[in] name of the unit to find
         * @return the position of a unit in the stratigraphic column
         */
        index_t get_index( const std::string& name ) const;

        /*!
         * @param[in] feature RockFeature used to find the corresponding unit in the column
         */
        const StratigraphicUnit* find_unit_from_rock_feature(
            const RockFeature& feature );
        /*!
         * @param[in] name Name of the RockFeature used to find the corresponding unit in the column
         */
        const StratigraphicUnit* find_unit_from_rock_feature_name(
            const std::string& name );
        /*!
         * @param[out] units Vector of StratigraphicUnit localized between top and base
         */
        void get_units_between(
            const StratigraphicUnit& top,
            const StratigraphicUnit& base,
            std::vector< const StratigraphicUnit* >& units );

        std::string name_;
        std::vector< const StratigraphicUnit* > layers_;
        STRATIGRAPHIC_PARADIGM type_;

    };
}
