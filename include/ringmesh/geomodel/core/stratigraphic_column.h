/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#pragma once

#include <ringmesh/geomodel/core/common.h>

#include <string>
#include <vector>

/*!
 * @file ringmesh/geomodel_tools/stratigraphic_column.h
 * @brief Declarations of a stratigraphic column, stratigraphic unit, rock
 * features
 * and so on.
 * @author Marie Sirvent, Pierre Anquez and Francois Bonneau
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( Interface );
    FORWARD_DECLARATION_DIMENSION_CLASS( Layer );

    ALIAS_3D( Interface );
    ALIAS_3D( Layer );
} // namespace RINGMesh

namespace RINGMesh
{
    // @todo To develop
    enum struct ROCKTYPE
    {
        NONE,
        MULTIPLE
    };

    /*!
     * @brief Manages the RockFeature, which contains a RockType and more
     * informations
     */
    class geomodel_core_api RockFeature
    {
    public:
        /*!
         * @brief Complete constructor of a RockFeature
         * @param[in] name Name of the feature
         * @param[in] type Rocktype
         */
        RockFeature( std::string name, ROCKTYPE type )
            : name_( std::move( name ) ), type_( type )
        {
        }
        /*!
         * @brief Simple constructor of RockFeature
         * @param[in] name Name of the feature
         */
        explicit RockFeature( std::string name )
            : RockFeature( std::move( name ), ROCKTYPE::NONE )
        {
        }

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

        const ROCKTYPE& get_rock_type() const
        {
            return type_;
        }

        void set_rock_type( ROCKTYPE type )
        {
            type_ = type;
        }

    private:
        std::string name_{};
        ROCKTYPE type_{ ROCKTYPE::NONE };
    };

    enum struct RELATION
    {
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
     * Each Unit has a name, two delimiting interfaces with two corresponding
     * relations, a layer, a RockFeature,
     * a minimum thickness and a maximum thickness. A StratigraphicColumn can be
     * a StratigraphicUnit.
     */
    class geomodel_core_api StratigraphicUnit
    {
        ringmesh_disable_copy_and_move( StratigraphicUnit );

    public:
        /*!
         * @brief Complete Constructor of StratigraphicUnit
         * @param[in] name Name of the unit
         * @param[in] interface_base Interface delimiting the bottom of the unit
         * @param[in] interface_top Interface delimiting the top of the unit
         * @param[in] layer Layer delimiting the volume of the unit
         * @param[in] relation_top Relation characterizing the upper contact of
         * the unit
         * @param[in] relation_base Relation characterizing the lower contact of
         * the unit
         * @param[in] rock RockFeature characterizing among others the rock type
         * associated with the unit
         * @param[in] min_thick Minimum thickness of the layer
         * @param[in] max_thick Maximum thickness of the layer
         */
        StratigraphicUnit( std::string name, RockFeature rock );

        virtual ~StratigraphicUnit() = default;

        virtual const std::string& get_name() const
        {
            return name_;
        }
        /*!
         * @param[out] out RockFeature returned after modification in the
         * function
         * @return the RockFeature of the unit,
         * @return if StratigraphicColumn, return a RockFeature with "multiple"
         * rocktype
         */
        void set_rock_feature( const RockFeature& rock_feature )
        {
            rock_ = rock_feature;
        }

        const RockFeature& get_rock_feature() const
        {
            return rock_;
        }

        virtual bool is_conformable_base() const = 0;
        virtual bool is_conformable_top() const = 0;
        virtual RELATION get_relation_base() const = 0;
        virtual RELATION get_relation_top() const = 0;
        virtual const Interface3D& get_interface_base() const = 0;
        virtual const Interface3D& get_interface_top() const = 0;
        virtual double get_min_thick() const = 0;
        virtual double get_max_thick() const = 0;

    protected:
        /*!
         * @brief Constructor of StratigraphicUnit
         */
        StratigraphicUnit();

    protected:
        std::string name_{};
        RockFeature rock_;
    };

    class geomodel_core_api UnsubdividedStratigraphicUnit
        : public StratigraphicUnit
    {
    public:
        UnsubdividedStratigraphicUnit( std::string name,
            const Interface3D& interface_base,
            const Interface3D& interface_top,
            const Layer3D& layer,
            RELATION relation_top,
            RELATION relation_base,
            RockFeature rock,
            double min_thick,
            double max_thick );

        bool is_conformable_base() const final
        {
            return ( relation_base_ == RELATION::CONFORMABLE );
        }

        bool is_conformable_top() const final
        {
            return ( relation_top_ == RELATION::CONFORMABLE );
        }

        RELATION get_relation_base() const final
        {
            return relation_base_;
        }

        RELATION get_relation_top() const final
        {
            return relation_top_;
        }

        const Interface3D& get_interface_base() const final
        {
            return *interface_base_;
        }

        const Interface3D& get_interface_top() const final
        {
            return *interface_top_;
        }

        double get_min_thick() const final
        {
            return min_thick_;
        }

        double get_max_thick() const final
        {
            return max_thick_;
        }

    private:
        const Interface3D* interface_top_;
        const Interface3D* interface_base_;
        const Layer3D* layer_;
        RELATION relation_top_;
        RELATION relation_base_;
        double min_thick_;
        double max_thick_;
    };

    class geomodel_core_api SubdividedStratigraphicUnit
        : public StratigraphicUnit
    {
    public:
        SubdividedStratigraphicUnit( std::string name,
            RockFeature rock,
            const std::vector< const StratigraphicUnit* >& sub_units )
            : StratigraphicUnit( std::move( name ), std::move( rock ) ),
              units_( sub_units )
        {
        }

        bool is_conformable_base() const final
        {
            return ( units_.back()->is_conformable_base() );
        }

        bool is_conformable_top() const final
        {
            return ( units_.front()->is_conformable_top() );
        }

        RELATION get_relation_base() const final
        {
            return units_.back()->get_relation_base();
        }

        RELATION get_relation_top() const final
        {
            return units_.front()->get_relation_top();
        }

        const Interface3D& get_interface_base() const final
        {
            return units_.back()->get_interface_base();
        }

        const Interface3D& get_interface_top() const final
        {
            return units_.front()->get_interface_top();
        }

        double get_min_thick() const final
        {
            double sum_min_thick = 0.;
            for( auto unit : units_ )
            {
                sum_min_thick += unit->get_min_thick();
            }
            return sum_min_thick;
        }

        double get_max_thick() const final
        {
            double sum_max_thick = 0.;
            for( auto unit : units_ )
            {
                sum_max_thick += unit->get_max_thick();
            }
            return sum_max_thick;
        }

    private:
        std::vector< const StratigraphicUnit* > units_{};
    };

    enum struct STRATIGRAPHIC_PARADIGM
    {
        CHRONOSTRATIGRAPHIC,
        BIOSTRATIGRAPHIC,
        LITHOSTRATIGRAPHIC,
        UNSPECIFIED
    };

    /*!
     * @brief A stratigraphic column is composed of several stratigraphic units
     */
    class geomodel_core_api StratigraphicColumn
    {
    public:
        /*!
         * @brief Complete constructor of StratigraphicColumn
         * @param[in] name Name of the stratigraphic column
         * @param[in] layers Vector of the StratigraphicUnit constituting the
         * StratigraphicColumn
         * @param[in] type Chronostratigraphic, Lithostratigraphic or
         * Biostratigraphic
         */
        StratigraphicColumn( std::string name,
            const std::vector< const StratigraphicUnit* >& units,
            STRATIGRAPHIC_PARADIGM type )
            : name_( std::move( name ) ), units_( units ), type_( type )
        {
        }

        /*!
         * @brief Simple Constructor of StratigraphicColumn
         * @param[in] name Name of the unit
         */
        StratigraphicColumn(
            std::string name, const STRATIGRAPHIC_PARADIGM type )
            : StratigraphicColumn( std::move( name ), {}, type )
        {
        }

        /*!
         * \name Stratigraphic Column edition
         * @{
         */

        /*!
         * @param[in] above Reference unit, the new unit will be added below it
         * @param[in] unit_to_add Unit you want to add to the column
         */
        void insert_unit_below( const StratigraphicUnit& above,
            const StratigraphicUnit& unit_to_add );
        /*!
         * @param[in] to_add Unit to add at the top of the column
         */
        void insert_top_unit( const StratigraphicUnit& to_add );
        /*!
         * @param[in] to_add Unit to add at the bottom of the column
         */
        void insert_base_unit( const StratigraphicUnit& to_add );

        void remove_unit( const StratigraphicUnit& unit );

        /*! @}
         * \name Get column units
         * @{
         */

        /*!
         * @return the top unit of the column
         */
        const StratigraphicUnit* get_top_unit() const
        {
            return units_.front();
        }

        /*!
         * @return the bottom unit of the column
         */
        const StratigraphicUnit* get_base_unit() const
        {
            return units_.back();
        }

        /*!
         * @param[in] unit Reference unit
         * @return the StratigraphicUnit which position in the
         * StratigraphicColumn
         * is just above the reference unit
         */
        const StratigraphicUnit* get_unit_above(
            const StratigraphicUnit& unit ) const;
        /*!
         * @param[in] unit Reference unit
         * @return the StratigraphicUnit which position in the
         * StratigraphicColumn
         * is just below the reference unit
         */
        const StratigraphicUnit* get_unit_below(
            const StratigraphicUnit& unit ) const;

        const StratigraphicUnit* get_unit( const index_t index ) const
        {
            return units_[index];
        }

        const StratigraphicUnit* get_unit( const std::string& name ) const;

        /*!
         * @return a vector of all the units of the column
         */
        const std::vector< const StratigraphicUnit* >& get_all_units() const
        {
            return units_;
        }

        /*! @}
         * \name Others
         * @{
         */

        STRATIGRAPHIC_PARADIGM get_paradigm() const
        {
            return type_;
        }

        /*!
         * @brief is_conformable_base for the Stratigraphic Column
         * @return true if the base of the last unit of the Stratigraphic Column
         * is conformable
         */
        bool is_conformable_base() const
        {
            return ( units_.back()->is_conformable_base() );
        }
        /*!
         * @brief is_conformable_top for the Stratigraphic Column
         * @return true if the top of the first unit of the Stratigraphic Column
         * is conformable
         */
        bool is_conformable_top() const
        {
            return ( units_.front()->is_conformable_top() );
        }

        /*!
         * @brief get_relation_base for the Stratigraphic Column
         * @return the relation of the base of the first unit of the
         * StratigraphicColumn
         */
        RELATION get_relation_base()
        {
            return ( units_.back()->get_relation_base() );
        }
        /*!
         * @brief get_relation_top for the Stratigraphic Column
         * @return the relation of the top of the first unit of the
         * StratigraphicColumn
         */
        RELATION get_relation_top()
        {
            return ( units_.front()->get_relation_top() );
        }

        /*!
         * @brief get_interface_base for the Stratigraphic Column
         * @return the base interface of the last unit in the Stratigraphic
         * Column
         */
        const Interface3D& get_interface_base() const
        {
            return ( units_.back()->get_interface_base() );
        }
        /*!
         * @brief get_interface_top for the Stratigraphic Column
         * @return the top interface of the first unit in the Stratigraphic
         * Column
         */
        const Interface3D& get_interface_top() const
        {
            return ( units_.front()->get_interface_top() );
        }

        /*!
         * @return the minimum thickness for the whole Column (i.e sum on units)
         */
        double get_column_min_thick() const;
        /*!
         * @return the maximum thickness for the whole Column (i.e sum on units)
         */
        double get_column_max_thick() const;

        const std::string& get_name() const
        {
            return name_;
        }

        /*! @}
         */

    private:
        /*!
         * @param[in] unit_name of the unit to find
         * @return the position of a unit in the stratigraphic column
         */
        index_t get_index( const std::string& unit_name ) const;

    private:
        std::string name_{};
        std::vector< const StratigraphicUnit* > units_{};
        STRATIGRAPHIC_PARADIGM type_{ STRATIGRAPHIC_PARADIGM::UNSPECIFIED };
    };
} // namespace RINGMesh
