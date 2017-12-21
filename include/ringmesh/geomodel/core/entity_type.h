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

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntityAccess );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntityAccess );
} // namespace RINGMesh

namespace RINGMesh
{
    /*
     * @brief Abstract class defining a Geomodel Entity Type
     * This class encapsulate a string which contains the name of the entity
     * type
     * It contains useful operator to compare and display the type
     * It is possible to do cast of an EntityType -> string
     */
    class geomodel_core_api EntityType
    {
    public:
        bool operator==( const EntityType& type2 ) const
        {
            return type_ == type2.type_;
        }
        bool operator!=( const EntityType& type2 ) const
        {
            return type_ != type2.type_;
        }
        friend std::ostream& operator<<(
            std::ostream& os, const EntityType& in )
        {
            os << in.type_;
            return os;
        }
        bool operator<( const EntityType& rhs ) const
        {
            return type_ < rhs.type_;
        }

        const std::string& string() const
        {
            return type_;
        }

    private:
        std::string type_{};

    protected:
        explicit EntityType( std::string type ) : type_( std::move( type ) ) {}
        EntityType() : EntityType( default_entity_type_string() ) {}

        void set_type( std::string type )
        {
            type_ = std::move( type );
        }

    private:
        std::string default_entity_type_string()
        {
            return "No_entity_type";
        }
    };

    /*!
     * @brief The MeshEntityType described the type of the meshed entities
     * There are 4 MeshEntityTypes corresponding to the 4 GeoModelMeshEntities:
     * Corner,
     * Line,
     * Surface,
     * Region
     */
    class geomodel_core_api MeshEntityType : public EntityType
    {
    public:
        explicit MeshEntityType( std::string type )
            : EntityType( std::move( type ) )
        {
        }
        MeshEntityType() = default;
    };

    /*!
     * @brief The GeologicalEntityType described the type of the Geological
     * entities
     * User can defined there own GeologicalEntityType even if there are some
     * already
     * defined (see geomodel_geological_entities.h
     * Contact,
     * Interface,
     * Layer
     */
    class geomodel_core_api GeologicalEntityType : public EntityType
    {
    public:
        explicit GeologicalEntityType( std::string type )
            : EntityType( std::move( type ) )
        {
        }
        GeologicalEntityType() = default;
    };

    /*!
     * @brief this is the MeshEntityType defined by default.
     * It is mainly used to test the validity of a created MeshEntityType
     */
    class geomodel_core_api ForbiddenMeshEntityType : public MeshEntityType
    {
    public:
        static ForbiddenMeshEntityType& type_name_static()
        {
            static ForbiddenMeshEntityType entity_type;
            return entity_type;
        }

    private:
        ForbiddenMeshEntityType() = default;
    };

    /*!
     * @brief this is the GeologicalEntityType defined by default.
     * It is mainly used to test the validity of a created GeologicalEntityType
     */
    class geomodel_core_api ForbiddenGeologicalEntityType
        : public GeologicalEntityType
    {
    public:
        static ForbiddenGeologicalEntityType& type_name_static()
        {
            static ForbiddenGeologicalEntityType entity_type;
            return entity_type;
        }

    private:
        ForbiddenGeologicalEntityType() = default;
    };

    /*!
     * @brief Unique identification of a GeoModelEntity in a GeoModel
     * It contains the EntityType and the index of the entity.
     * It is widely used in the code to easily access/modify/set a
     * GeoModelEntity
     */
    template < class Entity_type_template >
    struct gme_id
    {
        friend GeoModelMeshEntityAccess< 2 >;
        friend GeoModelGeologicalEntityAccess< 2 >;
        friend GeoModelMeshEntityAccess< 3 >;
        friend GeoModelGeologicalEntityAccess< 3 >;

    public:
        index_t index() const
        {
            return index_;
        }

        Entity_type_template type() const
        {
            return type_;
        }

        bool operator!=( const gme_id& rhs ) const
        {
            return type_ != rhs.type_ || index_ != rhs.index_;
        }

        bool operator==( const gme_id& rhs ) const
        {
            return type_ == rhs.type_ && index_ == rhs.index_;
        }

        friend std::ostream& operator<<( std::ostream& os, const gme_id& in )
        {
            os << in.type_ << " " << in.index_;
            return os;
        }

        bool operator<( const gme_id& rhs ) const
        {
            if( type_ != rhs.type_ )
            {
                /// @warning Is this now enough for EntityType = std::string?
                /// Did any code relied on that sorting? Maybe mine ... [JP]
                return type_ < rhs.type_;
            }
            if( index_ == NO_ID )
            {
                return true;
            }
            if( rhs.index_ == NO_ID )
            {
                return false;
            }
            return index_ < rhs.index_;
        }

    protected:
        gme_id() = default;

        gme_id( Entity_type_template entity_type, index_t index )
            : type_( std::move( entity_type ) ), index_( index )
        {
        }

    protected:
        Entity_type_template type_;
        /*!
         * Index of the GeoModelEntity in the GeoModel
         */
        index_t index_{ NO_ID };
    };

    /*!
     * @brief This template is a specialization of a gme_id to the
     * GeoModelGeologicalEntity
     */
    struct gmge_id : public gme_id< GeologicalEntityType >
    {
    public:
        gmge_id()
        {
            type_ = ForbiddenGeologicalEntityType::type_name_static();
        }

        gmge_id( GeologicalEntityType entity_type, index_t index )
            : gme_id< GeologicalEntityType >( std::move( entity_type ), index )
        {
        }

        bool is_defined() const
        {
            return type_ != ForbiddenGeologicalEntityType::type_name_static()
                   && index_ != NO_ID;
        }
    };
    /*!
     * @brief This template is a specialization of a gme_id to the
     * GeoModelMeshEntity
     */
    struct gmme_id : public gme_id< MeshEntityType >
    {
    public:
        gmme_id()
            : gme_id< MeshEntityType >(
                  ForbiddenMeshEntityType::type_name_static(), NO_ID )
        {
        }

        gmme_id( MeshEntityType entity_type, index_t index )
            : gme_id< MeshEntityType >( std::move( entity_type ), index )
        {
        }

        bool is_defined() const
        {
            return type_ != ForbiddenMeshEntityType::type_name_static()
                   && index_ != NO_ID;
        }
    };
} // namespace RINGMesh
