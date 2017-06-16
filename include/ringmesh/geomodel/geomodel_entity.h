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

/*!
 * @file Declaration of GeoModelEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/entity_type_manager.h>
#include <ringmesh/geomodel/geomodel_indexing_types.h>

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModel;
    template< index_t DIMENSION > class UniverseAccess;
    template< index_t DIMENSION > class UniverseAccessBase;
    template< index_t DIMENSION > class GeoModelBuilderTopology;
    template< index_t DIMENSION > class GeoModelBuilderRemoval;
}

namespace RINGMesh {
    /*!
     * @brief Abstract base class describing one entity of a GeoModel
     */
    template< index_t DIMENSION >
    class GeoModelEntity {
    ringmesh_disable_copy( GeoModelEntity );
    public:

        virtual ~GeoModelEntity() = default;

        virtual bool is_on_voi() const = 0;
        virtual bool is_valid() const = 0;

        const GeoModel< DIMENSION >& geomodel() const
        {
            return geomodel_;
        }
        const std::string& name() const
        {
            return name_;
        }
        index_t index() const
        {
            return id_;
        }

    protected:
        /*!
         * @details Client code should only create GeoModelEntities through
         * GeoModelBuilderTopology class.
         *
         * @param[in] geomodel Geomodel owning the Entity to create
         * @param[in] id Index of the entity in the corresponding vector in the geomodel
         * @param[in] name Name of the entity
         * @param[in] geological_feature Geological feature of the entity, none by default.
         */
        GeoModelEntity(
            const GeoModel< DIMENSION >& geomodel,
            index_t id,
            const std::string& name = "Unnamed" )
            : geomodel_( geomodel ), name_( name ), id_( id )
        {
        }

        void copy_name( const GeoModelEntity< DIMENSION >& from )
        {
            name_ = from.name_;
        }
        virtual bool is_index_valid() const = 0;

    protected:
        /// Reference to the GeoModel owning this entity
        const GeoModel< DIMENSION >& geomodel_;
        /// Name of the entity - default is "Unnamed"
        std::string name_;

        /// Index of the entity
        index_t id_;
    };

    template< index_t DIMENSION >
    class UniverseBase: public GeoModelEntity< DIMENSION > {
    public:
        friend class UniverseAccess< DIMENSION > ;
        friend class UniverseAccessBase< DIMENSION > ;

        UniverseBase( const GeoModel< DIMENSION >& geomodel );

        static const UniverseType universe_type_name()
        {
            return UniverseType();
        }

        virtual ~UniverseBase() = default;

        virtual bool is_valid() const override;
        virtual bool is_on_voi() const override
        {
            return true;
        }
        const UniverseType type_name() const
        {
            return universe_type_name();
        }

        index_t nb_boundaries() const
        {
            return static_cast< index_t >( universe_boundaries_.size() );
        }
        gmme_id boundary_gmme( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries() );
            return universe_boundaries_[i];
        }

        virtual bool is_identification_valid() const
        {
            return true;
        }

    protected:
        virtual bool is_index_valid() const override
        {
            return true;
        }

        void copy_universe_base( const UniverseBase& from )
        {
            universe_boundaries_ = from.universe_boundaries_;
        }

    private:
        std::vector< gmme_id > universe_boundaries_;
    };

    template< index_t DIMENSION >
    class Universe final: public UniverseBase< DIMENSION > {
    };

    template< >
    class Universe< 2 > final: public UniverseBase< 2 > {
        friend class UniverseAccessBase< 2 > ;
        friend class UniverseAccess< 2 > ;
    public:
        Universe( const GeoModel< 2 >& geomodel );

    private:
        void copy_universe( const Universe& from )
        {
            copy_universe_base( from );
        }
    };

    template< >
    class Universe< 3 > final: public UniverseBase< 3 > {
        friend class UniverseAccessBase< 3 > ;
        friend class UniverseAccess< 3 > ;
    public:
        Universe( const GeoModel< 3 >& geomodel );

        bool side( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries() );
            return universe_boundary_sides_[i];
        }
    private:
        void copy_universe( const Universe& from )
        {
            copy_universe_base( from );
            universe_boundary_sides_ = from.universe_boundary_sides_;
        }

    private:
        std::vector< bool > universe_boundary_sides_;
    };

    template< index_t DIMENSION >
    class UniverseAccessBase {
    ringmesh_disable_copy( UniverseAccessBase );
        friend class GeoModelBuilderTopology< DIMENSION > ;
        friend class GeoModelBuilderRemoval< DIMENSION > ;

    private:
        ~UniverseAccessBase() = default;

        std::vector< gmme_id >& modifiable_boundaries()
        {
            return universe_.universe_boundaries_;
        }

        void copy( const Universe< DIMENSION >& from )
        {
            universe_.copy_universe( from );
        }

    protected:
        UniverseAccessBase( Universe< DIMENSION >& universe )
            : universe_( universe )
        {
        }

    protected:
        Universe< DIMENSION >& universe_;
    };

    template< index_t DIMENSION >
    class UniverseAccess final: public UniverseAccessBase< DIMENSION > {
    };

    template< >
    class UniverseAccess< 2 > final: public UniverseAccessBase< 2 > {
        friend class GeoModelBuilderTopology< 2 > ;
        friend class GeoModelBuilderRemoval< 2 > ;
    private:
        UniverseAccess( Universe< 2 >& universe )
            : UniverseAccessBase( universe )
        {
        }
    };

    template< >
    class UniverseAccess< 3 > final: public UniverseAccessBase< 3 > {
        friend class GeoModelBuilderTopology< 3 > ;
        friend class GeoModelBuilderRemoval< 3 > ;
    private:
        UniverseAccess( Universe< 3 >& universe )
            : UniverseAccessBase( universe )
        {
        }

        std::vector< bool >& modifiable_sides()
        {
            return universe_.universe_boundary_sides_;
        }
    };

} // namespace
