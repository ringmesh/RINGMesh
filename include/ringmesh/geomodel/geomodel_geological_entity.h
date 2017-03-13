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
 * @file Declaration of GeoModelGeologicalEntity and all its children classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel_indexing_types.h>
#include <ringmesh/geomodel/geomodel_entity.h>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelMeshEntity ;
}

namespace RINGMesh {
    class RINGMESH_API GeoModelGeologicalEntity: public GeoModelEntity {
    ringmesh_disable_copy( GeoModelGeologicalEntity ) ;
    public:
        friend class GeoModelGeologicalEntityAccess ;

        static void initialize() ;

        virtual ~GeoModelGeologicalEntity()
        {
        }

        const gmge_t& gmge_id() const
        {
            return gmge_id_ ;
        }
        index_t index() const
        {
            return gmge_id().index() ;
        }
        const GeologicalEntityType& entity_type() const
        {
            return gmge_id().type() ;
        }

        virtual const MeshEntityType child_type_name() const = 0 ;
        virtual bool is_on_voi() const ;
        virtual bool is_connectivity_valid() const ;
        static const GeologicalEntityType type_name_static()
        {
            return DefaultGeologicalEntityType::default_entity_type() ;
        }
        virtual const GeologicalEntityType type_name() const
        {
            return type_name_static() ;
        }
        index_t nb_children() const
        {
            return static_cast< index_t >( children_.size() ) ;
        }
        const gmme_t& child_gmme( index_t x ) const
        {
            ringmesh_assert( x < nb_children() ) ;
            return children_[x] ;
        }
        const GeoModelMeshEntity& child( index_t x ) const ;

        virtual bool is_identification_valid() const ;

    protected:
        GeoModelGeologicalEntity(
            const GeoModel& geomodel,
            index_t id = NO_ID,
            const std::string& name = "unnamed",
            GEOL_FEATURE geological_feature = NO_GEOL )
            :
                GeoModelEntity( geomodel, name, geological_feature ),
                gmge_id_( type_name(), id )
        {
        }
        virtual void copy( const GeoModelGeologicalEntity& from )
        {
            GME::copy( from ) ;
            children_ = from.children_ ;
        }

        virtual bool is_index_valid() const ;

    protected:
        /// Entities constituting this one - see child_type( TYPE )
        std::vector< gmme_t > children_ ;
        gmge_t gmge_id_ ;
    } ;

    /// @todo Review: I am still not convinced that we always have to 
    /// derive the base class to define new entities. [JP]
    class RINGMESH_API Contact: public GeoModelGeologicalEntity {
    public:
        Contact( const GeoModel& geomodel )
            : GeoModelGeologicalEntity( geomodel )
        {
            gmge_id_.type() = type_name_static() ;
        }
        virtual ~Contact()
        {
        }

        static const GeologicalEntityType type_name_static()
        {
            return GeologicalEntityType( "Contact" ) ;
        }
        virtual const GeologicalEntityType type_name() const
        {
            return type_name_static() ;
        }
        virtual const MeshEntityType child_type_name() const ;

        virtual bool is_valid() const ;
    } ;

    class RINGMESH_API Interface: public GeoModelGeologicalEntity {
    public:
        Interface( const GeoModel& geomodel )
            : GeoModelGeologicalEntity( geomodel )
        {
            gmge_id_.type() = type_name_static() ;
        }
        virtual ~Interface()
        {
        }

        static const GeologicalEntityType type_name_static()
        {
            return GeologicalEntityType( "Interface" ) ;
        }
        virtual const GeologicalEntityType type_name() const
        {
            return type_name_static() ;
        }
        virtual const MeshEntityType child_type_name() const ;

        virtual bool is_valid() const ;
    } ;

    class RINGMESH_API Layer: public GeoModelGeologicalEntity {
    public:
        Layer( const GeoModel& geomodel )
            : GeoModelGeologicalEntity( geomodel )
        {
            gmge_id_.type() = type_name_static() ;
        }
        virtual ~Layer()
        {
        }

        static const GeologicalEntityType type_name_static()
        {
            return GeologicalEntityType( "Layer" ) ;
        }
        virtual const GeologicalEntityType type_name() const
        {
            return type_name_static() ;
        }
        virtual const MeshEntityType child_type_name() const ;

        virtual bool is_valid() const ;
    } ;

    class GeoModelGeologicalEntityAccess {
    ringmesh_disable_copy( GeoModelGeologicalEntityAccess ) ;
        friend class GeoModelBuilderTopology ;
        friend class GeoModelBuilderGeology ;
        friend class GeoModelBuilderInfo ;
        friend class GeoModelBuilderRemoval ;

    private:
        GeoModelGeologicalEntityAccess( GeoModelGeologicalEntity& gmge )
            : gmge_( gmge )
        {
        }

        std::string& modifiable_name()
        {
            return gmge_.name_ ;
        }

        index_t& modifiable_index()
        {
            return gmge_.gmge_id_.index() ;
        }

        GME::GEOL_FEATURE& modifiable_geol_feature()
        {
            return gmge_.geol_feature_ ;
        }

        std::vector< gmme_t >& modifiable_children()
        {
            return gmge_.children_ ;
        }

        static GeoModelGeologicalEntity* create_geological_entity(
            const GeologicalEntityType& type,
            const GeoModel& geomodel,
            index_t index_in_geomodel ) ;

        void copy( const GeoModelGeologicalEntity& from )
        {
            gmge_.copy( from ) ;
        }

    private:
        GeoModelGeologicalEntity& gmge_ ;
    } ;

}
