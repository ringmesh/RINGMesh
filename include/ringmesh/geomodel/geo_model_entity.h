/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#ifndef __RINGMESH_GEO_MODEL_ENTITY__
#define __RINGMESH_GEO_MODEL_ENTITY__

#include <ringmesh/basic/common.h>

#include <string>
#include <vector>
#include <map>

#include <ringmesh/geomodel/geomodel_indexing_types.h>


namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {
    /*!
     * @brief Abstract base class describing one entity of a GeoModel
     */
    class RINGMESH_API GeoModelEntity {
    public:
        friend class GeoModelEditor ;

        typedef std::string EntityType ;
        
        /*!
         * @brief Geological feature types for GeoModelEntity
         * @todo Remove this enum and find a nice way to add new types at runtime [JP]
         */
        enum GEOL_FEATURE {
            /// All geological features 
            ALL_GEOL,
            /// Default value - No geological feature defined
            NO_GEOL,
            /// Stratigraphic surface - an horizon
            STRATI,
            /// Unconformity
            UNCONFORMITY,
            /// A normal fault
            NORMAL_FAULT,
            /// A reverse fault 
            REVERSE_FAULT,
            /// An unspecified fault 
            FAULT,
            /// Volume Of Interest
            VOI
        } ;       

        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type(
            const std::vector< GEOL_FEATURE >& types ) ;
        static std::string geol_name( GEOL_FEATURE ) ;
        static bool is_fault( GEOL_FEATURE T )
        {
            return T == FAULT || T == REVERSE_FAULT || T == NORMAL_FAULT ;
        }
        static bool is_stratigraphic_limit( GEOL_FEATURE T )
        {
            return T == STRATI || T == UNCONFORMITY ;
        }

        static const EntityType type_name_static() ;
        
        virtual ~GeoModelEntity() {};
        virtual const EntityType type_name() const
        {
            return type_name_static() ;
        }
        virtual bool is_on_voi() const = 0 ;
        virtual bool is_valid() const = 0 ;
    
        
        const GeoModel& model() const
        {
            return model_ ;
        }
        const std::string& name() const
        {
            return name_ ;
        }
        const gme_t& gme_id() const
        {
            return id_ ;
        }
        index_t index() const
        {
            return gme_id().index ;
        }
        const EntityType& entity_type() const
        {
            return gme_id().type ;
        }
        bool has_geological_feature() const
        {
            return geol_feature_ != NO_GEOL ;
        }
        GEOL_FEATURE geological_feature() const
        {
            return geol_feature_ ;
        }
        /*!
         * @brief Basic checks on the minimum required information 
         */
        bool is_identification_valid() const ;

    protected:
        /*!
         * @details Client code should only create GeoModelEntities through
         * GeoModelEditor derived classes.
         *
         * @param[in] model Geomodel owning the Entity to create
         * @param[in] id Index of the entity in the corresponding vector in the model
         * @param[in] name Name of the entity
         * @param[in] geological_feature Geological feature of the entity, none by default.
         */
        GeoModelEntity(
            const GeoModel& model,
            index_t id,
            const std::string& name = "Unnamed",
            GEOL_FEATURE geological_feature = NO_GEOL )
            :
                model_( model ),
                id_( type_name_static(), id ),
                name_( name ),
                geol_feature_( geological_feature )
        {}
        GeoModelEntity& operator=( const GeoModelEntity& rhs )
        {
            id_ = rhs.id_;
            name_ = rhs.name_;
            geol_feature_ = rhs.geol_feature_ ;
            return *this ;
        }
        GeoModelEntity( const GeoModelEntity& in )
            :
            model_( in.model_ ),
            id_( in.id_ ),
            name_( in.name_ ),
            geol_feature_( in.geol_feature_ )
        {}

    protected:
        /// Reference to the GeoModel owning this entity
        const GeoModel& model_ ;
        /// Unique identifier of this GeoModelEntity in the model
        gme_t id_ ;
        /// Name of the entity - default is "Unnamed"
        std::string name_ ;
        /// Geological feature of this object - default is NO_GEOL
        GEOL_FEATURE geol_feature_ ;
    } ;

    typedef GeoModelEntity GME ;

    class RINGMESH_API Universe: public GeoModelEntity {
    public:       
        friend class GeoModelEditor ;

        Universe( const GeoModel& model ) ;        
        
        static const EntityType universe_type_name()
        {
            return "Universe" ;
        }
        
        virtual ~Universe() {};        
        virtual bool is_valid() const ;
        virtual bool is_on_voi() const
        {
            return true ;
        }
        virtual const EntityType type_name() const
        {
            return universe_type_name() ;
        }
        
        index_t nb_boundaries() const
        {
            return static_cast< index_t >( boundary_surfaces_.size() ) ;
        }
        gme_t boundary_gme( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries() ) ;
            return boundary_surfaces_[i] ;
        }        
        bool side( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries() ) ;
            return boundary_surface_sides_[i] ;
        }

    protected:
        Universe( const Universe& in )
            :
            GeoModelEntity( in ),
            boundary_surfaces_( in.boundary_surfaces_ ),
            boundary_surface_sides_( in.boundary_surface_sides_ )
        {}
        Universe& operator=( const Universe& rhs )
        {
            GME::operator=(rhs);
            boundary_surfaces_ = rhs.boundary_surfaces_ ;
            boundary_surface_sides_ = rhs.boundary_surface_sides_ ;
            return *this ;
        }

    private:
        std::vector< gme_t > boundary_surfaces_ ;
        std::vector< bool > boundary_surface_sides_ ;
    };
 

} // namespace

#endif
