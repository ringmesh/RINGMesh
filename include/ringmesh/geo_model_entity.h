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

#include <ringmesh/common.h>

#include <string>
#include <vector>
#include <map>

#include <ringmesh/mesh.h>

namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {


    /*!  
     * @brief Manages the type relationship between GeoModelEntities
     * One instance owned by the GeoModel.
     */
    class RINGMESH_API EntityRelationships {
    public:
        void register_relationship( const std::string& parent_type_name,
            const std::string& child_type_name )
        {
            register_child_type( parent_type_name, child_type_name ) ;
            register_parent_type( parent_type_name, child_type_name ) ;
        }
        void register_child_type( const std::string& parent_type_name, 
            const std::string& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name ;
        }

        void register_parent_type( const std::string& parent_type_name,
            const std::string& child_type_name )
        {
            child_to_parents_[child_type_name].insert( parent_type_name ) ;
        }

        const std::set< std::string >& parent_types( const std::string& child_type ) const 
        {
            std::map< std::string, std::set< std::string > >::const_iterator
                itr = child_to_parents_.find( child_type );
            return itr->second ;
        }
       /* index_t nb_parent_types( const std::string& child_type ) const
        {
            return parent_types( child_type ).size() ;
        }
        const std::string parent_type( const std::string child_type, index_t i ) const
        {
            const std::set< std::string >& parent_types
            return *parent_types( child_type ).begin()[i] ;
        } */ 

        const std::string& child_type( const std::string& parent_type ) const
        {
           std::map< std::string, std::string >::const_iterator
                itr = parent_to_child_.find( parent_type );
           return itr->second ;
        }

    private:
        std::map< std::string, std::string > parent_to_child_ ;
        std::map< std::string, std::set< std::string > > child_to_parents_ ;
    };

    /*!
     * @brief Generic class describing one entity of a GeoModel
     * 
     */
    class RINGMESH_API GeoModelEntity {
   // ringmesh_disable_copy( GeoModelEntity ) ;
        friend class GeoModelEditor ;
    public:
        /*!
         * @brief Geological feature types for GeoModelEntity
         * @todo Read all possible geological features used in RESQML.
         */
        enum GEOL_FEATURE {
            /// All geological features 
            ALL_GEOL,
            /// Default value - No geological feature defined
            NO_GEOL,
            /// Stratigraphic surface - an horizon
            STRATI,
            /*!
             * Unconformity
             * @todo: Distinguish between as erosive, baselap or intrusive
             * unconformities --GC
             */
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
  
        /*! 
         * @brief Unique identification of a GeoModelEntity in a GeoModel
         * @details Stores the TYPE of the entity and its index in the GeoModel.
         *          Default values are NO_TYPE and NO_ID
         * @todo Should we change this name? it looks like index_t but does not enforce
         *       the programming guidelines [JP]
         */
        struct gme_t {
            gme_t()
                : type(), index( NO_ID )
            {
            }
            gme_t( const std::string& entity_type, index_t id )
                : type( entity_type ), index( id )
            {
            }
            bool operator!=( const gme_t& rhs ) const
            {
                return type != rhs.type || index != rhs.index ;
            }
            bool operator==( const gme_t& rhs ) const
            {
                return type == rhs.type && index == rhs.index ;
            }
            /*!
             * @brief Sort GME identifiers
             * @details Compare first types, then compare indices, 
             *          beginning with NO_ID indices. 
             * @note In a sorted vector v of gme_t one can find the first surface with
             *       std::lower_bound( v.begin(), v.end(), gme_t( SURFACE, NO_ID ) ) ;
             */
            bool operator<( const gme_t& rhs ) const
            {
                if( type != rhs.type ) {
                    // Probleme ? est ce que la comparaison des strings suffit ?
                    // Ordre important pour qui ? Attention
                    return type < rhs.type ;
                } else {
                    if( index == NO_ID ) return true ;
                    if( rhs.index == NO_ID ) return false ;
                    return index < rhs.index ;
                }
            }
            friend std::ostream& operator<<( std::ostream& os, const gme_t& in )
            {
                os << in.type << " " << in.index ;
                return os ;    
            }

            bool is_defined() const ;
                
            /*!
             * TYPE of the GeoModelEntity
             */
            std::string type ;
            /*!
             * Index of the entity in the GeoModel
             */
            index_t index ;
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

        static const std::string type_name_ ;
        /*!
         * \name Key functions to access relationships between TYPEs 
         * @{
         */
        virtual const std::string& type_name() const = 0 ;

        /*!@}
         */

        virtual ~GeoModelEntity()
        {
        }

        GeoModelEntity& operator=( const GeoModelEntity& rhs )
        {
            ringmesh_assert( model() == rhs.model() ) ;
            id_ = rhs.id_;
            name_ = rhs.name_;
            geol_feature_ = rhs.geol_feature_ ;
            return *this ;
        }

        /*!@}
         * \name Validity checks
         * @{
         */

        /*!
         * @brief Global validity check of the GME
         */
        virtual bool is_valid() const = 0 ;

        /*!
         * @brief Basic checks on the minimum required information 
         * @details Required connectivity information depends on the TYPE.   
         */
        bool is_connectivity_valid() const ;

        /*!
         * \name Accessors to basic information
         * @{
         */

        const GeoModel& model() const
        {
            return model_ ;
        }
        bool has_name() const
        {
            return !name().empty() ;
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
        /*TYPE type() const
        {
            return gme_id().type ;
        }*/
        bool has_geological_feature() const
        {
            return geol_feature_ != NO_GEOL ;
        }
        GEOL_FEATURE geological_feature() const
        {
            return geol_feature_ ;
        }
        virtual bool is_on_voi() const = 0 ;

    protected:
        /*!
         * @brief Constructs a GeoModelEntity
         * Client code should only create GeoModelEntities through
         * GeoModelEditor derived classes.
         *
         * @param[in] model Constant reference to the parent model of this entity.
         * @param[in] id Index of the entity in the corresponding vector in the model
         * @param[in] name Name of the entity, empty by default.
         * @param[in] geological_feature Feature of the entity, none by default.
         */
        GeoModelEntity(
            const GeoModel& model,
            index_t id,
            const std::string& name = "",
            GEOL_FEATURE geological_feature = NO_GEOL )
            :
                model_( model ),
                id_( type_name_, id ),
                name_( name ),
                geol_feature_( geological_feature )
        {
        }

    protected:
        /// Reference to the GeoModel owning this entity
        const GeoModel& model_ ;

        /// Unique identifier of the GeoModelEntity in the model
        gme_t id_ ;

        /// Name of the entity - default is an empty string
        std::string name_ ;

        /// Geological feature of this object - default is NO_GEOL
        GEOL_FEATURE geol_feature_ ;
    } ;

    typedef GeoModelEntity GME ;

  

    class RINGMESH_API Universe : public GeoModelEntity
    {
    public:
        friend GeoModelEditor ;

        const std::string type_name_ = "Universe" ;
        Universe( const GeoModel& model ) ;
        virtual ~Universe() {};
        
        Universe& operator=(const Universe& rhs)
        {
            GME::operator=(rhs);
            boundary_surfaces_ = rhs.boundary_surfaces_ ;
            boundary_surface_sides_ = rhs.boundary_surface_sides_ ;
            return *this ;
        }

        bool is_valid() const
        {
            return false ;
        }
        bool is_on_voi() const
        {
            return true ;
        }
        const std::string& type_name() const
        {
            return type_name_ ;
        }
        index_t nb_boundaries() const
        {
            return static_cast< index_t >( boundary_surfaces_.size() ) ;
        }
        gme_t boundary_gme( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries ) ;
            return boundary_surfaces_[i] ;
        }        
        bool side( index_t i ) const
        {
            ringmesh_assert( i < nb_boundaries ) ;
            return boundary_surface_sides_[i] ;
        }

    private:
        std::vector< gme_t > boundary_surfaces_ ;
        std::vector< bool > boundary_surface_sides_ ;
    };
 

} // namespace

#endif
