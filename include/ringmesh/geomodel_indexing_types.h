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


#ifndef __RINGMESH_GEOMODEL_INDEXING_TYPES__
#define __RINGMESH_GEOMODEL_INDEXING_TYPES__

#include <string>

/*!
 * @brief Structures and classes used to index elements in a GeoModel,
 * in the meshes of its entities etc.
 * @note Goal: decrease dependancies between header files. 
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    typedef std::string EntityType ;

    /*! 
     * @brief Unique identification of a GeoModelEntity in a GeoModel
     * @todo Should we change this name? Looks like index_t but does not enforce
     *       the programming guidelines [JP]
     */
    struct gme_t {
        gme_t()
            : type( "No_entity_type" ), index( NO_ID )
            // Still not perfect  "No_entity_type" is also typed in in geo_model_entity.cpp
        {}
        gme_t( const EntityType& entity_type, index_t id )
            : type( entity_type ), index( id )
        {}
        bool operator!=(const gme_t& rhs) const
        {
            return type != rhs.type || index != rhs.index ;
        }
        bool operator==(const gme_t& rhs) const
        {
            return type == rhs.type && index == rhs.index ;
        }
        /*!
         * @details Compare first types, then compare indices,
         *          beginning with NO_ID indices.
         * @note In a sorted vector v of gme_t one can find the first surface with
         *       std::lower_bound( v.begin(), v.end(), gme_t( SURFACE, NO_ID ) ) ;
         */
        bool operator<(const gme_t& rhs) const
        {
            if( type != rhs.type ) {
                /// @warning Is this now enough for EntityType = std::string?  
                /// Did any code relied on that sorting? Maybe mine ... [JP]
                return type < rhs.type ;
            } else {
                if( index == NO_ID ) return true ;
                if( rhs.index == NO_ID ) return false ;
                return index < rhs.index ;
            }
        }
        friend std::ostream& operator<<(std::ostream& os, const gme_t& in)
        {
            os << in.type << " " << in.index ;
            return os ;
        }
        bool is_defined() const
        {
            /// @todo hard encoded default name to remove 
            return type != "No_entity_type" && index != NO_ID ;
        }

        EntityType type ;
        /*!
         * Index of the GeoModelEntity in the GeoModel
         */
        index_t index ;
    } ;

    /*!
     * @brief Vertex in a GeoModelEntity
     */
    struct GMEVertex {
        GMEVertex( gme_t t, index_t vertex_id_in )
            : gme_id( t ), v_id( vertex_id_in )
        {}
        GMEVertex()
            : gme_id(), v_id( NO_ID )
        {}
        bool operator<(const GMEVertex& rhs) const
        {
            if( gme_id != rhs.gme_id ) {
                return gme_id < rhs.gme_id ;
            } else {
                return v_id < rhs.v_id ;
            }
        }
        bool operator==(const GMEVertex& rhs) const
        {
            return gme_id == rhs.gme_id && v_id == rhs.v_id ;
        }
        bool is_defined() const
        {
            return gme_id.is_defined() && v_id != NO_ID ;
        }
        /// GeoModelEntity index in the GeoModel that owns it
        gme_t gme_id ;
        /// Index of the vertex in the GeoModelEntity
        index_t v_id ;
    } ;
}

#endif 