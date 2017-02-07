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

#ifndef __RINGMESH_ENTITY_TYPE_MANAGER__
#define __RINGMESH_ENTITY_TYPE_MANAGER__

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/algorithm.h>

#include <vector>

namespace RINGMesh {
    class GeoModel ;
    class Corner ;
    class Line ;
    class Surface ;
    class Region ;
    class EntityTypeManager ;
}

namespace RINGMesh {


    class RINGMESH_API EntityType {
    public:
        bool operator==(  const EntityType& type2 ) const
        {
            return type_ == type2.type_ ;
        }
        bool operator!=( const EntityType& type2 ) const
        {
            return type_ != type2.type_ ;
        }
        friend std::ostream& operator<<( std::ostream& os, const EntityType& in ) {
            os << in.type_  ;
            return os ;
        }

    private:
        std::string type_ ;
    protected:
        EntityType( const std::string& type )
            : type_( type )
        {
        }
        EntityType()
            : type_( default_entity_type_string() )
        {
        }
    private:
        const std::string default_entity_type_string()
        {
            return "No_entity_type" ;
        }
    } ;

    class RINGMESH_API DefaultEntityType : public EntityType {
    public:
        static DefaultEntityType& default_entity_type()
        {
            static DefaultEntityType default_entity_type ;
            return default_entity_type ;
        }
    private:
        DefaultEntityType() ;
    } ;




    class RINGMESH_API MeshEntityType: public EntityType {
    public:
        MeshEntityType( const std::string& type )
            : EntityType( type )
        {
        }
        MeshEntityType()
            : EntityType()
        {
        }

    } ;

    class RINGMESH_API GeologicalEntityType: public EntityType {
    public:
        GeologicalEntityType( const std::string& type )
            : EntityType( type )
        {
        }
        GeologicalEntityType()
            : EntityType()
        {
        }
    } ;

    /*!
     * @brief Unique identification of a GeoModelEntity in a GeoModel
     * @todo Should we change this name? Looks like index_t but does not enforce
     *       the programming guidelines [JP]
     */
    struct gme_t {
    public:
        index_t index() const
        {
            return index_ ;
        }
        const EntityType& type() const
        {
            return type_ ;
        }
        bool operator!=( const gme_t& rhs ) const
        {
            return type_ != rhs.type_ || index_ != rhs.index_ ;
        }
        bool operator==( const gme_t& rhs ) const
        {
            return type_ == rhs.type_ && index_ == rhs.index_ ;
        }
        /*!
         * @details Compare first types, then compare indices,
         *          beginning with NO_ID indices.
         * @note In a sorted vector v of gme_t one can find the first surface with
         *       std::lower_bound( v.begin(), v.end(), gme_t( SURFACE, NO_ID ) ) ;
         */

        friend std::ostream& operator<<( std::ostream& os, const gme_t& in )
        {
            os << in.type() << " " << in.index_ ;
            return os ;
        }
        bool is_defined() const
        {
            return type_ != DefaultEntityType::default_entity_type() && index_ != NO_ID ;
        }
    protected:
        gme_t()
            : type_( DefaultEntityType::default_entity_type()), index_( NO_ID )
        {
        }
        gme_t( const EntityType& entity_type, index_t index )
            : type_( entity_type ), index_( index )
        {
        }
    protected:
        EntityType type_ ;
        /*!
         * Index of the GeoModelEntity in the GeoModel
         */
        index_t index_ ;
    } ;

    struct gmge_t: public gme_t {
        gmge_t()
            : gme_t()
        {
        }
        gmge_t( const EntityType& entity_type, index_t index )
            : gme_t( entity_type, index )
        {
        }
    }
    ;

    struct gmme_t: public gme_t {
        gmme_t()
            : gme_t()
        {
        }
        gmme_t( const EntityType& entity_type, index_t index )
            : gme_t( entity_type, index )
        {
        }
    }
    ;

    class RINGMESH_API MeshEntityTypeManager {
    public:
        MeshEntityTypeManager(EntityTypeManager& type_manager) ;

        static bool is_corner( const MeshEntityType& type ) ;
        static bool is_line( const MeshEntityType& type ) ;
        static bool is_surface( const MeshEntityType& type ) ;
        static bool is_region( const MeshEntityType& type ) ;

        static const MeshEntityType& boundary_type( const MeshEntityType& type ) ;
        static const MeshEntityType& in_boundary_type( const MeshEntityType& type ) ;

        static const std::vector< MeshEntityType >& mesh_entity_types() ;
        static index_t nb_mesh_entity_types() ;


    private:
        typedef std::map< MeshEntityType, MeshEntityType > MeshEntityTypeMap ;
        struct MeshEntityTypeBoundaryMap {
            MeshEntityTypeBoundaryMap()
            {
                register_boundary< Corner, GeoModelEntity >() ;
                register_boundary< Line, Corner >() ;
                register_boundary< Surface, Line >() ;
                register_boundary< Region, Surface >() ;
            }
            template< typename TYPE, typename BOUNDARY >
            void register_boundary()
            {
                map.insert(
                    std::pair< MeshEntityType, MeshEntityType >(
                        TYPE::type_name_static(), BOUNDARY::type_name_static() ) ) ;
            }
            MeshEntityTypeMap map ;
        } ;

        struct MeshEntityTypeInBoundaryMap {
            MeshEntityTypeInBoundaryMap()
            {
                register_in_boundary< Corner, Line >() ;
                register_in_boundary< Line, Surface >() ;
                register_in_boundary< Surface, Region >() ;
                register_in_boundary< Region, GeoModelEntity >() ;
            }
            template< typename TYPE, typename IN_BOUNDARY >
            void register_in_boundary()
            {
                map.insert(
                    std::pair< MeshEntityType, MeshEntityType >(
                        TYPE::type_name_static(),
                        IN_BOUNDARY::type_name_static() ) ) ;
            }
            MeshEntityTypeMap map ;
        } ;
    private:
        const EntityTypeManager& type_manager_ ;
        const MeshEntityTypeBoundaryMap& boundary_relationships_ ;
        const MeshEntityTypeInBoundaryMap& in_boundary_relationships_ ;

    } ;

    class RINGMESH_API GeologicalTypeManager {
    public:
        GeologicalTypeManager(EntityTypeManager& type_manager) ;

        index_t nb_geological_entity_types() const ;
        const std::vector< GeologicalEntityType >& geological_entity_types() const ;
        const EntityType& geological_entity_type( index_t index ) const ;
        index_t geological_entity_type_index( const EntityType& type ) const ;
    private:
        std::vector< GeologicalEntityType > geological_entity_types_ ;
        const EntityTypeManager& type_manager_ ;
    private:
        void register_geological_entity_type(
            const GeologicalEntityType& geological_type_name )
        {
            if( find( geological_entity_types_, geological_type_name ) == NO_ID ) {
                geological_entity_types_.push_back( ( geological_type_name ) ) ;
            }
        }
    } ;
    class RINGMESH_API RelationshipManager {
    public:
        typedef std::map< GeologicalEntityType, MeshEntityType > GeologicalEntityToChild ;
        typedef std::map< MeshEntityType, std::set< GeologicalEntityType > > MeshEntityToParents ;

        RelationshipManager(EntityTypeManager& type_manager) ;
        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const ;
        index_t nb_parent_types( const MeshEntityType& child_type ) const ;
        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const ;
    private:
        MeshEntityToParents child_to_parents_ ;
        GeologicalEntityToChild parent_to_child_ ;

        const EntityTypeManager& type_manager_ ;
    private:
        void register_relationship(
            const GeologicalEntityType& parent_type_name,
            const MeshEntityType& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name ;
            child_to_parents_[child_type_name].insert( parent_type_name ) ;
        }

    } ;

    class RINGMESH_API EntityTypeManager {
    public:
        MeshEntityTypeManager mesh_entity_manager ;
        GeologicalTypeManager geological_entity_manager ;
        RelationshipManager relationship_manager ;
    } ;
}

#endif
