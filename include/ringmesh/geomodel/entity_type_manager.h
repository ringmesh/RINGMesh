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
#include <ringmesh/geomodel/entity_type.h>

#include <vector>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelMeshEntity ;
    class Corner ;
    class Line ;
    class Surface ;
    class Region ;
    class EntityTypeManager ;
}

namespace RINGMesh {

    class RINGMESH_API MeshEntityTypeManager {
    public:
        MeshEntityTypeManager() ;

        static bool is_corner( const MeshEntityType& type ) ;
        static bool is_line( const MeshEntityType& type ) ;
        static bool is_surface( const MeshEntityType& type ) ;
        static bool is_region( const MeshEntityType& type ) ;
        static bool is_valid_type( const MeshEntityType& type ) ;

        static const MeshEntityType& boundary_type( const MeshEntityType& type ) ;
        static const MeshEntityType& in_boundary_type( const MeshEntityType& type ) ;

        static const std::vector< MeshEntityType >& mesh_entity_types() ;
        static index_t nb_mesh_entity_types() ;

    private:
        typedef std::map< MeshEntityType, MeshEntityType > MeshEntityTypeMap ;
        struct MeshEntityTypeBoundaryMap {
            MeshEntityTypeBoundaryMap() ;
            void register_boundary(
                const MeshEntityType& type,
                const MeshEntityType& boundary )
            {
                map.insert(
                    std::pair< MeshEntityType, MeshEntityType >( type, boundary ) ) ;
            }
            MeshEntityTypeMap map ;
        } ;

        struct MeshEntityTypeInBoundaryMap {
            MeshEntityTypeInBoundaryMap() ;
            void register_in_boundary(
                const MeshEntityType& type,
                const MeshEntityType& in_boundary )
            {
                map.insert(
                    std::pair< MeshEntityType, MeshEntityType >( type,
                        in_boundary ) ) ;
            }
            MeshEntityTypeMap map ;
        } ;
    private:
        static MeshEntityTypeBoundaryMap boundary_relationships_ ;
        static MeshEntityTypeInBoundaryMap in_boundary_relationships_ ;

    } ;

    class RINGMESH_API GeologicalTypeManager {
        friend class GeoModelBuilderGeology ;
    public:
        index_t nb_geological_entity_types() const ;
        const std::vector< GeologicalEntityType >& geological_entity_types() const ;
        const GeologicalEntityType& geological_entity_type( index_t index ) const ;
        index_t geological_entity_type_index(
            const GeologicalEntityType& type ) const ;
    private:
        std::vector< GeologicalEntityType > geological_entity_types_ ;
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
        friend class GeoModelBuilderGeology ;
    public:
        typedef std::map< GeologicalEntityType, MeshEntityType > GeologicalEntityToChild ;
        typedef std::map< MeshEntityType, std::set< GeologicalEntityType > > MeshEntityToParents ;

        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const ;
        index_t nb_parent_types( const MeshEntityType& child_type ) const ;
        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const ;
    private:
        MeshEntityToParents child_to_parents_ ;
        GeologicalEntityToChild parent_to_child_ ;
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
