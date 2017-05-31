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
#include <ringmesh/basic/algorithm.h>
#include <ringmesh/geomodel/entity_type.h>

#include <deque>
#include <vector>

namespace RINGMesh {
    class GeoModel;
    template< index_t DIMENSION > class GeoModelMeshEntity;
    template< index_t DIMENSION > class Corner;
    template< index_t DIMENSION > class Line;
    template< index_t DIMENSION > class Surface;
    template< index_t DIMENSION > class Region;
    class EntityTypeManager;
}

namespace RINGMesh {

    /*!
     * @brief struct used to map the type of a Mesh Entity to the type of its boundary
     * "Corner" is boundary of "Line"
     * "Line" is boundary of "Surface"
     * "Surface" is boundary of "Region"
     */
    using MeshEntityTypeMap = std::map< MeshEntityType, MeshEntityType >;
    struct MeshEntityTypeBoundaryMap {
        MeshEntityTypeBoundaryMap();
        void register_boundary(
            const MeshEntityType& type,
            const MeshEntityType& boundary )
        {
            map.emplace( type, boundary );
        }
        MeshEntityTypeMap map;
    };

    /*!
     * @brief struct used to map the type of a Mesh Entity to the type of its incident entities
     * "Line" is incident entity of "Corner"
     * "Surface" is incident entity of "Line"
     * "Region" is incident entity of "Surface"
     */
    struct MeshEntityTypeIncidentEntityMap {
        MeshEntityTypeIncidentEntityMap();
        void register_incident_entity(
            const MeshEntityType& type,
            const MeshEntityType& incident_entity )
        {
            map.emplace( type, incident_entity );
        }
        MeshEntityTypeMap map;
    };

    /*!
     * @brief this class contains only static methods to manage the type of the
     * GeoModelMeshEntity. It gives access to the number of meshed entities of each
     * type and also their (in) boundary
     */
    class RINGMESH_API MeshEntityTypeManager {
    public:
        MeshEntityTypeManager();

        static bool is_corner( const MeshEntityType& type );
        static bool is_line( const MeshEntityType& type );
        static bool is_surface( const MeshEntityType& type );
        static bool is_region( const MeshEntityType& type );
        static bool is_valid_type( const MeshEntityType& type );

        static const MeshEntityType& boundary_type( const MeshEntityType& type );
        static const MeshEntityType& incident_entity_type( const MeshEntityType& type );

        static const std::vector< MeshEntityType >& mesh_entity_types();
        static index_t nb_mesh_entity_types();

    private:
        static MeshEntityTypeBoundaryMap boundary_relationships_;
        static MeshEntityTypeIncidentEntityMap incident_entity_relationships_;

    };

    /*!
     * @brief this class contains methods to manage the type of the
     * GeoModelGeologicalEntity. It gives access to the number of geological entities of each
     * type and also give the opportunity to create and manage new one.
     */
    class RINGMESH_API GeologicalTypeManager {
        friend class GeoModelBuilderGeology;
    public:
        index_t nb_geological_entity_types() const;
        const std::vector< GeologicalEntityType >& geological_entity_types() const;
        const GeologicalEntityType& geological_entity_type( index_t index ) const;
        index_t geological_entity_type_index(
            const GeologicalEntityType& type ) const;
        bool is_valid_type( const GeologicalEntityType& type ) const;

    private:
        std::vector< GeologicalEntityType > geological_entity_types_;
    private:
        void register_geological_entity_type(
            const GeologicalEntityType& geological_type_name )
        {
            if( find( geological_entity_types_, geological_type_name ) == NO_ID ) {
                geological_entity_types_.push_back( ( geological_type_name ) );
            }
        }
    };

    /*!
     * @brief this class contains methods to manage relations between Geological and
     * Mesh entities. For instance:
     * A "Contact" can be the parent of one or more "Line"
     * An "Interface" can the parent of one or more "Surface"
     * A "Layer" can be the parent of one or more "Region"
     *
     */
    class RINGMESH_API RelationshipManager {
        friend class GeoModelBuilderGeology;
        friend class GeoModelBuilderTopology;
    public:
        using GeologicalEntityToChild = std::map< GeologicalEntityType, MeshEntityType >;
        using MeshEntityToParents = std::map< MeshEntityType, std::set< GeologicalEntityType > >;

        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const;
        index_t nb_parent_types( const MeshEntityType& child_type ) const;
        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const;

        /*! @}
         * \name Access to the Boundary/Incident entity relations
         * @{
         */

        const gmme_id& boundary_gmme( index_t relation_id ) const
        {
            ringmesh_assert( relation_id < boundary_relationships_.size() );
            return boundary_relationships_[relation_id].boundary_id_;
        }
        const gmme_id& incident_entity_gmme( index_t relation_id ) const
        {
            ringmesh_assert( relation_id < boundary_relationships_.size() );
            return boundary_relationships_[relation_id].incident_entity_id_;
        }

        /*! @}
         * \name Access to the Parent / Child relations
         * @{
         */

        const gmge_id& parent_of_gmme( index_t relation_id ) const
        {
            ringmesh_assert( relation_id < parent_child_relationships_.size() );
            return parent_child_relationships_[relation_id].parent_id_;
        }
        const gmme_id& child_of_gmge( index_t relation_id ) const
        {
            ringmesh_assert( relation_id < parent_child_relationships_.size() );
            return parent_child_relationships_[relation_id].child_id_;
        }

    private:
        void register_geology_relationship(
            const GeologicalEntityType& parent_type_name,
            const MeshEntityType& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name;
            child_to_parents_[child_type_name].insert( parent_type_name );
        }

        /*! @}
         * \name Boundary Relationship manager
         * @{
         */

        index_t add_boundary_relationship(
            const gmme_id& incident_entity,
            const gmme_id& boundary )
        {
            index_t relationship_id =
                static_cast< index_t >( boundary_relationships_.size() );
            boundary_relationships_.emplace_back( incident_entity, boundary );
            return relationship_id;
        }

        index_t find_boundary_relationship(
            const gmme_id& incident_entity,
            const gmme_id& boundary )
        {
            return find( boundary_relationships_,
                BoundaryRelationship( incident_entity, boundary ) );
        }

        void set_boundary_to_boundary_relationship(
            index_t relationship_id,
            const gmme_id& boundary )
        {
            boundary_relationships_[relationship_id].boundary_id_ = boundary;
        }

        void set_incident_entity_to_boundary_relationship(

            index_t relationship_id,
            const gmme_id& incident_entity )
        {
            boundary_relationships_[relationship_id].incident_entity_id_ = incident_entity;
        }

        struct BoundaryRelationship {
            BoundaryRelationship(
                const gmme_id& incident_entity,
                const gmme_id& boundary )
                : incident_entity_id_( incident_entity ), boundary_id_( boundary )
            {
            }
            bool operator==( const BoundaryRelationship& rhs ) const
            {
                return incident_entity_id_ == rhs.incident_entity_id_
                    && boundary_id_ == rhs.boundary_id_;
            }
            gmme_id incident_entity_id_;
            gmme_id boundary_id_;
        };

        /*! @}
         * \name Parent/Child Relationship manager
         * @{
         */

        index_t add_parent_child_relationship(
            const gmge_id& parent,
            const gmme_id& child )
        {
            index_t relationship_id =
                static_cast< index_t >( parent_child_relationships_.size() );
            parent_child_relationships_.emplace_back( parent, child );
            return relationship_id;
        }

        index_t find_parent_child_relationship(
            const gmge_id& parent,
            const gmme_id& child )
        {
            return find( parent_child_relationships_,
                ParentChildRelationship( parent, child ) );
        }

        void set_parent_to_parent_child_relationship(
            index_t relationship_id,
            const gmge_id& parent )
        {
            parent_child_relationships_[relationship_id].parent_id_ = parent;
        }

        void set_child_to_parent_child_relationship(
            index_t relationship_id,
            const gmme_id& child )
        {
            parent_child_relationships_[relationship_id].child_id_ = child;
        }

        struct ParentChildRelationship {
            ParentChildRelationship( const gmge_id& parent, const gmme_id& child )
                : parent_id_( parent ), child_id_( child )
            {
            }
            bool operator==( const ParentChildRelationship& rhs ) const
            {
                return parent_id_ == rhs.parent_id_ && child_id_ == rhs.child_id_;
            }
            gmge_id parent_id_;
            gmme_id child_id_;
        };
    private:
        MeshEntityToParents child_to_parents_;
        GeologicalEntityToChild parent_to_child_;

        std::deque< BoundaryRelationship > boundary_relationships_;
        std::deque< ParentChildRelationship > parent_child_relationships_;

    };

    /*!
     * @brief Global entity manager which coulb be associated to a geomodel
     * to give access to different manager to handle the entity types
     */
    class RINGMESH_API EntityTypeManager {
    public:
        MeshEntityTypeManager mesh_entity_manager;
        GeologicalTypeManager geological_entity_manager;
        RelationshipManager relationship_manager;
    };
}
