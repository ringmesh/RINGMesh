/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/geomodel/core/entity_type_manager.h>

#include <deque>
#include <map>
#include <set>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/pimpl_impl.h>

#include <ringmesh/geomodel/core/entity_type_manager.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

namespace RINGMesh
{
    using MeshEntityTypeMap = std::map< MeshEntityType, MeshEntityType >;

    /*!
     * @brief struct used to map the type of a Mesh Entity to the type of its
     * boundary
     * "Corner" is boundary of "Line"
     * "Line" is boundary of "Surface"
     * "Surface" is boundary of "Region"
     */
    template < index_t DIMENSION >
    struct MeshEntityTypeBoundaryMap
    {
        MeshEntityTypeBoundaryMap();
        void register_boundary(
            const MeshEntityType& type, const MeshEntityType& boundary )
        {
            map.emplace( type, boundary );
        }
        MeshEntityTypeMap map;

    private:
        void initialize_base()
        {
            register_boundary( Corner< DIMENSION >::type_name_static(),
                ForbiddenMeshEntityType::type_name_static() );
            register_boundary( Line< DIMENSION >::type_name_static(),
                Corner< DIMENSION >::type_name_static() );
            register_boundary( Surface< DIMENSION >::type_name_static(),
                Line< DIMENSION >::type_name_static() );
        }
    };

    template <>
    MeshEntityTypeBoundaryMap< 3 >::MeshEntityTypeBoundaryMap()
    {
        initialize_base();
        register_boundary(
            Region3D::type_name_static(), Surface3D::type_name_static() );
    }

    template <>
    MeshEntityTypeBoundaryMap< 2 >::MeshEntityTypeBoundaryMap()
    {
        initialize_base();
    }

    /*!
     * @brief struct used to map the type of a Mesh Entity to the type of
     * its incident mesh entity
     * "Line" is incident of "Corner"
     * "Surface" is incident of "Line"
     * "Region" is incident of "Surface"
     */
    template < index_t DIMENSION >
    struct MeshEntityTypeIncidentEntityMap
    {
        MeshEntityTypeIncidentEntityMap();
        void register_incident_entity(
            const MeshEntityType& type, const MeshEntityType& incident_entity )
        {
            map.emplace( type, incident_entity );
        }
        MeshEntityTypeMap map;

    private:
        void initialize_base()
        {
            register_incident_entity( Corner< DIMENSION >::type_name_static(),
                Line< DIMENSION >::type_name_static() );
            register_incident_entity( Line< DIMENSION >::type_name_static(),
                Surface< DIMENSION >::type_name_static() );
        }
    };

    template <>
    MeshEntityTypeIncidentEntityMap< 3 >::MeshEntityTypeIncidentEntityMap()
    {
        initialize_base();
        register_incident_entity(
            Surface3D::type_name_static(), Region3D::type_name_static() );
        register_incident_entity( Region3D::type_name_static(),
            ForbiddenMeshEntityType::type_name_static() );
    }

    template <>
    MeshEntityTypeIncidentEntityMap< 2 >::MeshEntityTypeIncidentEntityMap()
    {
        initialize_base();
        register_incident_entity( Surface2D::type_name_static(),
            ForbiddenMeshEntityType::type_name_static() );
    }

    template < index_t DIMENSION >
    class MeshEntityTypes
    {
    public:
        MeshEntityTypes()
        {
            initialize_base();
        }

        index_t size() const
        {
            return static_cast< index_t >( mesh_entity_types_.size() );
        }
        const std::vector< MeshEntityType >& container() const
        {
            return mesh_entity_types_;
        }

    private:
        void initialize_base()
        {
            mesh_entity_types_.emplace_back( corner_type_name_static() );
            mesh_entity_types_.emplace_back( line_type_name_static() );
            mesh_entity_types_.emplace_back( surface_type_name_static() );
        }

    private:
        std::vector< MeshEntityType > mesh_entity_types_;
    };

    template <>
    MeshEntityTypes< 3 >::MeshEntityTypes()
    {
        initialize_base();
        mesh_entity_types_.emplace_back( region_type_name_static() );
    }

    template < index_t DIMENSION >
    class MeshEntityTypeManagerBase< DIMENSION >::Impl
    {
    public:
        const MeshEntityType& boundary_entity_type(
            const MeshEntityType& mesh_entity_type ) const
        {
            MeshEntityTypeMap::const_iterator itr{
                boundary_relationships_.map.find( mesh_entity_type )
            };
            ringmesh_assert( itr != boundary_relationships_.map.end() );
            return itr->second;
        }

        const MeshEntityType& incident_entity_type(
            const MeshEntityType& mesh_entity_type ) const
        {
            MeshEntityTypeMap::const_iterator itr{
                incident_entity_relationships_.map.find( mesh_entity_type )
            };
            ringmesh_assert( itr != incident_entity_relationships_.map.end() );
            return itr->second;
        }

        bool is_corner( const MeshEntityType& type ) const
        {
            return type == mesh_entity_types_.container()[0];
        }

        bool is_line( const MeshEntityType& type ) const
        {
            return type == mesh_entity_types_.container()[1];
        }

        bool is_surface( const MeshEntityType& type ) const
        {
            return type == mesh_entity_types_.container()[2];
        }

        bool is_region( const MeshEntityType& type ) const
        {
            return type == mesh_entity_types_.container()[3];
        }

        bool is_valid_type( const MeshEntityType& type ) const
        {
            return find( mesh_entity_types_.container(), type ) != NO_ID;
        }

        const std::vector< MeshEntityType >& mesh_entity_types() const
        {
            return mesh_entity_types_.container();
        }

    private:
        MeshEntityTypeBoundaryMap< DIMENSION > boundary_relationships_;
        MeshEntityTypeIncidentEntityMap< DIMENSION >
            incident_entity_relationships_;
        MeshEntityTypes< DIMENSION > mesh_entity_types_;
    };

    template < index_t DIMENSION >
    MeshEntityTypeManagerBase< DIMENSION >::MeshEntityTypeManagerBase()
    {
    }

    template < index_t DIMENSION >
    MeshEntityTypeManagerBase< DIMENSION >::~MeshEntityTypeManagerBase()
    {
    }

    template < index_t DIMENSION >
    const MeshEntityType&
        MeshEntityTypeManagerBase< DIMENSION >::boundary_entity_type(
            const MeshEntityType& mesh_entity_type ) const
    {
        return impl_->boundary_entity_type( mesh_entity_type );
    }

    template < index_t DIMENSION >
    const MeshEntityType&
        MeshEntityTypeManagerBase< DIMENSION >::incident_entity_type(
            const MeshEntityType& mesh_entity_type ) const
    {
        return impl_->incident_entity_type( mesh_entity_type );
    }

    template < index_t DIMENSION >
    const std::vector< MeshEntityType >&
        MeshEntityTypeManagerBase< DIMENSION >::mesh_entity_types() const
    {
        return impl_->mesh_entity_types();
    }

    template < index_t DIMENSION >
    bool MeshEntityTypeManagerBase< DIMENSION >::is_corner(
        const MeshEntityType& type ) const
    {
        return impl_->is_corner( type );
    }

    template < index_t DIMENSION >
    bool MeshEntityTypeManagerBase< DIMENSION >::is_line(
        const MeshEntityType& type ) const
    {
        return impl_->is_line( type );
    }

    bool MeshEntityTypeManager< 3 >::is_region(
        const MeshEntityType& type ) const
    {
        return impl_->is_region( type );
    }

    template < index_t DIMENSION >
    bool MeshEntityTypeManagerBase< DIMENSION >::is_surface(
        const MeshEntityType& type ) const
    {
        return impl_->is_surface( type );
    }

    template < index_t DIMENSION >
    bool MeshEntityTypeManagerBase< DIMENSION >::is_valid_type(
        const MeshEntityType& type ) const
    {
        return impl_->is_valid_type( type );
    }

    index_t GeologicalTypeManager::nb_geological_entity_types() const
    {
        return static_cast< index_t >( geological_entity_types_.size() );
    }

    const std::vector< GeologicalEntityType >&
        GeologicalTypeManager::geological_entity_types() const
    {
        return geological_entity_types_;
    }

    const GeologicalEntityType& GeologicalTypeManager::geological_entity_type(
        index_t index ) const
    {
        return geological_entity_types_.at( index );
    }

    index_t GeologicalTypeManager::geological_entity_type_index(
        const GeologicalEntityType& type ) const
    {
        return find( geological_entity_types_, type );
    }

    bool GeologicalTypeManager::is_valid_type(
        const GeologicalEntityType& type ) const
    {
        return contains( geological_entity_types_, type );
    }

    void GeologicalTypeManager::register_geological_entity_type(
        const GeologicalEntityType& geological_type_name )
    {
        if( !contains( geological_entity_types_, geological_type_name ) )
        {
            geological_entity_types_.push_back( geological_type_name );
        }
    }

    class RelationshipManager::Impl
    {
    public:
        void copy( const Impl& from )
        {
            child_to_parents_ = from.child_to_parents_;
            parent_to_child_ = from.parent_to_child_;
            boundary_relationships_ = from.boundary_relationships_;
            parent_child_relationships_ = from.parent_child_relationships_;
        }

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

        index_t add_parent_child_relationship(
            const gmge_id& parent, const gmme_id& child )
        {
            index_t relationship_id{ static_cast< index_t >(
                parent_child_relationships_.size() ) };
            parent_child_relationships_.emplace_back( parent, child );
            return relationship_id;
        }

        void set_parent_to_parent_child_relationship(
            index_t relationship_id, const gmge_id& parent )
        {
            parent_child_relationships_[relationship_id].parent_id_ = parent;
        }

        void set_child_to_parent_child_relationship(
            index_t relationship_id, const gmme_id& child )
        {
            parent_child_relationships_[relationship_id].child_id_ = child;
        }

        index_t find_parent_child_relationship(
            const gmge_id& parent, const gmme_id& child )
        {
            return find( parent_child_relationships_,
                ParentChildRelationship( parent, child ) );
        }

        index_t add_boundary_relationship(
            const gmme_id& incident_entity, const gmme_id& boundary )
        {
            index_t relationship_id{ static_cast< index_t >(
                boundary_relationships_.size() ) };
            boundary_relationships_.emplace_back( incident_entity, boundary );
            return relationship_id;
        }

        void set_boundary_to_boundary_relationship(
            index_t relationship_id, const gmme_id& boundary )
        {
            boundary_relationships_[relationship_id].boundary_id_ = boundary;
        }

        void set_incident_entity_to_boundary_relationship(

            index_t relationship_id, const gmme_id& incident_entity )
        {
            boundary_relationships_[relationship_id].incident_entity_id_ =
                incident_entity;
        }

        index_t find_boundary_relationship(
            const gmme_id& incident_entity, const gmme_id& boundary )
        {
            return find( boundary_relationships_,
                BoundaryRelationship( incident_entity, boundary ) );
        }

        void register_geology_relationship(
            const GeologicalEntityType& parent_type_name,
            const MeshEntityType& child_type_name )
        {
            parent_to_child_[parent_type_name] = child_type_name;
            child_to_parents_[child_type_name].insert( parent_type_name );
        }

        const MeshEntityType child_type(
            const GeologicalEntityType& parent_type ) const
        {
            GeologicalEntityToChild::const_iterator itr{ parent_to_child_.find(
                parent_type ) };
            if( itr == parent_to_child_.end() )
            {
                return ForbiddenMeshEntityType::type_name_static();
            }
            return itr->second;
        }

        std::vector< GeologicalEntityType > parent_types(
            const MeshEntityType& child_type ) const
        {
            MeshEntityToParents::const_iterator itr{ child_to_parents_.find(
                child_type ) };
            std::vector< GeologicalEntityType > result;
            if( itr != child_to_parents_.end() )
            {
                result.insert(
                    result.begin(), itr->second.begin(), itr->second.end() );
            }
            return result;
        }

    private:
        struct ParentChildRelationship
        {
            ParentChildRelationship( gmge_id parent, gmme_id child )
                : parent_id_( std::move( parent ) ),
                  child_id_( std::move( child ) )
            {
            }
            bool operator==( const ParentChildRelationship& rhs ) const
            {
                return parent_id_ == rhs.parent_id_
                       && child_id_ == rhs.child_id_;
            }
            gmge_id parent_id_;
            gmme_id child_id_;
        };

        struct BoundaryRelationship
        {
            BoundaryRelationship( gmme_id incident_entity, gmme_id boundary )
                : incident_entity_id_( std::move( incident_entity ) ),
                  boundary_id_( std::move( boundary ) )
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

    private:
        using GeologicalEntityToChild =
            std::map< GeologicalEntityType, MeshEntityType >;
        using MeshEntityToParents =
            std::map< MeshEntityType, std::set< GeologicalEntityType > >;
        MeshEntityToParents child_to_parents_;
        GeologicalEntityToChild parent_to_child_;

        std::deque< BoundaryRelationship > boundary_relationships_;
        std::deque< ParentChildRelationship > parent_child_relationships_;
    };

    RelationshipManager::RelationshipManager()
    {
    }
    RelationshipManager::~RelationshipManager()
    {
    }

    std::vector< GeologicalEntityType > RelationshipManager::parent_types(
        const MeshEntityType& child_type ) const
    {
        return impl_->parent_types( child_type );
    }

    void RelationshipManager::copy( const RelationshipManager& from )
    {
        impl_->copy( *from.impl_ );
    }

    const gmge_id& RelationshipManager::parent_of_gmme(
        index_t relation_id ) const
    {
        return impl_->parent_of_gmme( relation_id );
    }

    const gmme_id& RelationshipManager::child_of_gmge(
        index_t relation_id ) const
    {
        return impl_->child_of_gmge( relation_id );
    }

    const gmme_id& RelationshipManager::boundary_gmme(
        index_t relation_id ) const
    {
        return impl_->boundary_gmme( relation_id );
    }
    const gmme_id& RelationshipManager::incident_entity_gmme(
        index_t relation_id ) const
    {
        return impl_->incident_entity_gmme( relation_id );
    }

    index_t RelationshipManager::nb_parent_types(
        const MeshEntityType& child_type ) const
    {
        return static_cast< index_t >( parent_types( child_type ).size() );
    }

    const MeshEntityType RelationshipManager::child_type(
        const GeologicalEntityType& parent_type ) const
    {
        return impl_->child_type( parent_type );
    }

    index_t RelationshipManager::add_parent_child_relationship(
        const gmge_id& parent, const gmme_id& child )
    {
        return impl_->add_parent_child_relationship( parent, child );
    }

    void RelationshipManager::set_parent_to_parent_child_relationship(
        index_t relationship_id, const gmge_id& parent )
    {
        impl_->set_parent_to_parent_child_relationship(
            relationship_id, parent );
    }

    void RelationshipManager::set_child_to_parent_child_relationship(
        index_t relationship_id, const gmme_id& child )
    {
        impl_->set_child_to_parent_child_relationship( relationship_id, child );
    }

    index_t RelationshipManager::find_boundary_relationship(
        const gmme_id& incident_entity, const gmme_id& boundary )
    {
        return impl_->find_boundary_relationship( incident_entity, boundary );
    }

    index_t RelationshipManager::find_parent_child_relationship(
        const gmge_id& parent, const gmme_id& child )
    {
        return impl_->find_parent_child_relationship( parent, child );
    }

    index_t RelationshipManager::add_boundary_relationship(
        const gmme_id& incident_entity, const gmme_id& boundary )
    {
        return impl_->add_boundary_relationship( incident_entity, boundary );
    }

    void RelationshipManager::set_boundary_to_boundary_relationship(
        index_t relationship_id, const gmme_id& boundary )
    {
        impl_->set_boundary_to_boundary_relationship(
            relationship_id, boundary );
    }

    void RelationshipManager::set_incident_entity_to_boundary_relationship(
        index_t relationship_id, const gmme_id& incident_entity )
    {
        impl_->set_incident_entity_to_boundary_relationship(
            relationship_id, incident_entity );
    }

    void RelationshipManager::register_geology_relationship(
        const GeologicalEntityType& parent_type_name,
        const MeshEntityType& child_type_name )
    {
        impl_->register_geology_relationship(
            parent_type_name, child_type_name );
    }

    template class geomodel_core_api MeshEntityTypes< 2 >;
    template class geomodel_core_api MeshEntityTypeManagerBase< 2 >;
    template class geomodel_core_api MeshEntityTypeManager< 2 >;
    template struct geomodel_core_api MeshEntityTypeIncidentEntityMap< 2 >;
    template struct geomodel_core_api MeshEntityTypeBoundaryMap< 2 >;

    template class geomodel_core_api MeshEntityTypes< 3 >;
    template class geomodel_core_api MeshEntityTypeManagerBase< 3 >;
    template struct geomodel_core_api MeshEntityTypeIncidentEntityMap< 3 >;
    template struct geomodel_core_api MeshEntityTypeBoundaryMap< 3 >;
} // namespace RINGMesh
