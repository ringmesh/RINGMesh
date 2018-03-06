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

#include <stack>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/task_handler.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/builder/geomodel_builder_geology.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

/*!
 * @file ringmesh/geomodel/builder/geomodel_builder_geology.cpp
 * @brief Implementation of the classes to build GeoModel geological entities
 * @author Jeanne Pellerin
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel, index_t geomodel_point_id )
    {
        const auto& geomodel_vertices = geomodel.mesh.vertices;
        const auto& vertices =
            geomodel_vertices.gme_vertices( geomodel_point_id );
        for( const auto& vertex : vertices )
        {
            if( vertex.gmme.type() == Corner< DIMENSION >::type_name_static() )
            {
                return vertex.gmme;
            }
        }
        return gmme_id{};
    }

    /*!
     * @brief Returns true if the Line< DIMENSION > has exactly the given
     * vertices
     * @todo Reimplement using std::iterators
     */
    template < index_t DIMENSION >
    bool line_equal( const Line< DIMENSION >& line,
        const std::vector< index_t >& rhs_vertices )
    {
        if( line.nb_vertices() != rhs_vertices.size() )
        {
            return false;
        }
        const auto& geomodel_vertices = line.geomodel().mesh.vertices;
        bool equal = true;
        for( auto i : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id( line.gmme(), i ) )
            {
                equal = false;
                break;
            }
        }
        if( equal )
        {
            return true;
        }
        // If the order is the other one
        equal = true;
        for( auto i : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[i]
                != geomodel_vertices.geomodel_vertex_id(
                       line.gmme(), line.nb_vertices() - i - 1 ) )
            {
                equal = false;
                break;
            }
        }
        return equal;
    }

    /*!
     * @brief Reorders the line so that front() is a corner.
     * @note Closed line has front()==back().
     */
    template < index_t DIMENSION >
    void reorder_closed_line_vertices_to_start_at_corner(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< index_t >& line_vertices )
    {
        if( geomodel.nb_corners() == 0 )
        {
            // Maybe should throw an assertion, but I am not sure [JP]
            // this really may happen for sure, so no throw [RM]
            return;
        }
        if( line_vertices.empty() )
        {
            return;
        }
        for( auto i : range( 1, line_vertices.size() - 1 ) )
        {
            auto corner = find_corner( geomodel, line_vertices[i] );
            if( corner.is_defined() )
            {
                line_vertices.pop_back();
                std::rotate( line_vertices.begin(), line_vertices.begin() + i,
                    line_vertices.end() );
                line_vertices.push_back( line_vertices.front() );
                break;
            }
        }
    }

} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBuilderGeology< DIMENSION >::GeoModelBuilderGeology(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : builder_( builder ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::copy_geology(
        const GeoModel< DIMENSION >& from )
    {
        for( auto t : range( from.nb_geological_entity_types() ) )
        {
            builder_.geology.copy_geological_entity_topology(
                from, from.geological_entity_type( t ) );
        }
    }

    template < index_t DIMENSION >
    bool GeoModelBuilderGeology< DIMENSION >::create_geological_entities(
        const GeologicalEntityType& type, index_t nb_additional_entities )
    {
        find_or_create_geological_entity_type( type );
        auto& store = geomodel_access_.modifiable_geological_entities( type );
        auto old_size = static_cast< index_t >( store.size() );
        auto new_size = old_size + nb_additional_entities;
        store.reserve( new_size );
        for( auto i : range( old_size, new_size ) )
        {
            store.emplace_back( GeoModelGeologicalEntityAccess<
                DIMENSION >::create_geological_entity( type, geomodel_, i ) );
        }
        return true;
    }

    template < index_t DIMENSION >
    bool GeoModelBuilderGeology< DIMENSION >::
        check_if_children_to_parent_relation_already_exists(
            const gmge_id& parent, const gmme_id& children )
    {
        const auto& children_mesh_entity = geomodel_.mesh_entity( children );
        if( children_mesh_entity.has_parent( parent.type() ) )
        {
            ringmesh_assert(
                parent == children_mesh_entity.parent_gmge( parent.type() ) );
            return true;
        }
        return false;
    }

    template < index_t DIMENSION >
    bool GeoModelBuilderGeology< DIMENSION >::
        check_if_parent_to_children_relation_already_exists(
            const gmge_id& parent, const gmme_id& children )
    {
        const auto& parent_geological_entity =
            geomodel_.geological_entity( parent );

        for( auto children_id :
            range( parent_geological_entity.nb_children() ) )
        {
            if( parent_geological_entity.child_gmme( children_id ) == children )
            {
                return true;
            }
        }
        return false;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::set_mesh_entity_parent(
        const gmme_id& child_gmme, index_t id, const gmge_id& parent_gmge, bool update_parent )
    {
        /// No check on the validity of the index of the entity parents_
        /// NO_ID is used to flag entities to delete
        auto& mesh_entity =
            geomodel_access_.modifiable_mesh_entity( child_gmme );
        ringmesh_assert( id < mesh_entity.nb_parents() );
        GeoModelMeshEntityAccess< DIMENSION > gmme_access( mesh_entity );
        auto relationship_id = gmme_access.modifiable_parents()[id];
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        const auto old_parent_gmge = mesh_entity.parent_gmge( id );
        manager.set_parent_to_parent_child_relationship(
            relationship_id, parent_gmge );
        if( update_parent ) {
            // Add child in new parent relationships
            GeoModelGeologicalEntityAccess< DIMENSION > parent_access{
                geomodel_access_.modifiable_geological_entity( parent_gmge )
            };
            parent_access.modifiable_children().push_back( relationship_id );
            // Remove child in old parent relationships
            GeoModelGeologicalEntityAccess< DIMENSION > old_parent_access{
                geomodel_access_.modifiable_geological_entity( old_parent_gmge )
            };
            auto relationship_position = std::find(
                old_parent_access.modifiable_children().begin(),
                old_parent_access.modifiable_children().end(), relationship_id );
            old_parent_access.modifiable_children().erase( relationship_position );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::add_parent_children_relation(
        const gmge_id& parent, const gmme_id& children )
    {
        auto& parent_entity =
            geomodel_access_.modifiable_geological_entity( parent );
        const auto& parent_entity_types =
            geomodel_.entity_type_manager().relationship_manager.parent_types(
                children.type() );
        if( !contains( parent_entity_types, parent.type() ) )
        {
            throw RINGMeshException( "Entity",
                "Wrong Parent type in the parent children relation between ",
                parent, " and ", children );
        }

        if( check_if_children_to_parent_relation_already_exists(
                parent, children ) )
        {
            if( check_if_parent_to_children_relation_already_exists(
                    parent, children ) )
            {
                return;
            }
            else
            {
                ringmesh_assert_not_reached;
                return;
            }
        }

        auto& children_entity =
            geomodel_access_.modifiable_mesh_entity( children );
        const auto& children_type =
            geomodel_.entity_type_manager().relationship_manager.child_type(
                parent.type() );

        if( children_type != children.type() )
        {
            throw RINGMeshException( "Entity",
                "Wrong children type in the parent children relation between ",
                parent, " and ", children );
        }

        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        auto relation_id =
            manager.add_parent_child_relationship( parent, children );
        GeoModelGeologicalEntityAccess< DIMENSION > parent_access{
            parent_entity
        };
        parent_access.modifiable_children().push_back( relation_id );
        GeoModelMeshEntityAccess< DIMENSION > children_access{
            children_entity
        };
        children_access.modifiable_parents().push_back( relation_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::remove_parent_children_relation(
        const gmge_id& parent, const gmme_id& children )
    {
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        auto relation_id =
            manager.find_parent_child_relationship( parent, children );
        if( relation_id == NO_ID )
        {
            throw RINGMeshException( "Entity",
                "No parent children relation found between ", parent, " and ",
                children );
        }
        GeoModelGeologicalEntityAccess< DIMENSION > parent_access{
            geomodel_access_.modifiable_geological_entity( parent )
        };
        auto& childs = parent_access.modifiable_children();
        std::remove_if(
            childs.begin(), childs.end(), [relation_id]( index_t relation ) {
                return relation == relation_id;
            } );
        GeoModelMeshEntityAccess< DIMENSION > children_access{
            geomodel_access_.modifiable_mesh_entity( children )
        };
        auto& parents = children_access.modifiable_parents();
        std::remove_if(
            parents.begin(), parents.end(), [relation_id]( index_t relation ) {
                return relation == relation_id;
            } );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::set_geological_entity_child(
        const gmge_id& parent_gmge, index_t id, index_t child_id )
    {
        /// No check on the validity of the index of the entity child_index
        /// NO_ID is used to flag entities to delete
        auto& geol_entity =
            geomodel_access_.modifiable_geological_entity( parent_gmge );
        const auto& child_type =
            geomodel_.entity_type_manager().relationship_manager.child_type(
                parent_gmge.type() );
        gmme_id child{ child_type, child_id };
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access{ geol_entity };
        auto relationship_id = gmge_access.modifiable_children()[id];
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        manager.set_child_to_parent_child_relationship(
            relationship_id, child );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::delete_geological_entity(
        const GeologicalEntityType& type, index_t index )
    {
        geomodel_access_.modifiable_geological_entities( type )[index].reset();
    }

    template < index_t DIMENSION >
    gmge_id GeoModelBuilderGeology< DIMENSION >::create_geological_entity(
        const GeologicalEntityType& type )
    {
        auto index = find_or_create_geological_entity_type( type );
        auto id =
            static_cast< index_t >( geomodel_.nb_geological_entities( type ) );
        geomodel_access_.modifiable_geological_entities()[index].emplace_back(
            GeoModelGeologicalEntityAccess<
                DIMENSION >::create_geological_entity( type, geomodel_, id ) );
        return geomodel_access_.modifiable_geological_entities()[index]
            .back()
            ->gmge();
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderGeology< DIMENSION >::
        find_or_create_geological_entity_type(
            const GeologicalEntityType& type )
    {
        auto type_index =
            geomodel_.entity_type_manager()
                .geological_entity_manager.geological_entity_type_index( type );
        if( type_index == NO_ID )
        {
            type_index = create_geological_entity_type( type );
        }
        return type_index;
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderGeology< DIMENSION >::create_geological_entity_type(
        const GeologicalEntityType& type )
    {
        ringmesh_assert(
            GeoModelGeologicalEntityFactory< DIMENSION >::has_creator( type ) );

        geomodel_access_.modifiable_entity_type_manager()
            .geological_entity_manager.geological_entity_types_.push_back(
                type );
        geomodel_access_.modifiable_geological_entities().emplace_back();
        auto geol_entity = GeoModelGeologicalEntityFactory< DIMENSION >::create(
            type, geomodel_ );

        const auto child_type = geol_entity->child_type_name();
        auto& parentage = geomodel_access_.modifiable_entity_type_manager()
                              .relationship_manager;
        parentage.register_geology_relationship( type, child_type );

        return geomodel_.entity_type_manager()
                   .geological_entity_manager.nb_geological_entity_types()
               - 1;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::copy_geological_entity_topology(
        const GeoModel< DIMENSION >& from, const GeologicalEntityType& type )
    {
        create_geological_entities( type, from.nb_geological_entities( type ) );
        parallel_for( geomodel_.nb_geological_entities( type ),
            [&type, &from, this]( index_t i ) {
                gmge_id id{ type, i };
                GeoModelGeologicalEntityAccess< DIMENSION > gmge_access{
                    geomodel_access_.modifiable_geological_entity( id )
                };
                gmge_access.copy( from.geological_entity( id ) );
            } );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::build_contacts()
    {
        if( geomodel_.entity_type_manager()
                .geological_entity_manager.is_valid_type(
                    Contact< DIMENSION >::type_name_static() )
            && geomodel_.nb_geological_entities(
                   Contact< DIMENSION >::type_name_static() )
                   > 0 )
        {
            return;
        }

        std::vector< std::set< gmge_id > > interfaces;
        for( const auto& line : geomodel_.lines() )
        {
            std::set< gmge_id > cur_interfaces;
            for( auto j : range( line.nb_incident_entities() ) )
            {
                const auto& surface = line.incident_entity( j );
                auto parent_interface = surface.parent_gmge(
                    Interface< DIMENSION >::type_name_static() );
                cur_interfaces.insert( parent_interface );
            }
            gmge_id contact_id;
            for( auto j : range( interfaces.size() ) )
            {
                if( cur_interfaces.size() == interfaces[j].size()
                    && std::equal( cur_interfaces.begin(), cur_interfaces.end(),
                           interfaces[j].begin() ) )
                {
                    contact_id =
                        gmge_id{ Contact< DIMENSION >::type_name_static(), j };
                    break;
                }
            }
            if( !contact_id.is_defined() )
            {
                contact_id = create_geological_entity(
                    Contact< DIMENSION >::type_name_static() );
                ringmesh_assert( contact_id.index() == interfaces.size() );
                // Create a name for this contact
                std::string name = "contact";
                for( const auto& it : cur_interfaces )
                {
                    name += "_";
                    name += geomodel_.geological_entity( it ).name();
                }
                builder_.info.set_geological_entity_name( contact_id, name );
                interfaces.emplace_back( std::move( cur_interfaces ) );
            }
            add_parent_children_relation( contact_id, line.gmme() );
        }
    }
    template < index_t DIMENSION >
    void
        GeoModelBuilderGeology< DIMENSION >::set_geological_entity_geol_feature(
            const gmge_id& gmge_id,
            typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE
                geol_feature )
    {
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access{
            geomodel_access_.modifiable_geological_entity( gmge_id )
        };
        gmge_access.modifiable_geol_feature() = geol_feature;
    }

    template class geomodel_builder_api GeoModelBuilderGeology< 2 >;
    template class geomodel_builder_api GeoModelBuilderGeology< 3 >;

} // namespace RINGMesh
