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

#include <algorithm>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/task_handler.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/builder/geomodel_builder_topology.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

/*!
 * @file Implementation of GeoModel topological edition functions
 * @author Jeanne Pellerin
 * @author Pierre Anquez
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel, const vecn< DIMENSION >& point )
    {
        for( const auto& corner : geomodel.corners() )
        {
            if( corner.vertex( 0 ) == point )
            {
                return corner.gmme();
            }
        }
        return gmme_id();
    }

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
        return gmme_id();
    }

    /*!
     * @brief Returns true if the Line has exactly the given vertices
     *
     * @param[in] L the line to compare to
     * @param[in] rhs_vertices Vertices to compare to
     */
    template < index_t DIMENSION >
    bool line_equal( const Line< DIMENSION >& line,
        const std::vector< vecn< DIMENSION > >& rhs_vertices )
    {
        if( line.nb_vertices() != rhs_vertices.size() )
        {
            return false;
        }
        bool equal = true;
        for( auto v : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[v] != line.vertex( v ) )
            {
                equal = false;
                break;
            }
        }
        if( equal )
        {
            return true;
        }

        equal = true;
        for( auto v : range( line.nb_vertices() ) )
        {
            if( rhs_vertices[v] != line.vertex( line.nb_vertices() - v - 1 ) )
            {
                equal = false;
                break;
            }
        }
        return equal;
    }

    template < index_t DIMENSION >
    std::vector< index_t > get_sorted_incident_surfaces(
        const GeoModelMeshEntity< DIMENSION >& E )
    {
        std::vector< index_t > incident_surfaces;
        index_t nb{ E.nb_incident_entities() };
        incident_surfaces.resize( nb );
        for( auto i : range( nb ) )
        {
            incident_surfaces[i] = E.incident_entity_gmme( i ).index();
        }
        std::sort( incident_surfaces.begin(), incident_surfaces.end() );
        return incident_surfaces;
    }

    template < index_t DIMENSION >
    index_t add_to_set_children_of_geological_entities(
        const GeoModel< DIMENSION >& geomodel,
        const std::set< gmge_id >& geological_entities,
        std::set< gmme_id >& mesh_entities )
    {
        const auto old_number_of_entities = static_cast< index_t >(
            mesh_entities.size() + geological_entities.size() );
        for( const auto& cur_gmge_id : geological_entities )
        {
            const auto& cur_geol_entity =
                geomodel.geological_entity( cur_gmge_id );
            for( auto child_i : range( cur_geol_entity.nb_children() ) )
            {
                mesh_entities.insert( cur_geol_entity.child_gmme( child_i ) );
            }
        }
        const auto new_number_of_entities = static_cast< index_t >(
            mesh_entities.size() + geological_entities.size() );
        ringmesh_assert( new_number_of_entities >= old_number_of_entities );
        return new_number_of_entities - old_number_of_entities;
    }

    template < index_t DIMENSION >
    index_t add_to_set_geological_entities_which_have_no_child(
        const GeoModel< DIMENSION >& geomodel,
        const std::set< gmme_id >& mesh_entities,
        std::set< gmge_id >& geological_entities )
    {
        const auto old_number_of_entities = static_cast< index_t >(
            mesh_entities.size() + geological_entities.size() );
        const auto nb_geological_entity_types =
            geomodel.entity_type_manager()
                .geological_entity_manager.nb_geological_entity_types();
        for( auto geol_entity_type_i : range( nb_geological_entity_types ) )
        {
            const auto& geol_type =
                geomodel.entity_type_manager()
                    .geological_entity_manager.geological_entity_type(
                        geol_entity_type_i );

            for( auto geol_entity_i :
                range( geomodel.nb_geological_entities( geol_type ) ) )
            {
                bool no_child{ true };
                const auto& cur_geol_entity =
                    geomodel.geological_entity( geol_type, geol_entity_i );
                for( auto child_i : range( cur_geol_entity.nb_children() ) )
                {
                    if( mesh_entities.count(
                            cur_geol_entity.child_gmme( child_i ) )
                        == 0 )
                    {
                        no_child = false;
                        break;
                    }
                }
                if( no_child )
                {
                    geological_entities.insert( cur_geol_entity.gmge() );
                }
            }
        }
        const auto new_number_of_entities = static_cast< index_t >(
            mesh_entities.size() + geological_entities.size() );
        ringmesh_assert( new_number_of_entities >= old_number_of_entities );
        return new_number_of_entities - old_number_of_entities;
    }

    template < index_t DIMENSION >
    index_t add_to_set_mesh_entities_being_boundaries_of_no_mesh_entity(
        const GeoModel< DIMENSION >& geomodel,
        std::set< gmme_id >& mesh_entities )
    {
        const auto old_number_of_entities =
            static_cast< index_t >( mesh_entities.size() );
        for( const auto& mesh_type :
            geomodel.entity_type_manager()
                .mesh_entity_manager.mesh_entity_types() )
        {
            if( mesh_type == Region3D::type_name_static() )
            {
                // A Region3D cannot be a boundary.
                continue;
            }
            for( auto mesh_entity_i :
                range( geomodel.nb_mesh_entities( mesh_type ) ) )
            {
                bool no_incident{ true };
                const auto& cur_mesh_entity =
                    geomodel.mesh_entity( mesh_type, mesh_entity_i );
                for( auto incident_entity_i :
                    range( cur_mesh_entity.nb_incident_entities() ) )
                {
                    if( mesh_entities.count(
                            cur_mesh_entity.incident_entity_gmme(
                                incident_entity_i ) )
                        == 0 )
                    {
                        no_incident = false;
                        break;
                    }
                }
                if( no_incident )
                {
                    mesh_entities.insert( cur_mesh_entity.gmme() );
                }
            }
        }
        const auto new_number_of_entities =
            static_cast< index_t >( mesh_entities.size() );
        ringmesh_assert( new_number_of_entities >= old_number_of_entities );
        return new_number_of_entities - old_number_of_entities;
    }

} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBuilderTopologyBase< DIMENSION >::GeoModelBuilderTopologyBase(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : builder_( builder ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
    }

    template < index_t DIMENSION >
    void
        GeoModelBuilderTopologyBase< DIMENSION >::copy_all_mesh_entity_topology(
            const GeoModel< DIMENSION >& from )
    {
        copy_mesh_entity_topology< Corner >( from );
        copy_mesh_entity_topology< Line >( from );
        copy_mesh_entity_topology< Surface >( from );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::copy_topology(
        const GeoModel< DIMENSION >& from )
    {
        copy_all_mesh_entity_topology( from );

        geomodel_access_.modifiable_epsilon() = from.epsilon();
        geomodel_access_.modifiable_entity_type_manager()
            .relationship_manager.copy(
                from.entity_type_manager().relationship_manager );
    }

    template < index_t DIMENSION >
    bool GeoModelBuilderTopologyBase< DIMENSION >::get_dependent_entities(
        std::set< gmme_id >& mesh_entities,
        std::set< gmge_id >& geological_entities ) const
    {
        auto nb_added = add_to_set_children_of_geological_entities(
            geomodel_, geological_entities, mesh_entities );
        nb_added += add_to_set_geological_entities_which_have_no_child(
            geomodel_, mesh_entities, geological_entities );
        nb_added += add_to_set_mesh_entities_being_boundaries_of_no_mesh_entity(
            geomodel_, mesh_entities );

        // Recursive call till nothing is added
        if( nb_added > 0 )
        {
            return get_dependent_entities( mesh_entities, geological_entities );
        }
        return false;
    }

    template < index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_corner(
        const vecn< DIMENSION >& point, const MeshType& mesh_type )
    {
        gmme_id result{ find_corner( geomodel_, point ) };
        if( !result.is_defined() )
        {
            result = create_mesh_entity< Corner >( mesh_type );
            builder_.geometry.set_corner( result.index(), point );
        }
        return result;
    }

    template < index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_corner(
        index_t geomodel_point_id, const MeshType& mesh_type )
    {
        gmme_id result{ find_corner( geomodel_, geomodel_point_id ) };
        if( !result.is_defined() )
        {
            result = create_mesh_entity< Corner >( mesh_type );
            builder_.geometry.set_corner( result.index(), geomodel_point_id );
        }
        return result;
    }

    template < index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_line(
        const std::vector< vecn< DIMENSION > >& vertices,
        const MeshType& mesh_type )
    {
        gmme_id result;
        for( const auto& line : geomodel_.lines() )
        {
            if( line_equal( line, vertices ) )
            {
                result = line.gmme();
            }
        }
        if( !result.is_defined() )
        {
            result = create_mesh_entity< Line >( mesh_type );
            builder_.geometry.set_line( result.index(), vertices );

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_mesh_entity_boundary_relation(
                result, find_or_create_corner( vertices.front() ) );
            add_mesh_entity_boundary_relation(
                result, find_or_create_corner( vertices.back() ) );
        }
        return result;
    }

    template < index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        const gmme_id& first_corner,
        const gmme_id& second_corner,
        const MeshType& mesh_type )
    {
        for( const auto& line : geomodel_.lines() )
        {
            const auto& c0 = line.boundary_gmme( 0 );
            const auto& c1 = line.boundary_gmme( 1 );

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) )
            {
                std::vector< index_t > cur_adjacent_surfaces{
                    get_sorted_incident_surfaces( line )
                };
                if( cur_adjacent_surfaces.size()
                        == sorted_adjacent_surfaces.size()
                    && std::equal( cur_adjacent_surfaces.begin(),
                           cur_adjacent_surfaces.end(),
                           sorted_adjacent_surfaces.begin() ) )
                {
                    return line.gmme();
                }
            }
        }
        return create_mesh_entity< Line >( mesh_type );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::
        remove_mesh_entity_boundary_relation(
            const gmme_id& incident_entity, const gmme_id& boundary )
    {
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        index_t relation_id{ manager.find_boundary_relationship(
            incident_entity, boundary ) };
        if( relation_id == NO_ID )
        {
            throw RINGMeshException( "Entity",
                "No boundary relation found between ", boundary, " and ",
                incident_entity );
        }
        GeoModelMeshEntityAccess< DIMENSION > boundary_access(
            geomodel_access_.modifiable_mesh_entity( boundary ) );
        auto& incident_entities =
            boundary_access.modifiable_incident_entities();
        std::remove_if( incident_entities.begin(), incident_entities.end(),
            [relation_id](
                index_t relation ) { return relation == relation_id; } );
        GeoModelMeshEntityAccess< DIMENSION > incident_entity_access(
            geomodel_access_.modifiable_mesh_entity( incident_entity ) );
        auto& boundaries = incident_entity_access.modifiable_boundaries();
        std::remove_if( boundaries.begin(), boundaries.end(),
            [relation_id](
                index_t relation ) { return relation == relation_id; } );
    }

    template < index_t DIMENSION >
    template < template < index_t > class ENTITY >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entity(
        const MeshType& mesh_type )
    {
        const auto& entity_type = ENTITY< DIMENSION >::type_name_static();
        index_t nb_entities{ geomodel_.nb_mesh_entities( entity_type ) };
        index_t new_id{ nb_entities };
        geomodel_access_.modifiable_mesh_entities( entity_type )
            .emplace_back( GeoModelMeshEntityAccess< DIMENSION >::
                    template create_entity< ENTITY >(
                        geomodel_, new_id, mesh_type ) );
        return geomodel_access_.modifiable_mesh_entities( entity_type )
            .back()
            ->gmme();
    }

    template < index_t DIMENSION >
    template < template < index_t > class ENTITY >
    bool GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entities(
        index_t nb_additionnal_entities, const MeshType& type )
    {
        const auto& entity_type = ENTITY< DIMENSION >::type_name_static();
        auto& store = geomodel_access_.modifiable_mesh_entities( entity_type );
        index_t old_size{ static_cast< index_t >( store.size() ) };
        index_t new_size{ old_size + nb_additionnal_entities };
        store.reserve( new_size );
        for( auto i : range( old_size, new_size ) )
        {
            store.emplace_back( GeoModelMeshEntityAccess< DIMENSION >::
                    template create_entity< ENTITY >( geomodel_, i, type ) );
        }
        return true;
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderTopologyBase< DIMENSION >::
        check_if_boundary_incident_entity_relation_already_exists(
            const gmme_id& incident_entity, const gmme_id& boundary )
    {
        const auto& incident_entity_mesh_entity =
            geomodel_.mesh_entity( incident_entity );
        for( auto in_ent :
            range( incident_entity_mesh_entity.nb_incident_entities() ) )
        {
            if( incident_entity_mesh_entity.incident_entity_gmme( in_ent )
                == boundary )
            {
                GeoModelMeshEntityConstAccess< DIMENSION > entity_access(
                    incident_entity_mesh_entity );
                return entity_access.incident_entity_relation_ids()[in_ent];
            }
        }
        return NO_ID;
    }
    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::
        add_mesh_entity_boundary_relation(
            const gmme_id& incident_entity_id, const gmme_id& boundary_id )
    {
        const auto& incident_entity_type =
            geomodel_.entity_type_manager()
                .mesh_entity_manager.incident_entity_type( boundary_id.type() );
        if( incident_entity_id.type() != incident_entity_type )
        {
            throw RINGMeshException( "Entity",
                "Wrong incident entity type in the boundary relation between ",
                boundary_id, " and ", incident_entity_id );
        }
        const auto& boundary_type =
            geomodel_.entity_type_manager()
                .mesh_entity_manager.boundary_entity_type(
                    incident_entity_id.type() );
        if( boundary_id.type() != boundary_type )
        {
            throw RINGMeshException( "Entity",
                "Wrong boundary type in the boundary relation between ",
                boundary_id, " and ", incident_entity_id );
        }
        index_t relation_id{
            check_if_boundary_incident_entity_relation_already_exists(
                incident_entity_id, boundary_id )
        };
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        if( relation_id == NO_ID )
        {
            relation_id = manager.add_boundary_relationship(
                incident_entity_id, boundary_id );
        }
        auto& boundary_entity =
            geomodel_access_.modifiable_mesh_entity( boundary_id );
        GeoModelMeshEntityAccess< DIMENSION > boundary_access(
            boundary_entity );
        boundary_access.modifiable_incident_entities().push_back( relation_id );
        auto& incident_entity =
            geomodel_access_.modifiable_mesh_entity( incident_entity_id );
        GeoModelMeshEntityAccess< DIMENSION > incident_entity_access(
            incident_entity );
        incident_entity_access.modifiable_boundaries().push_back( relation_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::set_mesh_entity_boundary(
        const gmme_id& gmme,
        index_t current_local_boundary_id,
        index_t new_global_boundary_id )
    {
        ringmesh_assert( current_local_boundary_id
                         < geomodel_.mesh_entity( gmme ).nb_boundaries() );
        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        const auto& b_type =
            geomodel_.entity_type_manager()
                .mesh_entity_manager.boundary_entity_type( gmme.type() );
        gmme_id boundary( b_type, new_global_boundary_id );
        GeoModelMeshEntityAccess< DIMENSION > gme_access( mesh_entity );
        index_t relation_id{
            gme_access.modifiable_boundaries()[current_local_boundary_id]
        };
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        manager.set_boundary_to_boundary_relationship( relation_id, boundary );
    }

    template < index_t DIMENSION >
    template < template < index_t > class ENTITY >
    void GeoModelBuilderTopologyBase< DIMENSION >::copy_mesh_entity_topology(
        const GeoModel< DIMENSION >& from )
    {
        const auto& type = ENTITY< DIMENSION >::type_name_static();
        create_mesh_entities< ENTITY >( from.nb_mesh_entities( type ), "" );

        parallel_for( geomodel_.nb_mesh_entities( type ),
            [&type, &from, this]( index_t i ) {
                gmme_id id( type, i );
                GeoModelMeshEntityAccess< DIMENSION > gmme_access(
                    geomodel_access_.modifiable_mesh_entity( id ) );
                gmme_access.copy( from.mesh_entity( id ) );
            } );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::
        set_mesh_entity_incident_entity( const gmme_id& gmme,
            index_t current_local_incident_entity_id,
            index_t new_global_incident_entity_id )
    {
        /// No check on the validity of the index of the entity incident_entity
        /// NO_ID is used to flag entities to delete
        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        ringmesh_assert( current_local_incident_entity_id
                         < mesh_entity.nb_incident_entities() );
        const auto& in_ent_type =
            geomodel_.entity_type_manager()
                .mesh_entity_manager.incident_entity_type( gmme.type() );
        gmme_id incident_entity( in_ent_type, new_global_incident_entity_id );
        GeoModelMeshEntityAccess< DIMENSION > gme_access( mesh_entity );
        index_t relation_id{ gme_access.modifiable_incident_entities()
                                 [current_local_incident_entity_id] };
        auto& manager = geomodel_access_.modifiable_entity_type_manager()
                            .relationship_manager;
        manager.set_incident_entity_to_boundary_relationship(
            relation_id, incident_entity );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::delete_mesh_entity(
        const MeshEntityType& type, index_t index )
    {
        geomodel_access_.modifiable_mesh_entities( type )[index].reset();
    }

    template < index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entity(
        const MeshEntityType& entity_type, const MeshType& mesh_type )
    {
        const auto& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( entity_type ) )
        {
            return this->create_mesh_entity< Corner >( mesh_type );
        }
        if( manager.is_line( entity_type ) )
        {
            return create_mesh_entity< Line >( mesh_type );
        }
        if( manager.is_surface( entity_type ) )
        {
            return create_mesh_entity< Surface >( mesh_type );
        }
        ringmesh_assert_not_reached;
        return gmme_id();
    }

    template < index_t DIMENSION >
    bool GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entities(
        const MeshEntityType& entity_type,
        index_t nb_additional_entities,
        const MeshType& mesh_type )
    {
        const auto& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( entity_type ) )
        {
            return this->create_mesh_entities< Corner >(
                nb_additional_entities, mesh_type );
        }
        if( manager.is_line( entity_type ) )
        {
            return create_mesh_entities< Line >(
                nb_additional_entities, mesh_type );
        }
        if( manager.is_surface( entity_type ) )
        {
            return create_mesh_entities< Surface >(
                nb_additional_entities, mesh_type );
        }
        ringmesh_assert_not_reached;
        return false;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::set_line_corner_boundary(
        index_t incident_line_id,
        index_t current_local_boundary_corner_id,
        index_t new_global_boundary_corner_id )
    {
        set_mesh_entity_boundary(
            { Line< DIMENSION >::type_name_static(), incident_line_id },
            current_local_boundary_corner_id, new_global_boundary_corner_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::
        add_line_corner_boundary_relation(
            index_t incident_line_id, index_t boundary_corner_id )
    {
        add_mesh_entity_boundary_relation(
            { Line< DIMENSION >::type_name_static(), incident_line_id },
            { Corner< DIMENSION >::type_name_static(), boundary_corner_id } );
    }

    gmme_id GeoModelBuilderTopology< 3 >::create_mesh_entity(
        const MeshEntityType& entity_type, const MeshType& mesh_type )
    {
        const auto& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( entity_type ) )
        {
            return GeoModelBuilderTopologyBase3D::create_mesh_entity< Region >(
                mesh_type );
        }
        return GeoModelBuilderTopologyBase3D::create_mesh_entity(
            entity_type, mesh_type );
    }
    bool GeoModelBuilderTopology< 3 >::create_mesh_entities(
        const MeshEntityType& entity_type,
        index_t nb_additional_entities,
        const MeshType& mesh_type )
    {
        const auto& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( entity_type ) )
        {
            return GeoModelBuilderTopologyBase3D::
                create_mesh_entities< Region >(
                    nb_additional_entities, mesh_type );
        }
        return GeoModelBuilderTopologyBase3D::create_mesh_entities(
            entity_type, nb_additional_entities, mesh_type );
    }
    void GeoModelBuilderTopology< 3 >::copy_all_mesh_entity_topology(
        const GeoModel3D& from )
    {
        GeoModelBuilderTopologyBase3D::copy_all_mesh_entity_topology( from );
        copy_mesh_entity_topology< Region >( from );
    }

    void GeoModelBuilderTopology< 2 >::set_surface_line_boundary(
        index_t incident_surface_id,
        index_t current_local_boundary_line_id,
        index_t new_global_boundary_line_id,
        bool side )
    {
        GeoModelBuilderTopologyBase2D::set_mesh_entity_boundary(
            { Surface2D::type_name_static(), incident_surface_id },
            current_local_boundary_line_id, new_global_boundary_line_id );

        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            { Surface2D::type_name_static(), incident_surface_id } );
        GeoModelMeshEntityAccess2D gme_access( mesh_entity );
        gme_access.modifiable_sides()[current_local_boundary_line_id] = side;
    }

    void GeoModelBuilderTopology< 2 >::add_surface_line_boundary_relation(
        index_t incident_surface_id, index_t boundary_line_id, bool side )
    {
        GeoModelBuilderTopologyBase2D::add_mesh_entity_boundary_relation(
            { Surface2D::type_name_static(), incident_surface_id },
            { Line2D::type_name_static(), boundary_line_id } );

        auto& incident_entity = geomodel_access_.modifiable_mesh_entity(
            { Surface2D::type_name_static(), incident_surface_id } );
        GeoModelMeshEntityAccess2D incident_entity_access( incident_entity );
        incident_entity_access.modifiable_sides().push_back( side );
    }

    void GeoModelBuilderTopology< 3 >::set_surface_line_boundary(
        index_t incident_surface_id,
        index_t current_local_boundary_line_id,
        index_t new_boundary_gobal_line_id )
    {
        GeoModelBuilderTopologyBase3D::set_mesh_entity_boundary(
            { Surface3D::type_name_static(), incident_surface_id },
            current_local_boundary_line_id, new_boundary_gobal_line_id );
    }

    void GeoModelBuilderTopology< 3 >::add_surface_line_boundary_relation(
        index_t incident_surface_id, index_t boundary_line_id )
    {
        GeoModelBuilderTopologyBase3D::add_mesh_entity_boundary_relation(
            { Surface3D::type_name_static(), incident_surface_id },
            { Line3D::type_name_static(), boundary_line_id } );
    }

    void GeoModelBuilderTopology< 3 >::set_region_surface_boundary(
        index_t incident_region_id,
        index_t current_local_boundary_surface_id,
        index_t new_global_boundary_surface_id,
        bool side )
    {
        GeoModelBuilderTopologyBase3D::set_mesh_entity_boundary(
            { Region3D::type_name_static(), incident_region_id },
            current_local_boundary_surface_id, new_global_boundary_surface_id );

        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            { Region3D::type_name_static(), incident_region_id } );
        GeoModelMeshEntityAccess3D gme_access( mesh_entity );
        gme_access.modifiable_sides()[current_local_boundary_surface_id] = side;
    }

    void GeoModelBuilderTopology< 3 >::add_region_surface_boundary_relation(
        index_t incident_region_id, index_t boundary_surface_id, bool side )
    {
        GeoModelBuilderTopologyBase3D::add_mesh_entity_boundary_relation(
            { Region3D::type_name_static(), incident_region_id },
            { Surface3D::type_name_static(), boundary_surface_id } );

        auto& incident_entity = geomodel_access_.modifiable_mesh_entity(
            { Region3D::type_name_static(), incident_region_id } );
        GeoModelMeshEntityAccess3D incident_entity_access( incident_entity );
        incident_entity_access.modifiable_sides().push_back( side );
    }

    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 2 >::create_mesh_entity< Corner >(
            const MeshType& );
    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 2 >::create_mesh_entity< Line >(
            const MeshType& );
    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 2 >::create_mesh_entity< Surface >(
            const MeshType& );
    template class geomodel_builder_api GeoModelBuilderTopologyBase< 2 >;

    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 3 >::create_mesh_entity< Corner >(
            const MeshType& );
    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 3 >::create_mesh_entity< Line >(
            const MeshType& );
    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 3 >::create_mesh_entity< Surface >(
            const MeshType& );
    template gmme_id geomodel_builder_api
        GeoModelBuilderTopologyBase< 3 >::create_mesh_entity< Region >(
            const MeshType& );
    template class geomodel_builder_api GeoModelBuilderTopologyBase< 3 >;
} // namespace RINGMesh
