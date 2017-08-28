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

#include <ringmesh/geomodel/geomodel_builder_topology.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/geomodel_builder.h>

/*!
 * @file Implementation of GeoModel topological edition functions
 * @author Jeanne Pellerin
 * @author Pierre Anquez
 */

namespace {
    using namespace RINGMesh;

    template< index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel,
        const vecn< DIMENSION >& point )
    {
        for( const auto& corner : geomodel.corners() ) {
            if( corner.vertex( 0 ) == point ) {
                return corner.gmme();
            }
        }
        return gmme_id();
    }

    template< index_t DIMENSION >
    gmme_id find_corner(
        const GeoModel< DIMENSION >& geomodel,
        index_t geomodel_point_id )
    {
        const auto& geomodel_vertices = geomodel.mesh.vertices;
        const auto& vertices = geomodel_vertices.gme_vertices( geomodel_point_id );
        for( const auto& vertex : vertices ) {
            if( vertex.gmme.type() == Corner< DIMENSION >::type_name_static() ) {
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
    template< index_t DIMENSION >
    bool line_equal(
        const Line< DIMENSION >& line,
        const std::vector< vecn< DIMENSION > >& rhs_vertices )
    {
        if( line.nb_vertices() != rhs_vertices.size() ) {
            return false;
        }
        bool equal = true;
        for( auto v : range( line.nb_vertices() ) ) {
            if( rhs_vertices[v] != line.vertex( v ) ) {
                equal = false;
                break;
            }
        }
        if( equal ) return true;

        equal = true;
        for( auto v : range( line.nb_vertices() ) ) {
            if( rhs_vertices[v] != line.vertex( line.nb_vertices() - v - 1 ) ) {
                equal = false;
                break;
            }
        }
        return equal;
    }

    template< index_t DIMENSION >
    std::vector< index_t > get_sorted_incident_surfaces(
        const GeoModelMeshEntity< DIMENSION >& E )
    {
        std::vector< index_t > incident_surfaces;
        index_t nb { E.nb_incident_entities() };
        incident_surfaces.resize( nb );
        for( auto i : range( nb ) ) {
            incident_surfaces[i] = E.incident_entity_gmme( i ).index();
        }
        std::sort( incident_surfaces.begin(), incident_surfaces.end() );
        return incident_surfaces;
    }

    template< index_t DIMENSION >
    index_t add_to_set_children_of_geological_entities(
        const GeoModel< DIMENSION >& geomodel,
        const std::set< gmge_id >& geological_entities,
        std::set< gmme_id >& mesh_entities )
    {
        int nb_added =
            -1
                * static_cast< int >( mesh_entities.size()
                    + geological_entities.size() );
        for( const gmge_id& cur_gmge_id : geological_entities ) {
            const GeoModelGeologicalEntity< DIMENSION >& cur_geol_entity =
                geomodel.geological_entity( cur_gmge_id );
            for( index_t child_i = 0; child_i < cur_geol_entity.nb_children();
                ++child_i ) {
                mesh_entities.insert( cur_geol_entity.child_gmme( child_i ) );
            }
        }
        nb_added += static_cast< int >( mesh_entities.size()
            + geological_entities.size() );
        ringmesh_assert( nb_added >= 0 );
        return static_cast< index_t >( nb_added );
    }

    template< index_t DIMENSION >
    index_t add_to_set_geological_entities_which_have_no_child(
        const GeoModel< DIMENSION >& geomodel,
        const std::set< gmme_id >& mesh_entities,
        std::set< gmge_id >& geological_entities )
    {
        int nb_added =
            -1
                * static_cast< int >( mesh_entities.size()
                    + geological_entities.size() );
        const index_t nb_geological_entity_types =
            geomodel.entity_type_manager().geological_entity_manager.nb_geological_entity_types();
        for( index_t geol_entity_type_i = 0;
            geol_entity_type_i < nb_geological_entity_types; ++geol_entity_type_i ) {
            const GeologicalEntityType& geol_type =
                geomodel.entity_type_manager().geological_entity_manager.geological_entity_type(
                    geol_entity_type_i );

            for( index_t geol_entity_i = 0;
                geol_entity_i < geomodel.nb_geological_entities( geol_type );
                ++geol_entity_i ) {
                bool no_child = true;
                const GeoModelGeologicalEntity< DIMENSION >& cur_geol_entity =
                    geomodel.geological_entity( geol_type, geol_entity_i );
                for( index_t child_i = 0; child_i < cur_geol_entity.nb_children();
                    ++child_i ) {
                    if( mesh_entities.count( cur_geol_entity.child_gmme( child_i ) )
                        == 0 ) {
                        no_child = false;
                        break;
                    }
                }
                if( no_child ) {
                    geological_entities.insert( cur_geol_entity.gmge() );
                }
            }
        }
        nb_added += static_cast< int >( mesh_entities.size()
            + geological_entities.size() );
        ringmesh_assert( nb_added >= 0 );
        return static_cast< index_t >( nb_added );
    }

    template< index_t DIMENSION >
    index_t add_to_set_mesh_entities_being_boundaries_of_no_mesh_entity(
        const GeoModel< DIMENSION >& geomodel,
        std::set< gmme_id >& mesh_entities )
    {
        int nb_added { -1 * static_cast< int >( mesh_entities.size() ) };
        for( auto mesh_type : geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types() ) {
            ringmesh_assert( mesh_type != Region3D::type_name_static() );
            for( index_t mesh_entity_i : range(
                geomodel.nb_mesh_entities( mesh_type ) ) ) {
                bool no_incident = true;
                const GeoModelMeshEntity< DIMENSION >& cur_mesh_entity =
                    geomodel.mesh_entity( mesh_type, mesh_entity_i );
                for( index_t incident_entity_i : range(
                    cur_mesh_entity.nb_incident_entities() ) ) {
                    if( mesh_entities.count(
                        cur_mesh_entity.incident_entity_gmme( incident_entity_i ) )
                        == 0 ) {
                        no_incident = false;
                        break;
                    }
                }
                if( no_incident ) {
                    mesh_entities.insert( cur_mesh_entity.gmme() );
                }
            }
        }
        nb_added += static_cast< int >( mesh_entities.size() );
        ringmesh_assert( nb_added >= 0 );
        return static_cast< index_t >( nb_added );
    }

}

namespace RINGMesh {

    template< index_t DIMENSION >
    GeoModelBuilderTopologyBase< DIMENSION >::GeoModelBuilderTopologyBase(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::copy_all_mesh_entity_topology(
        const GeoModel< DIMENSION >& from )
    {

        copy_mesh_entity_topology< Corner >( from );
        copy_mesh_entity_topology< Line >( from );
        copy_mesh_entity_topology< Surface >( from );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::copy_topology(
        const GeoModel< DIMENSION >& from )
    {
        copy_all_mesh_entity_topology( from );

        UniverseAccess< DIMENSION > universe_access(
            geomodel_access_.modifiable_universe() );
        universe_access.copy( from.universe() );
        geomodel_access_.modifiable_epsilon() = from.epsilon();
        geomodel_access_.modifiable_entity_type_manager().relationship_manager =
            from.entity_type_manager().relationship_manager;
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderTopologyBase< DIMENSION >::get_dependent_entities(
        std::set< gmme_id >& mesh_entities,
        std::set< gmge_id >& geological_entities ) const
    {
        index_t nb_added = add_to_set_children_of_geological_entities( geomodel_,
            geological_entities, mesh_entities );
        nb_added += add_to_set_geological_entities_which_have_no_child( geomodel_,
            mesh_entities, geological_entities );
        nb_added += add_to_set_mesh_entities_being_boundaries_of_no_mesh_entity(
            geomodel_, mesh_entities );

        // Recursive call till nothing is added
        if( nb_added > 0 ) {
            return get_dependent_entities( mesh_entities, geological_entities );
        } else {
            return false;
        }
    }

    template< index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_corner(
        const vecn< DIMENSION >& point )
    {
        gmme_id result { find_corner( geomodel_, point ) };
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >();
            builder_.geometry.set_corner( result.index(), point );
        }
        return result;
    }

    template< index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_corner(
        index_t geomodel_point_id )
    {
        gmme_id result { find_corner( geomodel_, geomodel_point_id ) };
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >();
            builder_.geometry.set_corner( result.index(), geomodel_point_id );
        }
        return result;
    }

    template< index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_line(
        const std::vector< vecn< DIMENSION > >& vertices )
    {
        gmme_id result;
        for( const auto& line : geomodel_.lines() ) {
            if( line_equal( line, vertices ) ) {
                result = line.gmme();
            }
        }
        if( !result.is_defined() ) {
            result = create_mesh_entity< Line >();
            builder_.geometry.set_line( result.index(), vertices );

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_mesh_entity_boundary_relation( result,
                find_or_create_corner( vertices.front() ) );
            add_mesh_entity_boundary_relation( result,
                find_or_create_corner( vertices.back() ) );
        }
        return result;
    }

    template< index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        const gmme_id& first_corner,
        const gmme_id& second_corner )
    {
        for( const auto& line : geomodel_.lines() ) {
            gmme_id c0 { line.boundary_gmme( 0 ) };
            gmme_id c1 { line.boundary_gmme( 1 ) };

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) ) {
                std::vector< index_t > cur_adjacent_surfaces {
                    get_sorted_incident_surfaces( line ) };
                if( cur_adjacent_surfaces.size() == sorted_adjacent_surfaces.size()
                    && std::equal( cur_adjacent_surfaces.begin(),
                        cur_adjacent_surfaces.end(),
                        sorted_adjacent_surfaces.begin() ) ) {
                    return line.gmme();
                }
            }
        }
        return create_mesh_entity< Line >();
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::remove_mesh_entity_boundary_relation(
        const gmme_id& incident_entity,
        const gmme_id& boundary )
    {
        auto& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id { manager.find_boundary_relationship( incident_entity,
            boundary ) };
        if( relation_id == NO_ID ) {
            throw RINGMeshException( "Entity", "No boundary relation found between ",
                boundary, " and ", incident_entity );
        }
        GeoModelMeshEntityAccess< DIMENSION > boundary_access(
            geomodel_access_.modifiable_mesh_entity( boundary ) );
        auto& incident_entities = boundary_access.modifiable_incident_entities();
        std::remove_if( incident_entities.begin(), incident_entities.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
        GeoModelMeshEntityAccess< DIMENSION > incident_entity_access(
            geomodel_access_.modifiable_mesh_entity( incident_entity ) );
        auto& boundaries = incident_entity_access.modifiable_boundaries();
        std::remove_if( boundaries.begin(), boundaries.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
    }

    template< index_t DIMENSION >
    template< template< index_t > class ENTITY >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entity(
        const MeshType& mesh_type )
    {
        const auto& entity_type = ENTITY< DIMENSION >::type_name_static();
        index_t nb_entities { geomodel_.nb_mesh_entities( entity_type ) };
        index_t new_id { nb_entities };
        geomodel_access_.modifiable_mesh_entities( entity_type ).emplace_back(
            GeoModelMeshEntityAccess< DIMENSION >::template create_entity< ENTITY >(
                geomodel_, new_id, mesh_type ) );
        return geomodel_access_.modifiable_mesh_entities( entity_type ).back()->gmme();
    }

    template< index_t DIMENSION >
    template< template< index_t > class ENTITY >
    bool GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entities(
        index_t nb_additionnal_entities,
        const MeshType& type )
    {
        const auto& entity_type = ENTITY< DIMENSION >::type_name_static();
        auto& store = geomodel_access_.modifiable_mesh_entities( entity_type );
        index_t old_size { static_cast< index_t >( store.size() ) };
        index_t new_size { old_size + nb_additionnal_entities };
        store.reserve( new_size );
        for( auto i : range( old_size, new_size ) ) {
            store.emplace_back(
                GeoModelMeshEntityAccess< DIMENSION >::template create_entity< ENTITY >(
                    geomodel_, i, type ) );
        }
        return true;
    }

    template< index_t DIMENSION >
    index_t GeoModelBuilderTopologyBase< DIMENSION >::check_if_boundary_incident_entity_relation_already_exists(
        const gmme_id& incident_entity,
        const gmme_id& boundary )
    {
        const auto& incident_entity_mesh_entity = geomodel_.mesh_entity(
            incident_entity );
        for( auto in_ent : range(
            incident_entity_mesh_entity.nb_incident_entities() ) ) {
            if( incident_entity_mesh_entity.incident_entity_gmme( in_ent )
                == boundary ) {
                GeoModelMeshEntityConstAccess< DIMENSION > entity_access(
                    incident_entity_mesh_entity );
                return entity_access.incident_entity_relation_ids()[in_ent];
            }
        }
        return NO_ID;
    }
    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::add_mesh_entity_boundary_relation(
        const gmme_id& incident_entity_id,
        const gmme_id& boundary,
        bool side )
    {
        ringmesh_unused( side );
        const auto& incident_entity_type =
            geomodel_.entity_type_manager().mesh_entity_manager.incident_entity_type(
                boundary.type() );
        if( incident_entity_id.type() != incident_entity_type ) {
            throw RINGMeshException( "Entity",
                "Wrong incident entity type in the boundary relation between ",
                boundary, " and ", incident_entity_id );
        }
        const auto& boundary_type =
            geomodel_.entity_type_manager().mesh_entity_manager.boundary_entity_type(
                incident_entity_id.type() );
        if( boundary.type() != boundary_type ) {
            throw RINGMeshException( "Entity",
                "Wrong boundary type in the boundary relation between ", boundary,
                " and ", incident_entity_id );
        }
        index_t relation_id {
            check_if_boundary_incident_entity_relation_already_exists(
                incident_entity_id, boundary ) };
        auto& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        if( relation_id == NO_ID ) {
            relation_id = manager.add_boundary_relationship( incident_entity_id,
                boundary );
        }
        auto& boundary_entity = geomodel_access_.modifiable_mesh_entity( boundary );
        GeoModelMeshEntityAccess< DIMENSION > boundary_access( boundary_entity );
        boundary_access.modifiable_incident_entities().push_back( relation_id );
        auto& incident_entity = geomodel_access_.modifiable_mesh_entity(
            incident_entity_id );
        GeoModelMeshEntityAccess< DIMENSION > incident_entity_access(
            incident_entity );
        incident_entity_access.modifiable_boundaries().push_back( relation_id );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::set_mesh_entity_boundary(
        const gmme_id& gmme,
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_unused( side );
        ringmesh_assert( id < geomodel_.mesh_entity( gmme ).nb_boundaries() );
        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        const auto& b_type =
            geomodel_.entity_type_manager().mesh_entity_manager.boundary_entity_type(
                gmme.type() );
        gmme_id boundary( b_type, boundary_id );
        GeoModelMeshEntityAccess< DIMENSION > gme_access( mesh_entity );
        index_t relation_id { gme_access.modifiable_boundaries()[id] };
        auto& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        manager.set_boundary_to_boundary_relationship( relation_id, boundary );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::set_mesh_entity_incident_entity(
        const gmme_id& gmme,
        index_t id,
        index_t incident_entity_id )
    {
        /// No check on the validity of the index of the entity incident_entity
        /// NO_ID is used to flag entities to delete
        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        ringmesh_assert( id < mesh_entity.nb_incident_entities() );
        const auto& in_ent_type =
            geomodel_.entity_type_manager().mesh_entity_manager.incident_entity_type(
                gmme.type() );
        gmme_id incident_entity( in_ent_type, incident_entity_id );
        GeoModelMeshEntityAccess< DIMENSION > gme_access( mesh_entity );
        index_t relation_id { gme_access.modifiable_incident_entities()[id] };
        auto& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        manager.set_incident_entity_to_boundary_relationship( relation_id,
            incident_entity );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderTopologyBase< DIMENSION >::delete_mesh_entity(
        const MeshEntityType& type,
        index_t index )
    {
        geomodel_access_.modifiable_mesh_entities( type )[index].reset();
    }

    template< index_t DIMENSION >
    gmme_id GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entity(
        const MeshEntityType& type )
    {
        const auto& manager = geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( type ) ) {
            return this->create_mesh_entity< Corner >();
        } else if( manager.is_line( type ) ) {
            return create_mesh_entity< Line >();
        } else if( manager.is_surface( type ) ) {
            return create_mesh_entity< Surface >();
        } else {
            ringmesh_assert_not_reached;
            return gmme_id();
        }
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderTopologyBase< DIMENSION >::create_mesh_entities(
        const MeshEntityType& type,
        index_t nb_additional_entities )
    {
        const auto& manager = geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( type ) ) {
            return this->create_mesh_entities< Corner >( nb_additional_entities );
        } else if( manager.is_line( type ) ) {
            return create_mesh_entities< Line >( nb_additional_entities );
        } else if( manager.is_surface( type ) ) {
            return create_mesh_entities< Surface >( nb_additional_entities );
        } else {
            ringmesh_assert_not_reached;
            return false;
        }
    }

    void GeoModelBuilderTopology< 2 >::add_universe_boundary( index_t boundary_id,
    bool side )
    {
        gmme_id boundary( Line2D::type_name_static(), boundary_id );
        UniverseAccess2D universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries().push_back( boundary );
        universe_access.modifiable_sides().push_back( side );
    }

    void GeoModelBuilderTopology< 2 >::set_universe_boundary(
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.universe().nb_boundaries() );
        gmme_id boundary( Line2D::type_name_static(), boundary_id );
        UniverseAccess2D universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries()[id] = boundary;
        universe_access.modifiable_sides()[id] = side;
    }

    void GeoModelBuilderTopology< 3 >::add_universe_boundary( index_t boundary_id,
    bool side )
    {
        gmme_id boundary( Surface3D::type_name_static(), boundary_id );
        UniverseAccess3D universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries().push_back( boundary );
        universe_access.modifiable_sides().push_back( side );
    }

    void GeoModelBuilderTopology< 3 >::set_universe_boundary(
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.universe().nb_boundaries() );
        gmme_id boundary( Surface3D::type_name_static(), boundary_id );
        UniverseAccess3D universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries()[id] = boundary;
        universe_access.modifiable_sides()[id] = side;
    }

    void GeoModelBuilderTopology< 2 >::compute_universe()
    {
        if( geomodel_.universe().nb_boundaries() != 0 ) return;
        std::vector< bool > is_line_universe_boundary( geomodel_.nb_lines(), false );
        std::vector< bool > line_side( geomodel_.nb_lines() );
        for( const auto& surface : geomodel_.surfaces() ) {
            for( auto l : range( surface.nb_boundaries() ) ) {
                index_t line_id { surface.boundary_gmme( l ).index() };
                is_line_universe_boundary[line_id] =
                    !is_line_universe_boundary[line_id];
                line_side[line_id] = surface.side( l );
            }
        }

        for( const auto& line : geomodel_.lines() ) {
            if( is_line_universe_boundary[line.index()] ) {
                add_universe_boundary( line.index(), line_side[line.index()] );
            }
        }
    }

    void GeoModelBuilderTopology< 3 >::compute_universe()
    {
        if( geomodel_.universe().nb_boundaries() != 0 ) return;
        std::vector< bool > is_surface_universe_boundary( geomodel_.nb_surfaces(),
        false );
        std::vector< bool > surface_side( geomodel_.nb_surfaces() );
        for( const auto& region : geomodel_.regions() ) {
            for( auto s : range( region.nb_boundaries() ) ) {
                index_t surface_id { region.boundary_gmme( s ).index() };
                is_surface_universe_boundary[surface_id] =
                    !is_surface_universe_boundary[surface_id];
                surface_side[surface_id] = region.side( s );
            }
        }

        for( const auto& surface : geomodel_.surfaces() ) {
            if( is_surface_universe_boundary[surface.index()] ) {
                add_universe_boundary( surface.index(),
                    surface_side[surface.index()] );
            }
        }
    }

    gmme_id GeoModelBuilderTopology< 3 >::create_mesh_entity(
        const MeshEntityType& type )
    {
        const auto& manager = geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( type ) ) {
            return GeoModelBuilderTopologyBase3D::create_mesh_entity< Region >();
        } else {
            return GeoModelBuilderTopologyBase3D::create_mesh_entity( type );
        }
    }
    bool GeoModelBuilderTopology< 3 >::create_mesh_entities(
        const MeshEntityType& type,
        index_t nb_additional_entities )
    {
        const auto& manager = geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( type ) ) {
            return GeoModelBuilderTopologyBase3D::create_mesh_entities< Region >(
                nb_additional_entities );
        } else {
            return GeoModelBuilderTopologyBase3D::create_mesh_entities( type,
                nb_additional_entities );
        }
    }
    void GeoModelBuilderTopology< 3 >::copy_all_mesh_entity_topology(
        const GeoModel3D& from )
    {
        GeoModelBuilderTopologyBase3D::copy_all_mesh_entity_topology( from );
        copy_mesh_entity_topology < Region > ( from );
    }

    void GeoModelBuilderTopology< 2 >::set_mesh_entity_boundary(
        const gmme_id& gmme,
        index_t id,
        index_t boundary_id,
        bool side )
    {
        GeoModelBuilderTopologyBase2D::set_mesh_entity_boundary( gmme, id,
            boundary_id );

        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        GeoModelMeshEntityAccess2D gme_access( mesh_entity );
        if( gmme.type() == Surface2D::type_name_static() ) {
            gme_access.modifiable_sides()[id] = side;
        }
    }

    void GeoModelBuilderTopology< 2 >::add_mesh_entity_boundary_relation(
        const gmme_id& incident_entity_id,
        const gmme_id& boundary,
        bool side )
    {
        GeoModelBuilderTopologyBase2D::add_mesh_entity_boundary_relation(
            incident_entity_id, boundary );

        auto& incident_entity = geomodel_access_.modifiable_mesh_entity(
            incident_entity_id );
        GeoModelMeshEntityAccess2D incident_entity_access( incident_entity );
        if( incident_entity_id.type() == Surface2D::type_name_static() ) {
            incident_entity_access.modifiable_sides().push_back( side );
        }
    }

    void GeoModelBuilderTopology< 3 >::set_mesh_entity_boundary(
        const gmme_id& gmme,
        index_t id,
        index_t boundary_id,
        bool side )
    {
        GeoModelBuilderTopologyBase3D::set_mesh_entity_boundary( gmme, id,
            boundary_id );

        auto& mesh_entity = geomodel_access_.modifiable_mesh_entity( gmme );
        GeoModelMeshEntityAccess3D gme_access( mesh_entity );
        if( gmme.type() == Region3D::type_name_static() ) {
            gme_access.modifiable_sides()[id] = side;
        }
    }

    void GeoModelBuilderTopology< 3 >::add_mesh_entity_boundary_relation(
        const gmme_id& incident_entity_id,
        const gmme_id& boundary,
        bool side )
    {
        GeoModelBuilderTopologyBase3D::add_mesh_entity_boundary_relation(
            incident_entity_id, boundary );

        auto& incident_entity = geomodel_access_.modifiable_mesh_entity(
            incident_entity_id );
        GeoModelMeshEntityAccess3D incident_entity_access( incident_entity );
        if( incident_entity_id.type() == Region3D::type_name_static() ) {
            incident_entity_access.modifiable_sides().push_back( side );
        }
    }

    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 2 >::create_mesh_entity<
        Corner >( const MeshType& );
    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 2 >::create_mesh_entity<
        Line >( const MeshType& );
    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 2 >::create_mesh_entity<
        Surface >( const MeshType& );
    template class RINGMESH_API GeoModelBuilderTopologyBase< 2 > ;

    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 3 >::create_mesh_entity<
        Corner >( const MeshType& );
    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 3 >::create_mesh_entity<
        Line >( const MeshType& );
    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 3 >::create_mesh_entity<
        Surface >( const MeshType& );
    template gmme_id RINGMESH_API GeoModelBuilderTopologyBase< 3 >::create_mesh_entity<
        Region >( const MeshType& );
    template class RINGMESH_API GeoModelBuilderTopologyBase< 3 > ;
}
