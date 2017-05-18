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

#include <ringmesh/geomodel/geomodel_builder.h>

/*!
 * @file Implementation of GeoModel topological edition functions
 * @author Jeanne Pellerin
 * @author Pierre Anquez
 */

namespace {
    using namespace RINGMesh;

    gmme_id find_corner( const GeoModel& geomodel, const vec3& point )
    {
        for( index_t i = 0; i < geomodel.nb_corners(); ++i ) {
            if( geomodel.corner( i ).vertex( 0 ) == point ) {
                return gmme_id( Corner::type_name_static(), i );
            }
        }
        return gmme_id();
    }

    gmme_id find_corner( const GeoModel& geomodel, index_t geomodel_point_id )
    {
        const GeoModelMeshVertices& geomodel_vertices = geomodel.mesh.vertices;
        const std::vector< GMEVertex >& vertices = geomodel_vertices.gme_vertices(
            geomodel_point_id );
        for( const GMEVertex& vertex : vertices ) {
            if( vertex.gmme.type() == Corner::type_name_static() ) {
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
    bool line_equal( const Line& L, const std::vector< vec3 >& rhs_vertices )
    {
        if( L.nb_vertices() != rhs_vertices.size() ) {
            return false;
        }
        bool equal = true;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( i ) ) {
                equal = false;
                break;
            }
        }
        if( equal ) return true;

        equal = true;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( L.nb_vertices() - i - 1 ) ) {
                equal = false;
                break;
            }
        }
        return equal;
    }

    void get_sorted_incident_surfaces(
        const GeoModelMeshEntity& E,
        std::vector< index_t >& incident_surfaces )
    {
        index_t nb = E.nb_in_boundary();
        incident_surfaces.resize( nb );
        for( index_t i = 0; i < nb; ++i ) {
            incident_surfaces[i] = E.in_boundary_gmme( i ).index();
        }
        std::sort( incident_surfaces.begin(), incident_surfaces.end() );
    }
}

namespace RINGMesh {

    GeoModelBuilderTopology::GeoModelBuilderTopology(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    void GeoModelBuilderTopology::copy_topology( const GeoModel& from )
    {

        copy_mesh_entity_topology< Corner >( from );
        copy_mesh_entity_topology< Line >( from );
        copy_mesh_entity_topology< Surface >( from );
        copy_mesh_entity_topology< Region >( from );

        UniverseAccess universe_access( geomodel_access_.modifiable_universe() );
        universe_access.copy( from.universe() );
        geomodel_access_.modifiable_epsilon() = from.epsilon();
        geomodel_access_.modifiable_entity_type_manager().relationship_manager =
            from.entity_type_manager().relationship_manager;
    }

    bool GeoModelBuilderTopology::get_dependent_entities(
        std::set< gmme_id >& mesh_entities,
        std::set< gmge_id >& geological_entities ) const
    {
        std::size_t input_geological_size = geological_entities.size();
        std::size_t input_mesh_size = mesh_entities.size();

        // Add children of geological entities
        for( const gmge_id& cur : geological_entities ) {
            const GeoModelGeologicalEntity& E = geomodel_.geological_entity( cur );
            for( index_t j = 0; j < E.nb_children(); ++j ) {
                mesh_entities.insert( E.child_gmme( j ) );
            }
        }
        // Add geological entities which have no child
        index_t nb_geological_entity_types =
            geomodel_.entity_type_manager().geological_entity_manager.nb_geological_entity_types();
        for( index_t i = 0; i < nb_geological_entity_types; ++i ) {
            const GeologicalEntityType& type =
                geomodel_.entity_type_manager().geological_entity_manager.geological_entity_type(
                    i );

            for( index_t j = 0; j < geomodel_.nb_geological_entities( type ); ++j ) {
                bool no_child = true;
                const GeoModelGeologicalEntity& E = geomodel_.geological_entity(
                    type, j );
                for( index_t k = 0; k < E.nb_children(); ++k ) {
                    if( mesh_entities.count( E.child_gmme( k ) ) == 0 ) {
                        no_child = false;
                        break;
                    }
                }
                if( no_child ) {
                    geological_entities.insert( E.gmge() );
                }
            }
        }
        // Add mesh entities that are in the boundary of no mesh entity
        for( index_t i = 0; i < MeshEntityTypeManager::nb_mesh_entity_types();
            ++i ) {
            const MeshEntityType& type =
                MeshEntityTypeManager::mesh_entity_types()[i];
            for( index_t j = 0; j < geomodel_.nb_mesh_entities( type ); ++j ) {
                bool no_incident = true;
                const GeoModelMeshEntity& E = geomodel_.mesh_entity( type, j );
                for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
                    if( mesh_entities.count( E.in_boundary_gmme( k ) ) == 0 ) {
                        no_incident = false;
                        break;
                    }
                }
                if( no_incident ) {
                    mesh_entities.insert( E.gmme() );
                }
            }
        }
        // Recursive call till nothing is added
        if( mesh_entities.size() != input_mesh_size
            || geological_entities.size() != input_geological_size ) {
            return get_dependent_entities( mesh_entities, geological_entities );
        } else {
            return false;
        }
    }

    gmme_id GeoModelBuilderTopology::find_or_create_corner( const vec3& point )
    {
        gmme_id result = find_corner( geomodel_, point );
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >();
            builder_.geometry.set_corner( result.index(), point );
        }
        return result;
    }

    gmme_id GeoModelBuilderTopology::find_or_create_corner(
        index_t geomodel_point_id )
    {
        gmme_id result = find_corner( geomodel_, geomodel_point_id );
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >();
            builder_.geometry.set_corner( result.index(), geomodel_point_id );
        }
        return result;
    }

    gmme_id GeoModelBuilderTopology::find_or_create_line(
        const std::vector< vec3 >& vertices )
    {
        gmme_id result;
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            if( line_equal( geomodel_.line( i ), vertices ) ) {
                result = geomodel_.line( i ).gmme();
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

    gmme_id GeoModelBuilderTopology::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        const gmme_id& first_corner,
        const gmme_id& second_corner )
    {
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            const Line& line = geomodel_.line( i );
            gmme_id c0 = line.boundary_gmme( 0 );
            gmme_id c1 = line.boundary_gmme( 1 );

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) ) {
                std::vector< index_t > cur_adjacent_surfaces;
                get_sorted_incident_surfaces( line, cur_adjacent_surfaces );
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

    void GeoModelBuilderTopology::compute_universe()
    {
        if( geomodel_.universe().nb_boundaries() != 0 ) return;
        std::vector< bool > is_surface_universe_boundary( geomodel_.nb_surfaces(),
            false );
        std::vector< bool > surface_side( geomodel_.nb_surfaces() );
        for( index_t r = 0; r < geomodel_.nb_regions(); r++ ) {
            const Region& region = geomodel_.region( r );
            for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                index_t surface_id = region.boundary_gmme( s ).index();
                is_surface_universe_boundary[surface_id] =
                    !is_surface_universe_boundary[surface_id];
                surface_side[surface_id] = region.side( s );
            }
        }

        for( index_t s = 0; s < geomodel_.nb_surfaces(); s++ ) {
            if( !is_surface_universe_boundary[s] ) continue;
            add_universe_boundary( s, surface_side[s] );
        }
    }

    void GeoModelBuilderTopology::complete_entity_connectivity()
    {
        // Order is important
        complete_mesh_entity_connectivity( Line::type_name_static() );
        complete_mesh_entity_connectivity( Corner::type_name_static() );
        complete_mesh_entity_connectivity( Surface::type_name_static() );
        complete_mesh_entity_connectivity( Region::type_name_static() );

        // Geological entities
        for( index_t i = 0; i < geomodel_.nb_geological_entity_types(); i++ ) {
            const GeologicalEntityType& type = geomodel_.geological_entity_type( i );
            if( geomodel_.nb_geological_entities( type ) > 0 ) {
                if( geomodel_.geological_entity( type, 0 ).nb_children() == 0 ) {
                    builder_.geology.fill_geological_entities_children( type );
                }
            }
        }
    }

//    void GeoModelBuilderTopology::remove_mesh_entity_boundary_relation(
//        const gmme_id& in_boundary,
//        const gmme_id& boundary )
//    {
//
//    }

    void GeoModelBuilderTopology::add_mesh_entity_boundary_relation(
        const gmme_id& in_boundary,
        const gmme_id& boundary,
        bool side )
    {
        GeoModelMeshEntity& boundary_entity =
            geomodel_access_.modifiable_mesh_entity( boundary );
        const MeshEntityType& in_boundary_type =
            geomodel_.entity_type_manager().mesh_entity_manager.in_boundary_type(
                boundary.type() );
        if( in_boundary.type() != in_boundary_type ) {
            std::ostringstream message;
            message << "Wrong in boundary type in the boundary relation between "
                << boundary << " and " << in_boundary;
            throw RINGMeshException( "Entity", message.str() );
        }
        GeoModelMeshEntity& in_boundary_entity =
            geomodel_access_.modifiable_mesh_entity( in_boundary );
        const MeshEntityType& boundary_type =
            geomodel_.entity_type_manager().mesh_entity_manager.boundary_type(
                in_boundary.type() );
        if( boundary.type() != boundary_type ) {
            std::ostringstream message;
            message << "Wrong boundary type in the boundary relation between "
                << boundary << " and " << in_boundary;
            throw RINGMeshException( "Entity", message.str() );
        }
        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id = manager.add_boundary_relationship( in_boundary,
            boundary );
        GeoModelMeshEntityAccess boundary_access( boundary_entity );
        boundary_access.modifiable_in_boundaries().push_back( relation_id );
        GeoModelMeshEntityAccess in_boundary_access( in_boundary_entity );
        in_boundary_access.modifiable_boundaries().push_back( relation_id );

        if( in_boundary.type() == Region::type_name_static() ) {
            in_boundary_access.modifiable_sides().push_back( side );
        }
    }

    void GeoModelBuilderTopology::set_mesh_entity_boundary(
        const gmme_id& gmme,
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.mesh_entity( gmme ).nb_boundaries() );
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            gmme );
        const MeshEntityType& b_type =
            geomodel_.entity_type_manager().mesh_entity_manager.boundary_type(
                gmme.type() );
        gmme_id boundary( b_type, boundary_id );
        GeoModelMeshEntityAccess gme_access( mesh_entity );
        index_t relation_id = gme_access.modifiable_boundaries()[id];
        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        manager.set_boundary_to_boundary_relationship( relation_id, boundary );

        if( gmme.type() == Region::type_name_static() ) {
            gme_access.modifiable_sides()[id] = side;
        }
    }

    void GeoModelBuilderTopology::add_universe_boundary(
        index_t boundary_id,
        bool side )
    {
        gmme_id boundary( Surface::type_name_static(), boundary_id );
        UniverseAccess universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries().push_back( boundary );
        universe_access.modifiable_sides().push_back( side );
    }

    void GeoModelBuilderTopology::set_universe_boundary(
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.universe().nb_boundaries() );
        gmme_id boundary( Surface::type_name_static(), boundary_id );
        UniverseAccess universe_access( geomodel_access_.modifiable_universe() );
        universe_access.modifiable_boundaries()[id] = boundary;
        universe_access.modifiable_sides()[id] = side;
    }

    void GeoModelBuilderTopology::set_mesh_entity_in_boundary(
        const gmme_id& gmme,
        index_t id,
        index_t in_boundary_id )
    {
        /// No check on the validity of the index of the entity in_boundary
        /// NO_ID is used to flag entities to delete
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            gmme );
        ringmesh_assert( id < mesh_entity.nb_in_boundary() );
        const MeshEntityType& in_b_type =
            geomodel_.entity_type_manager().mesh_entity_manager.in_boundary_type(
                gmme.type() );
        gmme_id in_boundary( in_b_type, in_boundary_id );
        GeoModelMeshEntityAccess gme_access( mesh_entity );
        index_t relation_id = gme_access.modifiable_in_boundaries()[id];
        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        manager.set_in_boundary_to_boundary_relationship( relation_id, in_boundary );
    }

    void GeoModelBuilderTopology::delete_mesh_entity(
        const MeshEntityType& type,
        index_t index )
    {
        geomodel_access_.modifiable_mesh_entities( type )[index].reset();
    }

    void GeoModelBuilderTopology::complete_mesh_entity_connectivity(
        const MeshEntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) > 0 ) {
            const GeoModelMeshEntity& E = geomodel_.mesh_entity( type, 0 );
            if( E.nb_parents() == 0 ) {
                builder_.geology.fill_mesh_entities_parent( type );
            }
        }
    }
}
