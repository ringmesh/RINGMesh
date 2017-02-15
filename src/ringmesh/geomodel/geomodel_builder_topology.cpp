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
    using namespace RINGMesh ;

    gme_t find_corner( const GeoModel& geomodel, const vec3& point )
    {
        for( index_t i = 0; i < geomodel.nb_corners(); ++i ) {
            if( geomodel.corner( i ).vertex( 0 ) == point ) {
                return gme_t( Corner::type_name_static(), i ) ;
            }
        }
        return gme_t() ;
    }

    gme_t find_corner( const GeoModel& geomodel, index_t geomodel_point_id )
    {
        const GeoModelMeshVertices& geomodel_vertices = geomodel.mesh.vertices ;
        std::vector< GMEVertex > vertices ;
        geomodel_vertices.gme_vertices( geomodel_point_id, vertices ) ;
        for( index_t i = 0; i < vertices.size(); ++i ) {
            if( vertices[i].gme_id.type == Corner::type_name_static() ) {
                return vertices[i].gme_id ;
            }
        }
        return gme_t() ;
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
            return false ;
        }
        bool equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( i ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) return true ;

        equal = true ;
        for( index_t i = 0; i < L.nb_vertices(); i++ ) {
            if( rhs_vertices[i] != L.vertex( L.nb_vertices() - i - 1 ) ) {
                equal = false ;
                break ;
            }
        }
        return equal ;
    }

    void get_sorted_incident_surfaces(
        const GeoModelMeshEntity& E,
        std::vector< index_t >& incident_surfaces )
    {
        index_t nb = E.nb_in_boundary() ;
        incident_surfaces.resize( nb ) ;
        for( index_t i = 0; i < nb; ++i ) {
            incident_surfaces[i] = E.in_boundary_gme( i ).index ;
        }
        std::sort( incident_surfaces.begin(), incident_surfaces.end() ) ;
    }
}

namespace RINGMesh {

    GeoModelBuilderTopology::GeoModelBuilderTopology(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    bool GeoModelBuilderTopology::create_geological_entities(
        const EntityType& type,
        index_t nb_additional_entities )
    {
        find_or_create_geological_entity_type( type ) ;
        std::vector< GeoModelGeologicalEntity* >& store =
            geomodel_access_.modifiable_geological_entities( type ) ;
        index_t old_size = static_cast< index_t >( store.size() ) ;
        index_t new_size = old_size + nb_additional_entities ;
        store.resize( new_size, nil ) ;
        for( index_t i = old_size; i < new_size; i++ ) {
            ringmesh_assert( store[i] == nil ) ;
            store[i] = GeoModelGeologicalEntityAccess::create_geological_entity(
                type, geomodel_, i ) ;
        }
        return true ;
    }

    index_t GeoModelBuilderTopology::create_geological_entity_type(
        const EntityType& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory::has_creator( type ) ) ;

        geomodel_access_.modifiable_entity_type_manager().geological_entity_types_.push_back(
            type ) ;
        geomodel_access_.modifiable_geological_entities().push_back(
            std::vector< GeoModelGeologicalEntity* >() ) ;
        GeoModelGeologicalEntity * E = GeoModelGeologicalEntityFactory::create_object(
            type, geomodel_ ) ;

        const EntityType child_type = E->child_type_name() ;
        delete E ;
        EntityTypeManager& parentage =
            geomodel_access_.modifiable_entity_type_manager() ;
        parentage.register_relationship( type, child_type ) ;

        return geomodel_.entity_type_manager().nb_geological_entity_types() - 1 ;
    }

    void GeoModelBuilderTopology::copy_topology( const GeoModel& from )
    {

        copy_mesh_entity_topology< Corner >( from ) ;
        copy_mesh_entity_topology< Line >( from ) ;
        copy_mesh_entity_topology< Surface >( from ) ;
        copy_mesh_entity_topology< Region >( from ) ;

        for( index_t t = 0; t < from.nb_geological_entity_types(); t++ ) {
            copy_geological_entity_topology( from,
                from.geological_entity_type( t ) ) ;
        }

        UniverseAccess universe_access( geomodel_access_.modifiable_universe() ) ;
        universe_access.copy( from.universe() ) ;
        geomodel_access_.modifiable_epsilon() = from.epsilon() ;

    }

    bool GeoModelBuilderTopology::get_dependent_entities(
        std::set< gme_t >& in ) const
    {
        std::size_t input_size = in.size() ;

        // Add children of geological entities
        for( std::set< gme_t >::iterator it( in.begin() ); it != in.end(); ++it ) {
            gme_t cur = *it ;
            if( geomodel_.entity_type_manager().is_geological_entity_type(
                cur.type ) ) {
                const GeoModelGeologicalEntity& E = geomodel_.geological_entity(
                    cur ) ;
                for( index_t j = 0; j < E.nb_children(); ++j ) {
                    in.insert( E.child_gme( j ) ) ;
                }
            }
        }
        // Add geological entities which have no child
        index_t nb_geological_entity_types =
            geomodel_.entity_type_manager().nb_geological_entity_types() ;
        for( index_t i = 0; i < nb_geological_entity_types; ++i ) {
            const EntityType& type =
                geomodel_.entity_type_manager().geological_entity_type( i ) ;

            for( index_t j = 0; j < geomodel_.nb_geological_entities( type ); ++j ) {
                bool no_child = true ;
                const GeoModelGeologicalEntity& E = geomodel_.geological_entity(
                    type, j ) ;
                for( index_t k = 0; k < E.nb_children(); ++k ) {
                    if( in.count( E.child_gme( k ) ) == 0 ) {
                        no_child = false ;
                        break ;
                    }
                }
                if( no_child ) {
                    in.insert( E.gme_id() ) ;
                }
            }
        }
        // Add mesh entities that are in the boundary of no mesh entity
        for( index_t i = 0; i < EntityTypeManager::nb_mesh_entity_types(); ++i ) {
            const EntityType& type = EntityTypeManager::mesh_entity_types()[i] ;
            const EntityType& in_boundary_type = EntityTypeManager::in_boundary_type(
                type ) ;
            if( !EntityTypeManager::is_mesh_entity_type( in_boundary_type ) ) {
                continue ;
            } else {
                for( index_t j = 0; j < geomodel_.nb_mesh_entities( type ); ++j ) {
                    bool no_incident = true ;
                    const GeoModelMeshEntity& E = geomodel_.mesh_entity( type, j ) ;
                    for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
                        if( in.count( E.in_boundary_gme( k ) ) == 0 ) {
                            no_incident = false ;
                            break ;
                        }
                    }
                    if( no_incident ) {
                        in.insert( E.gme_id() ) ;
                    }
                }
            }
        }
        // Recursive call till nothing is added
        if( in.size() != input_size ) {
            return get_dependent_entities( in ) ;
        } else {
            return false ;
        }
    }

    gme_t GeoModelBuilderTopology::create_geological_entity( const EntityType& type )
    {
        index_t index = find_or_create_geological_entity_type( type ) ;
        index_t id =
            static_cast< index_t >( geomodel_.nb_geological_entities( type ) ) ;
        GeoModelGeologicalEntity* E =
            GeoModelGeologicalEntityAccess::create_geological_entity( type,
                geomodel_, id ) ;
        geomodel_access_.modifiable_geological_entities()[index].push_back( E ) ;
        return E->gme_id() ;
    }

    gme_t GeoModelBuilderTopology::find_or_create_corner( const vec3& point )
    {
        gme_t result = find_corner( geomodel_, point ) ;
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >() ;
            builder_.geometry.set_corner( result.index, point ) ;
        }
        return result ;
    }

    gme_t GeoModelBuilderTopology::find_or_create_corner( index_t geomodel_point_id )
    {
        gme_t result = find_corner( geomodel_, geomodel_point_id ) ;
        if( !result.is_defined() ) {
            result = create_mesh_entity< Corner >() ;
            builder_.geometry.set_corner( result.index, geomodel_point_id ) ;
        }
        return result ;
    }

    /*!
     * @brief Finds or creates a line
     * @param[in] vertices Coordinates of the vertices of the line
     * @return Index of the Line
     */
    gme_t GeoModelBuilderTopology::find_or_create_line(
        const std::vector< vec3 >& vertices )
    {
        gme_t result ;
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            if( line_equal( geomodel_.line( i ), vertices ) ) {
                result = geomodel_.line( i ).gme_id() ;
            }
        }
        if( !result.is_defined() ) {
            result = create_mesh_entity< Line >() ;
            builder_.geometry.set_line( result.index, vertices ) ;

            // Finds the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            add_mesh_entity_boundary( result,
                find_or_create_corner( vertices.front() ).index ) ;
            add_mesh_entity_boundary( result,
                find_or_create_corner( vertices.back() ).index ) ;
        }
        return result ;
    }

    /*!
     * @brief Finds or creates a line knowing its topological adjacencies
     */
    gme_t GeoModelBuilderTopology::find_or_create_line(
        const std::vector< index_t >& sorted_adjacent_surfaces,
        const gme_t& first_corner,
        const gme_t& second_corner )
    {
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            const Line& line = geomodel_.line( i ) ;
            gme_t c0 = line.boundary_gme( 0 ) ;
            gme_t c1 = line.boundary_gme( 1 ) ;

            if( ( c0 == first_corner && c1 == second_corner )
                || ( c0 == second_corner && c1 == first_corner ) ) {
                std::vector< index_t > cur_adjacent_surfaces ;
                get_sorted_incident_surfaces( line, cur_adjacent_surfaces ) ;
                if( cur_adjacent_surfaces.size() == sorted_adjacent_surfaces.size()
                    && std::equal( cur_adjacent_surfaces.begin(),
                        cur_adjacent_surfaces.end(),
                        sorted_adjacent_surfaces.begin() ) ) {
                    return line.gme_id() ;
                }
            }
        }
        return create_mesh_entity< Line >() ;
    }

    void GeoModelBuilderTopology::compute_universe()
    {
        if( geomodel_.universe().nb_boundaries() != 0 ) return ;
        std::vector< bool > is_surface_universe_boundary( geomodel_.nb_surfaces(),
            false ) ;
        std::vector< bool > surface_side( geomodel_.nb_surfaces() ) ;
        for( index_t r = 0; r < geomodel_.nb_regions(); r++ ) {
            const Region& region = geomodel_.region( r ) ;
            for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                index_t surface_id = region.boundary_gme( s ).index ;
                is_surface_universe_boundary[surface_id] =
                    !is_surface_universe_boundary[surface_id] ;
                surface_side[surface_id] = region.side( s ) ;
            }
        }

        for( index_t s = 0; s < geomodel_.nb_surfaces(); s++ ) {
            if( !is_surface_universe_boundary[s] ) continue ;
            add_universe_boundary( s, surface_side[s] ) ;
        }
    }

    void GeoModelBuilderTopology::complete_entity_connectivity()
    {
        // Order is important
        complete_mesh_entity_connectivity( Line::type_name_static() ) ;
        complete_mesh_entity_connectivity( Corner::type_name_static() ) ;
        complete_mesh_entity_connectivity( Surface::type_name_static() ) ;
        complete_mesh_entity_connectivity( Region::type_name_static() ) ;

        // Geological entities
        for( index_t i = 0; i < geomodel_.nb_geological_entity_types(); i++ ) {
            const EntityType& type = geomodel_.geological_entity_type( i ) ;
            if( geomodel_.nb_geological_entities( type ) > 0 ) {
                if( geomodel_.geological_entity( type, 0 ).nb_children() == 0 ) {
                    builder_.geology.fill_geological_entities_children( type ) ;
                }
            }
        }
    }

    void GeoModelBuilderTopology::fill_mesh_entities_boundaries(
        const EntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& b_type = geomodel_.entity_type_manager().boundary_type(
            type ) ;
        if( EntityTypeManager::is_defined_type( b_type ) ) {
            for( index_t i = 0; i < geomodel_.nb_mesh_entities( b_type ); ++i ) {
                const GeoModelMeshEntity& b =
                    geomodel_access_.modifiable_mesh_entity( gme_t( b_type, i ) ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_mesh_entity_boundary( b.in_boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelBuilderTopology::fill_mesh_entities_in_boundaries(
        const EntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& in_b_type =
            geomodel_.entity_type_manager().in_boundary_type( type ) ;
        if( EntityTypeManager::is_defined_type( in_b_type ) ) {
            for( index_t i = 0; i < geomodel_.nb_mesh_entities( in_b_type ); ++i ) {
                const GeoModelMeshEntity& in_b =
                    geomodel_access_.modifiable_mesh_entity(
                        gme_t( in_b_type, i ) ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_mesh_entity_in_boundary( in_b.boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelBuilderTopology::add_mesh_entity_boundary(
        const gme_t& gme_id,
        index_t boundary_id,
        bool side )
    {
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            gme_id ) ;
        const EntityType& b_type = geomodel_.entity_type_manager().boundary_type(
            gme_id.type ) ;
        gme_t boundary( b_type, boundary_id ) ;
        GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
        gme_access.modifiable_boundaries().push_back( boundary ) ;

        if( gme_id.type == Region::type_name_static() ) {
            gme_access.modifiable_sides().push_back( side ) ;
        }
    }

    void GeoModelBuilderTopology::set_mesh_entity_boundary(
        const gme_t& gme_id,
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.mesh_entity( gme_id ).nb_boundaries() ) ;
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            gme_id ) ;
        const EntityType& b_type = geomodel_.entity_type_manager().boundary_type(
            gme_id.type ) ;
        gme_t boundary( b_type, boundary_id ) ;
        GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
        gme_access.modifiable_boundaries()[id] = boundary ;

        if( gme_id.type == Region::type_name_static() ) {
            gme_access.modifiable_sides()[id] = side ;
        }
    }

    void GeoModelBuilderTopology::add_universe_boundary(
        index_t boundary_id,
        bool side )
    {
        gme_t boundary( Surface::type_name_static(), boundary_id ) ;
        UniverseAccess universe_access( geomodel_access_.modifiable_universe() ) ;
        universe_access.modifiable_boundaries().push_back( boundary ) ;
        universe_access.modifiable_sides().push_back( side ) ;
    }

    void GeoModelBuilderTopology::set_universe_boundary(
        index_t id,
        index_t boundary_id,
        bool side )
    {
        ringmesh_assert( id < geomodel_.universe().nb_boundaries() ) ;
        gme_t boundary( Surface::type_name_static(), boundary_id ) ;
        UniverseAccess universe_access( geomodel_access_.modifiable_universe() ) ;
        universe_access.modifiable_boundaries()[id] = boundary ;
        universe_access.modifiable_sides()[id] = side ;
    }

    void GeoModelBuilderTopology::add_mesh_entity_in_boundary(
        const gme_t& t,
        index_t in_boundary_id )
    {
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            t ) ;
        const EntityType& in_b_type =
            geomodel_.entity_type_manager().in_boundary_type( t.type ) ;
        gme_t in_boundary( in_b_type, in_boundary_id ) ;
        GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
        gme_access.modifiable_in_boundaries().push_back( in_boundary ) ;
    }

    void GeoModelBuilderTopology::set_mesh_entity_in_boundary(
        const gme_t& gme_id,
        index_t id,
        index_t in_boundary_id )
    {
        /// No check on the validity of the index of the entity in_boundary
        /// NO_ID is used to flag entities to delete
        GeoModelMeshEntity& mesh_entity = geomodel_access_.modifiable_mesh_entity(
            gme_id ) ;
        ringmesh_assert( id < mesh_entity.nb_in_boundary() ) ;
        const EntityType& in_b_type =
            geomodel_.entity_type_manager().in_boundary_type( gme_id.type ) ;
        gme_t in_boundary( in_b_type, in_boundary_id ) ;
        GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
        gme_access.modifiable_in_boundaries()[id] = in_boundary ;
    }

    void GeoModelBuilderTopology::delete_mesh_entity(
        const EntityType& type,
        index_t index )
    {
        std::vector< GeoModelMeshEntity* >& store =
            geomodel_access_.modifiable_mesh_entities( type ) ;
        delete store[index] ;
        store[index] = nil ;
    }
    void GeoModelBuilderTopology::delete_geological_entity(
        const EntityType& type,
        index_t index )
    {
        std::vector< GeoModelGeologicalEntity* >& store =
            geomodel_access_.modifiable_geological_entities( type ) ;
        delete store[index] ;
        store[index] = nil ;
    }

    index_t GeoModelBuilderTopology::find_or_create_geological_entity_type(
        const EntityType& type )
    {
        index_t type_index =
            geomodel_.entity_type_manager().geological_entity_type_index( type ) ;
        if( type_index == NO_ID ) {
            type_index = create_geological_entity_type( type ) ;
        }
        return type_index ;
    }

    void GeoModelBuilderTopology::complete_mesh_entity_connectivity(
        const EntityType& type )
    {
        if( geomodel_.nb_mesh_entities( type ) > 0 ) {
            const GeoModelMeshEntity& E = geomodel_.mesh_entity( type, 0 ) ;
            if( E.nb_boundaries() == 0 ) {
                fill_mesh_entities_boundaries( type ) ;
            }
            if( E.nb_in_boundary() == 0 ) {
                fill_mesh_entities_in_boundaries( type ) ;
            }
            if( E.nb_parents() == 0 ) {
                builder_.geology.fill_mesh_entities_parent( type ) ;
            }
        }
    }

    void GeoModelBuilderTopology::copy_geological_entity_topology(
        const GeoModel& from,
        const EntityType& type )
    {
        create_geological_entities( type, from.nb_geological_entities( type ) ) ;

        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < geomodel_.nb_geological_entities( type ); ++e ) {
            gme_t id( type, e ) ;
            GeoModelGeologicalEntityAccess gmge_access(
                geomodel_access_.modifiable_geological_entity( id ) ) ;
            gmge_access.copy( from.geological_entity( id ) ) ;
        }
    }

}
