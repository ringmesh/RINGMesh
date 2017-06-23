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

/*!
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geomodel/geomodel_mesh_entity.h>

#include <stack>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>
#include <ringmesh/geomodel/geomodel_validity.h>

namespace {
    using namespace RINGMesh;

    typedef GeoModelMeshEntity GMME;
    typedef GeoModelGeologicalEntity GMGE;

    /*!
     * @brief Checks that the geomodel vertex indices of @param E
     *       are in a valid range
     */
    bool check_range_model_vertex_ids( const GMME& E )
    {
        const GeoModelMeshVertices& geomodel_vertices = E.geomodel().mesh.vertices;
        /// Check that the stored geomodel vertex indices are in a valid range
        gmme_id id = E.gmme();
        for( index_t i = 0; i < E.nb_vertices(); ++i ) {
            if( geomodel_vertices.geomodel_vertex_id( id, i ) == NO_ID
                && geomodel_vertices.geomodel_vertex_id( id, i )
                    >= E.geomodel().mesh.vertices.nb() ) {
                Logger::warn( "GeoModelEntity", "Invalid geomodel vertex index in ",
                    id );
                return false;
            }
        }
        return true;
    }

    /*!
     * @brief Computes and returns the surface connected components
     * @details In Debug mode, the connected components are saved into 
     * an Attribute on surface polygons.
     */
    index_t compute_nb_surface_connected_components( const Surface& surface )
    {
        const index_t NO_COMPONENT = index_t( -1 );
        GEO::Attribute< index_t > component( surface.polygon_attribute_manager(),
            "component" );
        component.fill( NO_COMPONENT );
        index_t nb_components = 0;
        for( index_t polygon = 0; polygon < surface.nb_mesh_elements(); polygon++ ) {
            if( component[polygon] == NO_COMPONENT ) {
                std::stack< index_t > S;
                S.push( polygon );
                component[polygon] = nb_components;
                do {
                    index_t cur_polygon = S.top();
                    S.pop();
                    for( index_t edge = 0;
                        edge < surface.nb_mesh_element_vertices( cur_polygon );
                        edge++ ) {
                        index_t adj_polygon = surface.polygon_adjacent_index(
                            cur_polygon, edge );
                        if( adj_polygon != NO_ID
                            && component[adj_polygon] == NO_COMPONENT ) {
                            S.push( adj_polygon );
                            component[adj_polygon] = nb_components;
                        }
                    }
                } while( !S.empty() );
                nb_components++;
            }
        }
#ifndef RINGMESH_DEBUG
        component.destroy();
#endif
        return nb_components;
    }

    /*!
     * @brief Computes and returns the region connected components
     * @details In Debug mode, the connected components are saved into 
     * an Attribute on region cells.
     */
    index_t compute_nb_volume_connected_components( const Region& M )
    {
        const index_t NO_COMPONENT = index_t( -1 );
        GEO::Attribute< index_t > component( M.cell_attribute_manager(),
            "component" );
        component.fill( NO_COMPONENT );
        index_t nb_components = 0;
        for( index_t cell = 0; cell < M.nb_mesh_elements(); cell++ ) {
            if( component[cell] == NO_COMPONENT ) {
                std::stack< index_t > S;
                S.push( cell );
                component[cell] = nb_components;
                do {
                    index_t cur_cell = S.top();
                    S.pop();
                    for( index_t facet = 0; facet < M.nb_cell_facets( cur_cell );
                        facet++ ) {
                        index_t adj_cell = M.cell_adjacent_index( cur_cell, facet );
                        if( adj_cell != NO_ID
                            && component[adj_cell] == NO_COMPONENT ) {
                            S.push( adj_cell );
                            component[adj_cell] = nb_components;
                        }
                    }
                } while( !S.empty() );
                nb_components++;
            }
        }
#ifndef RINGMESH_DEBUG
        component.destroy();
#endif
        return nb_components;
    }

    /*!
     * @brief Count the number of times each vertex is in an edge or polygon
     *
     * @param[in] gmme The GeoModelMeshEntity
     * @param[out] nb Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an mesh_element collection of 
     *      the GeoModelMeshEntity edge or polygon of the mesh.
     */
    void count_vertex_occurences(
        const GeoModelMeshEntity& E,
        std::vector< index_t >& nb )
    {
        nb.resize( E.nb_vertices(), 0 );
        for( index_t mesh_element_index = 0;
            mesh_element_index < E.nb_mesh_elements(); ++mesh_element_index ) {
            for( index_t vertex = 0;
                vertex < E.nb_mesh_element_vertices( mesh_element_index );
                ++vertex ) {
                ++nb[E.mesh_element_vertex_index( mesh_element_index, vertex )];
            }
        }
    }

    index_t count_nb_isolated_vertices( const GeoModelMeshEntity& mesh )
    {
        std::vector< index_t > nb;
        count_vertex_occurences( mesh, nb );
        return static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) );
    }

    bool check_mesh_entity_vertices_are_different(
        std::vector< index_t >& vertices,
        std::vector< index_t >& vertices_global )
    {
        ringmesh_assert(
            std::count( vertices.begin(), vertices.end(), NO_ID ) == 0 );
        ringmesh_assert(
            std::count( vertices_global.begin(), vertices_global.end(), NO_ID )
                == 0 );
        // 0 is the default value of the geomodel_vertex_id
        // If we have only 0 either this is a degenerate polygons, but most certainly
        // geomodel vertex ids are not good
        ringmesh_assert(
            static_cast< index_t >( std::count( vertices_global.begin(),
                vertices_global.end(), 0 ) ) != vertices_global.size() );

        std::sort( vertices.begin(), vertices.end() );
        std::sort( vertices_global.begin(), vertices_global.end() );
        return std::unique( vertices.begin(), vertices.end() ) != vertices.end()
            || std::unique( vertices_global.begin(), vertices_global.end() )
                != vertices_global.end();
    }

    /*!
     * @brief Returns true if the surface polygon is incident twice to the same vertex
     */
    bool polygon_is_degenerate( const Surface& S, const gmme_id& id, index_t p )
    {
        index_t nb_polygon_vertices = S.nb_mesh_element_vertices( p );
        std::vector< index_t > corners( nb_polygon_vertices, NO_ID );
        std::vector< index_t > corners_global( nb_polygon_vertices, NO_ID );
        index_t v = 0;
        const GeoModelMeshVertices& geomodel_vertices = S.geomodel().mesh.vertices;
        for( index_t c = 0; c < S.nb_mesh_element_vertices( p ); ++c ) {
            index_t polygon_vertex_index = S.mesh_element_vertex_index( p, c );
            corners[v] = polygon_vertex_index;
            corners_global[v] = geomodel_vertices.geomodel_vertex_id( id, p, v );
            v++;
        }
        double area = S.mesh_element_size( p );
        return check_mesh_entity_vertices_are_different( corners, corners_global )
            || area < S.geomodel().epsilon2();
    }

    /*!
     * @brief Returns true if the region cell is incident twice to the same vertex
     * or if the cell volume is negative or inferior to epsilon
     */
    bool cell_is_degenerate( const Region& region, index_t cell_index )
    {
        index_t nb_vertices_in_cell = region.nb_mesh_element_vertices( cell_index );
        std::vector< index_t > vertices( nb_vertices_in_cell, NO_ID );
        std::vector< index_t > vertices_global( nb_vertices_in_cell, NO_ID );
        gmme_id id = region.gmme();
        const GeoModelMeshVertices& geomodel_vertices =
            region.geomodel().mesh.vertices;
        for( index_t v = 0; v < nb_vertices_in_cell; v++ ) {
            vertices[v] = region.mesh_element_vertex_index( cell_index, v );
            vertices_global[v] = geomodel_vertices.geomodel_vertex_id( id,
                cell_index, v );
        }
        double volume = region.mesh_element_size( cell_index );
        return check_mesh_entity_vertices_are_different( vertices, vertices_global )
            || volume < region.geomodel().epsilon3();
    }
}
/******************************************************************************/
namespace RINGMesh {
    bool GeoModelMeshEntity::is_inside_border( const GeoModelMeshEntity& rhs ) const
    {
        // Find out if this surface is twice in the incident_entity vector
        gmme_id rhs_id = rhs.gmme();
        const RelationshipManager& manager =
            geomodel().entity_type_manager().relationship_manager;
        return std::count_if( incident_entities_.begin(), incident_entities_.end(),
            [&rhs_id, &manager](index_t i) {return manager.incident_entity_gmme( i ) == rhs_id;} )
            > 1;
    }

    bool GeoModelMeshEntity::has_inside_border() const
    {
        for( index_t i = 0; i < nb_boundaries(); ++i ) {
            if( boundary( i ).is_inside_border( *this ) ) {
                return true;
            }
        }
        return false;
    }

    GeoModelMeshEntity::~GeoModelMeshEntity()
    {
#ifdef RINGMESH_DEBUG
        ringmesh_assert( mesh_ != nullptr );
        mesh_->print_mesh_bounded_attributes();
#endif
    }

    void GeoModelMeshEntity::unbind_vertex_mapping_attribute() const
    {
        GeoModel& modifiable_model = const_cast< GeoModel& >( geomodel() );
        modifiable_model.mesh.vertices.unbind_geomodel_vertex_map( gmme() );
    }

    void GeoModelMeshEntity::bind_vertex_mapping_attribute() const
    {
        GeoModel& modifiable_model = const_cast< GeoModel& >( geomodel() );
        modifiable_model.mesh.vertices.bind_geomodel_vertex_map( gmme() );
    }

    bool GeoModelMeshEntity::are_geomodel_vertex_indices_valid() const
    {
        bool valid = true;
        // For all vertices
        // Check that the global vertex has an index backward to 
        // the vertex of this entity
        const GeoModelMeshVertices& geomodel_vertices = geomodel().mesh.vertices;
        gmme_id id = gmme();
        for( index_t v = 0; v < nb_vertices(); ++v ) {
            index_t geomodel_v = geomodel_vertices.geomodel_vertex_id( id, v );

            if( geomodel_v == NO_ID ) {
                Logger::warn( "GeoModelEntity", id, " vertex ", v,
                    " is not mapped to the related global geomodel vertex indices." );
                valid = false;
            }

            std::vector< index_t > backward_vertices =
                geomodel_vertices.mesh_entity_vertex_id( id, geomodel_v );
            bool found_in_backward = false;
            for( index_t bv : backward_vertices ) {
                if( bv == v ) {
                    found_in_backward = true;
                }
            }
            if( !found_in_backward ) {
                Logger::warn( "GeoModelEntity", "Error in mapping of ", id,
                    " vertex ", v,
                    " to the related global geomodel vertex indices." );
                valid = false;
            }
        }
        return valid;
    }

    bool GeoModelMeshEntity::is_index_valid() const
    {
        return index() < geomodel().nb_mesh_entities( type_name() );
    }

    bool GeoModelMeshEntity::is_boundary_connectivity_valid() const
    {
        const MeshEntityTypeManager& family =
            geomodel().entity_type_manager().mesh_entity_manager;
        const MeshEntityType entity_type = type_name();
        const MeshEntityType& boundary_type = family.boundary_type( entity_type );

        bool valid = true;
        gmme_id id = gmme();
        if( family.is_valid_type( boundary_type ) ) {
            for( index_t i = 0; i < nb_boundaries(); ++i ) {
                const GeoModelMeshEntity& E = boundary( i );
                bool found = false;
                index_t j = 0;
                while( !found && j < E.nb_incident_entities() ) {
                    if( E.incident_entity_gmme( j ) == id ) {
                        found = true;
                    }
                    j++;
                }
                if( !found ) {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency boundary-incident_entity between ", id, " and ",
                        E.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    bool GeoModelMeshEntity::is_incident_entity_connectivity_valid() const
    {
        const MeshEntityTypeManager& family =
            geomodel().entity_type_manager().mesh_entity_manager;
        const MeshEntityType entity_type = type_name();
        const MeshEntityType& incident_entity_type = family.incident_entity_type(
            entity_type );

        bool valid = true;
        gmme_id id = gmme();
        if( family.is_valid_type( incident_entity_type ) ) {
            if( nb_incident_entities() == 0 ) {
                Logger::warn( "GeoModelEntity", id,
                    " is in the boundary of no entity " );
                valid = false;
            }
            for( index_t i = 0; i < nb_incident_entities(); ++i ) {
                const GeoModelMeshEntity& E = incident_entity( i );
                bool found = false;
                index_t j = 0;
                while( !found && j < E.nb_boundaries() ) {
                    if( E.boundary_gmme( j ) == id ) {
                        found = true;
                    }
                    j++;
                }
                if( !found ) {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency incident_entity-boundary between ", id, " and ",
                        E.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    bool GeoModelMeshEntity::is_parent_connectivity_valid() const
    {
        const RelationshipManager& family =
            geomodel().entity_type_manager().relationship_manager;
        const MeshEntityType entity_type = type_name();

        bool valid = true;
        const std::vector< GeologicalEntityType > parent_types = family.parent_types(
            entity_type );
        gmme_id id = gmme();
        for( const GeologicalEntityType& parent_type : parent_types ) {
            index_t nb_parent_entities_in_geomodel =
                geomodel_.nb_geological_entities( parent_type );
            if( nb_parent_entities_in_geomodel == 0 ) {
                continue;
            } else {
                // There must be one and only one parent of that type in this entity
                // And this parent must have this entity in its children
                index_t nb_found_parents = 0;
                for( index_t i = 0; i < nb_parents(); ++i ) {
                    const GeoModelGeologicalEntity& E = parent( i );
                    if( E.type_name() == parent_type ) {
                        nb_found_parents++;

                        // The parent must have this entity in its children
                        bool found = false;
                        index_t j = 0;
                        while( !found && j < E.nb_children() ) {
                            if( E.child_gmme( j ) == id ) {
                                found = true;
                            }
                            j++;
                        }
                        if( !found ) {
                            Logger::warn( "GeoModelEntity",
                                "Inconsistency parent-child between ", id, " and ",
                                E.gmge() );
                            valid = false;
                        }
                    }
                }
                if( nb_found_parents != 1 ) {
                    Logger::warn( "GeoModelEntity", id, " has ", nb_found_parents,
                        " geological parent entity of type ", parent_type,
                        " (expected one)" );
                    valid = false;
                }
            }
        }
        return valid;
    }

    bool GeoModelMeshEntity::is_connectivity_valid() const
    {
        return is_boundary_connectivity_valid()
            && is_incident_entity_connectivity_valid();
    }

    const GeoModelGeologicalEntity& GeoModelMeshEntity::parent(
        index_t parent_index ) const
    {
        gmge_id parent = parent_gmge( parent_index );
        ringmesh_assert( parent.is_defined() );
        return geomodel().geological_entity( parent );
    }

    const GeoModelGeologicalEntity& GeoModelMeshEntity::parent(
        const GeologicalEntityType& parent_type_name ) const
    {
        gmge_id id = parent_gmge( parent_type_name );
        ringmesh_assert( id.is_defined() );
        return geomodel().geological_entity( id );
    }

    gmge_id GeoModelMeshEntity::parent_gmge(
        const GeologicalEntityType& parent_type_name ) const
    {
        return defined_parent_gmge( parent_type_name );
    }

    gmge_id GeoModelMeshEntity::could_be_undefined_parent_gmge(
        const GeologicalEntityType& parent_type_name ) const
    {
        for( index_t i = 0; i < nb_parents(); ++i ) {
            if( parent_gmge( i ).type() == parent_type_name ) {
                return parent_gmge( i );
            }
        }
        return gmge_id( ForbiddenGeologicalEntityType::type_name_static(), NO_ID );
    }

    gmge_id GeoModelMeshEntity::defined_parent_gmge(
        const GeologicalEntityType& parent_type_name ) const
    {
        const gmge_id parent_gmge = could_be_undefined_parent_gmge( parent_type_name );
        ringmesh_assert(parent_gmge.is_defined());
        return parent_gmge;
    }

    const gmme_id& GeoModelMeshEntity::boundary_gmme( index_t x ) const
    {
        ringmesh_assert( x < nb_boundaries() );
        return geomodel().entity_type_manager().relationship_manager.boundary_gmme(
            boundaries_[x] );
    }

    const GeoModelMeshEntity& GeoModelMeshEntity::boundary( index_t x ) const
    {
        return geomodel().mesh_entity( boundary_gmme( x ) );
    }

    const GeoModelMeshEntity& GeoModelMeshEntity::incident_entity( index_t x ) const
    {
        return geomodel().mesh_entity( incident_entity_gmme( x ) );
    }

    const gmme_id& GeoModelMeshEntity::incident_entity_gmme( index_t x ) const
    {
        ringmesh_assert( x < nb_incident_entities() );
        return geomodel().entity_type_manager().relationship_manager.incident_entity_gmme(
            incident_entities_[x] );
    }

    const gmge_id& GeoModelMeshEntity::parent_gmge( index_t id ) const
    {
        ringmesh_assert( id < nb_parents() );
        return geomodel().entity_type_manager().relationship_manager.parent_of_gmme(
            parents_[id] );
    }
    /**************************************************************/

    bool Corner::is_on_voi() const
    {
        // True if one of the incident surface define the universe
        for( index_t i = 0; i < nb_incident_entities(); ++i ) {
            if( incident_entity( i ).is_on_voi() ) {
                return true;
            }
        }
        return false;

    }

    bool Corner::is_mesh_valid() const
    {
        bool valid = true;
        if( nb_vertices() != 1 ) {
            Logger::err( "GeoModelEntity", "Corner ", index(), " mesh has ",
                mesh0d_->nb_vertices(), " vertices " );
            valid = false;
        }
        return valid;
    }

    const Line& Corner::incident_entity( index_t x ) const
    {
        return static_cast< const Line& >( GeoModelMeshEntity::incident_entity( x ) );
    }

    /***************************************************************/

    bool Line::is_mesh_valid() const
    {
        bool valid = true;

        // Check that the GEO::Mesh has the expected entities
        if( nb_vertices() < 2 ) {
            Logger::err( "GeoModelEntity", "Line ", index(), " has ",
                mesh1d_->nb_vertices(), " vertices " );
            valid = false;
        }
        if( mesh1d_->nb_edges() == 0 ) {
            Logger::err( "GeoModelEntity", "Line ", index(), " mesh has ",
                mesh1d_->nb_edges(), " edges " );
            valid = false;
        }

        // Model indices must be valid
        valid = check_range_model_vertex_ids( *this ) && valid;

        if( nb_vertices() > 1 ) {
            // Count the number of edges in which each vertex is
            std::vector< index_t > nb;
            count_vertex_occurences( *this, nb );
            index_t nb0 = 0;
            index_t nb1 = 0;
            index_t nb2 = 0;
            for( index_t i : nb ) {
                if( i == 0 )
                    ++nb0;
                else if( i == 1 )
                    ++nb1;
                else if( i == 2 ) ++nb2;
            }

            // Vertices at extremitites must be in only one edge
            if( nb.front() != 1 || nb.back() != 1 ) {
                Logger::err( "GeoModelEntity", "Invalid extremity points in ",
                    gmme() );
                valid = false;
            }
            // No isolated vertices are allowed
            if( nb0 > 0 ) {
                Logger::warn( "GeoModelEntity", nb0, " isolated vertices in ",
                    gmme() );
                valid = false;
            }
            // Only the two extremities are in only 1 edge 
            // One connected component condition
            if( nb1 != 2 ) {
                Logger::warn( "GeoModelEntity",
                    "More than one connected component for ", gmme() );
                valid = false;
            }
            // All the others must be in 2 edges and 2 edges only
            // Manifold condition
            if( nb2 != nb.size() - 2 ) {
                Logger::warn( "GeoModelEntity", "Non-manifold entity", gmme() );
                valid = false;
            }
        }

        // No zero edge length
        index_t nb_degenerated = 0;
        for( index_t e = 0; e < nb_mesh_elements(); ++e ) {
            double l = length(
                mesh_element_vertex( e, 1 ) - mesh_element_vertex( e, 0 ) );
            if( l < geomodel().epsilon() ) {
                nb_degenerated++;
            }
        }
        if( nb_degenerated > 0 ) {
            Logger::warn( "GeoModelEntity", nb_degenerated, " degenerated edges in ",
                gmme() );
            valid = false;
        }
        if( !is_first_corner_first_vertex() ) {
            Logger::warn( "GeoModelEntity", "First and last vertex of Line", index(),
                " do not match first and second Corner respectively" );
            valid = false;
        }
        return valid;
    }

    bool Line::is_connectivity_valid() const
    {
        bool line_valid = GeoModelMeshEntity::is_connectivity_valid();

        // A Line must have 2 corners - they are identical if the Line is closed
        if( nb_boundaries() != 2 ) {
            Logger::warn( "Connectivity", gmme(), " does not have 2 corners" );
            line_valid = false;
        }
        return line_valid;
    }

    bool Line::is_first_corner_first_vertex() const
    {
        if( nb_boundaries() != 2 || nb_vertices() < 2 ) {
            return false;
        } else {
            // Geometric comparison - not great at all
            return boundary( 0 ).vertex( 0 ) == vertex( 0 )
                && boundary( 1 ).vertex( 0 ) == vertex( nb_vertices() - 1 );
        }
    }

    bool Line::is_on_voi() const
    {
        // True if one of the incident surface define the universe
        for( index_t i = 0; i < nb_incident_entities(); ++i ) {
            if( incident_entity( i ).is_on_voi() ) {
                return true;
            }
        }
        return false;
    }

    const Surface& Line::incident_entity( index_t x ) const
    {
        return static_cast< const Surface& >( GeoModelMeshEntity::incident_entity( x ) );
    }

    const Corner& Line::boundary( index_t x ) const
    {
        return static_cast< const Corner& >( GeoModelMeshEntity::boundary( x ) );
    }

    /********************************************************************/

    bool Surface::is_mesh_valid() const
    {
        bool valid = true;
        gmme_id id = gmme();
        // Check that the GEO::Mesh has the expected entities
        // at least 3 vertices and one polygon.
        if( nb_vertices() < 3 ) {
            Logger::warn( "GeoModelEntity", id, " has less than 3 vertices " );
            valid = false;
        }
        if( mesh2d_->nb_polygons() == 0 ) {
            Logger::warn( "GeoModelEntity", id, " has no polygons " );
            valid = false;
        }

        // No isolated vertices
        index_t nb_isolated_vertices = count_nb_isolated_vertices( *this );
        if( nb_isolated_vertices > 0 ) {
            Logger::warn( "GeoModelEntity", id, " mesh has ", nb_isolated_vertices,
                " isolated vertices " );
            valid = false;
        }

        // No zero area polygon
        // No polygon incident to the same vertex check local and global indices
        index_t nb_degenerate = 0;
        for( index_t p = 0; p < mesh2d_->nb_polygons(); p++ ) {
            if( polygon_is_degenerate( *this, id, p ) ) {
                nb_degenerate++;
            }
        }
        if( nb_degenerate != 0 ) {
            Logger::warn( "GeoModelEntity", id, " mesh has ", nb_degenerate,
                " degenerate polygons " );
            valid = false;
        }

        // No duplicated polygon
        GEO::vector< index_t > colocated;
        // GEO::mesh_detect_duplicated_facets( mesh_, colocated ) ; // not implemented yet 
        index_t nb_duplicated_p = 0;
        for( index_t p = 0; p < colocated.size(); ++p ) {
            if( colocated[p] != p ) {
                nb_duplicated_p++;
            }
        }
        if( nb_duplicated_p > 0 ) {
            Logger::warn( "GeoModelEntity", id, " mesh has ", nb_duplicated_p,
                " duplicated polygons " );
            valid = false;
        }

        // One connected component  
        index_t cc = compute_nb_surface_connected_components( *this );
        if( cc != 1 ) {
            Logger::warn( "GeoModelEntity", id, " mesh has ", cc,
                " connected components " );
            valid = false;
        }
        return valid;
    }

    bool Surface::is_on_voi() const
    {
        for( index_t i = 0; i < geomodel().universe().nb_boundaries(); ++i ) {
            if( geomodel().universe().boundary_gmme( i ) == gmme() ) {
                return true;
            }
        }
        return false;
    }

    const Line& Surface::boundary( index_t x ) const
    {
        return static_cast< const Line& >( GeoModelMeshEntity::boundary( x ) );
    }

    const Region& Surface::incident_entity( index_t x ) const
    {
        return static_cast< const Region& >( GeoModelMeshEntity::incident_entity( x ) );
    }

    /********************************************************************/

    const Surface& Region::boundary( index_t x ) const
    {
        return static_cast< const Surface& >( GeoModelMeshEntity::boundary( x ) );
    }

    bool Region::is_on_voi() const
    {
        return false;
    }

    bool Region::is_connectivity_valid() const
    {
        if( nb_boundaries() != sides_.size() ) {
            Logger::err( "GeoModelEntity", gmme(), " boundary sides are invalid " );
            return false;
        }
        bool region_valid = GeoModelMeshEntity::is_connectivity_valid();
        if( nb_boundaries() == 0 ) {
            Logger::warn( "Connectivity", gmme(), " has no boundaries " );
            region_valid = false;
        }
        return region_valid;
    }

    bool Region::is_mesh_valid() const
    {
        if( !is_meshed() ) {
            return true;
        } else {
            bool valid = true;
            // Check that the GEO::Mesh has the expected entities
            // at least 4 vertices and one cell.
            if( mesh3d_->nb_vertices() < 4 ) {
                Logger::warn( "GeoModelEntity", gmme(),
                    " has less than 4 vertices " );
                valid = false;
            }

            // No isolated vertices
            index_t nb_isolated_vertices = count_nb_isolated_vertices( *this );
            if( nb_isolated_vertices > 0 ) {
                Logger::warn( "GeoModelEntity", gmme(), " mesh has ",
                    nb_isolated_vertices, " isolated vertices " );
                valid = false;
            }

            // No cell with negative volume
            // No cell incident to the same vertex check local and global indices
            index_t nb_degenerate = 0;
            for( index_t c = 0; c < mesh3d_->nb_cells(); c++ ) {
                if( cell_is_degenerate( *this, c ) ) {
                    nb_degenerate++;
                }
            }
            if( nb_degenerate != 0 ) {
                Logger::warn( "GeoModelEntity", gmme(), " mesh has ", nb_degenerate,
                    " degenerate cells " );
                valid = false;
            }

            // One connected component
            index_t cc = compute_nb_volume_connected_components( *this );
            if( cc != 1 ) {
                Logger::warn( "GeoModelEntity", gmme(), " mesh has ", cc,
                    " connected components " );
                valid = false;
            }
            return valid;
        }
    }

    void Region::compute_region_volumes_per_cell_type(
        double& tet_volume,
        double& pyramid_volume,
        double& prism_volume,
        double& hex_volume,
        double& poly_volume ) const
    {
        for( index_t c = 0; c < nb_mesh_elements(); c++ ) {
            index_t nb_vertices = nb_mesh_element_vertices( c );
            double volume = mesh3d_->cell_volume( c );
            switch( nb_vertices ) {
                case 4:
                    tet_volume += volume;
                    break;
                case 5:
                    pyramid_volume += volume;
                    break;
                case 6:
                    prism_volume += volume;
                    break;
                case 8:
                    hex_volume += volume;
                    break;
                default:
                    poly_volume += volume;
                    break;
            }
        }
    }

    index_t Region::cells_around_vertex(
        index_t vertex_id,
        std::vector< index_t >& result,
        index_t cell_hint ) const
    {
        result.resize( 0 );

        if( cell_hint == NO_ID ) {
            return 0;
        }

        // Flag the visited cells
        std::vector< index_t > visited;
        visited.reserve( 10 );

        // Stack of the adjacent cells
        std::stack< index_t > S;
        S.push( cell_hint );
        visited.push_back( cell_hint );

        do {
            index_t c = S.top();
            S.pop();

            bool cell_includes_vertex = false;
            for( index_t v = 0; v < nb_mesh_element_vertices( c ); v++ ) {
                if( mesh_element_vertex_index( c, v ) == vertex_id ) {
                    result.push_back( c );
                    cell_includes_vertex = true;
                    break;
                }
            }
            if( !cell_includes_vertex ) {
                continue;
            }

            for( index_t f = 0; f < nb_cell_facets( c ); f++ ) {
                for( index_t v = 0; v < nb_cell_facet_vertices( c, f ); v++ ) {
                    index_t vertex = cell_facet_vertex_index( c, f, v );
                    if( vertex == vertex_id ) {
                        index_t adj_P = cell_adjacent_index( c, f );

                        if( adj_P != NO_ID ) {
                            if( !contains( visited, adj_P ) ) {
                                S.push( adj_P );
                                visited.push_back( adj_P );
                            }
                        }
                        break;
                    }
                }
            }
        } while( !S.empty() );

        return static_cast< index_t >( result.size() );
    }

    void GeoModelMeshEntityAccess::change_mesh_data_structure( const MeshType type )
    {
        if( gmme_.mesh_->type_name() != type ) {
            gmme_.unbind_vertex_mapping_attribute();
            gmme_.change_mesh_data_structure( type );
            gmme_.bind_vertex_mapping_attribute();
        }
    }

    void Corner::change_mesh_data_structure( const MeshType type )
    {
        std::unique_ptr< PointSetMesh > new_mesh = PointSetMesh::create_mesh( type );
        std::unique_ptr< PointSetMeshBuilder > builder =
            PointSetMeshBuilder::create_builder( *new_mesh );
        builder->copy( *mesh0d_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    void Line::change_mesh_data_structure( const MeshType type )
    {
        std::unique_ptr< LineMesh > new_mesh = LineMesh::create_mesh( type );
        std::unique_ptr< LineMeshBuilder > builder = LineMeshBuilder::create_builder(
            *new_mesh );
        builder->copy( *mesh1d_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    void Surface::change_mesh_data_structure( const MeshType type )
    {
        std::unique_ptr< SurfaceMesh > new_mesh = SurfaceMesh::create_mesh( type );
        std::unique_ptr< SurfaceMeshBuilder > builder =
            SurfaceMeshBuilder::create_builder( *new_mesh );
        builder->copy( *mesh2d_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    void Region::change_mesh_data_structure( const MeshType type )
    {
        std::unique_ptr< VolumeMesh > new_mesh = VolumeMesh::create_mesh( type );
        std::unique_ptr< VolumeMeshBuilder > builder =
            VolumeMeshBuilder::create_builder( *new_mesh );
        builder->copy( *mesh3d_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }
}
