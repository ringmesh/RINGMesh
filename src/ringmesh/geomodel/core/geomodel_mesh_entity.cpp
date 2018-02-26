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

/*!
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella
 */

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/entity_type.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

namespace
{
    using namespace RINGMesh;

    /*!
     * @brief Checks that the geomodel vertex indices of @param E
     *       are in a valid range
     */
    template < index_t DIMENSION >
    bool check_range_model_vertex_ids(
        const GeoModelMeshEntity< DIMENSION >& E )
    {
        const auto& geomodel_vertices = E.geomodel().mesh.vertices;
        /// Check that the stored geomodel vertex indices are in a valid range
        auto id = E.gmme();
        for( auto i : range( E.nb_vertices() ) )
        {
            if( geomodel_vertices.geomodel_vertex_id( id, i ) == NO_ID
                && geomodel_vertices.geomodel_vertex_id( id, i )
                       >= E.geomodel().mesh.vertices.nb() )
            {
                Logger::warn(
                    "GeoModelEntity", "Invalid geomodel vertex index in ", id );
                return false;
            }
        }
        return true;
    }

    /*!
     * @brief Computes and returns the surface connected components
     */
    template < index_t DIMENSION >
    index_t compute_nb_surface_connected_components(
        const SurfaceBase< DIMENSION >& surface )
    {
        index_t nb_connected_components;
        std::tie( nb_connected_components, std::ignore ) =
            surface.mesh().connected_components();
        return nb_connected_components;
    }

    /*!
     * @brief Computes and returns the region connected components
     */
    template < index_t DIMENSION >
    index_t compute_nb_volume_connected_components(
        const Region< DIMENSION >& region )
    {
        index_t nb_connected_components;
        std::tie( nb_connected_components, std::ignore ) =
            region.mesh().connected_components();
        return nb_connected_components;
    }

    /*!
     * @brief Count the number of times each vertex is in an edge or polygon
     *
     * @param[in] gmme The GeoModelMeshEntity
     * @return Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an mesh_element collection of
     *      the GeoModelMeshEntity edge or polygon of the mesh.
     */
    template < index_t DIMENSION >
    std::vector< index_t > count_vertex_occurences(
        const GeoModelMeshEntity< DIMENSION >& E )
    {
        std::vector< index_t > nb( E.nb_vertices(), 0 );
        for( auto mesh_element_index : range( E.nb_mesh_elements() ) )
        {
            for( auto vertex :
                range( E.nb_mesh_element_vertices( mesh_element_index ) ) )
            {
                ++nb[E.mesh_element_vertex_index(
                    { mesh_element_index, vertex } )];
            }
        }
        return nb;
    }

    template < index_t DIMENSION >
    bool polygon_is_degenerate(
        const SurfaceBase< DIMENSION >& surface, index_t polygon_id )
    {
        if( surface.mesh_element_size( polygon_id )
            < surface.geomodel().epsilon2() )
        {
            return true;
        }

        const auto& mesh = surface.mesh();
        auto min_length = surface.geomodel().epsilon();
        for( auto c : range( mesh.nb_polygon_vertices( polygon_id ) ) )
        {
            if( mesh.polygon_edge_length( { polygon_id, c } ) < min_length )
            {
                return false;
            }
        }
        return false;
    }

    template < index_t DIMENSION >
    bool cell_is_degenerate(
        const Region< DIMENSION >& region, index_t cell_id )
    {
        if( region.mesh_element_size( cell_id ) < region.geomodel().epsilon3() )
        {
            return true;
        }

        const auto& mesh = region.mesh();
        auto min_length = region.geomodel().epsilon();
        for( auto e : range( region.nb_cell_edges( cell_id ) ) )
        {
            if( mesh.cell_edge_length( cell_id, e ) < min_length )
            {
                return true;
            }
        }
        return false;
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::is_inside_border(
        const GeoModelMeshEntity& rhs ) const
    {
        // Find out if this surface is twice in the incident_entity vector
        auto rhs_id = rhs.gmme();
        const auto& manager =
            this->geomodel().entity_type_manager().relationship_manager;
        return std::count_if( incident_entities_.begin(),
                   incident_entities_.end(),
                   [&rhs_id, &manager]( index_t i ) {
                       return manager.incident_entity_gmme( i ) == rhs_id;
                   } )
               > 1;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::has_inside_border() const
    {
        for( auto i : range( nb_boundaries() ) )
        {
            if( boundary( i ).is_inside_border( *this ) )
            {
                return true;
            }
        }
        return false;
    }

    template < index_t DIMENSION >
    GeoModelMeshEntity< DIMENSION >::~GeoModelMeshEntity()
    {
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >&
        GeoModelMeshEntity< DIMENSION >::vertex_nn_search() const
    {
        return mesh_->vertex_nn_search();
    }

    template < index_t DIMENSION >
    GEO::AttributesManager&
        GeoModelMeshEntity< DIMENSION >::vertex_attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }

    template < index_t DIMENSION >
    void
        GeoModelMeshEntity< DIMENSION >::unbind_vertex_mapping_attribute() const
    {
        this->geomodel().mesh.vertices.unbind_geomodel_vertex_map( gmme() );
    }

    template < index_t DIMENSION >
    void GeoModelMeshEntity< DIMENSION >::bind_vertex_mapping_attribute() const
    {
        this->geomodel().mesh.vertices.bind_geomodel_vertex_map( gmme() );
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::are_geomodel_vertex_indices_valid()
        const
    {
        bool valid{ true };
        // For all vertices
        // Check that the global vertex has an index backward to
        // the vertex of this entity
        const auto& geomodel_vertices = this->geomodel().mesh.vertices;
        auto id = gmme();
        for( auto v : range( nb_vertices() ) )
        {
            index_t geomodel_v{ geomodel_vertices.geomodel_vertex_id( id, v ) };

            if( geomodel_v == NO_ID )
            {
                Logger::warn( "GeoModelEntity", id, " vertex ", v,
                    " is not mapped to the related global geomodel vertex "
                    "indices." );
                valid = false;
            }

            auto backward_vertices =
                geomodel_vertices.mesh_entity_vertex_id( id, geomodel_v );
            bool found_in_backward{ false };
            for( auto bv : backward_vertices )
            {
                if( bv == v )
                {
                    found_in_backward = true;
                }
            }
            if( !found_in_backward )
            {
                Logger::warn( "GeoModelEntity", "Error in mapping of ", id,
                    " vertex ", v,
                    " to the related global geomodel vertex indices." );
                valid = false;
            }
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::is_index_valid() const
    {
        return this->index() < this->geomodel().nb_mesh_entities( type_name() );
    }

    template < index_t DIMENSION >
    gmme_id GeoModelMeshEntity< DIMENSION >::gmme() const
    {
        return gmme_id( type_name(), this->index() );
    }

    template < index_t DIMENSION >
    MeshEntityType GeoModelMeshEntity< DIMENSION >::mesh_entity_type() const
    {
        return gmme().type();
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::is_boundary_connectivity_valid() const
    {
        const auto& family =
            this->geomodel().entity_type_manager().mesh_entity_manager;
        const auto entity_type = type_name();
        const auto& boundary_type = family.boundary_entity_type( entity_type );

        bool valid{ true };
        auto id = gmme();
        if( family.is_valid_type( boundary_type ) )
        {
            for( auto i : range( nb_boundaries() ) )
            {
                const auto& E = boundary( i );
                bool found{ false };
                index_t j{ 0 };
                while( !found && j < E.nb_incident_entities() )
                {
                    if( E.incident_entity_gmme( j ) == id )
                    {
                        found = true;
                    }
                    j++;
                }
                if( !found )
                {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency boundary-incident_entity between ", id,
                        " and ", E.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool
        GeoModelMeshEntity< DIMENSION >::is_incident_entity_connectivity_valid()
            const
    {
        const auto& family =
            this->geomodel().entity_type_manager().mesh_entity_manager;
        const auto entity_type = type_name();
        const auto& incident_entity_type =
            family.incident_entity_type( entity_type );

        bool valid{ true };
        auto id = gmme();
        if( family.is_valid_type( incident_entity_type ) )
        {
            if( nb_incident_entities() == 0 )
            {
                Logger::warn(
                    "GeoModelEntity", id, " is in the boundary of no entity " );
                valid = false;
            }
            for( auto i : range( nb_incident_entities() ) )
            {
                const auto& E = incident_entity( i );
                bool found{ false };
                index_t j{ 0 };
                while( !found && j < E.nb_boundaries() )
                {
                    if( E.boundary_gmme( j ) == id )
                    {
                        found = true;
                    }
                    j++;
                }
                if( !found )
                {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency incident_entity-boundary between ", id,
                        " and ", E.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::is_parent_connectivity_valid() const
    {
        const auto& family =
            this->geomodel().entity_type_manager().relationship_manager;
        const auto entity_type = type_name();

        bool valid{ true };
        const auto parent_types = family.parent_types( entity_type );
        auto id = gmme();
        for( const auto& parent_type : parent_types )
        {
            index_t nb_parent_entities_in_geomodel{
                this->geomodel_.nb_geological_entities( parent_type )
            };
            if( nb_parent_entities_in_geomodel == 0 )
            {
                continue;
            }
            // There must be one and only one parent of that type in this entity
            // And this parent must have this entity in its children
            index_t nb_found_parents{ 0 };
            for( auto i : range( nb_parents() ) )
            {
                const auto& E = parent( i );
                if( E.type_name() == parent_type )
                {
                    nb_found_parents++;

                    // The parent must have this entity in its children
                    bool found{ false };
                    index_t j{ 0 };
                    while( !found && j < E.nb_children() )
                    {
                        if( E.child_gmme( j ) == id )
                        {
                            found = true;
                        }
                        j++;
                    }
                    if( !found )
                    {
                        Logger::warn( "GeoModelEntity",
                            "Inconsistency parent-child between ", id, " and ",
                            E.gmge() );
                        valid = false;
                    }
                }
            }
            if( nb_found_parents != 1 )
            {
                Logger::warn( "GeoModelEntity", id, " has ", nb_found_parents,
                    " geological parent entity of type ", parent_type,
                    " (expected one)" );
                valid = false;
            }
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::is_connectivity_valid() const
    {
        return is_boundary_connectivity_valid()
               && is_incident_entity_connectivity_valid();
    }

    template < index_t DIMENSION >
    const GeoModelGeologicalEntity< DIMENSION >&
        GeoModelMeshEntity< DIMENSION >::parent( index_t parent_index ) const
    {
        auto parent = parent_gmge( parent_index );
        ringmesh_assert( parent.is_defined() );
        return this->geomodel().geological_entity( parent );
    }

    template < index_t DIMENSION >
    const GeoModelGeologicalEntity< DIMENSION >&
        GeoModelMeshEntity< DIMENSION >::parent(
            const GeologicalEntityType& parent_type ) const
    {
        auto id = parent_gmge( parent_type );
        ringmesh_assert( id.is_defined() );
        return this->geomodel().geological_entity( id );
    }

    template < index_t DIMENSION >
    void GeoModelMeshEntity< DIMENSION >::save(
        const std::string& filename ) const
    {
        mesh_->save_mesh( filename );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEntity< DIMENSION >::nb_vertices() const
    {
        return mesh_->nb_vertices();
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& GeoModelMeshEntity< DIMENSION >::vertex(
        index_t vertex_index ) const
    {
        return mesh_->vertex( vertex_index );
    }

    template < index_t DIMENSION >
    const gmge_id GeoModelMeshEntity< DIMENSION >::parent_gmge(
        const GeologicalEntityType& parent_type ) const
    {
        return defined_parent_gmge( parent_type );
    }

    template < index_t DIMENSION >
    gmge_id GeoModelMeshEntity< DIMENSION >::could_be_undefined_parent_gmge(
        const GeologicalEntityType& parent_type ) const
    {
        for( auto i : range( nb_parents() ) )
        {
            if( parent_gmge( i ).type() == parent_type )
            {
                return parent_gmge( i );
            }
        }
        return gmge_id(
            ForbiddenGeologicalEntityType::type_name_static(), NO_ID );
    }

    template < index_t DIMENSION >
    gmge_id GeoModelMeshEntity< DIMENSION >::defined_parent_gmge(
        const GeologicalEntityType& parent_type ) const
    {
        const auto parent_gmge = could_be_undefined_parent_gmge( parent_type );
        ringmesh_assert( parent_gmge.is_defined() );
        return parent_gmge;
    }

    template < index_t DIMENSION >
    const gmme_id& GeoModelMeshEntity< DIMENSION >::boundary_gmme(
        index_t x ) const
    {
        ringmesh_assert( x < nb_boundaries() );
        return this->geomodel()
            .entity_type_manager()
            .relationship_manager.boundary_gmme( boundaries_[x] );
    }

    template < index_t DIMENSION >
    const GeoModelMeshEntity< DIMENSION >&
        GeoModelMeshEntity< DIMENSION >::boundary( index_t x ) const
    {
        return this->geomodel().mesh_entity( boundary_gmme( x ) );
    }

    template < index_t DIMENSION >
    const GeoModelMeshEntity< DIMENSION >&
        GeoModelMeshEntity< DIMENSION >::incident_entity( index_t x ) const
    {
        return this->geomodel().mesh_entity( incident_entity_gmme( x ) );
    }

    template < index_t DIMENSION >
    const gmme_id& GeoModelMeshEntity< DIMENSION >::incident_entity_gmme(
        index_t x ) const
    {
        ringmesh_assert( x < this->nb_incident_entities() );
        return this->geomodel()
            .entity_type_manager()
            .relationship_manager.incident_entity_gmme( incident_entities_[x] );
    }

    template < index_t DIMENSION >
    const gmge_id& GeoModelMeshEntity< DIMENSION >::parent_gmge(
        index_t id ) const
    {
        ringmesh_assert( id < nb_parents() );
        return this->geomodel()
            .entity_type_manager()
            .relationship_manager.parent_of_gmme( parents_[id] );
    }

    template < index_t DIMENSION >
    bool GeoModelMeshEntity< DIMENSION >::has_parent(
        const GeologicalEntityType& parent_type ) const
    {
        return could_be_undefined_parent_gmge( parent_type ).is_defined();
    }

    /**************************************************************/

    template < index_t DIMENSION >
    Corner< DIMENSION >::Corner( const GeoModel< DIMENSION >& geomodel,
        index_t id,
        const MeshType& type )
        : GeoModelMeshEntity< DIMENSION >( geomodel, id )
    {
        update_mesh_storage_type(
            PointSetMesh< DIMENSION >::create_mesh( type ) );
    }

    template < index_t DIMENSION >
    void Corner< DIMENSION >::update_mesh_storage_type(
        std::unique_ptr< PointSetMesh< DIMENSION > > mesh )
    {
        point_set_mesh_ = std::move( mesh );
        GeoModelMeshEntity< DIMENSION >::set_mesh( point_set_mesh_ );
    }

    template < index_t DIMENSION >
    bool Corner< DIMENSION >::is_on_voi() const
    {
        // True if one of the incident lines defines the universe
        for( auto i : range( this->nb_incident_entities() ) )
        {
            if( incident_entity( i ).is_on_voi() )
            {
                return true;
            }
        }
        return false;
    }

    template < index_t DIMENSION >
    index_t Corner< DIMENSION >::nb_mesh_elements() const
    {
        return 1;
    }

    template < index_t DIMENSION >
    index_t Corner< DIMENSION >::nb_mesh_element_vertices(
        index_t mesh_element ) const
    {
        ringmesh_unused( mesh_element );
        index_t nb_vertices = point_set_mesh_->nb_vertices();
        ringmesh_assert( nb_vertices < 2 );
        return nb_vertices;
    }

    template < index_t DIMENSION >
    bool Corner< DIMENSION >::is_mesh_valid() const
    {
        bool valid{ true };
        if( this->nb_vertices() != 1 )
        {
            Logger::err( "GeoModelEntity", this->gmme(), " mesh has ",
                point_set_mesh_->nb_vertices(), " vertices " );
            valid = false;
        }
        if( !point_set_mesh_->is_mesh_valid() )
        {
            Logger::err( "GeoModelEntity", this->gmme(), " mesh is invalid" );
            valid = false;
        }
        return valid;
    }

    template < index_t DIMENSION >
    const Line< DIMENSION >& Corner< DIMENSION >::incident_entity(
        index_t x ) const
    {
        return static_cast< const Line< DIMENSION >& >(
            GeoModelMeshEntity< DIMENSION >::incident_entity( x ) );
    }

    template < index_t DIMENSION >
    MeshEntityType Corner< DIMENSION >::type_name_static()
    {
        return MeshEntityType( "Corner" );
    }

    template < index_t DIMENSION >
    MeshEntityType Corner< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    /***************************************************************/

    template < index_t DIMENSION >
    Line< DIMENSION >::Line( const GeoModel< DIMENSION >& geomodel,
        index_t id,
        const MeshType& type )
        : GeoModelMeshEntity< DIMENSION >( geomodel, id )
    {
        update_mesh_storage_type( LineMesh< DIMENSION >::create_mesh( type ) );
    }

    template < index_t DIMENSION >
    void Line< DIMENSION >::update_mesh_storage_type(
        std::unique_ptr< LineMesh< DIMENSION > > mesh )
    {
        line_mesh_ = std::move( mesh );
        GeoModelMeshEntity< DIMENSION >::set_mesh( line_mesh_ );
    }

    template < index_t DIMENSION >
    bool Line< DIMENSION >::is_mesh_valid() const
    {
        bool valid{ true };

        if( !line_mesh_->is_mesh_valid() )
        {
            Logger::err( "GeoModelEntity", this->gmme(), " mesh is invalid" );
            valid = false;
        }

        // Model indices must be valid
        valid = check_range_model_vertex_ids( *this ) && valid;

        if( this->nb_vertices() > 1 )
        {
            // Count the number of edges in which each vertex is
            auto nb = count_vertex_occurences( *this );
            index_t nb0{ 0 };
            index_t nb1{ 0 };
            index_t nb2{ 0 };
            for( auto i : nb )
            {
                if( i == 0 )
                {
                    ++nb0;
                }
                else if( i == 1 )
                {
                    ++nb1;
                }
                else if( i == 2 )
                {
                    ++nb2;
                }
            }

            // Vertices at extremitites must be in only one edge
            if( nb.front() != 1 || nb.back() != 1 )
            {
                Logger::err( "GeoModelEntity", "Invalid extremity points in ",
                    this->gmme() );
                valid = false;
            }
            // No isolated vertices are allowed
            if( nb0 > 0 )
            {
                Logger::warn( "GeoModelEntity", nb0, " isolated vertices in ",
                    this->gmme() );
                valid = false;
            }
            // Only the two extremities are in only 1 edge
            // One connected component condition
            if( nb1 != 2 )
            {
                Logger::warn( "GeoModelEntity",
                    "More than one connected component for ", this->gmme() );
                valid = false;
            }
            // All the others must be in 2 edges and 2 edges only
            // Manifold condition
            if( nb2 != nb.size() - 2 )
            {
                Logger::warn(
                    "GeoModelEntity", "Non-manifold entity", this->gmme() );
                valid = false;
            }
        }

        // No zero edge length
        index_t nb_degenerated{ 0 };
        for( auto e : range( nb_mesh_elements() ) )
        {
            double l = ( this->mesh_element_vertex( { e, 1 } )
                         - this->mesh_element_vertex( { e, 0 } ) )
                           .length();
            if( l < this->geomodel().epsilon() )
            {
                nb_degenerated++;
            }
        }
        if( nb_degenerated > 0 )
        {
            Logger::warn( "GeoModelEntity", nb_degenerated,
                " degenerated edges in ", this->gmme() );
            valid = false;
        }
        if( !is_first_corner_first_vertex() )
        {
            Logger::warn( "GeoModelEntity", "First and last vertex of Line",
                this->index(),
                " do not match first and second Corner respectively" );
            valid = false;
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool Line< DIMENSION >::is_connectivity_valid() const
    {
        bool line_valid{
            GeoModelMeshEntity< DIMENSION >::is_connectivity_valid()
        };

        // A Line must have 2 corners - they are identical if the Line is closed
        if( this->nb_boundaries() != 2 )
        {
            Logger::warn(
                "Connectivity", this->gmme(), " does not have 2 corners" );
            line_valid = false;
        }
        return line_valid;
    }

    template < index_t DIMENSION >
    const LineAABBTree< DIMENSION >& Line< DIMENSION >::edge_aabb() const
    {
        return line_mesh_->edge_aabb();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& Line< DIMENSION >::edge_nn_search() const
    {
        return line_mesh_->edge_nn_search();
    }

    template < index_t DIMENSION >
    index_t Line< DIMENSION >::nb_mesh_elements() const
    {
        return line_mesh_->nb_edges();
    }

    template < index_t DIMENSION >
    double Line< DIMENSION >::mesh_element_size( index_t edge_index ) const
    {
        ringmesh_assert( edge_index < nb_mesh_elements() );
        return line_mesh_->edge_length( edge_index );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > Line< DIMENSION >::mesh_element_barycenter(
        index_t edge_index ) const
    {
        ringmesh_assert( edge_index < nb_mesh_elements() );
        return line_mesh_->edge_barycenter( edge_index );
    }

    template < index_t DIMENSION >
    bool Line< DIMENSION >::is_first_corner_first_vertex() const
    {
        if( this->nb_boundaries() != 2 || this->nb_vertices() < 2 )
        {
            return false;
        }
        // Geometric comparison - not great at all
        return this->boundary( 0 ).vertex( 0 ) == this->vertex( 0 )
               && this->boundary( 1 ).vertex( 0 )
                      == this->vertex( this->nb_vertices() - 1 );
    }

    template < index_t DIMENSION >
    const Corner< DIMENSION >& Line< DIMENSION >::boundary( index_t x ) const
    {
        return static_cast< const Corner< DIMENSION >& >(
            GeoModelMeshEntity< DIMENSION >::boundary( x ) );
    }

    template < index_t DIMENSION >
    index_t Line< DIMENSION >::mesh_element_vertex_index(
        const ElementLocalVertex& element_local_vertex ) const
    {
        ringmesh_assert( element_local_vertex.element_id < nb_mesh_elements() );
        ringmesh_assert( element_local_vertex.local_vertex_id < 2 );
        return line_mesh_->edge_vertex( element_local_vertex );
    }

    template < index_t DIMENSION >
    MeshEntityType Line< DIMENSION >::type_name_static()
    {
        return MeshEntityType( "Line" );
    }

    template < index_t DIMENSION >
    MeshEntityType Line< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template <>
    bool geomodel_core_api Line< 2 >::is_on_voi() const
    {
        ringmesh_assert( this->nb_incident_entities() == 1
                         || this->nb_incident_entities() == 2 );
        return this->nb_incident_entities() == 1;
    }

    template <>
    bool geomodel_core_api Line< 3 >::is_on_voi() const
    {
        // True if one of the incident surfaces defines the universe
        for( auto i : range( this->nb_incident_entities() ) )
        {
            if( this->incident_entity( i ).is_on_voi() )
            {
                return true;
            }
        }
        return false;
    }

    /********************************************************************/

    template < index_t DIMENSION >
    SurfaceBase< DIMENSION >::SurfaceBase(
        const GeoModel< DIMENSION >& geomodel,
        index_t id,
        const MeshType& type )
        : GeoModelMeshEntity< DIMENSION >( geomodel, id )
    {
        update_mesh_storage_type(
            SurfaceMesh< DIMENSION >::create_mesh( type ) );
    }

    template < index_t DIMENSION >
    void SurfaceBase< DIMENSION >::update_mesh_storage_type(
        std::unique_ptr< SurfaceMesh< DIMENSION > > mesh )
    {
        surface_mesh_ = std::move( mesh );
        GeoModelMeshEntity< DIMENSION >::set_mesh( surface_mesh_ );
    }

    template < index_t DIMENSION >
    index_t SurfaceBase< DIMENSION >::nb_mesh_elements() const
    {
        return surface_mesh_->nb_polygons();
    }

    template < index_t DIMENSION >
    bool SurfaceBase< DIMENSION >::is_simplicial() const
    {
        return surface_mesh_->polygons_are_simplicies();
    }

    template < index_t DIMENSION >
    const SurfaceAABBTree< DIMENSION >&
        SurfaceBase< DIMENSION >::polygon_aabb() const
    {
        return surface_mesh_->polygon_aabb();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >&
        SurfaceBase< DIMENSION >::polygon_nn_search() const
    {
        return surface_mesh_->polygon_nn_search();
    }

    template < index_t DIMENSION >
    GEO::AttributesManager&
        SurfaceBase< DIMENSION >::polygon_attribute_manager() const
    {
        return surface_mesh_->polygon_attribute_manager();
    }

    template < index_t DIMENSION >
    index_t SurfaceBase< DIMENSION >::nb_mesh_element_vertices(
        index_t polygon_index ) const
    {
        ringmesh_assert( polygon_index < nb_mesh_elements() );
        return surface_mesh_->nb_polygon_vertices( polygon_index );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > SurfaceBase< DIMENSION >::mesh_element_barycenter(
        index_t polygon_index ) const
    {
        ringmesh_assert( polygon_index < nb_mesh_elements() );
        return surface_mesh_->polygon_barycenter( polygon_index );
    }

    template < index_t DIMENSION >
    double SurfaceBase< DIMENSION >::mesh_element_size(
        index_t polygon_index ) const
    {
        ringmesh_assert( polygon_index < nb_mesh_elements() );
        return surface_mesh_->polygon_area( polygon_index );
    }

    template < index_t DIMENSION >
    index_t SurfaceBase< DIMENSION >::polygon_adjacent_index(
        const PolygonLocalEdge& polygon_local_edge ) const
    {
        ringmesh_assert( polygon_local_edge.polygon_id < nb_mesh_elements() );
        ringmesh_assert(
            polygon_local_edge.local_edge_id
            < nb_mesh_element_vertices( polygon_local_edge.polygon_id ) );
        return surface_mesh_->polygon_adjacent( polygon_local_edge );
    }

    template < index_t DIMENSION >
    index_t SurfaceBase< DIMENSION >::mesh_element_vertex_index(
        const ElementLocalVertex& element_local_vertex ) const
    {
        ringmesh_assert( element_local_vertex.element_id < nb_mesh_elements() );
        ringmesh_assert(
            element_local_vertex.local_vertex_id
            < nb_mesh_element_vertices( element_local_vertex.element_id ) );
        return surface_mesh_->polygon_vertex( element_local_vertex );
    }

    template < index_t DIMENSION >
    bool SurfaceBase< DIMENSION >::is_mesh_valid_base() const
    {
        bool valid{ true };
        auto id = this->gmme();

        if( !surface_mesh_->is_mesh_valid() )
        {
            Logger::err( "GeoModelEntity", this->gmme(), " mesh is invalid" );
            valid = false;
        }

        // No zero area polygon
        // No polygon incident to the same vertex check local and global indices
        index_t nb_degenerate{ 0 };
        for( auto p : range( surface_mesh_->nb_polygons() ) )
        {
            if( polygon_is_degenerate( *this, p ) )
            {
                nb_degenerate++;
            }
        }
        if( nb_degenerate != 0 )
        {
            Logger::warn( "GeoModelEntity", id, " mesh has ", nb_degenerate,
                " degenerate polygons " );
            valid = false;
        }

        // One connected component
        index_t cc{ compute_nb_surface_connected_components( *this ) };
        if( cc != 1 )
        {
            Logger::warn( "GeoModelEntity", id, " mesh has ", cc,
                " connected components " );
            valid = false;
        }
        return valid;
    }

    template < index_t DIMENSION >
    const Line< DIMENSION >& SurfaceBase< DIMENSION >::boundary(
        index_t x ) const
    {
        return static_cast< const Line< DIMENSION >& >(
            GeoModelMeshEntity< DIMENSION >::boundary( x ) );
    }

    template < index_t DIMENSION >
    MeshEntityType SurfaceBase< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template < index_t DIMENSION >
    MeshEntityType SurfaceBase< DIMENSION >::type_name_static()
    {
        return MeshEntityType( "Surface" );
    }

    bool Surface< 2 >::is_on_voi() const
    {
        return false;
    }

    bool Surface< 2 >::is_mesh_valid() const
    {
        if( !is_meshed() )
        {
            return true;
        }
        return is_mesh_valid_base();
    }

    bool Surface< 3 >::is_on_voi() const
    {
        ringmesh_assert( this->nb_incident_entities() == 1
                         || this->nb_incident_entities() == 2 );
        return this->nb_incident_entities() == 1;
    }

    bool Surface< 3 >::is_mesh_valid() const
    {
        return is_mesh_valid_base();
    }

    const Region3D& Surface< 3 >::incident_entity( index_t x ) const
    {
        return static_cast< const Region3D& >(
            GeoModelMeshEntity3D::incident_entity( x ) );
    }

    bool Surface< 2 >::is_meshed() const
    {
        return mesh().nb_polygons() > 0;
    }

    /********************************************************************/

    template < index_t DIMENSION >
    Region< DIMENSION >::Region( const GeoModel< DIMENSION >& geomodel,
        index_t id,
        const MeshType& type )
        : GeoModelMeshEntity< DIMENSION >( geomodel, id )
    {
        update_mesh_storage_type(
            VolumeMesh< DIMENSION >::create_mesh( type ) );
    }

    template < index_t DIMENSION >
    void Region< DIMENSION >::update_mesh_storage_type(
        std::unique_ptr< VolumeMesh< DIMENSION > > mesh )
    {
        volume_mesh_ = std::move( mesh );
        GeoModelMeshEntity< DIMENSION >::set_mesh( volume_mesh_ );
    }

    template < index_t DIMENSION >
    const Surface< DIMENSION >& Region< DIMENSION >::boundary( index_t x ) const
    {
        return static_cast< const Surface< DIMENSION >& >(
            GeoModelMeshEntity< DIMENSION >::boundary( x ) );
    }

    template < index_t DIMENSION >
    bool Region< DIMENSION >::is_on_voi() const
    {
        return false;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::mesh_element_vertex_index(
        const ElementLocalVertex& element_local_vertex ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert(
                element_local_vertex.element_id < nb_mesh_elements() );
            ringmesh_assert(
                element_local_vertex.local_vertex_id
                < nb_mesh_element_vertices( element_local_vertex.element_id ) );
            return volume_mesh_->cell_vertex( element_local_vertex );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    bool Region< DIMENSION >::is_connectivity_valid() const
    {
        if( this->nb_boundaries() != sides_.size() )
        {
            Logger::err( "GeoModelEntity", this->gmme(),
                " boundary sides are invalid " );
            return false;
        }
        bool region_valid{
            GeoModelMeshEntity< DIMENSION >::is_connectivity_valid()
        };
        if( this->nb_boundaries() == 0 )
        {
            Logger::warn( "Connectivity", this->gmme(), " has no boundaries " );
            region_valid = false;
        }
        return region_valid;
    }

    template < index_t DIMENSION >
    bool Region< DIMENSION >::is_meshed() const
    {
        return volume_mesh_->nb_cells() > 0;
    }

    template < index_t DIMENSION >
    bool Region< DIMENSION >::is_simplicial() const
    {
        return volume_mesh_->cells_are_simplicies();
    }

    template < index_t DIMENSION >
    const VolumeAABBTree< DIMENSION >& Region< DIMENSION >::cell_aabb() const
    {
        return volume_mesh_->cell_aabb();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& Region< DIMENSION >::cell_nn_search() const
    {
        return volume_mesh_->cell_nn_search();
    }

    template < index_t DIMENSION >
    GEO::AttributesManager& Region< DIMENSION >::cell_attribute_manager() const
    {
        return volume_mesh_->cell_attribute_manager();
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::nb_mesh_elements() const
    {
        return volume_mesh_->nb_cells();
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::nb_mesh_element_vertices(
        index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->nb_cell_vertices( cell_index );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    CellType Region< DIMENSION >::cell_type( index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->cell_type( cell_index );
        }
        ringmesh_assert_not_reached;
        return CellType::UNDEFINED;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::nb_cell_edges( index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->nb_cell_edges( cell_index );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::nb_cell_facets( index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->nb_cell_facets( cell_index );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::nb_cell_facet_vertices(
        index_t cell_index, index_t facet_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
            return volume_mesh_->nb_cell_facet_vertices(
                { cell_index, facet_index } );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::cell_edge_vertex_index(
        index_t cell_index, index_t edge_index, index_t vertex_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            ringmesh_assert( edge_index < nb_cell_edges( cell_index ) );
            ringmesh_assert(
                vertex_index < nb_mesh_element_vertices( cell_index ) );
            return volume_mesh_->cell_edge_vertex(
                cell_index, edge_index, vertex_index );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::cell_facet_vertex_index(
        index_t cell_index, index_t facet_index, index_t vertex_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
            ringmesh_assert(
                vertex_index < nb_mesh_element_vertices( cell_index ) );
            return volume_mesh_->cell_facet_vertex(
                { cell_index, facet_index }, vertex_index );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    index_t Region< DIMENSION >::cell_adjacent_index(
        index_t cell_index, index_t facet_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            ringmesh_assert( facet_index < nb_cell_facets( cell_index ) );
            return volume_mesh_->cell_adjacent( { cell_index, facet_index } );
        }
        ringmesh_assert_not_reached;
        return NO_ID;
    }

    template < index_t DIMENSION >
    double Region< DIMENSION >::mesh_element_size( index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->cell_volume( cell_index );
        }
        ringmesh_assert_not_reached;
        return 0;
    }

    template < index_t DIMENSION >
    double Region< DIMENSION >::size() const
    {
        double result = 0.;
        for( auto i : range( this->nb_boundaries() ) )
        {
            const Surface< DIMENSION >& surface = boundary( i );
            for( auto t : range( surface.nb_mesh_elements() ) )
            {
                const vecn< DIMENSION >& p0 =
                    surface.mesh_element_vertex( ElementLocalVertex( t, 0 ) );
                for( auto v :
                    range( 1, surface.nb_mesh_element_vertices( t ) - 1 ) )
                {
                    double cur_volume =
                        ( dot(
                            p0, cross( surface.mesh_element_vertex(
                                           ElementLocalVertex( t, v ) ),
                                    surface.mesh_element_vertex(
                                        ElementLocalVertex( t, v + 1 ) ) ) ) )
                        / 6.;
                    side( i ) ? result -= cur_volume : result += cur_volume;
                }
            }
        }
        return std::fabs( result );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > Region< DIMENSION >::mesh_element_barycenter(
        index_t cell_index ) const
    {
        if( is_meshed() )
        {
            ringmesh_assert( cell_index < nb_mesh_elements() );
            return volume_mesh_->cell_barycenter( cell_index );
        }
        ringmesh_assert_not_reached;
        return vecn< DIMENSION >();
    }

    template < index_t DIMENSION >
    std::vector< index_t > Region< DIMENSION >::cells_around_vertex(
        index_t vertex_id, index_t cell_hint ) const
    {
        return volume_mesh_->cells_around_vertex( vertex_id, cell_hint );
    }

    template < index_t DIMENSION >
    MeshEntityType Region< DIMENSION >::type_name_static()
    {
        return MeshEntityType( "Region" );
    }

    template < index_t DIMENSION >
    MeshEntityType Region< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template < index_t DIMENSION >
    bool Region< DIMENSION >::is_mesh_valid() const
    {
        if( !is_meshed() )
        {
            return true;
        }
        bool valid{ true };

        if( !volume_mesh_->is_mesh_valid() )
        {
            Logger::err( "GeoModelEntity", this->gmme(), " mesh is invalid" );
            valid = false;
        }

        // No cell with negative volume
        // No cell incident to the same vertex check local and global indices
        index_t nb_degenerate{ 0 };
        for( auto c : range( volume_mesh_->nb_cells() ) )
        {
            if( cell_is_degenerate( *this, c ) )
            {
                nb_degenerate++;
            }
        }
        if( nb_degenerate != 0 )
        {
            Logger::warn( "GeoModelEntity", this->gmme(), " mesh has ",
                nb_degenerate, " degenerate cells " );
            valid = false;
        }

        // One connected component
        index_t cc{ compute_nb_volume_connected_components( *this ) };
        if( cc != 1 )
        {
            Logger::warn( "GeoModelEntity", this->gmme(), " mesh has ", cc,
                " connected components " );
            valid = false;
        }
        return valid;
    }

    template < index_t DIMENSION >
    void Corner< DIMENSION >::change_mesh_data_structure( const MeshType& type )
    {
        auto new_mesh = PointSetMesh< DIMENSION >::create_mesh( type );
        auto builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *new_mesh );
        builder->copy( *point_set_mesh_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    template < index_t DIMENSION >
    void Line< DIMENSION >::change_mesh_data_structure( const MeshType& type )
    {
        auto new_mesh = LineMesh< DIMENSION >::create_mesh( type );
        auto builder =
            LineMeshBuilder< DIMENSION >::create_builder( *new_mesh );
        builder->copy( *line_mesh_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    template < index_t DIMENSION >
    void SurfaceBase< DIMENSION >::change_mesh_data_structure(
        const MeshType& type )
    {
        auto new_mesh = SurfaceMesh< DIMENSION >::create_mesh( type );
        auto builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *new_mesh );
        builder->copy( *surface_mesh_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    template < index_t DIMENSION >
    ElementLocalVertex
        Region< DIMENSION >::find_cell_from_colocated_vertex_if_any(
            const vecn< DIMENSION >& vertex_vec ) const
    {
        ElementLocalVertex cell_local_vertex;
        volume_mesh_->find_cell_from_colocated_vertex_within_distance_if_any(
            vertex_vec, this->geomodel_.epsilon(), cell_local_vertex.element_id,
            cell_local_vertex.local_vertex_id );
        return cell_local_vertex;
    }

    template < index_t DIMENSION >
    void Region< DIMENSION >::change_mesh_data_structure( const MeshType& type )
    {
        auto new_mesh = VolumeMesh< DIMENSION >::create_mesh( type );
        auto builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *new_mesh );
        builder->copy( *volume_mesh_, true );
        update_mesh_storage_type( std::move( new_mesh ) );
    }

    MeshEntityType corner_type_name_static()
    {
        return Corner3D::type_name_static();
    }

    MeshEntityType line_type_name_static()
    {
        return Line3D::type_name_static();
    }

    MeshEntityType surface_type_name_static()
    {
        return Surface3D::type_name_static();
    }

    MeshEntityType region_type_name_static()
    {
        return Region3D::type_name_static();
    }

    template class geomodel_core_api GeoModelMeshEntity< 2 >;
    template class geomodel_core_api Corner< 2 >;
    template class geomodel_core_api Line< 2 >;
    template class geomodel_core_api SurfaceBase< 2 >;

    template class geomodel_core_api GeoModelMeshEntity< 3 >;
    template class geomodel_core_api Corner< 3 >;
    template class geomodel_core_api Line< 3 >;
    template class geomodel_core_api SurfaceBase< 3 >;
    template class geomodel_core_api Region< 3 >;
} // namespace RINGMesh
