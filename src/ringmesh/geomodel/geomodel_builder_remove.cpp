/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/geomodel/geomodel_builder_remove.h>

#include <ringmesh/basic/algorithm.h>

#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @file Implementation of GeoModelEntity removal
 * @author Antoine Mazuyer
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBuilderRemovalBase< DIMENSION >::GeoModelBuilderRemovalBase(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : builder_( builder ),
          geomodel_( geomodel ),
          geomodel_access_( geomodel )
    {
        nb_mesh_entity_types_ = geomodel_.entity_type_manager()
                                    .mesh_entity_manager.nb_mesh_entity_types();
        nb_geological_entity_types_ = geomodel_.nb_geological_entity_types();
        nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::remove_mesh_entities(
        const std::set< gmme_id >& entities )
    {
        if( entities.empty() )
        {
            return;
        }
        initialize_for_removal( entities );
        do_delete_flagged_mesh_entities();
        geomodel_.mesh.vertices.clear();
        update_mesh_entity_connectivity();
        flag_geological_entities_without_children();
        do_delete_flagged_geological_entities();
        update_geological_entity_connectivity();
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderRemovalBase< DIMENSION >::geological_entity_type_to_index(
        const GeologicalEntityType& type ) const
    {
        return find( geological_entity_types_, type );
    }

    template < index_t DIMENSION >
    index_t GeoModelBuilderRemovalBase< DIMENSION >::mesh_entity_type_to_index( const MeshEntityType& type ) const
    {
        return find( mesh_entity_types_, type );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::remove_geological_entities(
        const std::set< gmge_id >& entities )
    {
        std::set< gmme_id > mesh_entities;
        for( const gmge_id& it : entities )
        {
            const GeoModelGeologicalEntity< DIMENSION >& cur_gmge =
                geomodel_.geological_entity( it );
            for( auto i : range( cur_gmge.nb_children() ) )
            {
                mesh_entities.insert( cur_gmge.child( i ).gmme() );
            }
        }
        remove_mesh_entities( mesh_entities );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::
        do_delete_flagged_geological_entities()
    {
        for( auto i : range( nb_geological_entity_types_ ) )
        {
            const GeologicalEntityType& entity_type =
                index_to_geological_entity_type( i );
            for( auto j :
                range( geomodel_.nb_geological_entities( entity_type ) ) )
            {
                if( old_2_new_geological_entity_[i][j] == NO_ID )
                {
                    builder_.geology.delete_geological_entity( entity_type, j );
                }
            }
            clear_null_geological_entities( i );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_mesh_entity(
        index_t type, index_t index )
    {
        const MeshEntityType& type_name = index_to_mesh_entity_type( type );
        gmme_id id( type_name, index );
        builder_.topology.delete_mesh_entity( type_name, index );
    }

    template < index_t DIMENSION >
    index_t
        GeoModelBuilderRemovalBase< DIMENSION >::geological_entity_type_index(
            const GeoModelGeologicalEntity< DIMENSION >& E ) const
    {
        const GeologicalEntityType& type = E.type_name();
        return geological_entity_type_to_index( type );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::set_mesh_entity_index(
        GeoModelMeshEntity< DIMENSION >& mesh_entity, index_t new_index_in_geomodel )
    {
        GeoModelMeshEntityAccess< DIMENSION > gmme_access(
            dynamic_cast< GeoModelMeshEntity< DIMENSION >& >( mesh_entity ) );
        gmme_access.modifiable_index() = new_index_in_geomodel;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::do_delete_flagged_mesh_entities()
    {
        for( auto i : range( nb_mesh_entity_types_ ) )
        {
            for( auto j : range( nb_initial_mesh_entities_[i] ) )
            {
                if( mesh_entity_to_erase_[i][j] )
                {
                    const auto& type_name = index_to_mesh_entity_type( i );
                    for( auto p :
                        range( geomodel_.mesh_entity( type_name, j )
                                   .nb_parents() ) )
                    {
                        auto parent = geomodel_.mesh_entity( type_name, j )
                                .parent_gmge( p );
                        nb_childs_[geological_entity_type_to_index(
                            parent.type() )][parent.index()]--;
                    }

                    delete_mesh_entity( i, j );
                }
            }
            clear_null_mesh_entities( i );
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::clear_null_mesh_entities(
        index_t type )
    {
        const auto& type_name = index_to_mesh_entity_type( type );
        auto& store = geomodel_access_.modifiable_mesh_entities( type_name );
        store.erase(
            std::remove( store.begin(), store.end(),
                static_cast< std::
                        unique_ptr< GeoModelMeshEntity< DIMENSION > > >(
                    nullptr ) ),
            store.end() );

        // QC
        ringmesh_assert( geomodel_.nb_mesh_entities( type_name )
                         == nb_initial_mesh_entities_[type]
                                - nb_removed_mesh_entities_[type] );
    }


    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::clear_null_geological_entities(
        index_t type )
    {
        const auto& type_name = index_to_geological_entity_type( type );
        auto& store = geomodel_access_.modifiable_geological_entities( type_name );
        store.erase(
            std::remove( store.begin(), store.end(),
                static_cast< std::
                        unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >(
                    nullptr ) ),
            store.end() );

        // QC
        ringmesh_assert( geomodel_.nb_geological_entities( type_name )
                         == nb_initial_geological_entities_[type]
                                - nb_removed_geological_entities_[type] );
    }

    template< index_t DIMENSION >
    template < typename TEST, typename THINGS_TO_DELETE >
    void GeoModelBuilderRemovalBase< DIMENSION >::remove_invalid_values(
        std::vector< THINGS_TO_DELETE >& vector, const TEST& test )
    {
        auto new_end = std::remove_if( vector.begin(), vector.end(), test );
        if( new_end == vector.begin() )
        {
            // Clear instead of erase, because the behavior would be
            // undefined.
            vector.clear();
        }
        else if( new_end < vector.end() )
        {
            vector.erase( new_end, vector.end() );
        }
    }

    template< index_t DIMENSION >
    const MeshEntityType GeoModelBuilderRemovalBase< DIMENSION >::children_type(
        const GeologicalEntityType& type ) const
    {
        const auto& family =
            geomodel_.entity_type_manager().relationship_manager;
        return family.child_type( type );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::set_boundary_side(
        Region3D& region, index_t boundary_index, bool new_side )
    {
        ringmesh_assert( boundary_index < region.nb_boundaries() );
        GeoModelMeshEntityAccess< DIMENSION > region_access(
            geomodel_access_.modifiable_mesh_entity( region.gmme() ) );
        region_access.modifiable_sides()[boundary_index] = new_side;
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_invalid_children(
        GeoModelGeologicalEntity< DIMENSION >& E )
    {
        if( E.nb_children() == 0 )
        {
            return;
        }
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access{ E };
        const auto& manager =
            E.geomodel().entity_type_manager().relationship_manager;
        const auto& child_type = children_type( E.entity_type() );
        gmme_id invalid_child{ child_type, NO_ID };
        remove_invalid_values( gmge_access.modifiable_children(),
            [&invalid_child, &manager]( index_t i ) {
                return manager.child_of_gmge( i ) == invalid_child;
            } );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_invalid_boundaries(
        GeoModelMeshEntity< DIMENSION >& E )
    {
        const auto& b_type = boundary_entity_type( E.mesh_entity_type() );
        gmme_id invalid{ b_type, NO_ID };
        if( !geomodel_.entity_type_manager()
                 .mesh_entity_manager.is_valid_type( b_type ) )
        {
            return;
        }
        GeoModelMeshEntityAccess< DIMENSION > gmme_access{ E };
        const auto& manager =
            E.geomodel().entity_type_manager().relationship_manager;
        remove_invalid_values( gmme_access.modifiable_boundaries(),
            [&invalid, &manager]( index_t i ) {
                return manager.boundary_gmme( i ) == invalid;
            } );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_invalid_incident_entity(
        GeoModelMeshEntity< DIMENSION >& E )
    {
        const auto& in_ent_type = incident_entity_type( E.mesh_entity_type() );
        gmme_id invalid{ in_ent_type, NO_ID };
        if( !geomodel_.entity_type_manager()
                 .mesh_entity_manager.is_valid_type( in_ent_type ) )
        {
            return;
        }
        GeoModelMeshEntityAccess< DIMENSION > gmme_access{ E };
        const auto& manager =
            E.geomodel().entity_type_manager().relationship_manager;
        remove_invalid_values( gmme_access.modifiable_incident_entities(),
            [&invalid, &manager]( index_t i ) {
                return manager.incident_entity_gmme( i ) == invalid;
            } );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_invalid_parents(
        GeoModelMeshEntity< DIMENSION >& E )
    {
        GeoModelMeshEntityAccess< DIMENSION > gmme_access{ E };
        const auto& manager =
            E.geomodel().entity_type_manager().relationship_manager;
        remove_invalid_values(
            gmme_access.modifiable_parents(), [&manager]( index_t i ) {
                return manager.parent_of_gmme( i ).index() == NO_ID;
            } );
    }

    template< index_t DIMENSION >
    const MeshEntityType& GeoModelBuilderRemovalBase< DIMENSION >::incident_entity_type(
        const MeshEntityType& type ) const
    {
        const auto& family =
            geomodel_.entity_type_manager().mesh_entity_manager;
        return family.incident_entity_type( type );
    }

    template< index_t DIMENSION >
    const MeshEntityType& GeoModelBuilderRemovalBase< DIMENSION >::boundary_entity_type(
        const MeshEntityType& type ) const
    {
        const auto& family =
            geomodel_.entity_type_manager().mesh_entity_manager;
        return family.boundary_entity_type( type );
    }

    template< index_t DIMENSION >
    index_t GeoModelBuilderRemovalBase< DIMENSION >::children_type_index(
        const GeologicalEntityType& type ) const
    {
        const auto& child_type = children_type( type );
        return mesh_entity_type_to_index( child_type );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::fill_to_erase_vectors(
        const std::set< gmme_id >& mesh_entities_to_remove )
    {
        for( const gmme_id& cur : mesh_entities_to_remove )
        {
            index_t type_index = mesh_entity_type_to_index( cur.type() );
            mesh_entity_to_erase_[type_index][cur.index()] = true;
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::check_if_entities_are_meshed(
        const std::set< gmme_id >& mesh_entities_to_remove )
    {
        for( const gmme_id& it : mesh_entities_to_remove )
        {
            if( !geomodel_.entity_type_manager()
                     .mesh_entity_manager.is_valid_type( it.type() ) )
            {
                throw RINGMeshException( "REMOVE", "You try to remove a "
                                                   "Geological Entity "
                                                   "using mesh removal." );
            }
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::set_geological_entity_index(
        GeoModelGeologicalEntity< DIMENSION >& geological_entity,
        index_t new_index_in_geomodel )
    {
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access(
            dynamic_cast< GeoModelGeologicalEntity< DIMENSION >& >( geological_entity ) );
        gmge_access.modifiable_index() = new_index_in_geomodel;
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_index(
        GeoModelMeshEntity< DIMENSION >& mesh_entity )
    {
        index_t old_id = mesh_entity.index();
        index_t type = mesh_entity_type_index( mesh_entity );
        index_t new_id = old_2_new_mesh_entity_[type][old_id];
        ringmesh_assert( new_id != NO_ID );
        set_mesh_entity_index( mesh_entity, new_id );
    }

    template < index_t DIMENSION >
    void
        GeoModelBuilderRemovalBase< DIMENSION >::update_geological_entity_index(
            GeoModelGeologicalEntity< DIMENSION >& geological_entity )
    {
        index_t old_id = geological_entity.index();
        index_t type = geological_entity_type_index( geological_entity );
        index_t new_id = old_2_new_geological_entity_[type][old_id];
        ringmesh_assert( new_id != NO_ID );
        set_geological_entity_index( geological_entity, new_id );
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_boundaries(
        GeoModelMeshEntity< DIMENSION >& mesh_entity )
    {
        index_t type_index = boundary_type_index( mesh_entity.mesh_entity_type() );
        if( type_index == NO_ID )
        {
            return;
        }
        for( auto i : range( mesh_entity.nb_boundaries() ) )
        {
            index_t old_boundary = mesh_entity.boundary_gmme( i ).index();
            index_t new_boundary =
                old_2_new_mesh_entity_[type_index][old_boundary];
            builder_.topology.set_mesh_entity_boundary(
                mesh_entity.gmme(), i, new_boundary );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::
        update_mesh_entity_incident_entity( GeoModelMeshEntity< DIMENSION >& mesh_entity )
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const MeshEntityType& incident_entity_type =
            manager.incident_entity_type( mesh_entity.mesh_entity_type() );
        bool valid_type = manager.is_valid_type( incident_entity_type );
        if( !valid_type )
        {
            return;
        }
        index_t incident_entity_type_index =
            mesh_entity_type_to_index( incident_entity_type );
        for( auto i : range( mesh_entity.nb_incident_entities() ) )
        {
            index_t old_id = mesh_entity.incident_entity_gmme( i ).index();
            index_t new_id =
                old_2_new_mesh_entity_[incident_entity_type_index][old_id];
            builder_.topology.set_mesh_entity_incident_entity(
                mesh_entity.gmme(), i, new_id );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_parents(
        GeoModelMeshEntity< DIMENSION >& mesh_entity )
    {
        gmme_id id = mesh_entity.gmme();
        for( auto p : range( mesh_entity.nb_parents() ) )
        {
            const GeologicalEntityType& parent_type = mesh_entity.parent_gmge( p ).type();
            index_t parent_type_index =
                geological_entity_type_to_index( parent_type );

            index_t old_id = mesh_entity.parent_gmge( p ).index();
            index_t new_id =
                old_2_new_geological_entity_[parent_type_index][old_id];
            builder_.geology.set_mesh_entity_parent(
                id, p, gmge_id( parent_type, new_id ) );
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::
        update_geological_entity_children(
            GeoModelGeologicalEntity< DIMENSION >& geological_entity )
    {
        if( geological_entity.nb_children() > 0 )
        {
            index_t child_type = children_type_index( geological_entity.entity_type() );
            for( auto i : range( geological_entity.nb_children() ) )
            {
                index_t old_id = geological_entity.child_gmme( i ).index();
                index_t new_id = old_2_new_mesh_entity_[child_type][old_id];
                builder_.geology.set_geological_entity_child(
                    geological_entity.gmge(), i, new_id );
            }
        }
    }

    template < index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::fill_nb_children_vector()
    {
        for( auto i : range( nb_childs_.size() ) )
        {
            for( auto j : range( nb_childs_[i].size() ) )
            {
                nb_childs_[i][j] =
                    geomodel_
                        .geological_entity(
                            index_to_geological_entity_type( i ), j )
                        .nb_children();
            }
        }
    }

    template < index_t DIMENSION >
    GeoModelBuilderRemoval< DIMENSION >::GeoModelBuilderRemoval(
        GeoModelBuilder< DIMENSION >& builder, GeoModel< DIMENSION >& geomodel )
        : GeoModelBuilderRemovalBase< DIMENSION >( builder, geomodel )
    {
    }

    GeoModelBuilderRemoval< 3 >::GeoModelBuilderRemoval(
        GeoModelBuilder3D& builder, GeoModel3D& geomodel )
        : GeoModelBuilderRemovalBase< 3 >( builder, geomodel )
    {
    }

    void GeoModelBuilderRemoval< 3 >::update_mesh_entity(
        GeoModelMeshEntity3D& ME )
    {
        GeoModelBuilderRemovalBase3D::update_mesh_entity( ME );

        if( ME.mesh_entity_type() == Region3D::type_name_static() )
        {
            auto& R = dynamic_cast< Region3D& >( ME );
            update_region_boundary_signs( R );
            delete_invalid_signs( R );
        }
    }

    void GeoModelBuilderRemoval< 3 >::set_boundary_side(
        Region3D& R, index_t boundary_index, bool new_side )
    {
        ringmesh_assert( boundary_index < R.nb_boundaries() );
        GeoModelMeshEntityAccess3D region_access{
            geomodel_access_.modifiable_mesh_entity( R.gmme() ) };
        region_access.modifiable_sides()[boundary_index] = new_side;
    }

    void GeoModelBuilderRemoval< 3 >::update_region_boundary_signs( Region3D& R )
    {
        const auto& surface_type = boundary_entity_type( R.mesh_entity_type() );
        gmme_id invalid_value{ surface_type, NO_ID };

        index_t offset{ 0 };
        for( index_t i = 0; i + offset < R.nb_boundaries(); ++i )
        {
            if( R.boundary_gmme( i ) == invalid_value )
            {
                offset++;
            }
            else
            {
                bool new_side = R.side( i + offset );
                set_boundary_side( R, i, new_side );
            }
        }
    }

    void GeoModelBuilderRemoval< 3 >::delete_invalid_signs( Region3D& R )
    {
        GeoModelMeshEntityAccess3D region_access{
            geomodel_access_.modifiable_mesh_entity( R.gmme() ) };
        region_access.modifiable_sides().resize( R.nb_boundaries() );
    }

    template class RINGMESH_API GeoModelBuilderRemovalBase< 2 >;
    template class RINGMESH_API GeoModelBuilderRemoval< 2 >;

    template class RINGMESH_API GeoModelBuilderRemovalBase< 3 >;
} // namespace RINGMesh
