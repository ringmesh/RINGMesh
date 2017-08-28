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

#include <ringmesh/geomodel/geomodel_builder_remove.h>

#include <ringmesh/geomodel/geomodel_builder.h>

/*!
 * @file Implementation of GeoModelEntity removal
 * @author Antoine Mazuyer
 * @author Pierre Anquez
 */

namespace RINGMesh {

    template< index_t DIMENSION >
    GeoModelBuilderRemovalBase< DIMENSION >::GeoModelBuilderRemovalBase(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
        nb_mesh_entity_types_ =
            geomodel_.entity_type_manager().mesh_entity_manager.nb_mesh_entity_types();
        nb_geological_entity_types_ = geomodel_.nb_geological_entity_types();
        nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_;
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::remove_mesh_entities(
        const std::set< gmme_id >& entities )
    {
        if( entities.empty() ) {
            return;
        } else {
            initialize_for_removal( entities );
            do_delete_flagged_mesh_entities();
            geomodel_.mesh.vertices.clear();
            update_mesh_entity_connectivity();
            flag_geological_entities_without_children();
            do_delete_flagged_geological_entities();
            update_geological_entity_connectivity();
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::remove_geological_entities(
        const std::set< gmge_id >& entities )
    {
        std::set< gmme_id > mesh_entities;
        for( const gmge_id& it : entities ) {
            const GeoModelGeologicalEntity< DIMENSION >& cur_gmge =
                geomodel_.geological_entity( it );
            for( index_t i : range( cur_gmge.nb_children() ) ) {
                mesh_entities.insert( cur_gmge.child( i ).gmme() );
            }
        }
        remove_mesh_entities( mesh_entities );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::do_delete_flagged_geological_entities()
    {
        for( index_t i : range( nb_geological_entity_types_ ) ) {
            const GeologicalEntityType& entity_type =
                index_to_geological_entity_type( i );
            for( index_t j : range( geomodel_.nb_geological_entities( entity_type ) ) ) {
                if( old_2_new_geological_entity_[i][j] == NO_ID ) {
                    builder_.geology.delete_geological_entity( entity_type, j );
                }
            }
            clear_null_geological_entities( i );
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::delete_mesh_entity(
        index_t type,
        index_t index )
    {
        const MeshEntityType& type_name = index_to_mesh_entity_type( type );
        gmme_id id( type_name, index );
        builder_.topology.delete_mesh_entity( type_name, index );
    }

    template< index_t DIMENSION >
    index_t GeoModelBuilderRemovalBase< DIMENSION >::geological_entity_type_index(
        const GeoModelGeologicalEntity< DIMENSION >& E ) const
    {
        const GeologicalEntityType& type = E.type_name();
        return geological_entity_type_to_index( type );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::set_mesh_entity_index(
        GeoModelMeshEntity< DIMENSION >& E,
        index_t new_index_in_geomodel )
    {
        GeoModelMeshEntityAccess< DIMENSION > gmme_access(
            dynamic_cast< GeoModelMeshEntity< DIMENSION >& >( E ) );
        gmme_access.modifiable_index() = new_index_in_geomodel;
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::set_geological_entity_index(
        GeoModelGeologicalEntity< DIMENSION >& E,
        index_t new_index_in_geomodel )
    {
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access(
            dynamic_cast< GeoModelGeologicalEntity< DIMENSION >& >( E ) );
        gmge_access.modifiable_index() = new_index_in_geomodel;
        return;
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_index(
        GeoModelMeshEntity< DIMENSION >& ME )
    {
        index_t old_id = ME.index();
        index_t type = mesh_entity_type_index( ME );
        index_t new_id = old_2_new_mesh_entity_[type][old_id];
        ringmesh_assert( new_id != NO_ID );
        set_mesh_entity_index( ME, new_id );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_geological_entity_index(
        GeoModelGeologicalEntity< DIMENSION >& GE )
    {
        index_t old_id = GE.index();
        index_t type = geological_entity_type_index( GE );
        index_t new_id = old_2_new_geological_entity_[type][old_id];
        ringmesh_assert( new_id != NO_ID );
        set_geological_entity_index( GE, new_id );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_boundaries(
        GeoModelMeshEntity< DIMENSION >& ME )
    {
        index_t type_index = boundary_type_index( ME.mesh_entity_type() );
        if( type_index == NO_ID ) {
            return;
        }
        for( index_t i : range( ME.nb_boundaries() ) ) {
            index_t old_boundary = ME.boundary_gmme( i ).index();
            index_t new_boundary = old_2_new_mesh_entity_[type_index][old_boundary];
            builder_.topology.set_mesh_entity_boundary( ME.gmme(), i, new_boundary );
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_incident_entity(
        GeoModelMeshEntity< DIMENSION >& E )
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            geomodel_.entity_type_manager().mesh_entity_manager;
        const MeshEntityType& incident_entity_type = manager.incident_entity_type(
            E.mesh_entity_type() );
        bool valid_type = manager.is_valid_type( incident_entity_type );
        if( !valid_type ) {
            return;
        }
        index_t incident_entity_type_index = mesh_entity_type_to_index(
            incident_entity_type );
        for( index_t i : range( E.nb_incident_entities() ) ) {
            index_t old_id = E.incident_entity_gmme( i ).index();
            index_t new_id =
                old_2_new_mesh_entity_[incident_entity_type_index][old_id];
            builder_.topology.set_mesh_entity_incident_entity( E.gmme(), i, new_id );
        }
    }
    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_mesh_entity_parents(
        GeoModelMeshEntity< DIMENSION >& E )
    {
        gmme_id id = E.gmme();
        for( index_t p : range( E.nb_parents() ) ) {
            const GeologicalEntityType& parent_type = E.parent_gmge( p ).type();
            index_t parent_type_index = geological_entity_type_to_index(
                parent_type );

            index_t old_id = E.parent_gmge( p ).index();
            index_t new_id = old_2_new_geological_entity_[parent_type_index][old_id];
            builder_.geology.set_mesh_entity_parent( id, p,
                gmge_id( parent_type, new_id ) );
        }
    }
    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::update_geological_entity_children(
        GeoModelGeologicalEntity< DIMENSION >& E )
    {
        if( E.nb_children() > 0 ) {
            index_t child_type = children_type_index( E.entity_type() );
            for( index_t i : range( E.nb_children() ) ) {
                index_t old_id = E.child_gmme( i ).index();
                index_t new_id = old_2_new_mesh_entity_[child_type][old_id];
                builder_.geology.set_geological_entity_child( E.gmge(), i, new_id );
            }
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderRemovalBase< DIMENSION >::fill_nb_children_vector()
    {
        for( index_t i : range( nb_childs_.size() ) ) {
            for( index_t j : range( nb_childs_[i].size() ) ) {
                nb_childs_[i][j] = geomodel_.geological_entity(
                    index_to_geological_entity_type( i ), j ).nb_children();
            }
        }
    }
    template< index_t DIMENSION >
    GeoModelBuilderRemoval< DIMENSION >::GeoModelBuilderRemoval(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : GeoModelBuilderRemovalBase< DIMENSION >( builder, geomodel )
    {
    }

    GeoModelBuilderRemoval< 3 >::GeoModelBuilderRemoval(
        GeoModelBuilder3D& builder,
        GeoModel3D& geomodel )
        : GeoModelBuilderRemovalBase< 3 >( builder, geomodel )
    {
    }

    void GeoModelBuilderRemoval< 3 >::update_mesh_entity(
        GeoModelMeshEntity3D& ME )
    {
        GeoModelBuilderRemovalBase3D::update_mesh_entity( ME );

        if( ME.mesh_entity_type() == Region3D::type_name_static() ) {
            Region3D& R = dynamic_cast< Region3D& >( ME );
            update_region_boundary_signs( R );
            delete_invalid_signs( R );
        }
    }

    template class RINGMESH_API GeoModelBuilderRemovalBase< 2 > ;
    template class RINGMESH_API GeoModelBuilderRemoval< 2 > ;

    template class RINGMESH_API GeoModelBuilderRemovalBase< 3 > ;
} // namespace RINGMesh

