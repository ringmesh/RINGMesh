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

    GeoModelBuilderRemoval::GeoModelBuilderRemoval(
        GeoModelBuilder& builder,
        GeoModel& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
        nb_mesh_entity_types_ = EntityTypeManager::nb_mesh_entity_types() ;
        nb_geological_entity_types_ = geomodel_.nb_geological_entity_types() ;
        nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_ ;
    }

    void GeoModelBuilderRemoval::remove_mesh_entities(
        const std::set< gme_t >& entities )
    {
        if( entities.empty() ) {
            return ;
        } else {
            initialize_for_removal( entities ) ;
            do_delete_flagged_mesh_entities() ;
            geomodel_.mesh.vertices.clear() ;
            update_mesh_entity_connectivity() ;
            flag_geological_entities_without_children() ;
            do_delete_flagged_geological_entities() ;
            update_geological_entity_connectivity() ;
            update_universe() ;
        }
    }

    void GeoModelBuilderRemoval::remove_geological_entities(
        const std::set< gme_t >& entities )
    {
        check_if_entities_are_not_meshed_entities( entities ) ;
        std::set< gme_t > mesh_entities ;
        for( const gme_t& it : entities ) {
            const GeoModelGeologicalEntity& cur_gmge = geomodel_.geological_entity(
                it ) ;
            for( index_t i = 0; i < cur_gmge.nb_children(); i++ ) {
                mesh_entities.insert( cur_gmge.child( i ).gme_id() ) ;
            }
        }
        remove_mesh_entities( mesh_entities ) ;
    }

    void GeoModelBuilderRemoval::do_delete_flagged_geological_entities()
    {
        for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
            const EntityType& entity_type = index_to_geological_entity_type( i ) ;
            for( index_t j = 0; j < geomodel_.nb_geological_entities( entity_type );
                ++j ) {
                if( old_2_new_geological_entity_[i][j] == NO_ID ) {
                    builder_.geology.delete_geological_entity( entity_type, j ) ;
                }
            }
            clear_null_geological_entities( i ) ;
        }
    }

    void GeoModelBuilderRemoval::delete_mesh_entity( index_t type, index_t index )
    {
        const EntityType& type_name = index_to_mesh_entity_type( type ) ;
        gme_t id( type_name, index ) ;
        builder_.topology.delete_mesh_entity( type_name, index ) ;
    }

    index_t GeoModelBuilderRemoval::geological_entity_type_index(
        const GeoModelGeologicalEntity& E ) const
    {
        const EntityType& type = E.type_name() ;
        return geological_entity_type_to_index( type ) ;
    }

    void GeoModelBuilderRemoval::set_entity_index(
        GeoModelEntity& E,
        index_t new_index_in_geomodel )
    {
        if( EntityTypeManager::is_mesh_entity_type( E.type_name() ) ) {
            GeoModelMeshEntityAccess gmme_access(
                dynamic_cast< GeoModelMeshEntity& >( E ) ) ;
            gmme_access.modifiable_index() = new_index_in_geomodel ;
            return ;
        } else if( geomodel_.entity_type_manager().is_geological_entity_type(
            E.type_name() ) ) {
            GeoModelGeologicalEntityAccess gmge_access(
                dynamic_cast< GeoModelGeologicalEntity& >( E ) ) ;
            gmge_access.modifiable_index() = new_index_in_geomodel ;
            return ;
        } else {
            ringmesh_assert_not_reached ;
        }

    }

    void GeoModelBuilderRemoval::update_mesh_entity_index( GeoModelMeshEntity& ME )
    {
        index_t old_id = ME.index() ;
        index_t type = mesh_entity_type_index( ME ) ;
        index_t new_id = old_2_new_mesh_entity_[type][old_id] ;
        ringmesh_assert( new_id != NO_ID ) ;
        set_entity_index( ME, new_id ) ;
    }

    void GeoModelBuilderRemoval::update_geological_entity_index(
        GeoModelGeologicalEntity& GE )
    {
        index_t old_id = GE.index() ;
        index_t type = geological_entity_type_index( GE ) ;
        index_t new_id = old_2_new_geological_entity_[type][old_id] ;
        ringmesh_assert( new_id != NO_ID ) ;
        set_entity_index( GE, new_id ) ;
    }

    void GeoModelBuilderRemoval::update_mesh_entity_boundaries(
        GeoModelMeshEntity& ME )
    {
        index_t type_index = boundary_type_index( ME.entity_type() ) ;
        if( type_index == NO_ID ) {
            return ;
        }
        for( index_t i = 0; i < ME.nb_boundaries(); ++i ) {
            index_t old_boundary = ME.boundary_gme( i ).index ;
            index_t new_boundary = old_2_new_mesh_entity_[type_index][old_boundary] ;
            builder_.topology.set_mesh_entity_boundary( ME.gme_id(), i,
                new_boundary ) ;
        }
    }

    void GeoModelBuilderRemoval::update_mesh_entity_in_boundary(
        GeoModelMeshEntity& E )
    {
        const EntityType& in_boundary_type = EntityTypeManager::in_boundary_type(
            E.entity_type() ) ;
        bool valid_type = EntityTypeManager::is_defined_type( in_boundary_type ) ;
        if( !valid_type ) {
            return ;
        }
        index_t in_boundary_type_index = mesh_entity_type_to_index(
            in_boundary_type ) ;
        for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
            index_t old_id = E.in_boundary_gme( i ).index ;
            index_t new_id = old_2_new_mesh_entity_[in_boundary_type_index][old_id] ;
            builder_.topology.set_mesh_entity_in_boundary( E.gme_id(), i, new_id ) ;
        }
    }
    void GeoModelBuilderRemoval::update_mesh_entity_parents( GeoModelMeshEntity& E )
    {
        for( index_t p = 0; p < E.nb_parents(); ++p ) {
            const EntityType& parent_type = E.parent_gme( p ).type ;
            index_t parent_type_index = geological_entity_type_to_index(
                parent_type ) ;

            index_t old_id = E.parent_gme( p ).index ;
            index_t new_id = old_2_new_geological_entity_[parent_type_index][old_id] ;
            builder_.geology.set_mesh_entity_parent( E.gme_id(), p,
                gme_t( parent_type, new_id ) ) ;
        }
    }
    void GeoModelBuilderRemoval::update_geological_entity_children(
        GeoModelGeologicalEntity& E )
    {
        if( E.nb_children() > 0 ) {
            index_t child_type = children_type_index( E.entity_type() ) ;
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                index_t old_id = E.child_gme( i ).index ;
                index_t new_id = old_2_new_mesh_entity_[child_type][old_id] ;
                builder_.geology.set_geological_entity_child( E.gme_id(), i,
                    new_id ) ;
            }
        }
    }
    void GeoModelBuilderRemoval::update_universe_sided_boundaries( Universe& U )
    {
        index_t b_type_index = mesh_entity_type_to_index(
            Surface::type_name_static() ) ;
        index_t side_offset = 0 ;
        for( index_t i = 0; i < U.nb_boundaries(); ++i ) {
            index_t old_id = U.boundary_gme( i ).index ;
            index_t new_id = old_2_new_mesh_entity_[b_type_index][old_id] ;

            bool new_side = false ;
            // Mechanism to update the sides is not the same than to update
            // the boundary indices -- annoying
            if( new_id == NO_ID ) {
                side_offset++ ;
            } else if( i + side_offset < U.nb_boundaries() ) {
                // After that we do not care the values will be dropped
                new_side = U.side( i + side_offset ) ;
            }
            builder_.topology.set_universe_boundary( i, new_id, new_side ) ;
        }
    }

    void GeoModelBuilderRemoval::delete_invalid_parents( GeoModelMeshEntity& E )
    {
        //  Cannot use remove directly, do it by hand like the signs
        index_t offset = 0 ;
        index_t new_size = 0 ;
        for( index_t i = 0; i + offset < E.nb_parents(); ++i ) {
            if( E.parent( i ).index() == NO_ID ) {
                offset++ ;
            } else {
                gme_t new_id = E.parent_gme( i + offset ) ;
                builder_.geology.set_mesh_entity_parent( E.gme_id(), i, new_id ) ;
            }
            new_size = i + 1 ; /// @todo Check that this is the correct size
        }
        GeoModelMeshEntityAccess gmme_access( E ) ;
        gmme_access.modifiable_parents().resize( new_size ) ;
    }

    void GeoModelBuilderRemoval::fill_nb_children_vector()
    {
        for( index_t i = 0; i < nb_childs_.size(); i++ ) {
            for( index_t j = 0; j < nb_childs_[i].size(); j++ ) {
                nb_childs_[i][j] = geomodel_.geological_entity(
                    index_to_geological_entity_type( i ), j ).nb_children() ;
            }
        }
    }

    void GeoModelBuilderRemoval::delete_invalid_children(
        GeoModelGeologicalEntity& E )
    {
        if( E.nb_children() == 0 ) {
            return ;
        } else {
            const EntityType& child_type = children_type( E.entity_type() ) ;
            gme_t invalid_child( child_type, NO_ID ) ;
            GeoModelGeologicalEntityAccess gmge_access( E ) ;
            remove_invalid_values( gmge_access.modifiable_children(),
                invalid_child ) ;
        }
    }

}

