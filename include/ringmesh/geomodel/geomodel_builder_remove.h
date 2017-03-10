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

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

/*!
 * @brief Builder tools to remove entities from a GeoModel.
 * @author Antoine Mazuyer
 * @author Pierre Anquez
 */

namespace RINGMesh {
    class GeoModelBuilder ;
}

namespace RINGMesh {

    /*!
     * @brief Builder tools to remove entities from a GeoModel
     */
    class RINGMESH_API GeoModelBuilderRemoval {
    ringmesh_disable_copy( GeoModelBuilderRemoval ) ;
        friend class GeoModelBuilder ;

    public:

        /*!
         * @brief Remove a list of mesh entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_mesh_entities( const std::set< gme_t >& entities ) ;

        /*!
         * @brief Remove a list of geological entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_geological_entities( const std::set< gme_t >& entities ) ;

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gme_t >& entities_to_remove ) ;

    private:
        GeoModelBuilderRemoval( GeoModelBuilder& builder, GeoModel& geomodel ) ;

        // ---  High level functions ----------
        void initialize_for_removal(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            nb_mesh_entity_types_ = EntityTypeManager::nb_mesh_entity_types() ;
            nb_geological_entity_types_ = geomodel_.nb_geological_entity_types() ;
            nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_ ;
            nb_removed_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            nb_removed_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            fill_entity_type_to_index_map() ;
            fill_nb_initial_entities() ;
            initialize_costly_storage() ;
            fill_nb_children_vector() ;

            check_if_entities_are_meshed( mesh_entities_to_remove ) ;
            fill_to_erase_vectors( mesh_entities_to_remove ) ;
            fill_removed_entities_and_mapping() ;
        }
        void do_delete_flagged_mesh_entities()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_mesh_entities_[i]; ++j ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        const EntityType& type_name = index_to_mesh_entity_type(
                            i ) ;
                        for( index_t p = 0;
                            p < geomodel_.mesh_entity( type_name, j ).nb_parents();
                            p++ ) {
                            gme_t parent =
                                geomodel_.mesh_entity( type_name, j ).parent_gme(
                                    p ) ;
                            nb_childs_[geological_entity_type_to_index( parent.type )][parent.index]-- ;
                        }

                        delete_mesh_entity( i, j ) ;
                    }
                }
                clear_null_mesh_entities( i ) ;
            }
        }
        void do_delete_flagged_geological_entities() ;

        void check_if_entities_are_meshed(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                if( !RINGMesh::EntityTypeManager::is_mesh_entity_type( it->type ) ) {
                    throw RINGMeshException( "REMOVE",
                        "You try to remove a Geological Entity using mesh removal." ) ;
                }
            }
        }

        void check_if_entities_are_not_meshed_entities(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                if( RINGMesh::EntityTypeManager::is_mesh_entity_type( it->type ) ) {
                    throw RINGMeshException( "REMOVE",
                        "You try to remove a Mesh Entity using geological removal" ) ;
                }
            }
        }
        void initialize_costly_storage()
        {
            mesh_entity_to_erase_.resize( nb_mesh_entity_types_ ) ;

            old_2_new_mesh_entity_.resize( nb_mesh_entity_types_ ) ;
            old_2_new_geological_entity_.resize( nb_geological_entity_types_ ) ;
            nb_childs_.resize( nb_geological_entity_types_ ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                index_t size = geomodel_.nb_mesh_entities(
                    index_to_mesh_entity_type( i ) ) ;
                mesh_entity_to_erase_[i].resize( size, false ) ;
                old_2_new_mesh_entity_[i].resize( size, 0 ) ;
            }

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                index_t size = geomodel_.nb_geological_entities(
                    index_to_geological_entity_type( i ) ) ;
                old_2_new_geological_entity_[i].resize( size, 0 ) ;

                nb_childs_[i].resize( size, 0 ) ;

            }

        }
        void delete_mesh_entity( index_t type, index_t index ) ;

        void clear_null_mesh_entities( index_t type )
        {
            const EntityType& type_name = index_to_mesh_entity_type( type ) ;
            std::vector< GeoModelMeshEntity* >& store =
                geomodel_access_.modifiable_mesh_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelMeshEntity* >( nil ) ), store.end() ) ;

            // QC
            ringmesh_assert( geomodel_.nb_mesh_entities( type_name )
                == nb_initial_mesh_entities_[type] - nb_removed_mesh_entities_[type] ) ;
        }

        void clear_null_geological_entities( index_t type )
        {
            const EntityType& type_name = index_to_geological_entity_type( type ) ;
            std::vector< GeoModelGeologicalEntity* >& store =
                geomodel_access_.modifiable_geological_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelGeologicalEntity* >( nil ) ),
                store.end() ) ;

            // QC
            ringmesh_assert( geomodel_.nb_geological_entities( type_name )
                == nb_initial_geological_entities_[type] - nb_removed_geological_entities_[type] ) ;
        }
        void update_mesh_entity_connectivity()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < geomodel_.nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = geomodel_access_.modifiable_mesh_entity(
                        new_id ) ;
                    update_mesh_entity_index( ME ) ;
                    ringmesh_assert( new_id == ME.gme_id() ) ;
                    update_mesh_entity_boundaries( ME ) ;
                    delete_invalid_boundaries( ME ) ;

                    update_mesh_entity_in_boundary( ME ) ;
                    delete_invalid_in_boundary( ME ) ;

                    if( ME.entity_type() == Region::type_name_static() ) {
                        Region& R = dynamic_cast< Region& >( ME ) ;
                        update_region_boundary_signs( R ) ;
                        delete_invalid_signs( R ) ;
                    }
                }
            }
        }

        void update_geological_entity_connectivity()
        {

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_geological_entity_type(
                    i ) ;
                for( index_t j = 0;
                    j < geomodel_.nb_geological_entities( entity_type ); ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelGeologicalEntity& GE =
                        geomodel_access_.modifiable_geological_entity( new_id ) ;
                    update_geological_entity_index( GE ) ;
                    update_geological_entity_children( GE ) ;
                    delete_invalid_children( GE ) ;
                }
            }

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < geomodel_.nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = geomodel_access_.modifiable_mesh_entity(
                        new_id ) ;
                    update_mesh_entity_parents( ME ) ;
                    delete_invalid_parents( ME ) ;
                }
            }
        }

        void update_universe()
        {
            Universe& U = geomodel_access_.modifiable_universe() ;
            update_universe_sided_boundaries( U ) ;
            delete_invalid_universe_sided_boundaries( U ) ;
        }

        //        void remove_dependencies()
        //        {
        //            std::set< gme_t > new_gmme_to_remove ;
        //            for( index_t me = 0;
        //                me < geomodel().nb_mesh_entities( starting_dependency_ ); me++ ) {
        //                const GeoModelMeshEntity& cur_gmme = geomodel().mesh_entity(
        //                    starting_dependency_, me ) ;
        //                if( cur_gmme.in_boundary( 0 ).index() == NO_ID
        //                    && cur_gmme.nb_in_boundary() == 1 ) {
        //                    new_gmme_to_remove.insert( cur_gmme.gme_id() ) ;
        //                }
        //            }
        //            if( starting_dependency_ != Corner::type_name_static() ) {
        //                starting_dependency_ = EntityTypeManager::boundary_type(
        //                    starting_dependency_ ) ;
        //                remove_mesh_entities_with_dependencies( new_gmme_to_remove ) ;
        //            }
        //        }

        //------  Initialization -------
        void fill_removed_entities_and_mapping()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_mesh_entities_[i]; ++j ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        nb_removed_mesh_entities_[i]++ ;
                        old_2_new_mesh_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_mesh_entity_[i][j] = j
                            - nb_removed_mesh_entities_[i] ;
                    }
                }
            }
        }
        void fill_to_erase_vectors(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                gme_t cur = *it ;

                index_t type_index = mesh_entity_type_to_index( cur.type ) ;
                mesh_entity_to_erase_[type_index][cur.index] = true ;
            }
        }

        void fill_nb_children_vector() ;

        void fill_nb_initial_entities()
        {
            nb_initial_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& type = index_to_mesh_entity_type( i ) ;
                nb_initial_mesh_entities_[i] = geomodel_.nb_mesh_entities( type ) ;
            }

            nb_initial_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& type = index_to_geological_entity_type( i ) ;
                nb_initial_geological_entities_[i] =
                    geomodel_.nb_geological_entities( type ) ;
            }
        }
        void fill_entity_type_to_index_map()
        {
            const EntityTypeManager& manager = geomodel_.entity_type_manager() ;
            mesh_entity_types_.insert( mesh_entity_types_.end(),
                manager.mesh_entity_types().begin(),
                manager.mesh_entity_types().end() ) ;

            if( nb_geological_entity_types_ != 0 ) {
                geological_entity_types_.insert( geological_entity_types_.end(),
                    manager.geological_entity_types().begin(),
                    manager.geological_entity_types().end() ) ;
            }

        }

        // ---- Easier access to relationships between EntityTypes
        index_t mesh_entity_type_index( const GeoModelMeshEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return mesh_entity_type_to_index( type ) ;
        }
        index_t geological_entity_type_index(
            const GeoModelGeologicalEntity& E ) const ;
        index_t children_type_index( const EntityType& type ) const
        {
            const EntityType& child_type = children_type( type ) ;
            return mesh_entity_type_to_index( child_type ) ;
        }
        const EntityType children_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.child_type( type ) ;
        }
        index_t boundary_type_index( const EntityType& type ) const
        {
            const EntityType& b_type = boundary_type( type ) ;
            if( !EntityTypeManager::is_defined_type( b_type ) ) {
                return NO_ID ;
            } else {
                return mesh_entity_type_to_index( b_type ) ;
            }
        }
        const EntityType& boundary_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.boundary_type( type ) ;
        }
        index_t in_boundary_type_index( const EntityType& type ) const
        {
            const EntityType& in_b_type = in_boundary_type( type ) ;
            if( !EntityTypeManager::is_defined_type( in_b_type ) ) {
                return NO_ID ;
            } else {
                return mesh_entity_type_to_index( in_b_type ) ;
            }
        }
        const EntityType& in_boundary_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.in_boundary_type( type ) ;
        }
        bool is_mesh_entity( index_t i ) const
        {
            return i < nb_mesh_entity_types_ ;
        }
        bool is_geological_entity( index_t i ) const
        {
            return !is_mesh_entity( i ) ;
        }

        // ----  Update connectivity functions  ------

        void flag_geological_entities_without_children()
        {
            for( index_t i = 0; i < nb_childs_.size(); i++ ) {
                for( index_t j = 0; j < nb_childs_[i].size(); j++ ) {
                    if( nb_childs_[i][j] == 0 ) {
                        nb_removed_geological_entities_[i]++ ;
                        old_2_new_geological_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_geological_entity_[i][j] = j
                            - nb_removed_geological_entities_[i] ;
                    }

                }
            }
        }

        void set_entity_index( GeoModelEntity& E, index_t new_index_in_geomodel ) ;

        void update_mesh_entity_index( GeoModelMeshEntity& ME ) ;
        void update_geological_entity_index( GeoModelGeologicalEntity& GE ) ;
        void update_mesh_entity_boundaries( GeoModelMeshEntity& ME ) ;

        void set_boundary_side( Region& R, index_t boundary_index, bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() ) ;
            GeoModelMeshEntityAccess region_access(
                geomodel_access_.modifiable_mesh_entity( R.gme_id() ) ) ;
            region_access.modifiable_sides()[boundary_index] = new_side ;
        }

        void update_region_boundary_signs( Region& R )
        {
            const EntityType& surface_type = boundary_type( R.entity_type() ) ;
            gme_t invalid_value( surface_type, NO_ID ) ;

            index_t offset = 0 ;
            for( index_t i = 0; i + offset < R.nb_boundaries(); ++i ) {
                if( R.boundary_gme( i ) == invalid_value ) {
                    offset++ ;
                } else {
                    bool new_side = R.side( i + offset ) ;
                    set_boundary_side( R, i, new_side ) ;
                }
            }
        }
        void update_mesh_entity_in_boundary( GeoModelMeshEntity& E ) ;
        void update_mesh_entity_parents( GeoModelMeshEntity& E ) ;
        void update_geological_entity_children( GeoModelGeologicalEntity& E ) ;
        void update_universe_sided_boundaries( Universe& U ) ;

        // --- Deletion of some values the GeoModel storage
        void remove_invalid_values(
            std::vector< gme_t >& vector,
            const gme_t& invalid_value )
        {
            std::vector< gme_t >::iterator new_end = std::remove( vector.begin(),
                vector.end(), invalid_value ) ;
            if( new_end == vector.begin() ) {
                // Clear instead of erase, because the behavior would be undefined.
                vector.clear() ;
            } else if( new_end < vector.end() ) {
                vector.erase( new_end, vector.end() ) ;
            }
        }
        void delete_invalid_children( GeoModelGeologicalEntity& E ) ;
        void delete_invalid_boundaries( GeoModelMeshEntity& E )
        {
            const EntityType& b_type = boundary_type( E.entity_type() ) ;
            gme_t invalid( b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( b_type ) ) {
                return ;
            } else {
                GeoModelMeshEntityAccess gmme_access( E ) ;
                remove_invalid_values( gmme_access.modifiable_boundaries(),
                    invalid ) ;
            }
        }
        void delete_invalid_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_b_type = in_boundary_type( E.entity_type() ) ;
            gme_t invalid( in_b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( in_b_type ) ) {
                return ;
            } else {
                GeoModelMeshEntityAccess gmme_access( E ) ;
                remove_invalid_values( gmme_access.modifiable_in_boundaries(),
                    invalid ) ;
            }
        }
        void delete_invalid_parents( GeoModelMeshEntity& E ) ;
        void delete_invalid_signs( Region& R )
        {
            GeoModelMeshEntityAccess region_access(
                geomodel_access_.modifiable_mesh_entity( R.gme_id() ) ) ;
            region_access.modifiable_sides().resize( R.nb_boundaries() ) ;
        }
        void delete_invalid_universe_sided_boundaries( Universe& U )
        {
            const EntityType& b_type = Surface::type_name_static() ;
            gme_t invalid( b_type, NO_ID ) ;
            UniverseAccess universe_access( U ) ;
            remove_invalid_values( universe_access.modifiable_boundaries(),
                invalid ) ;
            universe_access.modifiable_sides().resize( U.nb_boundaries() ) ;
        }

        index_t mesh_entity_type_to_index( const EntityType& type ) const
        {
            return find( mesh_entity_types_, type ) ;
        }

        index_t geological_entity_type_to_index( const EntityType& type ) const
        {
            return find( geological_entity_types_, type ) ;
        }
        const EntityType& index_to_mesh_entity_type( index_t index ) const
        {
            return mesh_entity_types_.at( index ) ;
        }

        const EntityType& index_to_geological_entity_type( index_t index ) const
        {
            return geological_entity_types_.at( index ) ;
        }

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

        index_t nb_entity_types_ ;
        index_t nb_geological_entity_types_ ;
        index_t nb_mesh_entity_types_ ;

        std::vector< index_t > nb_initial_mesh_entities_ ;
        std::vector< index_t > nb_initial_geological_entities_ ;

        std::vector< index_t > nb_removed_mesh_entities_ ;
        std::vector< index_t > nb_removed_geological_entities_ ;

        /*! For each type of entity, store a vector of where the
         * entities to remove from the geomodel are flagged with NO_ID. */
        std::vector< std::vector< bool > > mesh_entity_to_erase_ ;
        /*! Stores the mapping table between indices for each type of
         *  element before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_mesh_entity_ ;
        std::vector< std::vector< index_t > > nb_childs_ ;

        std::vector< std::vector< index_t > > old_2_new_geological_entity_ ;

        //std::map< EntityType, index_t > entity_type_to_index_ ;
        //std::map< index_t, EntityType > index_to_entity_type_ ;
//        std::vector< EntityType > all_entity_types_ ;

        std::vector< EntityType > mesh_entity_types_ ;
        std::vector< EntityType > geological_entity_types_ ;

        EntityType starting_dependency_ ;
    } ;

}
