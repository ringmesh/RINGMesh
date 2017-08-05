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
    template< index_t DIMENSION > class GeoModelBuilderBase;
    template< index_t DIMENSION > class GeoModelBuilder;
}

namespace RINGMesh {

    /*!
     * @brief Builder tools to remove entities from a GeoModel
     */
    template< index_t DIMENSION >
    class GeoModelBuilderRemovalBase {
    ringmesh_disable_copy( GeoModelBuilderRemovalBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelBuilderRemovalBase() = default;

        /*!
         * @brief Remove a list of mesh entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_mesh_entities( const std::set< gmme_id >& entities );

        /*!
         * @brief Remove a list of geological entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_geological_entities( const std::set< gmge_id >& entities );

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gmme_id >& entities_to_remove );

    protected:
        GeoModelBuilderRemovalBase(
            GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

        virtual void update_mesh_entity( GeoModelMeshEntity< DIMENSION >& ME )
        {
            update_mesh_entity_index( ME );
            update_mesh_entity_boundaries( ME );
            delete_invalid_boundaries( ME );

            update_mesh_entity_incident_entity( ME );
            delete_invalid_incident_entity( ME );
        }

        // ---  High level functions ----------
        void initialize_for_removal(
            const std::set< gmme_id >& mesh_entities_to_remove )
        {
            nb_mesh_entity_types_ =
                geomodel_.entity_type_manager().mesh_entity_manager.nb_mesh_entity_types();
            nb_geological_entity_types_ = geomodel_.nb_geological_entity_types();
            nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_;
            nb_removed_mesh_entities_.resize( nb_mesh_entity_types_, 0 );
            nb_removed_geological_entities_.resize( nb_geological_entity_types_, 0 );
            fill_entity_type_to_index_map();
            fill_nb_initial_entities();
            initialize_costly_storage();
            fill_nb_children_vector();

            check_if_entities_are_meshed( mesh_entities_to_remove );
            fill_to_erase_vectors( mesh_entities_to_remove );
            fill_removed_entities_and_mapping();
        }
        void do_delete_flagged_mesh_entities()
        {
            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                for( index_t j : range( nb_initial_mesh_entities_[i] ) ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        const MeshEntityType& type_name = index_to_mesh_entity_type(
                            i );
                        for( index_t p : range(
                            geomodel_.mesh_entity( type_name, j ).nb_parents() ) ) {
                            gmge_id parent =
                                geomodel_.mesh_entity( type_name, j ).parent_gmge(
                                    p );
                            nb_childs_[geological_entity_type_to_index(
                                parent.type() )][parent.index()]--;
                        }

                        delete_mesh_entity( i, j );
                    }
                }
                clear_null_mesh_entities( i );
            }
        }
        void do_delete_flagged_geological_entities();

        void check_if_entities_are_meshed(
            const std::set< gmme_id >& mesh_entities_to_remove )
        {
            for( const gmme_id& it : mesh_entities_to_remove ) {
                if( !geomodel_.entity_type_manager().mesh_entity_manager.is_valid_type(
                    it.type() ) ) {
                    throw RINGMeshException( "REMOVE",
                        "You try to remove a Geological Entity using mesh removal." );
                }
            }

        }

        void initialize_costly_storage()
        {
            mesh_entity_to_erase_.resize( nb_mesh_entity_types_ );

            old_2_new_mesh_entity_.resize( nb_mesh_entity_types_ );
            old_2_new_geological_entity_.resize( nb_geological_entity_types_ );
            nb_childs_.resize( nb_geological_entity_types_ );
            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                index_t size = geomodel_.nb_mesh_entities(
                    index_to_mesh_entity_type( i ) );
                mesh_entity_to_erase_[i].resize( size, false );
                old_2_new_mesh_entity_[i].resize( size, 0 );
            }

            for( index_t i : range( nb_geological_entity_types_ ) ) {
                index_t size = geomodel_.nb_geological_entities(
                    index_to_geological_entity_type( i ) );
                old_2_new_geological_entity_[i].resize( size, 0 );
                nb_childs_[i].resize( size, 0 );
            }

        }
        void delete_mesh_entity( index_t type, index_t index );

        void clear_null_mesh_entities( index_t type )
        {
            const MeshEntityType& type_name = index_to_mesh_entity_type( type );
            std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >& store =
                geomodel_access_.modifiable_mesh_entities( type_name );
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >( nullptr ) ),
                store.end() );

            // QC
            ringmesh_assert(
                geomodel_.nb_mesh_entities( type_name )
                    == nb_initial_mesh_entities_[type]
                        - nb_removed_mesh_entities_[type] );
        }

        void clear_null_geological_entities( index_t type )
        {
            const GeologicalEntityType& type_name = index_to_geological_entity_type(
                type );
            std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& store =
                geomodel_access_.modifiable_geological_entities( type_name );
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< std::unique_ptr<
                        GeoModelGeologicalEntity< DIMENSION > > >( nullptr ) ),
                store.end() );

            // QC
            ringmesh_assert(
                geomodel_.nb_geological_entities( type_name )
                    == nb_initial_geological_entities_[type]
                        - nb_removed_geological_entities_[type] );
        }
        void update_mesh_entity_connectivity()
        {
            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                const MeshEntityType& entity_type = index_to_mesh_entity_type( i );

                for( index_t j : range( geomodel_.nb_mesh_entities( entity_type ) ) ) {
                    gmme_id new_id( entity_type, j );
                    GeoModelMeshEntity< DIMENSION >& ME =
                        geomodel_access_.modifiable_mesh_entity( new_id );
                    update_mesh_entity( ME );
                }
            }
        }

        void update_geological_entity_connectivity()
        {

            for( index_t i : range( nb_geological_entity_types_ ) ) {
                const GeologicalEntityType& entity_type =
                    index_to_geological_entity_type( i );
                for( index_t j : range(
                    geomodel_.nb_geological_entities( entity_type ) ) ) {
                    gmge_id new_id( entity_type, j );
                    GeoModelGeologicalEntity< DIMENSION >& GE =
                        geomodel_access_.modifiable_geological_entity( new_id );
                    update_geological_entity_index( GE );
                    update_geological_entity_children( GE );
                    delete_invalid_children( GE );
                }
            }

            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                const MeshEntityType& entity_type = index_to_mesh_entity_type( i );

                for( index_t j : range( geomodel_.nb_mesh_entities( entity_type ) ) ) {
                    gmme_id new_id( entity_type, j );
                    GeoModelMeshEntity< DIMENSION >& ME =
                        geomodel_access_.modifiable_mesh_entity( new_id );
                    update_mesh_entity_parents( ME );
                    delete_invalid_parents( ME );
                }
            }
        }

        //        void remove_dependencies()
        //        {
        //            std::set< gme_id > new_gmme_to_remove ;
        //            for( index_t me = 0;
        //                me < geomodel().nb_mesh_entities( starting_dependency_ ); me++ ) {
        //                const GeoModelMeshEntity& cur_gmme = geomodel().mesh_entity(
        //                    starting_dependency_, me ) ;
        //                if( cur_gmme.incident_entity( 0 ).index() == NO_ID
        //                    && cur_gmme.nb_incident_entities() == 1 ) {
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
            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                for( index_t j : range( nb_initial_mesh_entities_[i] ) ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        nb_removed_mesh_entities_[i]++;
                        old_2_new_mesh_entity_[i][j] = NO_ID;
                    } else {
                        old_2_new_mesh_entity_[i][j] = j
                            - nb_removed_mesh_entities_[i];
                    }
                }
            }
        }
        void fill_to_erase_vectors(
            const std::set< gmme_id >& mesh_entities_to_remove )
        {
            for( const gmme_id& cur : mesh_entities_to_remove ) {
                index_t type_index = mesh_entity_type_to_index( cur.type() );
                mesh_entity_to_erase_[type_index][cur.index()] = true;
            }
        }

        void fill_nb_children_vector();

        void fill_nb_initial_entities()
        {
            nb_initial_mesh_entities_.resize( nb_mesh_entity_types_, 0 );
            for( index_t i : range( nb_mesh_entity_types_ ) ) {
                const MeshEntityType& type = index_to_mesh_entity_type( i );
                nb_initial_mesh_entities_[i] = geomodel_.nb_mesh_entities( type );
            }

            nb_initial_geological_entities_.resize( nb_geological_entity_types_, 0 );
            for( index_t i : range( nb_geological_entity_types_ ) ) {
                const GeologicalEntityType& type = index_to_geological_entity_type(
                    i );
                nb_initial_geological_entities_[i] =
                    geomodel_.nb_geological_entities( type );
            }
        }
        void fill_entity_type_to_index_map()
        {
            const EntityTypeManager< DIMENSION >& manager =
                geomodel_.entity_type_manager();
            mesh_entity_types_.insert( mesh_entity_types_.end(),
                manager.mesh_entity_manager.mesh_entity_types().begin(),
                manager.mesh_entity_manager.mesh_entity_types().end() );

            if( nb_geological_entity_types_ != 0 ) {
                geological_entity_types_.insert( geological_entity_types_.end(),
                    manager.geological_entity_manager.geological_entity_types().begin(),
                    manager.geological_entity_manager.geological_entity_types().end() );
            }

        }

        // ---- Easier access to relationships between EntityTypes
        index_t mesh_entity_type_index(
            const GeoModelMeshEntity< DIMENSION >& E ) const
        {
            const MeshEntityType& type = E.type_name();
            return mesh_entity_type_to_index( type );
        }
        index_t geological_entity_type_index(
            const GeoModelGeologicalEntity< DIMENSION >& E ) const;
        index_t children_type_index( const GeologicalEntityType& type ) const
        {
            const MeshEntityType& child_type = children_type( type );
            return mesh_entity_type_to_index( child_type );
        }
        const MeshEntityType children_type( const GeologicalEntityType& type ) const
        {
            const RelationshipManager& family =
                geomodel_.entity_type_manager().relationship_manager;
            return family.child_type( type );
        }
        index_t boundary_type_index( const MeshEntityType& type ) const
        {
            const MeshEntityType& b_type = boundary_entity_type( type );
            if( !geomodel_.entity_type_manager().mesh_entity_manager.is_valid_type(
                b_type ) ) {
                return NO_ID;
            } else {
                return mesh_entity_type_to_index( b_type );
            }
        }
        const MeshEntityType& boundary_entity_type(
            const MeshEntityType& type ) const
        {
            const MeshEntityTypeManager< DIMENSION >& family =
                geomodel_.entity_type_manager().mesh_entity_manager;
            return family.boundary_entity_type( type );
        }
        /// TODO unused function. To handle during removal refactoring BC.
        index_t incident_entity_type_to_index( const MeshEntityType& type ) const
        {
            const MeshEntityType& in_ent_type = incident_entity_type( type );
            if( !geomodel_.entity_type_manager().mesh_entity_manager.is_valid_type(
                in_ent_type ) ) {
                return NO_ID;
            } else {
                return mesh_entity_type_to_index( in_ent_type );
            }
        }
        const MeshEntityType& incident_entity_type(
            const MeshEntityType& type ) const
        {
            const MeshEntityTypeManager< DIMENSION >& family =
                geomodel_.entity_type_manager().mesh_entity_manager;
            return family.incident_entity_type( type );
        }
        bool is_mesh_entity( index_t i ) const
        {
            return i < nb_mesh_entity_types_;
        }
        bool is_geological_entity( index_t i ) const
        {
            return !is_mesh_entity( i );
        }

        // ----  Update connectivity functions  ------

        void flag_geological_entities_without_children()
        {
            for( index_t i : range( nb_childs_.size() ) ) {
                for( index_t j : range( nb_childs_[i].size() ) ) {
                    if( nb_childs_[i][j] == 0 ) {
                        nb_removed_geological_entities_[i]++;
                        old_2_new_geological_entity_[i][j] = NO_ID;
                    } else {
                        old_2_new_geological_entity_[i][j] = j
                            - nb_removed_geological_entities_[i];
                    }

                }
            }
        }

        void set_mesh_entity_index(
            GeoModelMeshEntity< DIMENSION >& E,
            index_t new_index_in_geomodel );
        void set_geological_entity_index(
            GeoModelGeologicalEntity< DIMENSION >& E,
            index_t new_index_in_geomodel );

        void update_mesh_entity_index( GeoModelMeshEntity< DIMENSION >& ME );
        void update_geological_entity_index(
            GeoModelGeologicalEntity< DIMENSION >& GE );
        void update_mesh_entity_boundaries( GeoModelMeshEntity< DIMENSION >& ME );

        void set_boundary_side( Region3D& R, index_t boundary_index, bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() );
            GeoModelMeshEntityAccess< DIMENSION > region_access(
                geomodel_access_.modifiable_mesh_entity( R.gmme() ) );
            region_access.modifiable_sides()[boundary_index] = new_side;
        }

        void update_mesh_entity_incident_entity(
            GeoModelMeshEntity< DIMENSION >& E );
        void update_mesh_entity_parents( GeoModelMeshEntity< DIMENSION >& E );
        void update_geological_entity_children(
            GeoModelGeologicalEntity< DIMENSION >& E );

        // --- Deletion of some values the GeoModel storage
        template< typename TEST, typename THINGS_TO_DELETE >
        void remove_invalid_values(
            std::vector< THINGS_TO_DELETE >& vector,
            const TEST& test )
        {
            auto new_end = std::remove_if( vector.begin(), vector.end(), test );
            if( new_end == vector.begin() ) {
                // Clear instead of erase, because the behavior would be undefined.
                vector.clear();
            } else if( new_end < vector.end() ) {
                vector.erase( new_end, vector.end() );
            }
        }

        void delete_invalid_children( GeoModelGeologicalEntity< DIMENSION >& E )
        {
            if( E.nb_children() == 0 ) {
                return;
            } else {
                GeoModelGeologicalEntityAccess< DIMENSION > gmge_access( E );
                const RelationshipManager& manager =
                    E.geomodel().entity_type_manager().relationship_manager;
                const MeshEntityType& child_type = children_type( E.entity_type() );
                gmme_id invalid_child( child_type, NO_ID );
                remove_invalid_values( gmge_access.modifiable_children(),
                    [&invalid_child, &manager](index_t i) {return manager.boundary_gmme( i ) == invalid_child;} );
            }
        }

        void delete_invalid_boundaries( GeoModelMeshEntity< DIMENSION >& E )
        {
            const MeshEntityType& b_type = boundary_entity_type(
                E.mesh_entity_type() );
            gmme_id invalid( b_type, NO_ID );
            if( !geomodel_.entity_type_manager().mesh_entity_manager.is_valid_type(
                b_type ) ) {
                return;
            } else {
                GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
                const RelationshipManager& manager =
                    E.geomodel().entity_type_manager().relationship_manager;
                remove_invalid_values( gmme_access.modifiable_boundaries(),
                    [&invalid, &manager](index_t i) {return manager.boundary_gmme( i ) == invalid;} );
            }
        }

        void delete_invalid_incident_entity( GeoModelMeshEntity< DIMENSION >& E )
        {
            const MeshEntityType& in_ent_type = incident_entity_type(
                E.mesh_entity_type() );
            gmme_id invalid( in_ent_type, NO_ID );
            if( !geomodel_.entity_type_manager().mesh_entity_manager.is_valid_type(
                in_ent_type ) ) {
                return;
            } else {
                GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
                const RelationshipManager& manager =
                    E.geomodel().entity_type_manager().relationship_manager;
                remove_invalid_values( gmme_access.modifiable_incident_entities(),
                    [&invalid, &manager](index_t i) {return manager.incident_entity_gmme( i ) == invalid;} );
            }
        }

        void delete_invalid_parents( GeoModelMeshEntity< DIMENSION >& E )
        {
            GeoModelMeshEntityAccess< DIMENSION > gmme_access( E );
            const RelationshipManager& manager =
                E.geomodel().entity_type_manager().relationship_manager;
            remove_invalid_values( gmme_access.modifiable_parents(),
                [ &manager](index_t i) {return manager.parent_of_gmme( i ).index() == NO_ID;} );
        }

        index_t mesh_entity_type_to_index( const MeshEntityType& type ) const
        {
            return find( mesh_entity_types_, type );
        }

        index_t geological_entity_type_to_index(
            const GeologicalEntityType& type ) const
        {
            return find( geological_entity_types_, type );
        }
        const MeshEntityType& index_to_mesh_entity_type( index_t index ) const
        {
            return mesh_entity_types_[index];
        }

        const GeologicalEntityType& index_to_geological_entity_type(
            index_t index ) const
        {
            return geological_entity_types_[index];
        }

    protected:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;

        index_t nb_entity_types_;
        index_t nb_geological_entity_types_;
        index_t nb_mesh_entity_types_;

        std::vector< index_t > nb_initial_mesh_entities_;
        std::vector< index_t > nb_initial_geological_entities_;

        std::vector< index_t > nb_removed_mesh_entities_;
        std::vector< index_t > nb_removed_geological_entities_;

        /*! For each type of entity, store a vector of where the
         * entities to remove from the geomodel are flagged with NO_ID. */
        std::vector< std::vector< bool > > mesh_entity_to_erase_;
        /*! Stores the mapping table between indices for each type of
         *  element before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_mesh_entity_;
        std::vector< std::vector< index_t > > nb_childs_;

        std::vector< std::vector< index_t > > old_2_new_geological_entity_;

        //std::map< EntityType, index_t > entity_type_to_index_ ;
        //std::map< index_t, EntityType > index_to_entity_type_ ;
//        std::vector< EntityType > all_entity_types_ ;

        std::vector< MeshEntityType > mesh_entity_types_;
        std::vector< GeologicalEntityType > geological_entity_types_;
    };

    CLASS_DIMENSION_ALIASES( GeoModelBuilderRemovalBase );

    template< index_t DIMENSION >
    class GeoModelBuilderRemoval: public GeoModelBuilderRemovalBase< DIMENSION > {
        friend class GeoModelBuilderBase< DIMENSION > ;
        friend class GeoModelBuilder< DIMENSION > ;
    private:
        GeoModelBuilderRemoval(
            GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );
        virtual ~GeoModelBuilderRemoval() = default;
    };

    template< >
    class GeoModelBuilderRemoval< 3 > : public GeoModelBuilderRemovalBase< 3 > {
        friend class GeoModelBuilderBase< 3 > ;
        friend class GeoModelBuilder< 3 > ;
        using GeoModelBuilder3D = GeoModelBuilder< 3 >;
    private:
        GeoModelBuilderRemoval(
            GeoModelBuilder3D& builder,
            GeoModel3D& geomodel );
        virtual ~GeoModelBuilderRemoval() = default;

        void update_mesh_entity( GeoModelMeshEntity3D& ME ) override;

        void set_boundary_side(
            Region3D& R,
            index_t boundary_index,
            bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() );
            GeoModelMeshEntityAccess3D region_access(
                geomodel_access_.modifiable_mesh_entity( R.gmme() ) );
            region_access.modifiable_sides()[boundary_index] = new_side;
        }

        void update_region_boundary_signs( Region3D& R )
        {
            const MeshEntityType& surface_type = boundary_entity_type(
                R.mesh_entity_type() );
            gmme_id invalid_value( surface_type, NO_ID );

            index_t offset = 0;
            for( index_t i = 0; i + offset < R.nb_boundaries(); ++i ) {
                if( R.boundary_gmme( i ) == invalid_value ) {
                    offset++;
                } else {
                    bool new_side = R.side( i + offset );
                    set_boundary_side( R, i, new_side );
                }
            }
        }

        void delete_invalid_signs( Region3D& R )
        {
            GeoModelMeshEntityAccess3D region_access(
                geomodel_access_.modifiable_mesh_entity( R.gmme() ) );
            region_access.modifiable_sides().resize( R.nb_boundaries() );
        }
    };
}
