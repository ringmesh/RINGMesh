/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geomodel/geo_model_editor.h>

#include <algorithm>
#include <vector> 
#include <map>
#include <set>

#include <ringmesh/geomodel/geo_model.h>

/*!
 * @file Implementation of the GeoModelEditor
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    typedef std::string EntityType ;

    GeoModelEditor::GeoModelEditor( GeoModel& model )
        : model_( model ), create_entity_allowed_( true )
    {
    }
    GeoModelEditor::~GeoModelEditor()
    {
    }

    GeoModelGeologicalEntity* GeoModelEditor::create_geological_entity(
        const EntityType& type,
        index_t index_in_geomodel )
    {
        GeoModelGeologicalEntity* E = GeoModelGeologicalEntityFactory::create_object(
            type, model() ) ;
        E->id_.index = index_in_geomodel ;
        return E ;
    }

    bool GeoModelEditor::create_mesh_entities(
        const EntityType& type,
        index_t nb_additional_entities )
    {
        if( EntityTypeManager::is_corner( type ) ) {
            return create_mesh_entities< Corner >( nb_additional_entities ) ;
        } else if( EntityTypeManager::is_line( type ) ) {
            return create_mesh_entities< Line >( nb_additional_entities ) ;
        } else if( EntityTypeManager::is_surface( type ) ) {
            return create_mesh_entities< Surface >( nb_additional_entities ) ;
        } else {
            // Must be regions.
            return create_mesh_entities< Region >( nb_additional_entities ) ;
        }
        ringmesh_assert_not_reached ;
    }

    void GeoModelEditor::set_model_name( const std::string& name )
    {
        model_.geomodel_name_ = name ;
    }

    /*!
     * @brief Creates and store an entity of the given type
     * @return The index of the created mesh entity
     */
    gme_t GeoModelEditor::create_geological_entity( const EntityType& type )
    {
        index_t index = find_or_create_geological_entity_type( type ) ;
        index_t id =
            static_cast< index_t >( model_.geological_entities_[index].size() ) ;
        GeoModelGeologicalEntity* E = create_geological_entity( type, id ) ;
        model_.geological_entities_[index].push_back( E ) ;
        return E->gme_id() ;
    }

    index_t GeoModelEditor::create_geological_entity_type( const EntityType& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory::has_creator( type ) ) ;
        entity_type_manager().geological_entity_types_.push_back( type ) ;
        model_.geological_entities_.push_back(
            std::vector< GeoModelGeologicalEntity* >() ) ;
        GeoModelGeologicalEntity* E = GeoModelGeologicalEntityFactory::create_object(
            type, model() ) ;

        const EntityType child_type = E->child_type_name() ;

        EntityTypeManager& parentage = entity_type_manager() ;
        parentage.register_relationship( type, child_type ) ;

        return entity_type_manager().nb_geological_entity_types() - 1 ;
    }

    bool GeoModelEditor::create_geological_entities(
        const EntityType& type,
        index_t nb_additional_entities )
    {
        find_or_create_geological_entity_type( type ) ;
        std::vector< GeoModelGeologicalEntity* >& store =
            modifiable_geological_entities( type ) ;
        index_t old_size = static_cast< index_t >( store.size() ) ;
        index_t new_size = old_size + nb_additional_entities ;
        store.resize( new_size, nil ) ;
        for( index_t i = old_size; i < new_size; i++ ) {
            ringmesh_assert( store[i] == nil ) ;
            store[i] = create_geological_entity( type, i ) ;
        }
        return true ;
    }

    void GeoModelEditor::complete_entity_connectivity()
    {
        // Order is important
        complete_mesh_entity_connectivity< Line >() ;
        complete_mesh_entity_connectivity< Corner >() ;
        complete_mesh_entity_connectivity< Surface >() ;
        complete_mesh_entity_connectivity< Region >() ;

        // Geological entities
        for( index_t i = 0; i < model_.nb_geological_entity_types(); i++ ) {
            const EntityType& type = model().geological_entity_type( i ) ;
            if( model_.nb_geological_entities( type ) > 0 ) {
                if( model_.geological_entity( type, 0 ).nb_children() == 0 ) {
                    fill_geological_entities_children( type ) ;
                }
            }
        }
    }

    template< typename ENTITY >
    void GeoModelEditor::complete_mesh_entity_connectivity()
    {
        const EntityType& type = ENTITY::type_name_static() ;
        if( model_.nb_mesh_entities( type ) > 0 ) {
            const GeoModelMeshEntity& E = model_.mesh_entity( type, 0 ) ;
            if( E.nb_boundaries() == 0 ) {
                fill_mesh_entities_boundaries( type ) ;
            }
            if( E.nb_in_boundary() == 0 ) {
                fill_mesh_entities_in_boundaries( type ) ;
            }
            if( E.nb_parents() == 0 ) {
                fill_mesh_entities_parent( type ) ;
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_boundaries( const EntityType& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& b_type = entity_type_manager().boundary_type( type ) ;
        if( EntityTypeManager::is_defined_type( b_type ) ) {
            for( index_t i = 0; i < model().nb_mesh_entities( b_type ); ++i ) {
                const GeoModelMeshEntity& b = mesh_entity( b_type, i ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_mesh_entity_boundary( b.in_boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_in_boundaries( const EntityType& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& in_b_type = entity_type_manager().in_boundary_type(
            type ) ;
        if( EntityTypeManager::is_defined_type( in_b_type ) ) {
            for( index_t i = 0; i < model().nb_mesh_entities( in_b_type ); ++i ) {
                const GeoModelMeshEntity& in_b = mesh_entity( in_b_type, i ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_mesh_entity_in_boundary( in_b.boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_parent( const EntityType& type )
    {
        const GeoModel& M = model() ;
        if( M.nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const std::vector< EntityType > parent_types(
            entity_type_manager().parent_types( type ) ) ;
        for( index_t i = 0; i < parent_types.size(); ++i ) {
            const EntityType& parent_type = parent_types[i] ;
            if( EntityTypeManager::is_defined_type( parent_type ) ) {
                for( index_t j = 0; j < M.nb_geological_entities( parent_type );
                    ++j ) {
                    const GeoModelGeologicalEntity& parent = geological_entity(
                        parent_type, j ) ;
                    for( index_t k = 0; k < parent.nb_children(); ++k ) {
                        add_mesh_entity_parent( parent.child_gme( k ),
                            parent.gme_id() ) ;
                    }
                }
            }
        }
    }

    void GeoModelEditor::fill_geological_entities_children( const EntityType& type )
    {
        if( model().nb_geological_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& c_type = geological_entity( type, 0 ).child_type_name() ;
        if( EntityTypeManager::is_defined_type( c_type ) ) {
            for( index_t i = 0; i < model().nb_mesh_entities( c_type ); ++i ) {
                const GeoModelMeshEntity& p = mesh_entity( c_type, i ) ;
                for( index_t j = 0; j < p.nb_parents(); j++ ) {
                    add_geological_entity_child( p.parent_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::complete_mesh_entities_geol_feature_from_first_parent(
        const EntityType& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) {
            return ;
        }
        const std::vector< EntityType > parents = entity_type_manager().parent_types(
            type ) ;
        if( parents.size() == 0 ) {
            return ;
        } else {
            for( index_t i = 0; i < model().nb_mesh_entities( type ); ++i ) {
                GeoModelMeshEntity& E = mesh_entity( type, i ) ;
                if( !E.has_geological_feature() ) {
                    if( E.nb_parents() > 0
                        && E.parent( 0 ).has_geological_feature() ) {
                        E.geol_feature_ = E.parent( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }

    void GeoModelEditor::complete_geological_entities_geol_feature_from_first_child(
        const EntityType& type )
    {
        if( model().nb_geological_entities( type ) == 0 ) {
            return ;
        }
        const EntityType& child_type = entity_type_manager().child_type( type ) ;
        if( EntityTypeManager::is_defined_type( child_type ) ) {
            for( index_t i = 0; i < model().nb_geological_entities( type ); ++i ) {
                GeoModelGeologicalEntity& p = geological_entity( type, i ) ;
                if( !p.has_geological_feature() ) {
                    if( p.nb_children() > 0
                        && p.child( 0 ).has_geological_feature() ) {
                        p.geol_feature_ = p.child( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }

    /*!
     * @brief Add to the vector the entities which cannot exist if
     *        an entity in the set does not exist.
     * @return True if at least one entity was added.
     */
    bool GeoModelEditor::get_dependent_entities( std::set< gme_t >& in ) const
    {
        std::size_t input_size = in.size() ;

        // Add children of geological entities
        for( std::set< gme_t >::iterator it( in.begin() ); it != in.end(); ++it ) {
            gme_t cur = *it ;
            if( entity_type_manager().is_geological_entity_type( cur.type ) ) {
                const GeoModelGeologicalEntity& E = model_.geological_entity( cur ) ;
                for( index_t j = 0; j < E.nb_children(); ++j ) {
                    in.insert( E.child_gme( j ) ) ;
                }
            }
        }
        // Add geological entities which have no child
        index_t nb_geological_entity_types =
            entity_type_manager().nb_geological_entity_types() ;
        for( index_t i = 0; i < nb_geological_entity_types; ++i ) {
            const EntityType& type = entity_type_manager().geological_entity_type(
                i ) ;

            for( index_t j = 0; j < model_.nb_geological_entities( type ); ++j ) {
                bool no_child = true ;
                const GeoModelGeologicalEntity& E = model_.geological_entity( type,
                    j ) ;
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
                for( index_t j = 0; j < model_.nb_mesh_entities( type ); ++j ) {
                    bool no_incident = true ;
                    const GeoModelMeshEntity& E = model_.mesh_entity( type, j ) ;
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
        // Recusive call till nothing is added
        if( in.size() != input_size ) {
            return get_dependent_entities( in ) ;
        } else {
            return false ;
        }
    }

    /*!
     * @brief Class in charge of removing entities from a GeoModel
     */
    class GeoModelEntityRemoval: public GeoModelEditor {
    public:
        typedef std::string EntityType ;
        typedef std::map< EntityType, index_t > TypeToIndex ;
        typedef std::map< index_t, EntityType > IndexToType ;

        GeoModelEntityRemoval( GeoModel& model )
            : GeoModelEditor( model )
        {
            nb_mesh_entity_types_ = EntityTypeManager::nb_mesh_entity_types() ;
            nb_geological_entity_types_ =
                GeoModelEditor::model().nb_geological_entity_types() ;
            nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_ ;
            nb_removed_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            nb_removed_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            fill_entity_type_to_index_map() ;
            fill_nb_initial_entities() ;
            initialize_costly_storage() ;
            fill_nb_children_vector() ;

        }
        /*!
         * @brief Removes the given entities from the model
         * @warning ONLY takes care of deleting these entities and update
         * all references ( gme indices ) all over the model.
         * The client MUST:
         *    - ensure that the provided set is consistent to ensure
         *      the GeoModel validity.
         *    - ensure that eventual new connections between remaining entities
         *      are set to ensure the GeoModel validity
         */
        void remove_mesh_entities( const std::set< gme_t >& entities )
        {
            initialize_for_removal( entities ) ;
            do_delete_flagged_mesh_entities() ;
            clear_model_mesh_vertices() ;
            update_mesh_entity_connectivity() ;
            flag_geological_entities_without_children() ;
            do_delete_flagged_geological_entities() ;
            update_geological_entity_connectivity() ;
            update_universe() ;

        }

//        void remove_mesh_entities_with_dependencies(
//            const std::set< gme_t >& entities )
//        {
//            remove_mesh_entities( entities ) ;
//            remove_dependencies() ;
//        }
// TODO it doesn't work ! Will be done with the new struct of
//      relations between GME

        void remove_geological_entities( const std::set< gme_t >& entities )
        {
            check_if_entities_are_not_meshed_entities( entities ) ;
            std::set< gme_t > mesh_entities ;
            for( std::set< gme_t >::const_iterator it = entities.begin();
                it != entities.end(); ++it ) {
                const GeoModelGeologicalEntity& cur_gmge = geological_entity( *it ) ;
                for( index_t i = 0; i < cur_gmge.nb_children(); i++ ) {
                    mesh_entities.insert( cur_gmge.child( i ).gme_id() ) ;
                }
            }
            remove_mesh_entities( mesh_entities ) ;
        }

//        void remove_geological_entities_with_dependencies(
//            const std::set< gme_t >& entities )
//        {
//            remove_geological_entities( entities ) ;
//            remove_dependencies() ;
//        }
// TODO it doesn't work ! Will be done with the new struct of
//      relations between GME

    private:
        // ---  High level functions ----------
        void initialize_for_removal(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
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
                            p < model().mesh_entity( type_name, j ).nb_parents();
                            p++ ) {
                            gme_t parent =
                                model().mesh_entity( type_name, j ).parent_gme( p ) ;
                            nb_childs_[geological_entity_type_to_index( parent.type )][parent.index]-- ;
                        }

                        delete_mesh_entity( i, j ) ;
                    }
                }
                clear_null_mesh_entities( i ) ;
            }
        }
        void do_delete_flagged_geological_entities()
        {
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_geological_entity_type(
                    i ) ;
                for( index_t j = 0;
                    j < model().nb_geological_entities( entity_type ); ++j ) {
                    if( old_2_new_geological_entity_[i][j] == NO_ID ) {
                        delete_geological_entity( entity_type, j ) ;
                    }
                }
                clear_null_geological_entities( i ) ;
            }
        }

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
        void clear_model_mesh_vertices()
        {
            model().mesh.vertices.clear() ;
        }
        void initialize_costly_storage()
        {
            mesh_entity_to_erase_.resize( nb_mesh_entity_types_ ) ;

            old_2_new_mesh_entity_.resize( nb_mesh_entity_types_ ) ;
            old_2_new_geological_entity_.resize( nb_geological_entity_types_ ) ;
            nb_childs_.resize( nb_geological_entity_types_ ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                index_t size = model().nb_mesh_entities(
                    index_to_mesh_entity_type( i ) ) ;
                mesh_entity_to_erase_[i].resize( size, false ) ;
                old_2_new_mesh_entity_[i].resize( size, 0 ) ;
            }

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                index_t size = model().nb_geological_entities(
                    index_to_geological_entity_type( i ) ) ;
                old_2_new_geological_entity_[i].resize( size, 0 ) ;

                nb_childs_[i].resize( size, 0 ) ;

            }

        }
        void delete_mesh_entity( index_t type, index_t index )
        {
            const EntityType& type_name = index_to_mesh_entity_type( type ) ;
            gme_t id( type_name, index ) ;
            GeoModelEditor::delete_mesh_entity( type_name, index ) ;
        }
        void clear_null_mesh_entities( index_t type )
        {
            const EntityType& type_name = index_to_mesh_entity_type( type ) ;
            std::vector< GeoModelMeshEntity* >& store = modifiable_mesh_entities(
                type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelMeshEntity* >( nil ) ), store.end() ) ;

            // QC
            ringmesh_assert( model().nb_mesh_entities( type_name )
                == nb_initial_mesh_entities_[type] - nb_removed_mesh_entities_[type] ) ;
        }

        void clear_null_geological_entities( index_t type )
        {
            const EntityType& type_name = index_to_geological_entity_type( type ) ;
            std::vector< GeoModelGeologicalEntity* >& store =
                modifiable_geological_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelGeologicalEntity* >( nil ) ),
                store.end() ) ;

            // QC
            ringmesh_assert( model().nb_geological_entities( type_name )
                == nb_initial_geological_entities_[type] - nb_removed_geological_entities_[type] ) ;
        }
        void update_mesh_entity_connectivity()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < model().nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = modifiable_mesh_entity( new_id ) ;
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
                    j < model().nb_geological_entities( entity_type ); ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelGeologicalEntity& GE = modifiable_geological_entity(
                        new_id ) ;
                    update_geological_entity_index( GE ) ;
                    update_geological_entity_children( GE ) ;
                    delete_invalid_children( GE ) ;
                }
            }

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < model().nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = modifiable_mesh_entity( new_id ) ;
                    update_mesh_entity_parents( ME ) ;
                    delete_invalid_parents( ME ) ;
                }
            }
        }

        void update_universe()
        {
            Universe& U = universe() ;
            update_universe_sided_boundaries( U ) ;
            delete_invalid_universe_sided_boundaries( U ) ;
        }

//        void remove_dependencies()
//        {
//            std::set< gme_t > new_gmme_to_remove ;
//            for( index_t me = 0;
//                me < model().nb_mesh_entities( starting_dependency_ ); me++ ) {
//                const GeoModelMeshEntity& cur_gmme = model().mesh_entity(
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

        void fill_nb_children_vector()
        {
            for( index_t i = 0; i < nb_childs_.size(); i++ ) {
                for( index_t j = 0; j < nb_childs_[i].size(); j++ ) {
                    nb_childs_[i][j] = model().geological_entity(
                        index_to_geological_entity_type( i ), j ).nb_children() ;
                }
            }
        }

        void fill_nb_initial_entities()
        {
            nb_initial_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& type = index_to_mesh_entity_type( i ) ;
                nb_initial_mesh_entities_[i] = model().nb_mesh_entities( type ) ;
            }

            nb_initial_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& type = index_to_geological_entity_type( i ) ;
                nb_initial_geological_entities_[i] = model().nb_geological_entities(
                    type ) ;
            }
        }
        void fill_entity_type_to_index_map()
        {
            const EntityTypeManager& manager = entity_type_manager() ;
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
            const GeoModelGeologicalEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return geological_entity_type_to_index( type ) ;
        }
        index_t children_type_index( const EntityType& type ) const
        {
            const EntityType& child_type = children_type( type ) ;
            return mesh_entity_type_to_index( child_type ) ;
        }
        const EntityType children_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = model().entity_type_manager() ;
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
            const EntityTypeManager& family = model().entity_type_manager() ;
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
            const EntityTypeManager& family = model().entity_type_manager() ;
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

        void update_mesh_entity_index( GeoModelMeshEntity& ME )
        {
            index_t old_id = ME.index() ;
            index_t type = mesh_entity_type_index( ME ) ;
            index_t new_id = old_2_new_mesh_entity_[type][old_id] ;
            ringmesh_assert( new_id != NO_ID ) ;
            set_entity_index( ME, new_id ) ;
        }

        void update_geological_entity_index( GeoModelGeologicalEntity& GE )
        {
            index_t old_id = GE.index() ;
            index_t type = geological_entity_type_index( GE ) ;
            index_t new_id = old_2_new_geological_entity_[type][old_id] ;
            ringmesh_assert( new_id != NO_ID ) ;
            set_entity_index( GE, new_id ) ;
        }

        void update_mesh_entity_boundaries( GeoModelMeshEntity& ME )
        {
            index_t type_index = boundary_type_index( ME.entity_type() ) ;
            if( type_index == NO_ID ) {
                return ;
            }
            for( index_t i = 0; i < ME.nb_boundaries(); ++i ) {
                index_t old_boundary = ME.boundary_gme( i ).index ;
                index_t new_boundary =
                    old_2_new_mesh_entity_[type_index][old_boundary] ;
                set_mesh_entity_boundary( ME.gme_id(), i, new_boundary ) ;
            }
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
                    set_boundary_sign( R, i, new_side ) ;
                }
            }
        }
        void update_mesh_entity_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_boundary_type = EntityTypeManager::in_boundary_type(
                E.entity_type() ) ;
            bool valid_type = EntityTypeManager::is_defined_type(
                in_boundary_type ) ;
            if( !valid_type ) {
                return ;
            }
            index_t in_boundary_type_index = mesh_entity_type_to_index(
                in_boundary_type ) ;
            for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
                index_t old_id = E.in_boundary_gme( i ).index ;
                index_t new_id =
                    old_2_new_mesh_entity_[in_boundary_type_index][old_id] ;
                set_mesh_entity_in_boundary( E.gme_id(), i, new_id ) ;
            }
        }
        void update_mesh_entity_parents( GeoModelMeshEntity& E )
        {
            for( index_t p = 0; p < E.nb_parents(); ++p ) {
                const EntityType& parent_type = E.parent_gme( p ).type ;
                index_t parent_type_index = geological_entity_type_to_index(
                    parent_type ) ;

                index_t old_id = E.parent_gme( p ).index ;
                index_t new_id =
                    old_2_new_geological_entity_[parent_type_index][old_id] ;
                set_mesh_entity_parent( E.gme_id(), p,
                    gme_t( parent_type, new_id ) ) ;
            }
        }
        void update_geological_entity_children( GeoModelGeologicalEntity& E )
        {
            if( E.nb_children() > 0 ) {
                index_t child_type = children_type_index( E.entity_type() ) ;
                for( index_t i = 0; i < E.nb_children(); ++i ) {
                    index_t old_id = E.child_gme( i ).index ;
                    index_t new_id = old_2_new_mesh_entity_[child_type][old_id] ;
                    set_geological_entity_child( E.gme_id(), i, new_id ) ;
                }
            }
        }
        void update_universe_sided_boundaries( Universe& U )
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
                set_universe_boundary( i, new_id, new_side ) ;
            }
        }

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
        void delete_invalid_children( GeoModelGeologicalEntity& E )
        {
            if( E.nb_children() == 0 ) {
                return ;
            } else {
                const EntityType& child_type = children_type( E.entity_type() ) ;
                gme_t invalid_child( child_type, NO_ID ) ;
                remove_invalid_values( modifiable_children( E ), invalid_child ) ;
            }
        }
        void delete_invalid_boundaries( GeoModelMeshEntity& E )
        {
            const EntityType& b_type = boundary_type( E.entity_type() ) ;
            gme_t invalid( b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( b_type ) ) {
                return ;
            } else {
                remove_invalid_values( modifiable_boundaries( E ), invalid ) ;
            }
        }
        void delete_invalid_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_b_type = in_boundary_type( E.entity_type() ) ;
            gme_t invalid( in_b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( in_b_type ) ) {
                return ;
            } else {
                remove_invalid_values( modifiable_in_boundaries( E ), invalid ) ;
            }
        }
        void delete_invalid_parents( GeoModelMeshEntity& E )
        {
            //  Cannot use remove directly, do it by hand like the signs
            index_t offset = 0 ;
            index_t new_size = 0 ;
            for( index_t i = 0; i + offset < E.nb_parents(); ++i ) {
                if( E.parent( i ).index() == NO_ID ) {
                    offset++ ;
                } else {
                    gme_t new_id = E.parent_gme( i + offset ) ;
                    set_mesh_entity_parent( E.gme_id(), i, new_id ) ;
                }
                new_size = i + 1 ; /// @todo Check that this is the correct size
            }
            modifiable_parents( E ).resize( new_size ) ;
        }
        void delete_invalid_signs( Region& R )
        {
            modifiable_sides( R ).resize( R.nb_boundaries() ) ;
        }
        void delete_invalid_universe_sided_boundaries( Universe& U )
        {
            const EntityType& b_type = Surface::type_name_static() ;
            gme_t invalid( b_type, NO_ID ) ;
            remove_invalid_values( modifiable_boundaries( U ), invalid ) ;
            modifiable_sides( U ).resize( U.nb_boundaries() ) ;
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
        index_t nb_entity_types_ ;
        index_t nb_geological_entity_types_ ;
        index_t nb_mesh_entity_types_ ;

        std::vector< index_t > nb_initial_mesh_entities_ ;
        std::vector< index_t > nb_initial_geological_entities_ ;

        std::vector< index_t > nb_removed_mesh_entities_ ;
        std::vector< index_t > nb_removed_geological_entities_ ;

        /*! For each type of entity, store a vector of where the
         * entities to remove from the model are flagged with NO_ID. */
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

    /*!
     * @brief Remove a list of mesh entities of the model
     * @details No check is done on the consistency of this removal
     *          The entities and all references to them are removed.
     *          All dependent entities should be in the set of entities to remove,
     *          with a prior call to get_dependent_entities function.
     *
     */
    void GeoModelEditor::remove_mesh_entities( const std::set< gme_t >& entities )
    {
        if( entities.empty() ) {
            return ;
        } else {
            GeoModelEntityRemoval remover( model() ) ;
            remover.remove_mesh_entities( entities ) ;
        }
    }

    /*!
     * @brief Remove a list of geological entities of the model
     * @details No check is done on the consistency of this removal
     *          The entities and all references to them are removed.
     *          All dependent entities should be in the set of entities to remove,
     *          with a prior call to get_dependent_entities function.
     *
     */
    void GeoModelEditor::remove_geological_entities(
        const std::set< gme_t >& entities )
    {
        if( entities.empty() ) {
            return ;
        } else {
            GeoModelEntityRemoval remover( model() ) ;
            remover.remove_geological_entities( entities ) ;
        }
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy all the model entities and their relationship ignoring their geometry
     *
     * @param[in] from Model to copy the information from
     */
    void GeoModelEditor::copy_macro_topology( const GeoModel& from )
    {
        copy_mesh_entity_topology< Corner >( from ) ;
        copy_mesh_entity_topology< Line >( from ) ;
        copy_mesh_entity_topology< Surface >( from ) ;
        copy_mesh_entity_topology< Region >( from ) ;

        for( index_t t = 0; t < from.nb_geological_entity_types(); t++ ) {
            copy_geological_entity_topology( from,
                from.geological_entity_type( t ) ) ;
        }

        model().universe_.copy( from.universe_ ) ;
        model().epsilon_ = from.epsilon() ;
    }

    template< typename ENTITY >
    void GeoModelEditor::copy_mesh_entity_topology( const GeoModel& from )
    {
        const EntityType& type = ENTITY::type_name_static() ;
        create_mesh_entities< ENTITY >( from.nb_mesh_entities( type ) ) ;

        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < model_.nb_mesh_entities( type ); ++e ) {
            gme_t id( type, e ) ;
            mesh_entity( id ).copy( from.mesh_entity( id ) ) ;
        }
    }

    void GeoModelEditor::copy_geological_entity_topology(
        const GeoModel& from,
        const EntityType& type )
    {
        create_geological_entities( type, from.nb_geological_entities( type ) ) ;

        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < model_.nb_geological_entities( type ); ++e ) {
            gme_t id( type, e ) ;
            geological_entity( id ).copy( from.geological_entity( id ) ) ;
        }
    }
}
