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

    void GeoModelEditor::remove_entities_and_dependencies(
        const std::set< gme_t >& entities_to_remove )
    {
        // Asserts to remove when implementation is completed
        ringmesh_assert(
            entities_to_remove.size() == 1
            && entities_to_remove.begin()->type == Region::type_name_static() ) ;

        // Copy because it is not logical to have in output the removed elements. BC
        /// @todo Review : youpiii what did the comment on the function just said ? [JP]
        std::set< gme_t > entities = entities_to_remove ;
        // TODO Handle the case of several objects in elements

        const GeoModelMeshEntity& reg = model_.mesh_entity( *( entities.begin() ) ) ;
        get_dependent_entities( entities ) ;

        // TODO Dirty duplication of code------------------------
        // We need to remove elements type by type since they are
        // stored in different vectors and since we use indices in these
        // vectors to identify them.
        // Initialize the vector
        std::vector< std::vector< index_t > > to_erase_by_type ;
        const index_t nb_type = 7 ; // = GME::NO_TYPE. TODO to rework
        to_erase_by_type.reserve( nb_type ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_mesh_entities( Corner::type_name_static() ), 0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_mesh_entities( Line::type_name_static() ), 0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_mesh_entities( Surface::type_name_static() ), 0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_mesh_entities( Region::type_name_static() ), 0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_geological_entities( Contact::type_name_static() ), 0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_geological_entities( Interface::type_name_static() ),
                0 ) ) ;
        to_erase_by_type.push_back(
            std::vector< index_t >(
                model_.nb_geological_entities( Layer::type_name_static() ), 0 ) ) ;

        // Flag the elements to erase
        for( std::set< gme_t >::const_iterator it = entities.begin();
            it != entities.end(); ++it ) {
            gme_t cur = *it ;
            ringmesh_assert( NO_ID != 0 ) ; // If one day NO_ID changes of value.
            if( cur.type == Corner::type_name_static() ) {
                to_erase_by_type[0][cur.index] = NO_ID ;
            } else if( cur.type == Line::type_name_static() ) {
                to_erase_by_type[1][cur.index] = NO_ID ;
            } else if( cur.type == Surface::type_name_static() ) {
                to_erase_by_type[2][cur.index] = NO_ID ;
            } else if( cur.type == Region::type_name_static() ) {
                to_erase_by_type[3][cur.index] = NO_ID ;
            } else if( cur.type == Contact::type_name_static() ) {
                to_erase_by_type[4][cur.index] = NO_ID ;
            } else if( cur.type == Interface::type_name_static() ) {
                to_erase_by_type[5][cur.index] = NO_ID ;
            } else if( cur.type == Layer::type_name_static() ) {
                to_erase_by_type[6][cur.index] = NO_ID ;
            }
        }

        // Number of elements deleted for each TYPE
        std::vector< index_t > nb_removed( to_erase_by_type.size(), 0 ) ;

        /// 1. Get the mapping between old indices of the elements
        ///    and new ones (when elements to remove will actually be removed)
        for( index_t i = 0; i < to_erase_by_type.size(); ++i ) {
            for( index_t j = 0; j < to_erase_by_type[i].size(); ++j ) {
                if( to_erase_by_type[i][j] == NO_ID ) {
                    nb_removed[i]++ ;
                } else {
                    to_erase_by_type[i][j] = j - nb_removed[i] ;
                }
            }
        }
        // TODO Dirty duplication of code--------------------------

        std::vector< gme_t > to_add_in_universe ;

        if( reg.type_name() == Region::type_name_static() ) {
            index_t nb_added = 0 ;
            for( index_t b_i = 0; b_i < reg.nb_boundaries(); ++b_i ) {
                if( !reg.boundary( b_i ).is_on_voi() ) {
                    to_add_in_universe.push_back( reg.boundary( b_i ).gme_id() ) ;
                    if( reg.boundary( b_i ).type_name()
                        == Corner::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[0][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[0][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Line::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[1][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[1][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Surface::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[2][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[2][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Region::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[3][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[3][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Contact::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[4][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[4][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Interface::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[5][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[5][reg.boundary( b_i ).index()] ;
                    } else if( reg.boundary( b_i ).type_name()
                        == Layer::type_name_static() ) {
                        ringmesh_assert(
                            to_erase_by_type[6][reg.boundary(
                                b_i ).index()] != NO_ID ) ;
                        to_add_in_universe[nb_added].index =
                            to_erase_by_type[6][reg.boundary( b_i ).index()] ;
                    }
                    ++nb_added ;
                }
            }
        }

        remove_entities( entities ) ;

        // Update Universe
        /// @todo You first need to clean the existing universe [JP]
        for( std::vector< gme_t >::const_iterator itr = to_add_in_universe.begin();
            itr != to_add_in_universe.end(); ++itr ) {
            add_universe_boundary( itr->index, true ) ;
        }

        //ringmesh_assert( model_.check_model_validity() ) ;
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
            nb_removed_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            nb_removed_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            initialize_maps() ;

            fill_nb_initial_entities() ;
            initialize_costly_storage() ;
        }
        /*!
         * @brief Removes the given Geological entities from the model
         * @warning ONLY takes care of deleting these entities and update
         * all references ( gme indices ) all over the model.
         * The client MUST:
         *    - ensure that the provided set is consistent to ensure
         *      the GeoModel validity.
         *    - ensure that eventual new connections between remaining entities
         *      are set to ensure the GeoModel validity
         */
        void remove_geological_entities( const std::set< gme_t >& entities )
        {
            if( verify_geological_entities( entities ) ) {
                throw RINGMeshException( "REMOVE",
                    "You try try to remove a mesh entity using geological entity method" ) ;
            }
            initialize_for_removal_geological_entities( entities ) ;
            do_delete_flagged_geological_entities() ;
            clear_model_mesh_vertices() ;
            update_entity_connectivity() ;
        }
        /*!
         * @brief Removes the given Geological entities from the model
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
            if( verify_geological_entities( entities ) ) {
                throw RINGMeshException( "REMOVE",
                    "You try try to remove a geological entity using mesh entity method" ) ;
            }
            initialize_for_removal_mesh_entities( entities ) ;
            do_delete_flagged_mesh_entities() ;
            clear_model_mesh_vertices() ;
            update_entity_connectivity() ;
        }

    private:
        // ---  High level functions ----------
        bool verify_geological_entities( const std::set< gme_t >& entities )
        {
            for( std::set< gme_t >::const_iterator it = entities.begin();
                it != entities.end(); ++it ) {
                gme_t cur = *it ;
                const EntityTypeManager& manager = model().entity_type_manager() ;
                if( !manager.is_geological_entity_type( cur.type ) ) {
                    return false ;
                }
            }
            return true ;
        }

        bool verify_mesh_entities( const std::set< gme_t >& entities )
        {
            for( std::set< gme_t >::const_iterator it = entities.begin();
                it != entities.end(); ++it ) {
                gme_t cur = *it ;
                const EntityTypeManager& manager = model().entity_type_manager() ;
                if( manager.is_geological_entity_type( cur.type ) ) {
                    return false ;
                }
            }
            return true ;
        }

        void initialize_for_removal_geological_entities(
            const std::set< gme_t >& entities_to_remove )
        {
            fill_geological_entities_to_erase_vector( entities_to_remove ) ;
            fill_removed_geological_entities_and_mapping() ;
        }

        void initialize_for_removal_mesh_entities(
            const std::set< gme_t >& entities_to_remove )
        {
            fill_mesh_entities_to_erase_vector( entities_to_remove ) ;
            fill_removed_mesh_entities_and_mapping() ;
        }
        void do_delete_flagged_geological_entities()
        {

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType entity_type =
                    index_to_geological_type_.find( i )->second ;
                for( index_t j = 0; j < nb_initial_geological_entities_[i]; ++j ) {
                    if( geological_entities_to_erase_[i][j] == NO_ID ) {
                        delete_geological_entity( entity_type, j ) ;
                    }
                }
                clear_null_geological_entities( entity_type ) ;
            }
        }

        void do_delete_flagged_mesh_entities()
        {

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType entity_type = index_to_mesh_type_.find( i )->second ;
                for( index_t j = 0; j < nb_initial_geological_entities_[i]; ++j ) {
                    if( mesh_entities_to_erase_[i][j] == NO_ID ) {
                        delete_mesh_entity( entity_type, j ) ;
                    }
                }
                clear_null_geological_entities( entity_type ) ;
            }
        }
        void clear_model_mesh_vertices()
        {
            model().mesh.vertices.clear() ;
        }
        void initialize_costly_storage()
        {
            mesh_entities_to_erase_.resize( nb_mesh_entity_types_ ) ;
            geological_entities_to_erase_.resize( nb_geological_entity_types_ ) ;
            old_2_new_mesh_entity_.resize( nb_mesh_entity_types_ ) ;
            old_2_new_geological_entity_.resize( nb_geological_entity_types_ ) ;

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                geological_entities_to_erase_[i].resize(
                    nb_initial_geological_entities_[i], 0 ) ;
                old_2_new_geological_entity_[i].resize( nb_geological_entity_types_,
                    0 ) ;
            }

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                mesh_entities_to_erase_[i].resize( nb_initial_mesh_entities_[i],
                    0 ) ;
                old_2_new_mesh_entity_[i].resize( nb_mesh_entity_types_, 0 ) ;
            }

        }

        void initialize_maps()
        {
            mesh_type_to_index_[Corner::type_name_static()] = 0 ;
            mesh_type_to_index_[Line::type_name_static()] = 1 ;
            mesh_type_to_index_[Surface::type_name_static()] = 2 ;
            mesh_type_to_index_[Region::type_name_static()] = 3 ;

            geological_type_to_index_[Contact::type_name_static()] = 0 ;
            geological_type_to_index_[Interface::type_name_static()] = 1 ;
            geological_type_to_index_[Layer::type_name_static()] = 2 ;

            index_to_mesh_type_[0] = Corner::type_name_static() ;
            index_to_mesh_type_[1] = Line::type_name_static() ;
            index_to_mesh_type_[2] = Surface::type_name_static() ;
            index_to_mesh_type_[3] = Region::type_name_static() ;

            index_to_geological_type_[0] = Contact::type_name_static() ;
            index_to_geological_type_[1] = Interface::type_name_static() ;
            index_to_geological_type_[2] = Layer::type_name_static() ;

        }

        void clear_null_geological_entities( const EntityType& type )
        {
            std::vector< GeoModelGeologicalEntity* >& store =
                modifiable_geological_entities( type ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GME* >( nil ) ), store.end() ) ;
        }
        void update_entity_connectivity()
        {

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& type = index_to_geological_type_.find( i )->second ;
                for( index_t j = 0; j < model().nb_geological_entities( type );
                    ++j ) {
                    gme_t new_id( type, j ) ;
                    GeoModelGeologicalEntity& GE = modifiable_geological_entity(
                        new_id ) ;
                    update_geological_entity_index( GE ) ;
                    ringmesh_assert( new_id == GE.gme_id() ) ;
                    update_entity_children( GE ) ;
                    delete_invalid_children( GE ) ;

                }
            }

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& type = index_to_mesh_type_.find( i )->second ;
                for( index_t j = 0; j < model().nb_mesh_entities( type ); ++j ) {
                    gme_t new_id( type, j ) ;
                    GeoModelMeshEntity& ME = modifiable_mesh_entity( new_id ) ;
                    update_mesh_entity_index( ME ) ;
                    ringmesh_assert( new_id == ME.gme_id() ) ;

                    update_entity_boundaries( ME ) ;
                    delete_invalid_boundaries( ME ) ;

                    update_entity_in_boundary( ME ) ;
                    delete_invalid_in_boundary( ME ) ;

                    update_entity_parents( ME ) ;
                    delete_invalid_parents( ME ) ;

                    if( ME.type_name() == Region::type_name_static() ) {
                        Region& R = dynamic_cast< Region& >( ME ) ;
                        update_region_boundary_signs( R ) ;
                        delete_invalid_signs( R ) ;
                    }
                }
            }

            Universe& U = universe() ;
            update_universe_sided_boundaries( U ) ;
            delete_invalid_universe_sided_boundaries( U ) ;
        }

        //------  Initialization -------
        void fill_removed_geological_entities_and_mapping()
        {
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_geological_entities_[i]; ++j ) {
                    if( geological_entities_to_erase_[i][j] == NO_ID ) {
                        nb_removed_geological_entities_[i]++ ;
                        old_2_new_geological_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_geological_entity_[i][j] = j
                            - nb_removed_geological_entities_[i] ;
                    }
                }
            }
        }

        void fill_removed_mesh_entities_and_mapping()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_mesh_entities_[i]; ++j ) {
                    if( mesh_entities_to_erase_[i][j] == NO_ID ) {
                        nb_removed_mesh_entities_[i]++ ;
                        old_2_new_mesh_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_mesh_entity_[i][j] = j
                            - nb_removed_mesh_entities_[i] ;
                    }
                }
            }
        }
        void fill_geological_entities_to_erase_vector(
            const std::set< gme_t >& entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it = entities_to_remove.begin();
                it != entities_to_remove.end(); ++it ) {
                gme_t cur = *it ;

                index_t type_index =
                    geological_type_to_index_.find( cur.type )->second ;
                geological_entities_to_erase_[type_index][cur.index] = NO_ID ;
            }
        }

        void fill_mesh_entities_to_erase_vector(
            const std::set< gme_t >& entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it = entities_to_remove.begin();
                it != entities_to_remove.end(); ++it ) {
                gme_t cur = *it ;

                index_t type_index = mesh_type_to_index_.find( cur.type )->second ;
                mesh_entities_to_erase_[type_index][cur.index] = NO_ID ;
            }
        }
        void fill_nb_initial_entities()
        {
            nb_initial_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& type = index_to_geological_type_[i] ;
                nb_initial_geological_entities_[i] = model().nb_geological_entities(
                    type ) ;
            }

            nb_initial_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& type = index_to_mesh_type_.find( i )->second ;
                nb_initial_mesh_entities_[i] = model().nb_mesh_entities( type ) ;
            }
        }

        // ---- Easier access to relationships between EntityTypes
        index_t mesh_entity_type_index( const GeoModelMeshEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return mesh_type_to_index_.find( type )->second ;
        }

        index_t geological_entity_type_index(
            const GeoModelGeologicalEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return geological_type_to_index_.find( type )->second ;
        }
        index_t children_type_index( const EntityType& type ) const
        {
            const EntityType& child_type = children_type( type ) ;
            return mesh_type_to_index_.find( child_type )->second ;
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
                return mesh_type_to_index_.find( b_type )->second ;
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
                return mesh_type_to_index_.find( in_b_type )->second ;
            }
        }
        const EntityType& in_boundary_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = model().entity_type_manager() ;
            return family.in_boundary_type( type ) ;
        }

        // ----  Update connectivity functions  ------
        void update_geological_entity_index( GeoModelGeologicalEntity& E )
        {
            index_t old_id = E.index() ;
            index_t type = geological_entity_type_index( E ) ;
            index_t new_id = old_2_new_geological_entity_[type][old_id] ;
            ringmesh_assert( new_id != NO_ID ) ;
            set_entity_index( E, new_id ) ;
        }

        void update_mesh_entity_index( GeoModelMeshEntity& E )
        {
            index_t old_id = E.index() ;
            index_t type = mesh_entity_type_index( E ) ;
            index_t new_id = old_2_new_mesh_entity_[type][old_id] ;
            ringmesh_assert( new_id != NO_ID ) ;
            set_entity_index( E, new_id ) ;
        }
        void update_entity_boundaries( GeoModelMeshEntity& E )
        {
            const EntityType& b_type = boundary_type( E.type_name() ) ;
            index_t type_index = mesh_type_to_index_.find( b_type )->second ;
            if( type_index == NO_ID ) {
                return ;
            }
            for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
                index_t old_boundary = E.boundary_gme( i ).index ;
                index_t new_boundary =
                    old_2_new_mesh_entity_[type_index][old_boundary] ;
                set_mesh_entity_boundary( E.gme_id(), i, new_boundary ) ;
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
        void update_entity_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_b_type = in_boundary_type( E.type_name() ) ;
            bool valid_type = EntityTypeManager::is_defined_type( in_b_type ) ;
            if( !valid_type ) {
                return ;
            }
            index_t in_boundary_type_index = mesh_type_to_index_[in_b_type] ;
            for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
                index_t old_id = E.in_boundary_gme( i ).index ;
                index_t new_id =
                    old_2_new_mesh_entity_[in_boundary_type_index][old_id] ;
                set_mesh_entity_in_boundary( E.gme_id(), i, new_id ) ;
            }
        }
        void update_entity_parents( GeoModelMeshEntity& E )
        {
            const EntityTypeManager& family_tree = model().entity_type_manager() ;

            const std::vector< EntityType > parents = family_tree.parent_types(
                E.entity_type() ) ;
            for( index_t i = 0; i < parents.size(); ++i ) {
                const EntityType& parent_type = parents[i] ;
                index_t parent_type_index = geological_type_to_index_.find(
                    parent_type )->second ;

                index_t p_id = E.parent_id( parent_type ) ;
                index_t old_id = E.parent_gme( p_id ).index ;
                index_t new_id =
                    old_2_new_geological_entity_[parent_type_index][old_id] ;
                set_mesh_entity_parent( E.gme_id(), p_id,
                    gme_t( parent_type, new_id ) ) ;
            }
        }
        void update_entity_children( GeoModelGeologicalEntity& E )
        {
            if( E.nb_children() > 0 ) {
                const EntityType& child_type = children_type( E.type_name() ) ;
                index_t child_type_index = geological_type_to_index_.find(
                    child_type )->second ;
                for( index_t i = 0; i < E.nb_children(); ++i ) {
                    index_t old_id = E.child_gme( i ).index ;
                    index_t new_id =
                        old_2_new_geological_entity_[child_type_index][old_id] ;
                    set_geological_entity_child( E.gme_id(), i, new_id ) ;
                }
            }
        }
        void update_universe_sided_boundaries( Universe& U )
        {
            index_t b_type_index = mesh_type_to_index_.find(
                Surface::type_name_static() )->second ;
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

                const EntityType& child_type = children_type( E.type_name() ) ;
                gme_t invalid_child( child_type, NO_ID ) ;
                remove_invalid_values( modifiable_children( E ), invalid_child ) ;
            }
        }
        void delete_invalid_boundaries( GeoModelMeshEntity& E )
        {

            const EntityType& b_type = boundary_type( E.type_name() ) ;
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

    private:

        /// Number of Geological Entity types
        index_t nb_geological_entity_types_ ;

        /// Number of Mesh Entity types
        index_t nb_mesh_entity_types_ ;

        /// Number of Geological Entities before the removal process
        std::vector< index_t > nb_initial_geological_entities_ ;

        /// Number of Mesh Entities before the removal process
        std::vector< index_t > nb_initial_mesh_entities_ ;

        /// Number of removed Mesh Entities per type
        std::vector< index_t > nb_removed_mesh_entities_ ;

        /// Number of removed Geological Entities per type
        std::vector< index_t > nb_removed_geological_entities_ ;

        /*! For each type of Geological entity, store a vector of where the
         * entities to remove from the model are flagged with NO_ID. */
        std::vector< std::vector< index_t > > geological_entities_to_erase_ ;

        /*! For each type of Mesh entity, store a vector of where the
         * entities to remove from the model are flagged with NO_ID. */
        std::vector< std::vector< index_t > > mesh_entities_to_erase_ ;

        /*! Stores the mapping table between indices for each type of
         *  geological entitiy before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_geological_entity_ ;

        /*! Stores the mapping table between indices for each type of
         *  mesh entitiy before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_mesh_entity_ ;

        /// map to convert Geological Type to index
        /// Contact -> 0
        /// Interface -> 1
        /// Layer -> 2
        TypeToIndex geological_type_to_index_ ;

        /// map to convert Mesh Type to index
        /// Corner -> 0
        /// Line -> 1
        /// Surface -> 2
        /// Region -> 3
        TypeToIndex mesh_type_to_index_ ;

        /// map to convert index to Geological Type
        /// 0 -> Contact
        /// 1 -> Interface
        /// 2 -> Layer
        IndexToType index_to_geological_type_ ;

        /// map to convert index to Mesh Type
        /// 0 -> Corner
        /// 1 -> Contact
        /// 2 -> Surface
        /// 3 -> Region
        IndexToType index_to_mesh_type_ ;
    } ;

    /*!
     * @brief Remove a list of entities of the model
     * @details No check is done on the consistency of this removal
     *          The entities and all references to them are removed.
     *          All dependent entities should be in the set of entities to remove,
     *          with a prior call to get_dependent_entities function.
     *
     * @warning NOT FULLY TESTED.
     *          The client is responsible to set the proper connectivity
     *          information between the remaining model entities.
     */
    void GeoModelEditor::remove_entities( const std::set< gme_t >& entities )
    {
        if( entities.empty() ) {
            return ;
        } else {
            throw RINGMeshException( "REMOVE",
                "Entity removal is not fully implemented" ) ;
//            GeoModelEntityRemoval remover( model() ) ;
//            remover.remove_entities( entities ) ;
        }
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model entities and their relationship ignoring their geometry
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
