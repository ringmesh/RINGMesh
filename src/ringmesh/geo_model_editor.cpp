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

#include <ringmesh/geo_model_editor.h>

#include <algorithm>

#include <ringmesh/geo_model.h>

/*!
 * @file Implementation of the GeoModelEditor
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    typedef GeoModelEntity::gme_t gme_t ;

    GeoModelGeologicalEntity* GeoModelEditor::new_geological_entity(
        const std::string& type, index_t id )
    {
        GeoModelGeologicalEntity* E = GeoModelGeologicalEntityFactory::create_object(
            type, model() ) ;
        E->id_.index = id ;
        return E ;
    }

    void GeoModelEditor::set_model_name( const std::string& name )
    {
        model_.name_ = name ;
    }

    /*!
     * @brief Creates a mesh entity of the given type and add it to the correct vector
     * The GeoModelMeshEntity is created from its type and its index
     *
     * @param[in] type Type of the mesh entity to create
     * @return The index of the created mesh entity
     */
    gme_t GeoModelEditor::create_mesh_entity( const std::string& type )
    {
        assert_entity_creation_allowed() ;
        if( type == Corner::type_name_ ) {
            Corner* corner = new Corner( model(), model_.corners_.size() ) ;
            model_.corners_.push_back( corner ) ;
            return corner->gme_id() ;
        } else if( type == Line::type_name_ ) {
            Line* line = new Line( model(), model_.lines_.size() ) ;
            model_.lines_.push_back( line ) ;
            return line->gme_id() ;
        } else if( type == Surface::type_name_ ) {
            Surface* surface = new Surface( model(), model_.surfaces_.size() ) ;
            model_.surfaces_.push_back( surface ) ;
            return surface->gme_id() ;
        } else if( type == Region::type_name_ ) {
            Region* region = new Region( model(), model_.regions_.size() ) ;
            model_.regions_.push_back( region ) ;
            return region->gme_id() ;
        } else {
            ringmesh_assert_not_reached ;
            return gme_t() ;
        }
    }

    /*!
     * @brief Creates a mesh entity of the given type and add it to the correct vector
     * The GeoModelMeshEntity is created from its type and its index
     *
     * @param[in] type Type of the mesh entity to create
     * @return The index of the created mesh entity
     */
    gme_t GeoModelEditor::create_geological_entity( const std::string& type )
    {
        index_t index = model_.geological_entity_type( type ) ;
        if( index == NO_ID ) {
            index = create_geological_entity_type( type ) ;
        }
        index_t id = static_cast< index_t >( model_.geological_entities_[index].size() ) ;
        GeoModelGeologicalEntity* E = new_geological_entity( type, id ) ;
        model_.geological_entities_[index].push_back( E ) ;
        return E->gme_id() ;
    }


    index_t GeoModelEditor::create_geological_entity_type( const std::string& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory::has_creator( type ) ) ;
        model_.geological_entity_types_.push_back( type ) ;
        GeoModelGeologicalEntity* E = GeoModelGeologicalEntityFactory::create_object(
                type, model() ) ;

        const std::string child_type = E->child_type_name() ;

        EntityRelationships& parentage = model().entity_relationships_ ; 
        parentage.register_relationship( type, child_type ) ; 
        
        return static_cast< index_t >( model_.geological_entity_types_.size() - 1 ) ;
    }

    index_t GeoModelEditor::create_mesh_entities( const std::string& type, index_t nb )
    {
        assert_entity_creation_allowed() ;
        std::vector< GeoModelMeshEntity* >& store = model_.modifiable_mesh_entities(
            type ) ;
        index_t old_size = static_cast< index_t >( store.size() ) ;
        index_t new_size = old_size + nb ;
        store.resize( new_size, nil ) ;
        if( type == Corner::type_name_ ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Corner( model(), i ) ;
            }
        } else if( type == Line::type_name_ ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Line( model(), i ) ;
            }
        } else if( type == Surface::type_name_ ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Surface( model(), i ) ;
            }
        } else if( type == Region::type_name_ ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Region( model(), i ) ;
            }
        } else {
            ringmesh_assert_not_reached ;
            return NO_ID ;
        }

        return old_size ;
    }


    index_t GeoModelEditor::create_geological_entities( const std::string& type, index_t nb )
    {
        assert_entity_creation_allowed() ;
        std::vector< GeoModelGeologicalEntity* >& store =
            model_.modifiable_geological_entities( type ) ;
        index_t old_size = static_cast< index_t >( store.size() ) ;
        index_t new_size = old_size + nb ;
        store.resize( new_size, nil ) ;

        for( index_t i = old_size; i < new_size; i++ ) {
            ringmesh_assert( store[i] == nil ) ;
            store[i] = new_geological_entity( type, i ) ;
        }
        return old_size ;
    }

    /*! @details For all 7 types of entities, check what information is available
     * for the first one and fill the entities of the same type accordingly
     * THIS MEANS that the all the entities of the same type have been initialized with
     * the same information.
     */
    void GeoModelEditor::complete_entity_connectivity()
    {
        // Order is important
        complete_mesh_entity_connectivity< Line >() ;
        complete_mesh_entity_connectivity< Corner >() ;
        complete_mesh_entity_connectivity< Surface >() ;
        complete_mesh_entity_connectivity< Region >() ;

        // Geological entities
        for( index_t i = 0; i < model_.nb_geological_entity_type(); i++ ) {
            const std::string& type = model().geological_entity_type( i ) ;
            if( model_.nb_geological_entities( type ) > 0 ) {
                if( model_.geological_entity( type, 0 ).nb_children() == 0 ) {
                    fill_geological_entities_children( type ) ;
                }
            }
        }
    }

    template< typename E >
    void GeoModelEditor::complete_mesh_entity_connectivity()
    {
        const std::string& type = E::type_name_ ;
        if( model_.nb_mesh_entities( type ) > 0 ) {
            if( model_.mesh_entity( type, 0 ).nb_boundaries() == 0 ) {
                fill_mesh_entities_boundaries( type ) ;
            }
            if( model_.mesh_entity( type, 0 ).nb_in_boundary() == 0 ) {
                fill_mesh_entities_in_boundaries( type ) ;
            }
            if( model_.mesh_entity( type, 0 ).nb_parents() == 0 ) {
                fill_mesh_entities_parent( type ) ;
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_boundaries( const std::string& type )
    {
        // We have a problem if this is called for regions
        // No way yet to know the surface orientation
        ringmesh_assert( type != Region::type_name_ ) ;

        if( model().nb_mesh_entities( type ) == 0 ) return ;

        const std::string& b_type = mesh_entity( type, 0 ).boundary_type() ;
        if( b_type != GME::type_name_ ) {
            for( index_t i = 0; i < model().nb_mesh_entities( b_type ); ++i ) {
                gme_t cur_gme( b_type, i ) ;
                const GeoModelMeshEntity& b = mesh_entity( cur_gme ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_mesh_entity_boundary( b.in_boundary_gme( j ), cur_gme ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_in_boundaries( const std::string& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) return ;

        const std::string& in_b_type = mesh_entity( type, 0 ).boundary_type() ;
        if( in_b_type != GME::type_name_ ) {
            for( index_t i = 0; i < model().nb_mesh_entities( in_b_type ); ++i ) {
                gme_t cur_gme( in_b_type, i ) ;
                const GeoModelMeshEntity& in_b = mesh_entity( cur_gme ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_mesh_entity_in_boundary( in_b.boundary_gme( j ), cur_gme ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_parent( const std::string& type )
    {
        const GeoModel& M = model() ;
        if( M.nb_mesh_entities( type ) == 0 ) return ;

        const GeoModelMeshEntity& first_mesh_entity = mesh_entity( type, 0 ) ;

        const std::set< std::string >& parent_types( entity_relationships().parent_types( type ) ) ;

        for( std::set< std::string >::const_iterator it = parent_types.begin(); it != parent_types.end(); ++it ) {
            const std::string& parent_type = *it ;
            if( parent_type != GME::type_name_ ) {
                for( index_t i = 0; i < M.nb_geological_entities( parent_type ); ++i ) {
                    const GeoModelGeologicalEntity& parent = geological_entity(
                        parent_type, i ) ;
                    for( index_t j = 0; j < parent.nb_children(); ++j ) {
                        add_mesh_entity_parent( parent.child_id( j ), parent.gme_id() ) ;
                    }
                }
            }
        }
    }

    void GeoModelEditor::fill_geological_entities_children( const std::string& type )
    {
        if( model().nb_geological_entities( type ) == 0 ) return ;        

        const std::string& c_type = geological_entity( type, 0 ).child_type_name() ;
        if( c_type != GME::type_name_ ) {
            for( index_t i = 0; i < model().nb_mesh_entities( c_type ); ++i ) {
                gme_t cur_gme = gme_t( c_type, i ) ;
                const GeoModelMeshEntity& p = mesh_entity( cur_gme ) ;
                for( index_t j = 0; j < p.nb_parents(); j++ ) {
                    add_geological_entity_child( p.parent_id( j ), cur_gme ) ;
                }
            }
        }
    }

    void GeoModelEditor::complete_mesh_entities_geol_feature_from_first_parent(
        const std::string& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) return ;

        const std::string& parent_type = *entity_relationships().parent_types( type ).begin() ;   // beurk
        if( parent_type != GME::type_name_ ) {
            for( index_t i = 0; i < model().nb_mesh_entities( type ); ++i ) {
                GeoModelMeshEntity& E = mesh_entity( type, i ) ;
                if( !E.has_geological_feature() ) {
                    if( E.nb_parents() > 0 && E.parent( 0 ).has_geological_feature() ) {
                        E.geol_feature_ = E.parent( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }

    void GeoModelEditor::complete_geological_entities_geol_feature_from_first_child(
        const std::string& type )
    {
        if( model().nb_geological_entities( type ) == 0 ) {
            return ;
        }
        const std::string& child_type = entity_relationships().child_type( type ) ;
        // @todo change failuer return value for entity relationship requests
        if( child_type != GME::type_name_ ) {
            for( index_t i = 0; i < model().nb_geological_entities( type ); ++i ) {
                GeoModelGeologicalEntity& p = geological_entity( child_type, i ) ;
                if( !p.has_geological_feature() ) {
                    if( p.nb_children() > 0 && p.child( 0 ).has_geological_feature() ) {
                        p.geol_feature_ = p.child( 0 ).geological_feature() ;
                    }
                }
            }
        }
    }


    /*!
     * @brief Add to the vector the entities which cannot exist if
     *        an entity in the set does not exist.
     * @details These entities are added to the set.
     *          Recursive call till nothing is added.
     *
     * @return True if at least one entity was added, otherwise false.
     */
    bool GeoModelEditor::get_dependent_entities( std::set< gme_t >& in ) const
    {
//        index_t input_size = static_cast< index_t >( in.size() ) ;
//
//        for( std::set< gme_t >::iterator it( in.begin() ); it != in.end(); ++it ) {
//            gme_t cur = *it ;
//            /// If an entity has children entities - add them
//            if( GME::child_allowed( cur.type ) ) {
//                const GME& E = model_.entity( cur ) ;
//                for( index_t j = 0; j < E.nb_children(); ++j ) {
//                    in.insert( E.child_id( j ) ) ;
//                }
//            }
//        }
//
//        /// If a parent has no children anymore - add it
//        for( index_t p = GME::CONTACT; p < GME::NO_TYPE; ++p ) {
//            const std::string& P = (GME::TYPE) p ;
//            for( index_t j = 0; j < model_.nb_entities( P ); ++j ) {
//                bool no_child = true ;
//                const GME& E = model_.entity( GME::gme_t( P, j ) ) ;
//                for( index_t k = 0; k < E.nb_children(); ++k ) {
//                    if( in.count( E.child_id( k ) ) == 0 ) {
//                        no_child = false ;
//                        break ;
//                    }
//                }
//                if( no_child ) {
//                    in.insert( E.gme_id() ) ;
//                }
//            }
//        }
//
//        /// If an entity is in the boundary of nothing - add it
//        for( index_t t = GME::CORNER; t < GME::REGION; ++t ) {
//            const std::string& T = (GME::TYPE) t ;
//            for( index_t j = 0; j < model_.nb_entities( T ); ++j ) {
//                bool no_incident = true ;
//                const GME& E = model_.entity( GME::gme_t( T, j ) ) ;
//                for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
//                    if( in.count( E.in_boundary_gme( k ) ) == 0 ) {
//                        no_incident = false ;
//                        break ;
//                    }
//                }
//                if( no_incident ) {
//                    in.insert( E.gme_id() ) ;
//                }
//            }
//        }
//
//        if( in.size() != input_size ) {
//            return get_dependent_entities( in ) ;
//        } else {
//            return false ;
//        }
        return false ;
    }

    /*!
     * @brief Remove a list of entities of the model
     * @details No check is done on the consistency of this removal
     *          The entities and all references to them are removed.
     *          All dependent entities should be in the set of entities to remove,
     *          with a prior call to get_dependent_entities function.
     *
     * @warning NOT TESTED.
     *          The client is responsible to set the proper connectivity
     *          information between the remaining model entities.
     *
     * @todo TEST IT
     */
    void GeoModelEditor::remove_entities( const std::set< gme_t >& entities )
    {
//        if( entities.empty() ) {
//            return ;
//        }
//
//        // We need to remove entities type by type since they are
//        // stored in different vectors and since we use indices in these
//        // vectors to identify them.
//        // Initialize the vector
//        std::vector< std::vector< index_t > > to_erase_by_type ;
//        to_erase_by_type.reserve( GME::NO_TYPE ) ;
//        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
//            to_erase_by_type.push_back(
//                std::vector< index_t >(
//                    model_.nb_entities( static_cast< const std::string& >( i ) ), 0 ) ) ;
//        }
//        // Flag the entities to erase
//        for( std::set< gme_t >::const_iterator it = entities.begin();
//            it != entities.end(); ++it ) {
//            gme_t cur = *it ;
//            if( cur.type < GME::NO_TYPE ) {
//                ringmesh_assert( NO_ID != 0 ) ; // If one day NO_ID changes of value.
//                to_erase_by_type[cur.type][cur.index] = NO_ID ;
//            }
//        }
//
//        delete_elements( to_erase_by_type ) ;
    }

    /*!
     * @brief ONGOING WORK. Removes properly some entities of the Boundary Model.
     *
     * @param[in] entities_to_remove: in input the entities the client wants to
     * remove, in output all the removed entities (dependencies of ).
     *
     * Calls get_dependent_entities on each entities of \p entities_to_remove.
     * Then do remove these entities and updates the universe.
     *
     * @pre Assert that the entity to remove is one region that has only one neighbor
     *
     * @todo Finish to implement it for any kind of GME. BC must continue this work.
     *
     * @todo Review : Error in the comments [JP]
     */
    void GeoModelEditor::remove_entities_and_dependencies(
        const std::set< GME::gme_t >& entities_to_remove )
    {
//        // Asserts to remove when implementation is completed
//        ringmesh_assert(
//            entities_to_remove.size() == 1
//                && entities_to_remove.begin()->type == GME::REGION ) ;
//
//        // Copy because it is not logical to have in output the removed entities. BC
//        /// @todo Review : youpiii what did the comment on the function just said ? [JP]
//        std::set< GME::gme_t > entities = entities_to_remove ;
//        // TODO Handle the case of several objects in entities
//
//        const GeoModelEntity& reg = model_.entity( *( entities.begin() ) ) ;
//        get_dependent_entities( entities ) ;
//
//        // TODO Dirty duplication of code------------------------
//        // We need to remove entities type by type since they are
//        // stored in different vectors and since we use indices in these
//        // vectors to identify them.
//        // Initialize the vector
//        std::vector< std::vector< index_t > > to_erase_by_type ;
//        to_erase_by_type.reserve( GME::NO_TYPE ) ;
//        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
//            to_erase_by_type.push_back(
//                std::vector< index_t >(
//                    model_.nb_entities( static_cast< const std::string& >( i ) ), 0 ) ) ;
//        }
//        // Flag the entities to erase
//        for( std::set< gme_t >::const_iterator it = entities.begin();
//            it != entities.end(); ++it ) {
//            gme_t cur = *it ;
//            if( cur.type < GME::NO_TYPE ) {
//                ringmesh_assert( NO_ID != 0 ) ; // If one day NO_ID changes of value.
//                to_erase_by_type[cur.type][cur.index] = NO_ID ;
//            }
//        }
//
//        // Number of entities deleted for each TYPE
//        std::vector< index_t > nb_removed( to_erase_by_type.size(), 0 ) ;
//
//        /// 1. Get the mapping between old indices of the entities
//        ///    and new ones (when entities to remove will actually be removed)
//        for( index_t i = 0; i < to_erase_by_type.size(); ++i ) {
//            for( index_t j = 0; j < to_erase_by_type[i].size(); ++j ) {
//                if( to_erase_by_type[i][j] == NO_ID ) {
//                    nb_removed[i]++ ;
//                } else {
//                    to_erase_by_type[i][j] = j - nb_removed[i] ;
//                }
//            }
//        }
//        // TODO Dirty duplication of code--------------------------
//
//        std::vector< GME::gme_t > to_add_in_universe ;
//
//        if( reg.type() == GME::REGION ) {
//            index_t nb_added = 0 ;
//            for( index_t b_i = 0; b_i < reg.nb_boundaries(); ++b_i ) {
//                if( !reg.boundary( b_i ).is_on_voi() ) {
//                    to_add_in_universe.push_back( reg.boundary( b_i ).gme_id() ) ;
//                    ringmesh_assert(
//                        to_erase_by_type[reg.boundary( b_i ).type()][reg.boundary(
//                            b_i ).index()] != NO_ID ) ;
//                    to_add_in_universe[nb_added].index =
//                        to_erase_by_type[reg.boundary( b_i ).type()][reg.boundary(
//                            b_i ).index()] ;
//                    ++nb_added ;
//                }
//            }
//        }
//
//        remove_entities( entities ) ;
//
//        // Update Universe
//        /// @todo You first need to clean the existing universe [JP]
//        for( std::vector< GME::gme_t >::const_iterator itr =
//            to_add_in_universe.begin(); itr != to_add_in_universe.end(); ++itr ) {
//            add_entity_boundary( gme_t( GME::REGION, NO_ID ),
//                gme_t( GME::SURFACE, itr->index ), true ) ;
//        }
//
//        //ringmesh_assert( model_.check_model_validity() ) ;
    }

    /*!
     * @brief Delete entities and remove all references to them in the model
     *
     * @param[in,out] to_erase For each type of entity T,
     *        store a vector of the size of model_.nb_entities(T) in which
     *        entities are flagged with NO_ID.
     *        In output it stores the mapping table between old and new indices
     *        for the entities.
     * @todo TEST IT
     */
    void GeoModelEditor::delete_elements(
        std::vector< std::vector< index_t > >& to_erase )
    {
        // Number of entities deleted for each TYPE
//        std::vector< index_t > nb_removed( to_erase.size(), 0 ) ;
//
//        /// 1. Get the mapping between old indices of the entities
//        ///    and new ones (when entities to remove will actually be removed)
//        for( index_t i = 0; i < to_erase.size(); ++i ) {
//            for( index_t j = 0; j < to_erase[i].size(); ++j ) {
//                if( to_erase[i][j] == NO_ID ) {
//                    nb_removed[i]++ ;
//                } else {
//                    to_erase[i][j] = j - nb_removed[i] ;
//                }
//            }
//        }
//
//        /// 2. Effectively delete the entities
//        for( index_t i = 0; i < to_erase.size(); ++i ) {
//            for( index_t j = 0; j < to_erase[i].size(); ++j ) {
//                if( to_erase[i][j] == NO_ID ) {
//                    GME::gme_t cur( static_cast< const std::string& >( i ), j ) ;
//                    delete model_.entity_ptr( cur ) ;
//                    // Set the current entity to nil
//                    model_.modifiable_entities( cur.type )[cur.index] = nil ;
//                }
//            }
//            std::vector< GME* >& store = model_.modifiable_entities(
//                static_cast< const std::string& >( i ) ) ;
//            store.erase(
//                std::remove( store.begin(), store.end(),
//                    static_cast< GME* >( nil ) ), store.end() ) ;
//        }
//
//        /// 3. Deal with the model vertices
//        model_.mesh.vertices.clear() ;
//
//        /// 4. Update all possible indices in remaining entities
//        for( index_t i = 0; i < to_erase.size(); ++i ) {
//            const std::string& T = static_cast< GME::TYPE >( i ) ;
//
//            // Update all indices stored by the GME of that type
//            ringmesh_assert(
//                model_.nb_entities( T ) == to_erase[i].size() - nb_removed[i] ) ;
//            for( index_t j = 0; j < model_.nb_entities( T ); ++j ) {
//                gme_t entity_id( T, j ) ;
//                GeoModelEntity& E = model_.modifiable_entity( entity_id ) ;
//
//                // Not the same than j - since we have erased some entities
//                index_t old_id = E.index() ;
//                ringmesh_assert( to_erase[i][old_id] != NO_ID ) ;
//
//                // id_
//                E.id_.index = to_erase[i][old_id] ;
//                // boundary_
//                if( E.nb_boundaries() > 0 ) {
//                    const std::string& B = GME::boundary_type( T ) ;
//                    ringmesh_assert( B < GME::NO_TYPE ) ;
//                    for( index_t k = 0; k < E.nb_boundaries(); ++k ) {
//                        set_entity_boundary( E.gme_id(), k,
//                            gme_t( B, to_erase[B][E.boundary_gme( k ).index] ) ) ;
//                    }
//                }
//                // in_boundary
//                if( E.nb_in_boundary() > 0 ) {
//                    const std::string& IB = GME::in_boundary_type( T ) ;
//                    ringmesh_assert( IB < GME::NO_TYPE ) ;
//                    for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
//                        set_entity_in_boundary( E.gme_id(), k,
//                            gme_t( IB,
//                                to_erase[IB][E.in_boundary_gme( k ).index] ) ) ;
//                    }
//                }
//                // parent_
//                if( E.has_parent() ) {
//                    const std::string& P = GME::parent_type( T ) ;
//                    ringmesh_assert( P < GME::NO_TYPE ) ;
//                    set_entity_parent( E.gme_id(),
//                        gme_t( P, to_erase[P][E.parent_id().index] ) ) ;
//                }
//                // children_
//                if( E.nb_children() > 0 ) {
//                    const std::string& C = GME::child_type( T ) ;
//                    ringmesh_assert( C < GME::NO_TYPE ) ;
//                    for( index_t k = 0; k < E.nb_children(); ++k ) {
//                        set_entity_child( E.gme_id(), k,
//                            gme_t( C, to_erase[C][E.child_id( k ).index] ) ) ;
//                    }
//                }
//                // Clean the vectors in the entity
//                erase_invalid_entity_references( entity_id ) ;
//            }
//        }
//
//        // Do not forget the universe
//        /// @todo Put the universe in the list of regions - so annoying to think of it each time
//        /// @todo BUG ? sides are lost for the universe ? not sure ... [JP]
//        {
//            GME::gme_t univers_id( GME::REGION, NO_ID ) ;
//            Region& U = dynamic_cast< Region& >( model_.modifiable_entity(
//                univers_id ) ) ;
//            for( index_t i = 0; i < U.nb_boundaries(); ++i ) {
//                set_entity_boundary( U.gme_id(), i,
//                    gme_t( GME::SURFACE,
//                        to_erase[GME::SURFACE][U.boundary_gme( i ).index] ) ) ;
//            }
//            erase_invalid_entity_references( univers_id ) ;
//        }
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model entities and their relationship ignoring their geometry
     *
     * @param[in] from Model to copy the information from
     */
    void GeoModelEditor::copy_macro_topology( const GeoModel& from )
    {
        assert_entity_creation_allowed() ;
        copy_mesh_entity_topology< Corner >( from ) ;
        copy_mesh_entity_topology< Line >( from ) ;
        copy_mesh_entity_topology< Surface >( from ) ;
        copy_mesh_entity_topology< Region >( from ) ;

        model().universe_ = from.universe_ ;
    }

    template< typename ENTITY >
    void GeoModelEditor::copy_mesh_entity_topology( const GeoModel& from )
    {
        const std::string& type = ENTITY::type_name_ ;
        std::vector< GeoModelMeshEntity* >& store = model_.modifiable_mesh_entities( type ) ;
        store.resize( from.nb_mesh_entities( type ), nil ) ;

        for( index_t e = 0; e < model_.nb_mesh_entities( type ); ++e ) {
            store[e] = new ENTITY( model(), e ) ;
            ringmesh_assert( store[ e ] != nil ) ;
        }
        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < model_.nb_mesh_entities( type ); ++e ) {
            GME::gme_t id( type, e ) ;
            GeoModelEntity& lhs = mesh_entity( id ) ;
            const GeoModelEntity& rhs = from.mesh_entity( id ) ;
            lhs = rhs ;
        }
    }

    /*!
     * @brief Remove invalid reference to entities
     *       boundary, in_boundary and children vectors
     *       Invalid entities have a NO_ID index
     *
     * @param[in] E the entity gme_t to process
     */
    void GeoModelEditor::erase_invalid_entity_references( const GME::gme_t& E_id )
    {
//        GME& E = entity( E_id ) ;
//        const std::string& T = E.type() ;
//        if( E.nb_children() > 0 ) {
//            gme_t invalid_child( E.child_type( T ), NO_ID ) ;
//            index_t nb_found = static_cast< index_t >( std::count(
//                E.children_.begin(), E.children_.end(), invalid_child ) ) ;
//            if( nb_found == E.children_.size() ) {
//                // Calling erase on all entities -> undefined behavior
//                E.children_.clear() ;
//            } else {
//                E.children_.erase(
//                    std::remove( E.children_.begin(), E.children_.end(),
//                        invalid_child ), E.children_.end() ) ;
//            }
//        }
//        if( E.nb_boundaries() > 0 ) {
//            gme_t invalid_boundary( E.boundary_type( T ), NO_ID ) ;
//            if( E.type() == GME::REGION ) {
//                Region& R = dynamic_cast< Region& >( E ) ;
//                // Change side values if necessary
//                index_t offset = 0 ;
//                for( index_t i = 0; i + offset < E.nb_boundaries(); ++i ) {
//                    if( E.boundaries_[i] == invalid_boundary ) {
//                        offset++ ;
//                    } else {
//                        R.sides_[i] = R.side( i + offset ) ;
//                    }
//                }
//            }
//
//            index_t end = static_cast< index_t >( std::remove( E.boundaries_.begin(),
//                E.boundaries_.end(), invalid_boundary ) - E.boundaries_.begin() ) ;
//            if( end == 0 ) {
//                E.boundaries_.clear() ;
//                if( E.type() == GME::REGION ) {
//                    Region& R = dynamic_cast< Region& >( E ) ;
//                    R.sides_.clear() ;
//                }
//            } else {
//                E.boundaries_.erase( E.boundaries_.begin() + end,
//                    E.boundaries_.end() ) ;
//                if( E.type() == GME::REGION ) {
//                    Region& R = dynamic_cast< Region& >( E ) ;
//                    R.sides_.erase( R.sides_.begin() + end, R.sides_.end() ) ;
//                }
//            }
//        }
//        if( E.nb_in_boundary() > 0 ) {
//            gme_t invalid_in_boundary( E.in_boundary_type( T ), NO_ID ) ;
//            index_t nb_found = static_cast< index_t >( std::count(
//                E.in_boundary_.begin(), E.in_boundary_.end(), invalid_in_boundary ) ) ;
//            if( nb_found == E.in_boundary_.size() ) {
//                E.in_boundary_.clear() ;
//            } else {
//                E.in_boundary_.erase(
//                    std::remove( E.in_boundary_.begin(), E.in_boundary_.end(),
//                        invalid_in_boundary ), E.in_boundary_.end() ) ;
//            }
//        }
    }
}
