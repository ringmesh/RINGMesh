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
#include <vector> 
#include <map>
#include <set>

#include <ringmesh/geo_model.h>

/*!
 * @file Implementation of the GeoModelEditor
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    typedef GeoModelEntity::gme_t gme_t ;

    GeoModelEditor::GeoModelEditor( GeoModel& model )
        : model_( model ), create_entity_allowed_( true )
    {
    }
    GeoModelEditor::~GeoModelEditor()
    {
    }

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
     *
     gme_t GeoModelEditor::create_mesh_entity( const std::string& type )
     {
     assert_entity_creation_allowed() ;
     if( type == Corner::type_name_static() ) {
     index_t id = static_cast< index_t >( model_.corners_.size() ) ;
     Corner* corner = new Corner( model(), id ) ;
     model_.corners_.push_back( corner ) ;
     return corner->gme_id() ;
     } else if( type == Line::type_name_static() ) {
     index_t id = static_cast< index_t >( model_.lines_.size() ) ;
     Line* line = new Line( model(), id ) ;
     model_.lines_.push_back( line ) ;
     return line->gme_id() ;
     } else if( type == Surface::type_name_static() ) {
     index_t id = static_cast< index_t >( model_.surfaces_.size() ) ;
     Surface* surface = new Surface( model(), id ) ;
     model_.surfaces_.push_back( surface ) ;
     return surface->gme_id() ;
     } else if( type == Region::type_name_static() ) {
     index_t id = static_cast< index_t >( model_.regions_.size() ) ;
     Region* region = new Region( model(), id ) ;
     model_.regions_.push_back( region ) ;
     return region->gme_id() ;
     } else {
     ringmesh_assert_not_reached ;
     return gme_t() ;
     }
     }*/

    /*!
     * @brief Creates a mesh entity of the given type and add it to the correct vector
     * The GeoModelMeshEntity is created from its type and its index
     *
     * @param[in] type Type of the mesh entity to create
     * @return The index of the created mesh entity
     */
    gme_t GeoModelEditor::create_geological_entity( const std::string& type )
    {
        bool type_exists = model_.is_geological_entity_type( type ) ;
        if( !type_exists ) {
            create_geological_entity_type( type ) ;
        }
        index_t index = model().geological_entity_type( type ) ;
        index_t id = static_cast<index_t>(model_.geological_entities_[index].size()) ;
        GeoModelGeologicalEntity* E = new_geological_entity( type, id ) ;
        model_.geological_entities_[index].push_back( E ) ;
        return E->gme_id() ;
    }


    index_t GeoModelEditor::create_geological_entity_type( const std::string& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory::has_creator( type ) ) ;
        model_.geological_entity_types_.push_back( type ) ;
        model_.geological_entities_.push_back( std::vector< GeoModelGeologicalEntity* >() ) ;
        GeoModelGeologicalEntity* E = GeoModelGeologicalEntityFactory::create_object(
            type, model() ) ;

        const std::string child_type = E->child_type_name() ;

        EntityRelationships& parentage = model().entity_relationships_ ;
        parentage.register_relationship( type, child_type ) ;

        return static_cast<index_t>(model_.geological_entity_types_.size() - 1) ;
    }

    index_t GeoModelEditor::create_mesh_entities( const std::string& type, index_t nb )
    {
        assert_entity_creation_allowed() ;
        std::vector< GeoModelMeshEntity* >& store = modifiable_mesh_entities(
            type ) ;
        index_t old_size = static_cast<index_t>(store.size()) ;
        index_t new_size = old_size + nb ;
        store.resize( new_size, nil ) ;
        if( type == Corner::type_name_static() ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Corner( model(), i ) ;
            }
        } else if( type == Line::type_name_static() ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Line( model(), i ) ;
            }
        } else if( type == Surface::type_name_static() ) {
            for( index_t i = old_size; i < new_size; i++ ) {
                ringmesh_assert( store[i] == nil ) ;
                store[i] = new Surface( model(), i ) ;
            }
        } else if( type == Region::type_name_static() ) {
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
        index_t index = model_.geological_entity_type( type ) ;
        if( index == NO_ID ) {
            index = create_geological_entity_type( type ) ;
        }
        std::vector< GeoModelGeologicalEntity* >& store =
            modifiable_geological_entities( type ) ;
        index_t old_size = static_cast<index_t>(store.size()) ;
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
        for( index_t i = 0; i < model_.nb_geological_entity_types(); i++ ) {
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
        const std::string& type = E::type_name_static() ;
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
        if( model().nb_mesh_entities( type ) == 0 ) {
            return ;
        }

        const std::string& b_type = model().entity_relationships().boundary_type( type ) ;
        if( b_type != GME::type_name_static() ) {
            for( index_t i = 0; i < model().nb_mesh_entities( b_type ); ++i ) {
                const GeoModelMeshEntity& b = mesh_entity( b_type, i ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_mesh_entity_boundary( b.in_boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_in_boundaries( const std::string& type )
    {
        if( model().nb_mesh_entities( type ) == 0 ) {
            return ;
        }

        const std::string& in_b_type = model().entity_relationships().in_boundary_type( type ) ;
        if( in_b_type != GME::type_name_static() ) {
            for( index_t i = 0; i < model().nb_mesh_entities( in_b_type ); ++i ) {
                const GeoModelMeshEntity& in_b = mesh_entity( in_b_type, i ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_mesh_entity_in_boundary( in_b.boundary_gme( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_mesh_entities_parent( const std::string& type )
    {
        const GeoModel& M = model() ;
        if( M.nb_mesh_entities( type ) == 0
            || entity_relationships().nb_parent_types( type ) == 0 ) {
            return ;
        }

        const std::set< std::string >& parent_types( entity_relationships().parent_types( type ) ) ;

        for( std::set< std::string >::const_iterator it = parent_types.begin(); it != parent_types.end(); ++it ) {
            const std::string& parent_type = *it ;
            if( parent_type != GME::type_name_static() ) {
                for( index_t i = 0; i < M.nb_geological_entities( parent_type ); ++i ) {
                    const GeoModelGeologicalEntity& parent = geological_entity(
                        parent_type, i ) ;
                    for( index_t j = 0; j < parent.nb_children(); ++j ) {
                        add_mesh_entity_parent( parent.child_gme( j ), parent.gme_id() ) ;
                    }
                }
            }
        }
    }

    void GeoModelEditor::fill_geological_entities_children( const std::string& type )
    {
        if( model().nb_geological_entities( type ) == 0 ) {
            return ;
        }

        const std::string& c_type = geological_entity( type, 0 ).child_type_name() ;
        if( c_type != GME::type_name_static() ) {
            for( index_t i = 0; i < model().nb_mesh_entities( c_type ); ++i ) {
                const GeoModelMeshEntity& p = mesh_entity( c_type, i ) ;
                for( index_t j = 0; j < p.nb_parents(); j++ ) {
                    add_geological_entity_child( p.parent_id( j ), i ) ;
                }
            }
        }
    }

    void GeoModelEditor::complete_mesh_entities_geol_feature_from_first_parent(
        const std::string& type )
    {
        if( model().nb_mesh_entities( type ) == 0
            || entity_relationships().nb_parent_types( type ) == 0 ) {
            return ;
        }

        const std::string& parent_type = *entity_relationships().parent_types( type ).begin() ;   // beurk
        if( parent_type != GME::type_name_static() ) {
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
        if( child_type != GME::type_name_static() ) {
            for( index_t i = 0; i < model().nb_geological_entities( type ); ++i ) {
                GeoModelGeologicalEntity& p = geological_entity( type, i ) ;
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
        index_t input_size = static_cast<index_t>(in.size()) ;

        for( std::set< gme_t >::iterator it( in.begin() ); it != in.end(); ++it ) {
            gme_t cur = *it ;
            if( model().is_geological_entity_type( cur.type ) ) {
                const GeoModelGeologicalEntity& E = model_.geological_entity( cur ) ;
                for( index_t j = 0; j < E.nb_children(); ++j ) {
                    in.insert( E.child_gme( j ) ) ;
                }
            }
        }

        index_t nb_geological_entity_types = model().nb_geological_entity_types() ;
        for( index_t i = 0; i < nb_geological_entity_types; ++i ) {
            const std::string& type = model().geological_entity_type( i ) ;

            for( index_t j = 0; j < model_.nb_geological_entities( type ); ++j ) {
                bool no_child = true ;
                const GeoModelGeologicalEntity& E = model_.geological_entity( type, j ) ;
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

        std::string mesh_entity_types[4] = {Corner::type_name_static(), Line::type_name_static(),
            Surface::type_name_static(), Region::type_name_static()} ;

        for( index_t i = 0; i < 4 ; ++i ) {
            const std::string& type = mesh_entity_types[i] ;
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

        if( in.size() != input_size ) {
            return get_dependent_entities( in ) ;
        } else {
            return false ;
        }
    }

    /*!
     * @brief Class in charge of removing entities from a GeoModel
     */
    class GeoModelEntityRemoval : public GeoModelEditor {        
    public:
        typedef std::map< std::string, index_t > TypeToIndex ;
        typedef std::map< index_t, std::string > IndexToType ;
        typedef std::string EntityType ;

        GeoModelEntityRemoval( GeoModel& model ) :
            GeoModelEditor( model )
        { 
            nb_geological_entity_types_ = GeoModelEditor::model().nb_geological_entity_types() ;
            nb_mesh_entity_types_ = GeoModel::nb_mesh_entity_types() ;
            nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_ ;

            // Allocate storage 
            nb_removed_entities_.resize( nb_entity_types_, 0 ) ;
            
            fill_entity_type_to_index_map() ;           

            fill_nb_initial_entities() ;
            initialize_costly_storage() ;
        }

        void remove_entities( const std::set< gme_t >& entities )
        {
            initialize_for_removal(entities) ;
            do_delete_flagged_entities() ;

            clear_model_mesh_vertices() ;

            update_entity_connectivity() ;
        }

    private:
        void initialize_costly_storage()
        {
            to_erase_.resize( nb_entity_types_ ) ;
            old_2_new_entity_.resize( nb_entity_types_ ) ;
            for( index_t i = 0; i < nb_entity_types_; ++i ) {
                index_t size = model().nb_entities( index_to_entity_type_[i] ) ;
                to_erase_[i].resize( size, 0 ) ;
                old_2_new_entity_[i].resize( size, 0 ) ;
            }
        }
        void initialize_for_removal( const std::set< gme_t >& entities_to_remove )
        {
            fill_to_erase_vectors( entities_to_remove ) ;
            fill_removed_entities_and_mapping() ;
        }

        void do_delete_flagged_entities()
        {
            for( index_t i = 0; i < nb_entity_types_; ++i ) {
                const EntityType& type = index_to_entity_type_[i] ;
                for( index_t j = 0; j < nb_initial_entities_[i]; ++j ) {
                    if( to_erase_[i][j] == NO_ID ) {
                        delete_entity( i, j ) ;
                    }
                }
                clear_null_entities( i ) ;
            }
        }

        void delete_entity( index_t type, index_t index )
        {
            const EntityType& type_name = index_to_entity_type_[type] ;
            gme_t id( type_name, index ) ;
            GeoModelEditor::delete_entity( type_name, index ) ;
        }
        void clear_null_entities( index_t type )
        {
            const EntityType& type_name = index_to_entity_type_[type] ;
            std::vector< GME* >& store = modifiable_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                static_cast<GME*>(nil) ), store.end() ) ;

            // QC
            ringmesh_assert( model().nb_entities( type_name )
                == nb_initial_entities_[type] - nb_removed_entities_[type] ) ;
        }

        /*!
         * Computes the mapping between old indices of the entities
         * and new ones, once some entities are actually be removed.
         *
         */
        void fill_removed_entities_and_mapping()
        {
            for( index_t i = 0; i < nb_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_entities_[i]; ++j ) {
                    if( to_erase_[i][j] == NO_ID ) {
                        nb_removed_entities_[i]++ ;
                        old_2_new_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_entity_[i][j] = j - nb_removed_entities_[i] ;
                    }
                }
            }
        }
        void fill_to_erase_vectors( const std::set< gme_t >& entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it = entities_to_remove.begin();
                it != entities_to_remove.end(); ++it ) {
                gme_t cur = *it ;

                index_t type_index = entity_type_to_index_[cur.type] ;
                to_erase_[type_index][cur.index] = NO_ID ;

            }
        }

        void fill_nb_initial_entities()
        {
            nb_initial_entities_.resize( nb_entity_types_, 0 ) ;
            for( index_t i = 0; i < nb_entity_types_; ++i ) {
                const EntityType& type = index_to_entity_type_[i] ;
                nb_initial_entities_[i] = model().nb_entities( type )  ;
            }
        }

        void fill_entity_type_to_index_map()
        {
            entity_type_to_index_[Corner::type_name_static()] = 0 ;
            entity_type_to_index_[Line::type_name_static()] = 1 ;
            entity_type_to_index_[Surface::type_name_static()] = 2 ;
            entity_type_to_index_[Region::type_name_static()] = 3 ;

            index_t counter = 4 ;
            for( index_t i = 0 ; i < model().nb_geological_entity_types(); ++i ) {
                entity_type_to_index_[model().geological_entity_type( i )] = counter ;
                counter++ ;
            }

            index_to_entity_type_[0] = Corner::type_name_static() ;
            index_to_entity_type_[1] = Line::type_name_static() ;
            index_to_entity_type_[2] = Surface::type_name_static() ;
            index_to_entity_type_[3] = Region::type_name_static() ;

            counter = 4 ;
            for( index_t i = 0 ; i < model().nb_geological_entity_types(); ++i ) {
                index_to_entity_type_[counter] = model().geological_entity_type( i ) ; ;
                counter++ ;
            }
        }

        void clear_model_mesh_vertices()
        {
            model().mesh.vertices.clear() ;
        }

        index_t entity_type_index( const GeoModelEntity& E ) const
        {
            const EntityType& type = E.type_name() ;            
            return entity_type_to_index_.find(type)->second ;
        }

        void update_entity_index( GeoModelEntity& E )
        {
            index_t old_id = E.index() ;
            index_t type = entity_type_index( E ) ;

            index_t new_id = old_2_new_entity_[type][old_id] ;

            ringmesh_assert( new_id != NO_ID ) ;
            set_entity_index( E.gme_id(), new_id ) ;
        }
        void update_entity_boundaries( GeoModelMeshEntity& E )
        {
            index_t type_index = boundary_type_index( E.entity_type() ) ;
            if( type_index == NO_ID ) {
                return ;
            }

            for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
                index_t old_boundary = E.boundary_gme( i ).index ;
                index_t new_boundary = old_2_new_entity_[type_index][old_boundary] ;
                set_mesh_entity_boundary( E.gme_id(), i, new_boundary ) ;
            }
        }

        void update_region_boundary_signs( Region& R )
        {
            const EntityType& surface_type = boundary_type( R.entity_type() ) ;
            gme_t invalid_value( surface_type, NO_ID ) ;

            index_t offset = 0 ;
            for( index_t i = 0; i + offset < R.nb_boundaries(); ++i ) {
                if( R.boundary_gme(i) == invalid_value ) {
                    offset++ ;
                } else {
                    bool new_side = R.side( i + offset ) ;
                    set_boundary_sign( R, i, new_side ) ;
                }
            }
        }


        void update_entity_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_boundary_type = EntityRelationships::in_boundary_type( E.entity_type() ) ;
            bool valid_type = EntityRelationships::is_valid_type( in_boundary_type ) ;
            if( !valid_type ) {
                return ;
            }
            index_t in_boundary_type_index = entity_type_to_index_[in_boundary_type] ;
            for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
                index_t old_id = E.in_boundary_gme( i ).index ;
                index_t new_id = old_2_new_entity_[in_boundary_type_index][old_id] ;
                set_mesh_entity_in_boundary( E.gme_id(), i, new_id ) ;
            }
        }

        void update_entity_parents( GeoModelMeshEntity& E )
        {
            const EntityRelationships& family_tree = model().entity_relationships() ;

            if( family_tree.nb_parent_types( E.type_name() ) == 0 ) {
                return ;
            } else {
                const std::set< EntityType >& parents = family_tree.parent_types( E.entity_type() ) ;

                for( std::set< EntityType >::const_iterator it( parents.begin() ); it != parents.end(); ++it ) {
                    const EntityType& parent_type = *it ;
                    index_t parent_type_index = entity_type_to_index_.find( parent_type )->second ;

                    index_t p_id = E.parent_id( parent_type ) ;
                    index_t old_id = E.parent( p_id ).index() ;
                    index_t new_id = old_2_new_entity_[parent_type_index][old_id] ;
                    set_mesh_entity_parent( E.gme_id(), p_id, gme_t( parent_type, new_id ) ) ;
                }
            }
        }

        index_t children_type_index( const EntityType& type ) const
        {
            const EntityType& child_type = children_type( type ) ;
            return entity_type_to_index_.find( child_type )->second ;
        }
        const EntityType& children_type( const EntityType& type ) const
        {
            const EntityRelationships& family = model().entity_relationships() ;
            return family.child_type( type ) ;
        }
        index_t boundary_type_index( const EntityType& type ) const
        {           
            const EntityType& b_type = boundary_type( type ) ;
            if( !EntityRelationships::is_valid_type( b_type ) ) {
                return NO_ID ;
            } else {
                return entity_type_to_index_.find(b_type)->second ;
            }
        }
        const EntityType& boundary_type( const EntityType& type ) const
        {
            const EntityRelationships& family = model().entity_relationships() ;
            return family.boundary_type( type ) ;
        }

        index_t in_boundary_type_index( const EntityType& type ) const
        {           
            const EntityType& in_b_type = in_boundary_type( type ) ;
            if( !EntityRelationships::is_valid_type( in_b_type ) ) {
                return NO_ID ;
            } else {
                return entity_type_to_index_.find(in_b_type)->second ;
            }
        }
        const EntityType& in_boundary_type( const EntityType& type ) const
        {
            const EntityRelationships& family = model().entity_relationships() ;
            return family.in_boundary_type( type ) ;
        }


        void update_entity_children( GeoModelGeologicalEntity& E )
        {
            if( E.nb_children() > 0 ) {
                index_t child_type = children_type_index( E.entity_type() ) ;
                for( index_t i = 0; i < E.nb_children(); ++i ) {
                    index_t old_id = E.child_gme( i ).index ;
                    index_t new_id = old_2_new_entity_[child_type][old_id] ;
                    set_geological_entity_child( E.gme_id(), i, new_id ) ;
                }
            }
        }

        void remove_invalid_values( std::vector< gme_t >& vector, const gme_t& invalid_value )
        {
            std::vector< gme_t >::iterator new_end = std::remove(
                vector.begin(), vector.end(), invalid_value ) ;
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

                remove_invalid_values( modifiable_children(E), invalid_child ) ;
            }
        }
        void delete_invalid_boundaries( GeoModelMeshEntity& E )
        {
            const EntityType& b_type = boundary_type( E.entity_type() ) ;
            gme_t invalid( b_type, NO_ID ) ;
            if( !EntityRelationships::is_valid_type( b_type ) ) {
                return ;
            } else {
                remove_invalid_values( modifiable_boundaries(E), invalid ) ;
            }
        }
        void delete_invalid_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_b_type = in_boundary_type( E.entity_type() ) ;
            gme_t invalid( in_b_type, NO_ID ) ;
            if( !EntityRelationships::is_valid_type( in_b_type ) ) {
                return ;
            } else {
                remove_invalid_values( modifiable_in_boundaries(E), invalid ) ;
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
                    gme_t new_id = E.parent( i + offset ).gme_id() ;
                    set_mesh_entity_parent( E.gme_id(), i, new_id ) ;
                }
                new_size = i +1 ; // je suis pas sure de la taille .. to check
            }
            modifiable_parents(E).resize( new_size ) ; 
        }
        // To call after delete invalid boundaries on the Region !!
        void delete_invalid_signs( Region& R )
        {
            modifiable_sides(R).resize( R.nb_boundaries() ) ;
        }

        bool is_mesh_entity( index_t i ) const
        {
            return i < nb_mesh_entity_types_ ;
        }
        bool is_geological_entity( index_t i ) const
        {
            return !is_mesh_entity( i ) ;
        }
        bool is_region_entity( index_t i ) const
        {
            return i == 3 ; // Magic number = bad 
        }

        void update_entity_connectivity()
        {
            for( index_t i = 0; i < nb_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_entity_type_[i] ;

                for( index_t j = 0; j < model().nb_entities( entity_type ); ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelEntity& E = modifiable_entity( new_id ) ;

                    update_entity_index( E ) ;

                    if( is_mesh_entity(i) ){
                        GeoModelMeshEntity& ME = dynamic_cast<GeoModelMeshEntity&>(E);

                        update_entity_boundaries( ME ) ;
                        delete_invalid_boundaries( ME ) ;

                        update_entity_in_boundary( ME ) ;
                        delete_invalid_in_boundary( ME ) ;

                        update_entity_parents( ME ) ;
                        delete_invalid_parents( ME ) ;
                    }
                    if( is_region_entity(i) ) {
                        Region& R = dynamic_cast<Region&>(E);
                        update_region_boundary_signs( R ) ;
                        delete_invalid_signs( R ) ;
                    }
                    if( is_geological_entity( i ) ) {
                        GeoModelGeologicalEntity& GE = dynamic_cast<GeoModelGeologicalEntity&>(E);                         
                        update_entity_children( GE ) ;
                        delete_invalid_children( GE ) ;
                    }
          
                }
                update_universe_connectivity() ;
            }
        }

        void update_universe_connectivity()
        {
            ringmesh_assert_not_reached ; // todo implement !!
        }



        private:
            index_t nb_entity_types_ ;
            index_t nb_geological_entity_types_ ;
            index_t nb_mesh_entity_types_ ;


            std::vector< index_t > nb_initial_entities_ ;
            std::vector< index_t > nb_removed_entities_ ;

            /*! For each type of entity, store a vector of where the
                entities to remove from the model are flagged with NO_ID.
                */
            std::vector< std::vector< index_t > > to_erase_ ;
            /*! Stores the mapping table between indices for each type of
             *  element before and after the removal of entities
             */
            std::vector< std::vector< index_t > > old_2_new_entity_ ;

            std::map< EntityType, index_t > entity_type_to_index_ ;
            std::map< index_t, EntityType > index_to_entity_type_ ;

        };
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
        */
    void GeoModelEditor::remove_entities( const std::set< gme_t >& entities )
    {
        if( entities.empty() ) {
            return ;
        } else {
            GeoModelEntityRemoval remover( model() ) ;
            remover.remove_entities( entities ) ;
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
        const std::string& type = ENTITY::type_name_static() ;
        std::vector< GeoModelMeshEntity* >& store = modifiable_mesh_entities( type ) ;
        store.resize( from.nb_mesh_entities( type ), nil ) ;

        for( index_t e = 0; e < model_.nb_mesh_entities( type ); ++e ) {
            store[e] = new ENTITY( model(), e ) ;
            ringmesh_assert( store[e] != nil ) ;
        }
        RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < model_.nb_mesh_entities( type ); ++e ) {
                GME::gme_t id( type, e ) ;
                GeoModelEntity& lhs = mesh_entity( id ) ;
                const GeoModelEntity& rhs = from.mesh_entity( id ) ;
                lhs = rhs ;
            }
    }
}
