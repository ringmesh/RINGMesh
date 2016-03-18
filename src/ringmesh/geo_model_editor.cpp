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

    typedef GeoModelElement::gme_t gme_t ;


    GeoModelElement& GeoModelEditor::element( const GME::gme_t& id ) const
    {
        return *model_.element_ptr( id ) ;
    }

    void GeoModelEditor::set_model_name( const std::string& name )
    {
        model_.name_ = name ;
    }

    /*!
     * @brief Creates a element of the given type and add it to the correct vector
     * The GeoModelElement is created from its type and its index
     *
     * @param[in] type Type of the element to create
     * @return The index of the created element
     */
    gme_t GeoModelEditor::create_element( GME::TYPE type )
    {   
        assert_element_creation_allowed() ;
        if( type >= GME::CORNER && type < GME::NO_TYPE ) {
            GME* E = new_element( type ) ;
            ringmesh_assert( E != nil ) ;

            model_.modifiable_elements( type ).push_back( E ) ;                
            return E->gme_id();
        } else {
            ringmesh_assert_not_reached;
            return gme_t() ;
        }
    }

    void GeoModelEditor::create_elements( GME::TYPE type, index_t nb )
    {
        assert_element_creation_allowed() ;
        if( type >= GME::NO_TYPE ) {
            return ;
        }
        std::vector< GME* >& store = model_.modifiable_elements( type ) ;
        index_t old_size = static_cast<index_t> ( store.size() ) ;
        index_t new_size = old_size + nb ;
        store.resize( new_size, nil ) ;
        for( index_t i = old_size; i < new_size; i++ ) {
            ringmesh_assert( store[i] == nil ) ;
            store[i] = new_element( type, i ) ;
        }
    }

    /*! @details For all 7 types of elements, check what information is available
     * for the first one and fill the elements of the same type accordingly
     * THIS MEANS that the all the elements of the same type have been initialized with
     * the same information.
     */
    void GeoModelEditor::complete_element_connectivity()
    {
        // Lines
        if( model_.nb_lines() > 0 ) {
            if( model_.line( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( GME::LINE ) ;
            }
            if( model_.line( 0 ).nb_in_boundary() == 0 ) {
                fill_elements_in_boundaries( GME::LINE ) ;
            }
            if( !model_.line( 0 ).parent_id().is_defined()
                && model_.nb_contacts() > 0 ) {
                fill_elements_parent( GME::LINE ) ;
            }
        }
        // Corners
        if( model_.nb_corners() > 0 && model_.corner( 0 ).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( GME::CORNER ) ;
        }
        // Surfaces - There MUST be at least one
        if( model_.surface( 0 ).nb_boundaries() == 0 ) {
            fill_elements_boundaries( GME::SURFACE ) ;
        }
        if( model_.surface( 0 ).nb_in_boundary() == 0 ) {
            fill_elements_in_boundaries( GME::SURFACE ) ;
        }
        if( !model_.surface( 0 ).parent_id().is_defined() ) {
            fill_elements_parent( GME::SURFACE ) ;
        }
        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( GME::REGION ) ;
            }
            if( !model_.region( 0 ).parent_id().is_defined()
                && model_.nb_layers() > 0 ) {
                fill_elements_parent( GME::REGION ) ;
            }
        }
        // Contacts
        if( model_.nb_contacts() > 0 && model_.contact( 0 ).nb_children() == 0 ) {
            fill_elements_children( GME::CONTACT ) ;
        }
        // Interfaces
        if( model_.nb_interfaces() > 0
            && model_.one_interface( 0 ).nb_children() == 0 ) {
            fill_elements_children( GME::INTERFACE ) ;
        }
        // Layers
        if( model_.nb_layers() > 0 && model_.layer( 0 ).nb_children() == 0 ) {
            fill_elements_children( GME::LAYER ) ;
        }
    }

    void GeoModelEditor::fill_elements_boundaries( GME::TYPE type ) const
    {
        // We have a problem if this is called for regions
        // No way yet to know the surface orientation
        ringmesh_assert( type != GME::REGION ) ;

        GME::TYPE b_type = GME::boundary_type( type ) ;
        if( b_type != GME::NO_TYPE ) {
            for( index_t i = 0; i < model().nb_elements( b_type ); ++i ) {
                const GME& b = model().element( gme_t( b_type, i ) ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_element_boundary( b.in_boundary_gme( j ),
                        gme_t( b_type, i ) ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_elements_in_boundaries( GME::TYPE type ) const
    {
        GME::TYPE in_b_type = GME::in_boundary_type( type ) ;
        if( in_b_type != GME::NO_TYPE ) {
            for( index_t i = 0; i < model().nb_elements( in_b_type ); ++i ) {
                const GME& in_b = element( gme_t( in_b_type, i ) ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_element_in_boundary( in_b.boundary_gme( j ),
                        gme_t( in_b_type, i ) ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_elements_parent( GME::TYPE type ) const
    {
        GME::TYPE p_type = GME::parent_type( type ) ;
        if( p_type != GME::NO_TYPE ) {
            for( index_t i = 0; i < model().nb_elements( p_type ); ++i ) {
                const GME& p = model().element( gme_t( p_type, i ) ) ;
                for( index_t j = 0; j < p.nb_children(); ++j ) {
                    set_element_parent( p.child_id( j ), gme_t( p_type, i ) ) ;
                }
            }
        }
    }

    void GeoModelEditor::fill_elements_children( GME::TYPE type ) const
    {
        GME::TYPE c_type = GME::child_type( type ) ;
        if( c_type != GME::NO_TYPE ) {
            for( index_t i = 0; i < model().nb_elements( c_type ); ++i ) {
                gme_t cur_child = gme_t( c_type, i ) ;
                const gme_t& parent = model().element( cur_child ).parent_id() ;
                if( parent.is_defined() ) {
                    add_element_child( parent, cur_child ) ;
                }
            }
        }
    }

    void GeoModelEditor::set_element(
        const GME::gme_t& id,
        GeoModelElement* E ) const
    {
        if( id.type < GME::NO_TYPE ) {
            model_.modifiable_elements( id.type )[id.index] = E ;
        } else {
            throw RINGMeshException( "GeoModelEditor", "GME type is not valid" ) ;
        }
    }

    /*!
     * @brief Add to the vector the elements which cannot exist if
     *        an element in the set does not exist.
     * @details These elements are added to the set.
     *          Recursive call till nothing is added.
     *
     * @return True if at least one element was added, otherwise false.
     */
    bool GeoModelEditor::get_dependent_elements( std::set< gme_t >& in ) const
    {
        index_t input_size = static_cast< index_t >( in.size() ) ;

        for( std::set< gme_t >::iterator it( in.begin() ); it != in.end(); ++it ) {
            gme_t cur = *it ;
            /// If an element has children elements - add them 
            if( GME::child_allowed( cur.type ) ) {
                const GME& E = model_.element( cur ) ;
                for( index_t j = 0; j < E.nb_children(); ++j ) {
                    in.insert( E.child_id( j ) ) ;
                }
            }
        }

        /// If a parent has no children anymore - add it 
        for( index_t p = GME::CONTACT; p < GME::NO_TYPE; ++p ) {
            GME::TYPE P = (GME::TYPE) p ;
            for( index_t j = 0; j < model_.nb_elements( P ); ++j ) {
                bool no_child = true ;
                const GME& E = model_.element( GME::gme_t( P, j ) ) ;
                for( index_t k = 0; k < E.nb_children(); ++k ) {
                    if( in.count( E.child_id( k ) ) == 0 ) {
                        no_child = false ;
                        break ;
                    }
                }
                if( no_child ) {
                    in.insert( E.gme_id() ) ;
                }
            }
        }

        /// If an element is in the boundary of nothing - add it
        for( index_t t = GME::CORNER; t < GME::REGION; ++t ) {
            GME::TYPE T = (GME::TYPE) t ;
            for( index_t j = 0; j < model_.nb_elements( T ); ++j ) {
                bool no_incident = true ;
                const GME& E = model_.element( GME::gme_t( T, j ) ) ;
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
            return get_dependent_elements( in ) ;
        } else {
            return false ;
        }
    }

    /*!
     * @brief Remove a list of elements of the model
     * @details No check is done on the consistency of this removal
     *          The elements and all references to them are removed.
     *          All dependent elements should be in the set of elements to remove,
     *          with a prior call to get_dependent_elements function.
     *
     * @warning NOT TESTED.
     *          The client is responsible to set the proper connectivity
     *          information between the remaining model elements.
     *
     * @todo TEST IT
     */
    void GeoModelEditor::remove_elements( const std::set< gme_t >& elements )
    {
        if( elements.empty() ) {
            return ;
        }

        // We need to remove elements type by type since they are 
        // stored in different vectors and since we use indices in these 
        // vectors to identify them.
        // Initialize the vector
        std::vector< std::vector< index_t > > to_erase_by_type ;
        to_erase_by_type.reserve( GME::NO_TYPE ) ;
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
            to_erase_by_type.push_back(
                std::vector< index_t >(
                    model_.nb_elements( static_cast< GME::TYPE >( i ) ), 0 ) ) ;
        }
        // Flag the elements to erase
        for( std::set< gme_t >::const_iterator it = elements.begin();
            it != elements.end(); ++it ) {
            gme_t cur = *it ;
            if( cur.type < GME::NO_TYPE ) {
                ringmesh_assert( NO_ID != 0 ) ; // If one day NO_ID changes of value.
                to_erase_by_type[cur.type][cur.index] = NO_ID ;
            }
        }

        delete_elements( to_erase_by_type ) ;

        // Re-initialize global access to elements
        model_.init_global_model_element_access() ;
    }

    /*!
     * @brief ONGOING WORK. Removes properly some elements of the Boundary Model.
     *
     * @param[in] elements_to_remove: in input the elements the client wants to
     * remove, in output all the removed elements (dependencies of ).
     *
     * Calls get_dependent_elements on each elements of \p elements_to_remove.
     * Then do remove these elements and updates the universe.
     *
     * @pre Assert that the element to remove is one region that has only one neighbor
     *
     * @todo Finish to implement it for any kind of BME. BC must continue this work.
     *
     * @todo Review : Error in the comments [JP]
     */
    void GeoModelEditor::remove_elements_and_dependencies(
        const std::set< GME::gme_t >& elements_to_remove )
    {
        // Asserts to remove when implementation is completed
        ringmesh_assert(
            elements_to_remove.size() == 1
                && elements_to_remove.begin()->type == GME::REGION ) ;

        // Copy because it is not logical to have in output the removed elements. BC
        /// @todo Review : youpiii what did the comment on the function just said ? [JP]
        std::set< GME::gme_t > elements = elements_to_remove ;
        // TODO Handle the case of several objects in elements

        const GeoModelElement& reg = model_.element( *( elements.begin() ) ) ;
        get_dependent_elements( elements ) ;

        // TODO Dirty duplication of code------------------------
        // We need to remove elements type by type since they are
        // stored in different vectors and since we use indices in these
        // vectors to identify them.
        // Initialize the vector
        std::vector< std::vector< index_t > > to_erase_by_type ;
        to_erase_by_type.reserve( GME::NO_TYPE ) ;
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; ++i ) {
            to_erase_by_type.push_back(
                std::vector< index_t >(
                    model_.nb_elements( static_cast< GME::TYPE >( i ) ), 0 ) ) ;
        }
        // Flag the elements to erase
        for( std::set< gme_t >::const_iterator it = elements.begin();
            it != elements.end(); ++it ) {
            gme_t cur = *it ;
            if( cur.type < GME::NO_TYPE ) {
                ringmesh_assert( NO_ID != 0 ) ; // If one day NO_ID changes of value.
                to_erase_by_type[cur.type][cur.index] = NO_ID ;
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

        std::vector< GME::gme_t > to_add_in_universe ;

        if( reg.type() == GME::REGION ) {
            index_t nb_added = 0 ;
            for( index_t b_i = 0; b_i < reg.nb_boundaries(); ++b_i ) {
                if( !reg.boundary( b_i ).is_on_voi() ) {
                    to_add_in_universe.push_back( reg.boundary( b_i ).gme_id() ) ;
                    ringmesh_assert(
                        to_erase_by_type[reg.boundary( b_i ).type()][reg.boundary(
                            b_i ).index()] != NO_ID ) ;
                    to_add_in_universe[nb_added].index =
                        to_erase_by_type[reg.boundary( b_i ).type()][reg.boundary(
                            b_i ).index()] ;
                    ++nb_added ;
                }
            }
        }

        remove_elements( elements ) ;

        // Update Universe
        /// @todo You first need to clean the existing universe [JP]
        for( std::vector< GME::gme_t >::const_iterator itr =
            to_add_in_universe.begin(); itr != to_add_in_universe.end(); ++itr ) {
            add_element_boundary( gme_t( GME::REGION, NO_ID ),
                gme_t( GME::SURFACE, itr->index ), true ) ;
        }

        //ringmesh_assert( model_.check_model_validity() ) ;
    }

    /*!
     * @brief Delete elements and remove all references to them in the model
     *
     * @param[in,out] to_erase For each type of element T,
     *        store a vector of the size of model_.nb_elements(T) in which
     *        elements are flagged with NO_ID.
     *        In output it stores the mapping table between old and new indices
     *        for the elements.
     * @todo TEST IT
     */
    void GeoModelEditor::delete_elements(
        std::vector< std::vector< index_t > >& to_erase )
    {
        // Number of elements deleted for each TYPE
        std::vector< index_t > nb_removed( to_erase.size(), 0 ) ;

        /// 1. Get the mapping between old indices of the elements
        ///    and new ones (when elements to remove will actually be removed)
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            for( index_t j = 0; j < to_erase[i].size(); ++j ) {
                if( to_erase[i][j] == NO_ID ) {
                    nb_removed[i]++ ;
                } else {
                    to_erase[i][j] = j - nb_removed[i] ;
                }
            }
        }

        /// 2. Effectively delete the elements 
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            for( index_t j = 0; j < to_erase[i].size(); ++j ) {
                if( to_erase[i][j] == NO_ID ) {
                    GME::gme_t cur( static_cast< GME::TYPE >( i ), j ) ;
                    delete model_.element_ptr( cur ) ;
                    // Set the current element to nil
                    model_.modifiable_elements( cur.type )[cur.index] = nil ;
                }
            }
            std::vector< GME* >& store = model_.modifiable_elements(
                static_cast< GME::TYPE >( i ) ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GME* >( nil ) ), store.end() ) ;
        }

        /// 3. Deal with the model vertices
        model_.mesh.vertices.clear() ;

        /// 4. Update all possible indices in remaining elements
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            GME::TYPE T = static_cast< GME::TYPE >( i ) ;

            // Update all indices stored by the BME of that type 
            ringmesh_assert(
                model_.nb_elements( T ) == to_erase[i].size() - nb_removed[i] ) ;
            for( index_t j = 0; j < model_.nb_elements( T ); ++j ) {
                GeoModelElement& E = model_.modifiable_element( gme_t( T, j ) ) ;

                // Not the same than j - since we have erased some elements
                index_t old_id = E.index() ;
                ringmesh_assert( to_erase[i][old_id] != NO_ID ) ;

                // id_ 
                E.id_.index = to_erase[i][old_id] ;
                // boundary_
                if( E.nb_boundaries() > 0 ) {
                    GME::TYPE B = GME::boundary_type( T ) ;
                    ringmesh_assert( B < GME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_boundaries(); ++k ) {
                        set_element_boundary( E.gme_id(), k,
                            gme_t( B, to_erase[B][E.boundary_gme( k ).index] ) ) ;
                    }
                }
                // in_boundary
                if( E.nb_in_boundary() > 0 ) {
                    GME::TYPE IB = GME::in_boundary_type( T ) ;
                    ringmesh_assert( IB < GME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
                        set_element_in_boundary( E.gme_id(), k,
                            gme_t( IB,
                                to_erase[IB][E.in_boundary_gme( k ).index] ) ) ;
                    }
                }
                // parent_
                if( E.has_parent() ) {
                    GME::TYPE P = GME::parent_type( T ) ;
                    ringmesh_assert( P < GME::NO_TYPE ) ;
                    set_element_parent( E.gme_id(),
                        gme_t( P, to_erase[P][E.parent_id().index] ) ) ;
                }
                // children_ 
                if( E.nb_children() > 0 ) {
                    GME::TYPE C = GME::child_type( T ) ;
                    ringmesh_assert( C < GME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_children(); ++k ) {
                        set_element_child( E.gme_id(), k,
                            gme_t( C, to_erase[C][E.child_id( k ).index] ) ) ;
                    }
                }
                // Clean the vectors in the element
                erase_invalid_element_references( E ) ;
            }
        }

        // Do not forget the universe
        /// @todo Put the universe in the list of regions - so annoying to think of it each time
        /// @todo BUG ? sides are lost for the universe ? not sure ... [JP]
        {
            Region& U = dynamic_cast< Region& >( model_.modifiable_element(
                GME::gme_t( GME::REGION, NO_ID ) ) ) ;
            for( index_t i = 0; i < U.nb_boundaries(); ++i ) {
                set_element_boundary( U.gme_id(), i,
                    gme_t( GME::SURFACE,
                        to_erase[GME::SURFACE][U.boundary_gme( i ).index] ) ) ;
            }
            erase_invalid_element_references( U ) ;
        }
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model elements and their relationship ignoring their geometry
     *
     * @param[in] from Model to copy the information from
     */
    void GeoModelEditor::copy_macro_topology( const GeoModel& from )
    {
        assert_element_creation_allowed() ;
        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            std::vector< GME* >& store = model_.modifiable_elements( T ) ;
            store.resize( from.nb_elements( T ), nil ) ;

            for( index_t e = 0; e < model_.nb_elements( T ); ++e ) {
                store[e] = new_element( T, e ) ;
                ringmesh_assert( store[ e ] != nil ) ;
            }
            RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < model_.nb_elements( T ); ++e ) {
                copy_element_topology( *store[e], from.element( gme_t( T, e ) ) ) ;
            }
        }
        copy_element_topology( model_.universe_, from.universe_ ) ;

        model_.nb_elements_per_type_ = from.nb_elements_per_type_ ;
    }

    void GeoModelEditor::copy_element_topology(
        GeoModelElement& lhs,
        const GeoModelElement& rhs )
    {
        lhs.name_ = rhs.name_ ;
        lhs.geol_feature_ = rhs.geol_feature_ ;
        lhs.boundaries_ = rhs.boundaries_ ;
        lhs.in_boundary_ = rhs.in_boundary_ ;
        lhs.parent_ = rhs.parent_ ;
        lhs.children_ = rhs.children_ ;

        if( lhs.type() == GME::REGION ) {
            Region& R_lhs = dynamic_cast< Region& >( lhs ) ;
            const Region& R_rhs = dynamic_cast< const Region& >( rhs ) ;
            R_lhs.sides_ = R_rhs.sides_ ;
        }
    }

    /*!
     * @brief Remove invalid reference to elements
     *       boundary, in_boundary and children vectors
     *       Invalid elements have a NO_ID index
     *
     * @param[in] E the element to process
     */
    void GeoModelEditor::erase_invalid_element_references( GeoModelElement& E )
    {
        GME::TYPE T = E.type() ;
        if( E.nb_children() > 0 ) {
            gme_t invalid_child( E.child_type( T ), NO_ID ) ;
            index_t nb_found = static_cast< index_t >( std::count(
                E.children_.begin(), E.children_.end(), invalid_child ) ) ;
            if( nb_found == E.children_.size() ) {
                // Calling erase on all elements -> undefined behavior
                E.children_.clear() ;
            } else {
                E.children_.erase(
                    std::remove( E.children_.begin(), E.children_.end(),
                        invalid_child ), E.children_.end() ) ;
            }
        }
        if( E.nb_boundaries() > 0 ) {
            gme_t invalid_boundary( E.boundary_type( T ), NO_ID ) ;
            if( E.type() == GME::REGION ) {
                Region& R = dynamic_cast< Region& >( E ) ;
                // Change side values if necessary
                index_t offset = 0 ;
                for( index_t i = 0; i + offset < E.nb_boundaries(); ++i ) {
                    if( E.boundaries_[i] == invalid_boundary ) {
                        offset++ ;
                    } else {
                        R.sides_[i] = R.side( i + offset ) ;
                    }
                }
            }

            index_t end = static_cast< index_t >( std::remove( E.boundaries_.begin(),
                E.boundaries_.end(), invalid_boundary ) - E.boundaries_.begin() ) ;
            if( end == 0 ) {
                E.boundaries_.clear() ;
                if( E.type() == GME::REGION ) {
                    Region& R = dynamic_cast< Region& >( E ) ;
                    R.sides_.clear() ;
                }
            } else {
                E.boundaries_.erase( E.boundaries_.begin() + end,
                    E.boundaries_.end() ) ;
                if( E.type() == GME::REGION ) {
                    Region& R = dynamic_cast< Region& >( E ) ;
                    R.sides_.erase( R.sides_.begin() + end, R.sides_.end() ) ;
                }
            }
        }
        if( E.nb_in_boundary() > 0 ) {
            gme_t invalid_in_boundary( E.in_boundary_type( T ), NO_ID ) ;
            index_t nb_found = static_cast< index_t >( std::count(
                E.in_boundary_.begin(), E.in_boundary_.end(), invalid_in_boundary ) ) ;
            if( nb_found == E.in_boundary_.size() ) {
                E.in_boundary_.clear() ;
            } else {
                E.in_boundary_.erase(
                    std::remove( E.in_boundary_.begin(), E.in_boundary_.end(),
                        invalid_in_boundary ), E.in_boundary_.end() ) ;
            }
        }
    }

    GME* GeoModelEditor::new_element( GME::TYPE T, index_t id )
    {
        assert_element_creation_allowed() ;
        if( T == GME::CORNER ) {
            return new Corner( model(), id ) ;
        } else if( T == GME::LINE ) {
            return new Line( model(), id ) ;
        } else if( T == GME::SURFACE ) {
            return new Surface( model(), id ) ;
        } else if( T == GME::REGION ) {
            return new Region( model(), id ) ;
        } else if( T > GME::REGION && T < GME::NO_TYPE ) {
            return new GeoModelElement( model(), T, id ) ;
        } else {
            return nil ;
        }
    }
    
    GME* GeoModelEditor::new_element( GME::TYPE T )
    {
        index_t id = model_.nb_elements( T ) ;
        return new_element( T, id ) ;
    }


}
