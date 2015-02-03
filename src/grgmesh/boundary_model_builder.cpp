/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
*/

/*! \author Jeanne Pellerin */


#include <grgmesh/boundary_model_builder.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>
#include <stack>

namespace GRGMesh {

     /*!
     * @brief Rebuild a model ???
     * \todo Comment rebuild()... 
     * \todo Valgrind finds errors !!!!!!!
     * @return 
     */
    bool BoundaryModelBuilder::rebuild()
    {
   /*     std::vector< index_t > sp_to_remove ;
        std::vector< index_t > new_sp_id( model_.nb_surfaces() ) ;
        index_t offset = 0 ;
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            if( model_.surfaces_[sp].nb_cells() == 0 ) {
                offset++ ;
                sp_to_remove.push_back( sp ) ;
            } else {
                model_.surfaces_[sp - offset] = model_.surfaces_[sp] ;
            }
            new_sp_id[sp] = sp - offset ;
        }
        if( offset == 0 ) return false ;
        model_.surfaces_.erase( model_.surfaces_.end() - offset,
            model_.surfaces_.end() ) ;

        offset = 0 ;
        std::vector< index_t > cp_to_remove ;
        std::vector< index_t > new_cp_id( model_.nb_lines() ) ;
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < line.nb_in_boundary(); sp++ ) {
                if( Utils::contains( sp_to_remove,
                    line.in_boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    line.in_boundary_[sp - nb_sp_removed] =
                        new_sp_id[line.in_boundary_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == line.nb_in_boundary() ) {
                    offset++ ;
                    cp_to_remove.push_back( cp ) ;
                    continue ;
                } else {
                    line.in_boundary_.erase(
                        line.in_boundary_.end() - nb_sp_removed,
                        line.in_boundary_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.lines_[cp - offset] = model_.lines_[cp] ;
            }
            new_cp_id[cp] = cp - offset ;
        }
        if( offset > 0 ) {
            model_.lines_.erase( model_.lines_.end() - offset,
                model_.lines_.end() ) ;
        }

        // Build the contacts
        // Update surfaces
        offset = 0 ;
        std::vector< index_t > s_to_remove ;
        std::vector< index_t > new_s_id( model_.nb_interfaces() ) ;
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            BoundaryModelElement& surface = model_.interfaces_[s] ;
            surface.boundaries_.clear() ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < surface.nb_children(); sp++ ) {
                if( Utils::contains( sp_to_remove, surface.child_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    surface.children_[sp - nb_sp_removed] =
                        new_sp_id[surface.children_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == surface.nb_children() ) {
                    offset++ ;
                    s_to_remove.push_back( s ) ;
                    continue ;
                } else {
                    surface.children_.erase( surface.children_.end() - nb_sp_removed,
                        surface.children_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.interfaces_[s - offset] = model_.interfaces_[s] ;
            }
            new_s_id[s] = s - offset ;
        }
        if( offset > 0 ) {
            model_.interfaces_.erase( model_.interfaces_.end() - offset,
                model_.interfaces_.end() ) ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            BoundaryModelElement& surface = model_.surfaces_[sp] ;
            surface.parent_ = new_s_id[surface.parent_] ;
            for( index_t cp = 0; cp < surface.nb_boundaries(); cp++ ) {
                surface.boundaries_[cp] =
                    new_cp_id[surface.boundaries_[cp]] ;
            }
        }
        model_.contacts_.clear() ;
        build_contacts() ;

        // Then finish the job (the order matters)
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].boundaries_.clear() ;
        }
        offset = 0 ;
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            BoundaryModelElement& region = model_.regions_[r] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < region.nb_boundaries(); sp++ ) {
                if( Utils::contains( sp_to_remove, region.boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    region.boundaries_[sp - nb_sp_removed] =
                        new_sp_id[region.boundaries_[sp]] ;
                    region.sides_[sp - nb_sp_removed] = region.sides_[sp] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == region.nb_boundaries() ) {
                    offset++ ;
                    continue ;
                } else {
                    region.sides_.erase( region.sides_.end() - nb_sp_removed,
                        region.sides_.end() ) ;
                    region.boundaries_.erase(
                        region.boundaries_.end() - nb_sp_removed,
                        region.boundaries_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.regions_[r - offset] = model_.regions_[r] ;
            }
        }
        if( offset > 0 ) {
            model_.regions_.erase( model_.regions_.end() - offset,
                model_.regions_.end() ) ;
        }
        end_surfaces() ;

        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            for( index_t j = 0; j < model_.contacts_[i].nb_in_boundary();
                ++j ) {
                index_t b = model_.contacts_[i].in_boundary_id( j ) ;
                add_interface_boundary( b, i ) ;
            }
        }

        end_lines() ;
        end_contacts() ;

        offset = 0 ;
        std::vector< index_t > co_to_remove ;
        std::vector< index_t > new_co_id( model_.nb_corners() ) ;
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            BoundaryModelElement& corner = model_.corners_[co] ;
            index_t nb_cp_removed = 0 ;
            for( index_t cp = 0; cp < corner.nb_in_boundary(); cp++ ) {
                if( Utils::contains( cp_to_remove, corner.in_boundary_id( cp ) ) ) {
                    nb_cp_removed++ ;
                } else {
                    corner.in_boundary_[cp - nb_cp_removed] =
                        new_cp_id[corner.in_boundary_[cp]] ;
                }
            }
            if( nb_cp_removed > 0 ) {
                if( nb_cp_removed == corner.nb_in_boundary() ) {
                    offset++ ;
                    co_to_remove.push_back( co ) ;
                    continue ;
                } else {
                    corner.in_boundary_.erase(
                        corner.in_boundary_.end() - nb_cp_removed,
                        corner.in_boundary_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.corners_[co - offset] = model_.corners_[co] ;
            }
            new_co_id[co] = co - offset ;
        }
        if( offset > 0 ) {
            model_.corners_.erase( model_.corners_.end() - offset,
                model_.corners_.end() ) ;
        }

        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            line.boundaries_[0] = new_co_id[line.boundaries_[0]] ;
            line.boundaries_[1] = new_co_id[line.boundaries_[1]] ;
        }

        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            BoundaryModelElement& contact = model_.contacts_[c] ;
            for( index_t co = 0; co < contact.nb_boundaries(); co++ ) {
                contact.boundaries_[co] = new_co_id[contact.boundaries_[co]] ;
            }
        }

        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            BoundaryModelElement& layer = model_.layers_[l] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < layer.nb_boundaries(); sp++ ) {
                if( Utils::contains( sp_to_remove, layer.boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    layer.boundaries_[sp - nb_sp_removed] =
                        new_sp_id[layer.boundaries_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                layer.boundaries_.erase( layer.boundaries_.end() - nb_sp_removed,
                    layer.boundaries_.end() ) ;
            }
        }

        update_all_ids() ;
        return true ;*/

        grgmesh_assert_not_reached ;
        /// \todo Implement the rebuild function for the BoundaryModel
        return false ;
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model elements and their relationship ignoring their geometry
     * 
     * @param[in] from Model to copy the information from
     */
    void BoundaryModelBuilder::copy_macro_topology( const BoundaryModel& from )
    {
        model_.name_ = from.name_ ;
        model_.corners_.resize( from.nb_corners(), Corner( &model_ ) ) ;
        model_.lines_.resize( from.nb_lines(), Line( &model_ ) ) ;
        model_.surfaces_.resize( from.nb_surfaces(), Surface( &model_ ) ) ;
        model_.regions_.resize( from.nb_regions(), BoundaryModelElement( &model_, BoundaryModelElement::REGION ) ) ;
        model_.layers_.resize( from.nb_layers(), BoundaryModelElement( &model_, BoundaryModelElement::LAYER ) ) ;
        model_.contacts_.resize( from.nb_contacts(), BoundaryModelElement( &model_, BoundaryModelElement::CONTACT ) ) ;
        model_.interfaces_.resize( from.nb_interfaces(), BoundaryModelElement( &model_, BoundaryModelElement::INTERFACE ) ) ;
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_corners(); i++ ) {
            model_.corners_[i].copy_macro_topology( from.corner( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_lines(); i++ ) {
            model_.lines_[i].copy_macro_topology( from.line( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_surfaces(); i++ ) {
            model_.surfaces_[i].copy_macro_topology( from.surface( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_layers(); i++ ) {
            model_.layers_[i].copy_macro_topology( from.layer( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_regions(); i++ ) {
            model_.regions_[i].copy_macro_topology( from.region( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_contacts(); i++ ) {
            model_.contacts_[i].copy_macro_topology( from.contact( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_interfaces(); i++ ) {
            model_.interfaces_[i].copy_macro_topology( from.one_interface( i ), model_ ) ;
        }
        model_.universe_.copy_macro_topology( from.universe_, model_ ) ;
    }

    /*!
     * @brief Update the indices stored by each element of the model \
     * according to its actual position in the corresponding vector in the model
     */
    void BoundaryModelBuilder::update_all_ids()
    {
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            model_.corners_[co].set_id( co ) ;
        }
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            model_.lines_[cp].set_id( cp ) ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].set_id( sp ) ;
        }
        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            model_.contacts_[c].set_id( c ) ;
        }
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            model_.interfaces_[s].set_id( s ) ;
        }
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            model_.regions_[r].set_id( r );
        }
        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            model_.layers_[l].set_id( l ) ;
        }
    }

    /*!
     * @brief Remove duplicates in the model vertices 
     * @details When reading the file the vertices are duplicated between the different Surfaces,
     * and new vertices are added for Corners. 
     * Compute the duplicates inside the vertices_ vector - Update the vertex vector - 
     * and update the reference to vertices for all model Corners and Surfaces
     * 
     * \todo Check the compute_unique_kdtree function - It is not correct the number of neighbors 
     * is predefined 
     */
    void BoundaryModelBuilder::make_vertices_unique()
    {
        MakeUnique unique( model_.vertices_ ) ;
        unique.unique( 5 ) ;
        model_.vertices_.resize(0) ;
        unique.unique_points( model_.vertices_ ) ;
        const std::vector< index_t >& old2new = unique.indices() ;

        for( index_t s = 0; s < model_.nb_surfaces(); s++ ) {
            Surface& surface = model_.surfaces_[s] ;
            for( index_t p = 0; p < surface.nb_vertices(); p++ ) {
                surface.set_vertex( p, old2new[surface.model_vertex_id(p)] ) ;
            }
        }
        for( index_t l = 0; l < model_.nb_lines(); l++ ) {
            Line& line = model_.lines_[l] ;
            for( index_t p = 0; p < line.nb_vertices(); p++ ) {
                line.set_vertex( p, old2new[line.model_vertex_id(p)] ) ;
            }
        }
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            Corner& corner = model_.corners_[co] ;
            corner.set_vertex( old2new[corner.model_vertex_id()] ) ;
        }       
    }


    /*!
     * @brief Creates a element of the given type and add it to the correct vector
     * The BoundaryModelElement is created from its type and its index
     *
     * @param[in] type Type of the element to create
     * @return The index of the created element
     */
    index_t BoundaryModelBuilder::create_element( BoundaryModelElement::TYPE type ) {
        index_t id = model_.nb_elements( type ) ;
        grgmesh_assert( id != NO_ID ) ;
        switch( type ) {
        case BoundaryModelElement::CORNER: 
            {
                model_.corners_.push_back( Corner( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::LINE:
            {
                model_.lines_.push_back( Line( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::SURFACE:
            {
                model_.surfaces_.push_back( Surface( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::REGION:
            {
                model_.regions_.push_back( BoundaryModelElement( &model_, BoundaryModelElement::REGION, id ) ) ;
                break ;
            }
        case BoundaryModelElement::CONTACT:
            {
                model_.contacts_.push_back( BoundaryModelElement( &model_,BoundaryModelElement::CONTACT, id ) ) ;
                break ;
            }
        case BoundaryModelElement::INTERFACE:
            {
                model_.interfaces_.push_back( BoundaryModelElement( &model_,BoundaryModelElement::INTERFACE, id ) ) ;
                break ;
            }
        case BoundaryModelElement::LAYER:
            {
                model_.layers_.push_back( BoundaryModelElement( &model_, BoundaryModelElement::LAYER, id ) ) ;
                break ;
            }
        default:   
            grgmesh_assert_not_reached ;
        }
        return id ;
    }

    /*!
     * @brief Get the index of the Corner at a given model point
     * @param[in] p_id Index of the point
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilder::find_corner( index_t p_id ) const 
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).model_vertex_id() == p_id ) return i ;
        }
        return NO_ID ;
    }
   
     /*!
     * @brief Create a corner at a vertex.
     * 
     * @param[in] index Index of the vertex in the model
     * @return Index of the Corner
     */
    index_t BoundaryModelBuilder::create_corner( index_t index )
    {
       index_t id = create_element( BoundaryModelElement::CORNER ) ; 
       set_corner( id, index ) ;       
       return id ;
    }

     /*!
     * @brief Find or create a corner at a vertex.
     * 
     * @param[in] index Index of the vertex in the model
     * @return Index of the Corner
     */
    index_t BoundaryModelBuilder::find_or_create_corner( index_t index )
    {
        index_t result = find_corner( index ) ;
        if( result != NO_ID ) return result ;
        else return create_corner( index ) ;
    }

    /*!
     * @brief Looks for a line in the model
     * 
     * @param[in] corner0 Starting corner index
     * @param[in] corner1 Ending corner index
     * @param[in] vertices Indices of the vertices on the line
     * @return NO_ID or the index of the Line
     */
    index_t BoundaryModelBuilder::find_line(
        const std::vector< index_t >& vertices ) const
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            if( model_.line(i).equal( vertices ) ) 
                return i ;
        }
        return NO_ID ;
    }

    /*! 
     * @brief Add a Line knowing only the indices of its points and set its boundary corners 
     * The corners are created if they do not exist.
     * Used in Geomodeling to convert a surface to a model
     */
    index_t BoundaryModelBuilder::create_line( const std::vector< index_t >& points ) {
        index_t id = create_element( BoundaryModelElement::LINE ) ;
        set_line( id, points ) ;
               
        // Find the indices of the corner at both extremities
        index_t c0 = find_or_create_corner( points.front() ) ;
        index_t c1 = find_or_create_corner( points.back() ) ;       
        add_element_boundary( BoundaryModelElement::LINE, id, c0 ) ;
        if( c1 != c0 ) add_element_boundary( BoundaryModelElement::LINE, id, c1 ) ;         

        return id ;
    }

     /*!
     * @brief Find or create a line 
     * 
     * @param[in] corner0 Starting corner index
     * @param[in] corner1 Ending corner index
     * @param[in] vertices Indices of the vertices on this Line
     * @return Index of the Line
     */
    index_t BoundaryModelBuilder::find_or_create_line(
        const std::vector< index_t >& vertices )
    {      
        index_t result = find_line( vertices ) ;
        if( result != NO_ID ) return result ;
        else return create_line( vertices ) ;
    }

  
    /*!
     * @brief Create a surface
     * 
     * @return Index of the Surface in the surfaces_ vector 
     */
    index_t BoundaryModelBuilder::create_surface()
    {      
        return create_element( BoundaryModelElement::SURFACE ) ; 
    }

    /*!
     * @brief Find a Contact
     * @param[in] interfaces Indices of the Interfaces determining the contact
     * @return NO_ID or index of the contact
     */
    index_t BoundaryModelBuilder::find_contact( const std::vector< index_t >& interfaces ) const
    {
        std::vector< const BoundaryModelElement* > comp( interfaces.size() ) ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            comp[i] = &model_.one_interface(interfaces[i]) ;
        }

        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            if( comp.size() == model_.contact(i).nb_in_boundary() ) {
                bool equal = true ;
                for( index_t j = 0; j < model_.contact(i).nb_in_boundary(); j++ ) {
                    if( comp[j] != &model_.contact(i).in_boundary( j ) ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Create a contact between the given Interfaces
     * The name of the contact is determined from the names of the interfaces.
     * 
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    index_t BoundaryModelBuilder::create_contact( const std::vector< index_t >& interfaces ) {        
        // Create a name for this contact
        std::string name = "contact_" ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            name += model_.interfaces_[interfaces[i]].name() ;
            name += "_" ;
        }
        
        index_t id = create_element( BoundaryModelElement::CONTACT ) ;
        set_element_name( BoundaryModelElement::CONTACT, id, name ) ;
        
        /*for( index_t i = 0; i < interfaces.size(); ++i ) {
            add_element_in_boundary( BoundaryModelElement::CONTACT, id, interfaces[i] ) ;
        }*/        
        return id ;
    }

    
    /*!
     * @brief Find or create a contact between given Interfaces
     * 
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    index_t BoundaryModelBuilder::find_or_create_contact(
        const std::vector< index_t >& interfaces )
    {
        index_t result = find_contact( interfaces ) ;
        if( result != NO_ID ) return result ;
        else return create_contact( interfaces ) ;        
    }

     /*!
     * @brief Get the index of an Interface from its name
     * 
     * @param[in] name Name of the Interface
     * @return Index of the interface in the model, NO_ID if not found.
     */
    index_t BoundaryModelBuilder::find_interface( const std::string& name ) const
    {
        for( index_t i = 0; i < model_.nb_interfaces(); ++i ) {
            if( model_.one_interface(i).name() == name ) {
                return i ;
            }
        }
        return NO_ID ;
    }


    
    /*!
     * @brief Create a new Interface 
     * 
     * @param[in] name Name of the interface     
     * @param[in] type Type of the interface
     * @return The Interface index.
     */
    index_t BoundaryModelBuilder::create_interface(
        const std::string& name,
        BoundaryModelElement::GEOL_FEATURE type )
    {
        index_t id = create_element( BoundaryModelElement::INTERFACE ) ;
        set_element_geol_feature( BoundaryModelElement::INTERFACE, id, type ) ;
        set_element_name( BoundaryModelElement::INTERFACE, id, name ) ;
        return id ;
    }

     /*
    * @brief Adds an empty region to the model
    *
    *  Used in Geomodeling to convert a surface to a model
    */
    index_t BoundaryModelBuilder::create_region() {
        return create_element( BoundaryModelElement::REGION ) ;
    }

    /*!
     * @brief Adds a new region to the model
     *
     * @param[in] name Name of the region
     * @param[in] boundaries Indices of the surfaces on the region boundary, plus indication on which
     *            side of the surface is the region
     * @return Index of the created region
     */
    index_t BoundaryModelBuilder::create_region(
        const std::string& name,
        const std::vector< std::pair< index_t, bool > >& boundaries )
    {
        index_t id = create_element( BoundaryModelElement::REGION ) ;
        set_element_name( BoundaryModelElement::REGION, id, name ) ;
        for( index_t i = 0; i < boundaries.size(); ++i ) {            
            add_element_boundary( BoundaryModelElement::REGION, id, boundaries[i].first, boundaries[i].second ) ;
        }
        return id ;
    }

    
    /*!
     * @brief Creates a new empty Layer with the given name 
     *
     * @param[in] name Name of the layer
     * @return The layer index
     */
    index_t BoundaryModelBuilder::create_layer(
        const std::string& name )
    {
        index_t id = create_element( BoundaryModelElement::LAYER ) ;
        set_element_name( BoundaryModelElement::LAYER, id, name ) ;
        return id ;
    }
    
    
    /*!
     * @brief Fill the model universe_
     *
     * @param[in] boundaries Indices of the surfaces on the model boundary
     * plus indication on which side of the surface is universe_ ( outside of the model )
     */
    void BoundaryModelBuilder::set_universe(
        const std::vector< std::pair< index_t, bool > >& boundaries )
    {
        model_.universe_.set_name( "Universe" ) ;
        model_.universe_.set_element_type( BoundaryModelElement::REGION ) ;
        model_.universe_.set_model( &model_ ) ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            grgmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary( boundaries[i].first,
                boundaries[i].second ) ;
            // If this surface have no type, set it at BoundaryModelElement::VOI
            model_.surfaces_[boundaries[i].first].set_geological_feature( BoundaryModelElement::VOI ) ;
        }
    }

    /*!
     * @brief Remove the universe from the vector of regions and update their indices
     *
     * @param[in] id Index of the Universe region in the regions_ vector
     */
    void BoundaryModelBuilder::remove_universe_from_regions( index_t id ) 
    {
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            index_t cur_id = model_.region(i).id() ;            
            if( i > id ) model_.regions_[i].set_id( cur_id-1 ) ;
        }
        model_.regions_.erase( model_.regions_.begin() + id ) ;
    }

   
    /*!
     * @brief Set the vertex for a Corner
     *
     * @param[in] corner_id Index of the corner
     * @param[in] vertex_id Index of the vertex in the model
     */
    void BoundaryModelBuilder::set_corner( index_t corner_id,  index_t vertex_id ) 
    {
        grgmesh_debug_assert( vertex_id < model_.nb_vertices() ) ;
        model_.corners_.at(corner_id).set_vertex( vertex_id ) ;
    }


    /*!
     * @brief Set one Line vertices
     *
     * @param[in] id Line index 
     * @param[in] vertices Indices of the vertices on the line
     */
    void BoundaryModelBuilder::set_line( 
        index_t id, 
        const std::vector< index_t >& vertices )
    {
        grgmesh_assert( id < model_.nb_lines() ) ;
        model_.lines_[id].set_vertices( vertices );
    }

   
    /*!
     * @brief Set the vertices and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] vertices Model indices of the vertices 
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     * @param[in] surface_adjacencies Adjacent facet (size of facet_ptr)
     */
    void BoundaryModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& vertices,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr,
        const std::vector< index_t >& surface_adjacencies )
    {        
        if( facets.size() == 0 ) return ;

        model_.surfaces_[surface_id].set_geometry( vertices, facets, facet_ptr ) ;

        if( surface_adjacencies.empty() )
            set_surface_adjacencies( surface_id ) ;
        else
            model_.surfaces_[surface_id].set_adjacent( surface_adjacencies ) ;
    }


    /*!
     * @brief Set the facets of a surface
     * @details If facet_adjacencies are not given they are computed
     *
     * @param[in] surface_id Index of the surface
     * @param[in] facets Indices of the model vertices defining the facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     * @param[in] surface_adjacencies Adjacent facet (size of facet_ptr)
     */
    void BoundaryModelBuilder::set_surface_geometry_bis(
        index_t surface_id,
        const std::vector< index_t >& corners,
        const std::vector< index_t >& facet_ptr,
        const std::vector< index_t >& corner_adjacent_facets )
    {
        if( corners.size() == 0 ) return ;

        model_.surfaces_[surface_id].set_geometry( corners, facet_ptr ) ;

        if( corner_adjacent_facets.empty() )
            set_surface_adjacencies( surface_id ) ;
        else
            model_.surfaces_[surface_id].set_adjacent( corner_adjacent_facets ) ;
    }

   
    /*!
     * @brief Compute and set the adjacencies between the facets 
     * @details The adjacent facet is given for each vertex of each facet for the edge
     * starting at this vertex.
     * If there is no neighbor inside the same Surface adjacent is set to NO_ADJACENT
     *
     * @param[in] surface_id Index of the surface
     */
    void BoundaryModelBuilder::set_surface_adjacencies( index_t surface_id ) 
    { 
        Surface& S = model_.surfaces_[surface_id] ;
        grgmesh_assert( S.nb_cells() > 0  ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( S.nb_cells()-1 ), Surface::NO_ADJACENT ) ;
              
        index_t nb_facets = S.nb_cells() ;
        index_t nb_vertices = S.nb_vertices() ;        
    
        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > toto ;
        toto.reserve( 10 ) ;
        std::vector< std::vector< index_t > > 
            vertex_to_facets( nb_vertices, toto ) ;              

        for( index_t f = 0; f < nb_facets; ++f )
        {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                vertex_to_facets[S.surf_vertex_id( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ){
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f )
        {
            for( index_t v = 0; v < S.nb_vertices_in_facet(f); v++ )
            {                
                index_t cur = S.surf_vertex_id(f, v) ;
                index_t prev = S.surf_vertex_id(f, S.prev_in_facet(f,v)) ;

                const std::vector< index_t >& f_prev = vertex_to_facets[prev] ;
                const std::vector< index_t >& f_cur = vertex_to_facets[cur] ;

                std::vector< index_t > inter(
                    std::min( f_prev.size(), f_cur.size() ) ) ;
                index_t end = std::set_intersection( f_prev.begin(),
                    f_prev.end(), f_cur.begin(), f_cur.end(), inter.begin() )
                    - inter.begin() ;

                if( end == 2 ) {
                    // There is one neighbor
                    index_t f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[ S.facet_begin(f) + S.prev_in_facet(f,v) ] = f2 ;
                } else {
                    grgmesh_debug_assert( end == 1 ) ;
                }
            }
        }
        S.set_adjacent( adjacent ) ;
    }

    /*! 
     * @brief Basic check of the validity of a BoundaryModelElement     
     * 
     * \todo Write meaningful message when the test fails ? 
     */
    bool BoundaryModelBuilder::check_basic_element_validity( const BoundaryModelElement& E ) const 
    {
        /// Verify that E points to the right BoundaryModel
        /// that its index and type are the right one.
        if( &E.model() != &model_ ) return false ;
        if( E.element_type() == BoundaryModelElement::NO_TYPE ) return false ;
        if( E.id() == NO_ID ) return false ;
        if( E.id() >= model_.nb_elements( E.element_type() ) ) return false ;
        if( !(model_.element(E.element_type(), E.id()) == E) ) return false ;


        /// Verify that the stored model vertex indices are in a valid range
        for( index_t i = 0; i < E.nb_vertices(); ++i ){            
            if( E.model_vertex_id(i) == NO_ID && 
                E.model_vertex_id(i) >= model_.nb_vertices() ) return false ;
        }
        return true ;
    }

    /*!
     * @brief Check that one element has the expected connectivity information
     * @details Requirements depend on the element type
     * See the static functions   ***_type( TYPE ) in class BME  
     */
    bool BoundaryModelBuilder::check_element_connectivity( const BoundaryModelElement& E ) const {
        BME::TYPE T = E.element_type() ;
        if( BME::boundary_type( T ) != BME::NO_TYPE ) {
            if( T != BME::SURFACE ) {
                // A closed surface - bubble might have no boundary
                // The others Line - and Region must have one
                if( E.nb_boundaries() == 0 ){
                    return false ;
                }
            }
        }
        // In_boundary
        if( BME::boundary_type( T ) != BME::NO_TYPE ) {
            if( E.nb_boundaries() == 0 ){
                return false ;
            }
        }
        // Parent - High level elements are not mandatory
        // But if the model has elements of the parent type, the element must have a parent
        if( BME::parent_type( T ) != BME::NO_TYPE ) {
            if( E.parent_id() == NO_ID && 
                model_.nb_elements( BME::parent_type(T) ) > 0 ){
                    return false ;
            }
        }
        // Children
        if( BME::child_type( T ) != BME::NO_TYPE ) {
            if( E.nb_children() == 0 ){
                return false ;
            }
        }
        return true ;
    }


    /*!
     * @brief Complete missing information in BoundaryModelElements
     * boundaries - in_boundary - parent - children 
     *
     * @details For all 7 types of elements, check what information is available
     * for the first one and fill the elements of the same type accordingly
     * THIS MEANS that the all the elements of the same type have bee initialized with
     * the the same information.
     */
    bool BoundaryModelBuilder::complete_element_connectivity() {             
        // Lines
        if( model_.nb_lines() > 0 ) {
            if( model_.line(0).nb_boundaries() == 0 ){
                fill_elements_boundaries( BoundaryModelElement::LINE ) ;
            }
            if( model_.line(0).nb_in_boundary() == 0 ){
                fill_elements_in_boundaries( BoundaryModelElement::LINE ) ;    
            }
            if( model_.line(0).parent_id() == NO_ID && model_.nb_contacts() > 0 ){
                fill_elements_parent( BoundaryModelElement::LINE ) ;
            }
        }
        // Corners
        if( model_.nb_corners() > 0 && model_.corner(0).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( BoundaryModelElement::CORNER ) ;
        }    
        // Surfaces - There MUST be at least one
        if( model_.surface(0).nb_boundaries() == 0 ) {            
            fill_elements_boundaries( BoundaryModelElement::SURFACE ) ;
        }        
        if( model_.surface(0).nb_in_boundary() == 0 ) {            
            fill_elements_in_boundaries( BoundaryModelElement::SURFACE ) ;    
        }
        if( model_.surface(0).parent_id() == NO_ID ) {            
            fill_elements_parent( BoundaryModelElement::SURFACE ) ;    
        }
        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region(0).nb_boundaries() == 0 ) {
                fill_elements_boundaries( BoundaryModelElement::REGION ) ;
            }
            if( model_.region(0).parent_id() == NO_ID && model_.nb_layers() > 0 ){
                fill_elements_parent( BoundaryModelElement::REGION ) ;
            }
        }
        // Contacts
        if( model_.nb_contacts() > 0 &&
            model_.contact(0).nb_children() == 0 ) {
            fill_elements_children( BoundaryModelElement::CONTACT ) ;
        }
        // Interfaces
        if( model_.nb_interfaces() > 0 &&
            model_.one_interface(0).nb_children() == 0 ) {
            fill_elements_children( BoundaryModelElement::INTERFACE ) ;
        }
        // Layers
        if( model_.nb_layers() > 0 &&
            model_.layer(0).nb_children() == 0 ) {
                fill_elements_children( BoundaryModelElement::LAYER ) ;
        }
        return true ;
    }

    /*!
     * @brief Fill the boundaries of all elements of the given type
     *
     * @details If the boundary elements do not have any in_boundary
     * information, nothing is done, and model construction will eventually fail.
     */
    void BoundaryModelBuilder::fill_elements_boundaries( BoundaryModelElement::TYPE type ) 
    {
        // We have a problem if this is called for regions
        // No way yet to know the surface orientation
        grgmesh_debug_assert( type != BME::REGION ) ;

        BME::TYPE b_type = BME::boundary_type( type ) ;
        if( b_type != BME::NO_TYPE ) {            
            for( index_t i = 0; i < model_.nb_elements( b_type ); ++i ) {
                const BME& b = model_.element( b_type, i ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_element_boundary( type, b.in_boundary_id(j), i ) ;
                }
            }
        }
    }
    /*!
     * @brief Fill the in_boundary vector of all elements of the given type
     *
     * @details If the in_boundary elements do not have any boundary
     * information, nothing is done, and model construction will eventually fail.
     */
    void BoundaryModelBuilder::fill_elements_in_boundaries( BoundaryModelElement::TYPE type ) 
    {
        BME::TYPE in_b_type = BME::in_boundary_type( type ) ;
        if( in_b_type != BME::NO_TYPE ) {            
            for( index_t i = 0; i < model_.nb_elements( in_b_type ); ++i ) {
                const BME& in_b = model_.element( in_b_type, i ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_element_in_boundary( type, in_b.boundary_id(j), i ) ;
                }
            }
        }        
    }
    /*!
     * @brief Fill the parent of all elements of the given type
     *
     * @details If the parents do not have any child 
     *  nothing is done, and model construction will eventually fail.
     */
    void BoundaryModelBuilder::fill_elements_parent( BoundaryModelElement::TYPE type ) 
    {
        BME::TYPE p_type = BME::parent_type( type ) ;
        if( p_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < model_.nb_elements( p_type ); ++i ) {
                const BME& p = model_.element( p_type, i ) ;
                for( index_t j = 0; j < p.nb_children(); ++j ) {
                    set_parent( type, p.child_id(j), i ) ;
                }
            }
        }              
    }
    /*!
     * @brief Fill the children of all elements of the given type
     *
     * @details If the children elements do not have any parent information
     * nothing is done, and model construction will eventually fail.
     */
    void BoundaryModelBuilder::fill_elements_children( BoundaryModelElement::TYPE type ) 
    {
        BME::TYPE c_type = BME::child_type( type ) ;
        if( c_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < model_.nb_elements( c_type ); ++i ) {
                index_t parent = model_.element( c_type, i ).parent_id() ;
                if( parent != NO_ID ) add_child( type, parent, i ) ;
            }
        }
    }                                                                          

    /*! 
     * @brief Fills the model nb_elements_per_type_ vector
     * @details See global element access with BoundaryModel::element( BME::TYPE, index_t )
     */
    void BoundaryModelBuilder::init_global_model_element_access() {        
        index_t count = 0 ;
        model_.nb_elements_per_type_.push_back( count ) ;
        for( index_t type = BME::CORNER; type < BME::NO_TYPE; type++ ) {
            count += model_.nb_elements( (BME::TYPE) type ) ;
            model_.nb_elements_per_type_.push_back( count ) ;
        }       
    }
    /*!
    * @brief Fills the model nb_facets_in_surfaces_ vector
    * @detail See global facet index accessor BoundaryModel::surface_facet
    */
    void BoundaryModelBuilder::init_global_model_facet_access() {
        model_.nb_facets_in_surfaces_.resize( model_.nb_surfaces()+1, 0 ) ;
        index_t count = 0 ;
        for( index_t i = 1; i < model_.nb_facets_in_surfaces_.size(); ++i ) {
            count += model_.surface( i-1 ).nb_cells() ;
            model_.nb_facets_in_surfaces_[i] = count ;
        }
    }

     /*!
     * @brief Last function to call when building a model
     *
     * @details check that the model is correct and has all required information
     * Calls the complete_element_connectivity function 
     * Fills nb_elements_per_type_ and nb_facets_in_surfaces_ vectors
     * 
     * @return False if the model is not valid and cannot be fixed
     * otherwise returns true.
     *
     *  \todo FULL CHECK AND FIX OF THE MODEL CORRECTNESS !!
        \todo Trade the end_something functions in the BoundaryModelBuilder
        functions for a smart end_model function
        that checks model validity and complete all missing parts
     */
    bool BoundaryModelBuilder::end_model()
    {
        // The name should exist
        if( model_.name() == "" ) set_model_name("model_default_name") ;
        // There must be at least 3 vertices
        if( model_.nb_vertices() == 0 ) return false ;
        // And at least one surface
        if( model_.nb_surfaces() == 0 ) return false ;
        
        // The Universe
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)

        init_global_model_element_access() ;
        init_global_model_facet_access() ;

        complete_element_connectivity() ;
        
        /// 1. Check that all the elements of the BoundaryModel have 
        ///    the required attributes.       
        for( index_t i = 0; i < model_.nb_elements( BME::ALL_TYPES ); ++i ){
            const BME& E = model_.element( BME::ALL_TYPES, i ) ;
            if( !check_basic_element_validity( E ) ){
                return false ;
            }
            if( !check_element_connectivity( E ) ){
                return false ;
            }
        }
   
        /// 2. \todo Check the consistency of connectivity relationships between the elements


        /// 3. \todo Check the geometrical consistency of the topological relationships
        
        

#ifdef GRGMESH_DEBUG
        std::cout << "Model " << model_.name() <<" has " << std::endl 
            << std::setw(10) << std::left << model_.nb_vertices()   << " vertices "   << std::endl 
            << std::setw(10) << std::left << model_.nb_facets()   << " facets "   << std::endl  
            << std::setw(10) << std::left << model_.nb_regions()  << " regions "  << std::endl
            << std::setw(10) << std::left << model_.nb_surfaces() << " surfaces " << std::endl
            << std::setw(10) << std::left << model_.nb_lines()    << " lines "    << std::endl 
            << std::setw(10) << std::left << model_.nb_corners()  << " corners "  << std::endl
            << std::endl ;
#endif
        return true ;
    }


    /*!
     * @brief Load and build a BoundaryModel from a Gocad .ml file
     * 
     *  @details This is pretty tricky because of the annoying not well adapted file format. 
     * The correspondance between Gocad::Model3D elements and BoundaryModel elements is :
     * Gocad TSurf  <-> BoundaryModel Interface
     * Gocad TFace  <-> BoundaryModel Surface
     * Gocad Region <-> BoundaryModel Region
     * Gocad Layer  <-> BoundaryModel Layer
     *
     * @param[in] in Input .ml file stream
     */
    void BoundaryModelBuilderGocad::load_ml_file( std::istream& in )
    {
        // Clear the model_ 
        //  Not sure this is actually useful
        model_.clear() ;

        time_t start_load, end_load ;
        time( &start_load ) ;        

        // Count the number of TSurf - Interface
        index_t nb_tsurf = 0 ;
        // Count the number of TFace - Surface
        index_t nb_tface = 0 ;

        // Counters identifying the currently read TSurf or TFace
        index_t tsurf_count = 0 ;
        index_t tface_count = 0 ;

        index_t current_nb_tfaces = 0 ;
        index_t nb_tface_in_prev_tsurf = 0 ;

        // The file contains 2 parts and is read in 2 steps
        // 1. Read model info (true)
        // 2. Read surfaces geometry (false)
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and betweeen TSurfs
        
        // Indices of the vertices of the currently built TSurf in the model
        std::vector< index_t > tsurf_vertex_ptr ;        
        // Where the vertices of a TFace start in the vertices of the TSurf (offest)
        std::vector< index_t > tface_vertex_start ;

        // Indices of vertices in facets (triangles) of the currently built TFace
        std::vector< index_t > tface_facets ;
        // Starting and ending indices of each facet triangle in the tface_facets vector
        std::vector< index_t > tface_facets_ptr ;
        tface_facets_ptr.push_back( 0 ) ;
        
        // Intermediate information for contact parts building
        std::vector< Border > borders_to_build ;

        // Surfaces for which the KeyFacet orientation should be changed
        // because it does not match triangle orientations.
        std::vector< index_t > change_key_facet ;

        // Begin reading the file 
        InputStream lis( in ) ;

        while( !lis.eof() ) {
            lis.get_line() ;
            std::string keyword ;
            lis >> keyword ;

            if( read_model ) {
                if( std::string( keyword, 0, 5 ) == "name:" ) {                    
                    set_model_name( std::string( keyword, 5 ) );
                }
                else if( keyword == "TSURF" ) {
                    /// 1. Read the TSurf information and create 
                    /// the corresponding Interface from its name
                    std::string temp_str ;
                    std::stringstream tsurf_name ;
                    lis >> temp_str ;
                    tsurf_name << temp_str ;
                    while( !lis.eol() ) {
                        lis >> temp_str ;
                        tsurf_name << "_" << temp_str ;
                    }
                    create_interface( tsurf_name.str() ) ;
                    nb_tsurf++ ;
                } else if( keyword == "TFACE" ) {
                    /// 2. Read the TFace info and build the 
                    /// corresponding Surface from its parent Interface, its type, and 
                    /// its key facet - from which + and - side are determined
                    index_t id ;
                    lis >> id ;
                    std::string type ;
                    lis >> type ;

                    std::string temp_str ;
                    std::stringstream tsurf_name ;
                    lis >> temp_str ;
                    tsurf_name << temp_str ;
                    while( !lis.eol() ) {
                        lis >> temp_str ;
                        tsurf_name << "_" << temp_str ;
                    }
                    // Get the key facet that give the orientation of the surface part
                    // Triangles in Gocad clockwise
                    vec3 p0, p1, p2 ;
                    lis.get_line() ;
                    lis >> p0 ;
                    lis.get_line() ;
                    lis >> p1 ;
                    lis.get_line() ;
                    lis >> p2 ;

                    create_surface( tsurf_name.str(), type,
                        Surface::KeyFacet( p0, p1, p2 ) ) ;
                    nb_tface++ ;

                } else if( keyword == "REGION" ) {
                    /// 3. Read Region information and create them from their name,
                    /// the surfaces on their boundary                    
                    index_t id ;
                    std::string name ;
                    lis >> id >> name ;

                    std::vector< std::pair< index_t, bool > > region_boundaries ;
                    bool end_region = false ;

                    while( !end_region ) {
                        lis.get_line() ;
                        for( index_t i = 0; i < 5; ++i ) {
                            signed_index_t tface_id ;
                            lis >> tface_id ;
                            if( tface_id == 0 ) {
                                end_region = true ;
                                break ;
                            } else {
                                // Correction because ids begin at 1 in the file
                                tface_id = tface_id > 0 ? tface_id-1 : -tface_id-1 ;                                            
                                region_boundaries.push_back(
                                    std::pair< index_t, bool >( tface_id, false ) ) ;
                            }
                        }
                    }
                    // The Universe is not a regular region, it is the regions
                    // outside the volume of interest 
                    if( name != "Universe" ) {
                        create_region( name, region_boundaries ) ;
                    }
                    else {
                        set_universe( region_boundaries ) ;
                    }
                } else if( keyword == "LAYER" ) {
                    /// 4. Build the volumetric layers from their name and 
                    /// the regions in them
                    std::string name ;
                    lis >> name ;
                    index_t layer_id = create_layer( name ) ;
                    bool end_layer = false ;
                    while( !end_layer ) {
                        lis.get_line() ;
                        for( index_t i = 0; i < 5; ++i ) {
                            index_t region_id ;
                            lis >> region_id ;
                            if( region_id == 0 ) {
                                end_layer = true ;
                                break ;
                            } else {
                                region_id -= nb_tface+1 ; // Remove Universe region
                                // Correction because ids begin at 1 in the file
                                add_child( BoundaryModelElement::LAYER, layer_id, region_id-1 ) ;
                            }
                        }
                    }
                } else if( keyword == "END" ) {
                    // End of the high level information on the model
                    // Switch to reading the geometry of the model surfaces
                    read_model = false ;
                    continue ;
                }
            } else {
                if( keyword == "GOCAD" ) {
                    // This is the beginning of a new TSurf = Interface
                    tsurf_count++ ;
                }
                if( keyword == "ZPOSITIVE" ) {
                    std::string positive ;
                    lis >> positive ;
                    if( positive == "Elevation" ) z_sign = 1 ;
                    else if( positive == "Depth" ) z_sign = -1 ;
                } else if( keyword == "END" ) {
                    // This the END of a TSurf
                    if( tsurf_count > 0 ) {
                        // End the last TFace - Surface of this TSurf
                        set_surface_geometry( 
                            tface_count-1,
                            std::vector< index_t >(
                                tsurf_vertex_ptr.begin() + tface_vertex_start.back(),
                                tsurf_vertex_ptr.end() ),
                            tface_facets,
                            tface_facets_ptr ) ;

                        if( !check_key_facet_orientation( tface_count-1 ) ){
                            change_key_facet.push_back( tface_count-1 ) ;
                        }

                        tface_facets.clear() ;
                        tface_facets_ptr.clear() ;
                        tface_facets_ptr.push_back( 0 ) ;

                        // End this TSurf - Interface
                        nb_tface_in_prev_tsurf += tface_vertex_start.size() ;
                        tsurf_vertex_ptr.clear() ;
                        tface_vertex_start.clear() ;
                    }
                } else if( keyword == "TFACE" ) {
                    // Beginning of a new TFace - Surface
                    if( tface_vertex_start.size() > 0 ) {
                        // End the previous TFace - Surface

                        set_surface_geometry( 
                            tface_count-1,std::vector< index_t >(
                                tsurf_vertex_ptr.begin() + tface_vertex_start.back(),
                                tsurf_vertex_ptr.end() ),
                            tface_facets,
                            tface_facets_ptr ) ;

                        if( !check_key_facet_orientation( tface_count-1 ) ) {
                            change_key_facet.push_back( tface_count-1 ) ;
                        }
                        
                        tface_facets.clear() ;
                        tface_facets_ptr.clear() ;
                        tface_facets_ptr.push_back( 0 ) ;
                    }
                    // Register where begin the new TFace vertices
                    tface_vertex_start.push_back( tsurf_vertex_ptr.size() ) ;

                    tface_count++ ;
                }
                /// 4. Read the surface vertices and facets (only triangles in Gocad Model3d files)
                else if( keyword == "VRTX" || keyword == "PVRTX" ) {
                    index_t id ;
                    vec3 p ;
                    lis >> id >> p ;
                    p.z *= z_sign ;
                    tsurf_vertex_ptr.push_back( add_vertex( p ) ) ;
                } else if( keyword == "PATOM" || keyword == "ATOM" ) {
                    index_t id ;
                    index_t v_id ;
                    lis >> id >> v_id ;
                    tsurf_vertex_ptr.push_back( tsurf_vertex_ptr[v_id - 1] ) ;
                } else if( keyword == "TRGL" ) {
                    // Ids of the vertices of each triangle in the TSurf
                    index_t p1, p2, p3 ;
                    lis >> p1 >> p2 >> p3 ;
                    // Change to ids in the TFace
                    p1 += -tface_vertex_start.back()-1 ;
                    p2 += -tface_vertex_start.back()-1 ;
                    p3 += -tface_vertex_start.back()-1 ;

                    tface_facets.push_back( p1 ) ;
                    tface_facets.push_back( p2 ) ;
                    tface_facets.push_back( p3 ) ;
                    tface_facets_ptr.push_back( tface_facets.size() ) ;
                }
                /// 5. Build the corners from their position and the surface parts
                ///    containing them
                else if( keyword == "BSTONE" ) {
                    index_t v_id ;
                    lis >> v_id ;
                    // correction to start at 0
                    v_id-- ;

                    // Get the TFace
                    index_t part_id = tface_vertex_start.size() - 1 ;
                    for( index_t i = 0; i < tface_vertex_start.size(); ++i ) {
                        if( v_id < tface_vertex_start[i] ) {
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    part_id += nb_tface_in_prev_tsurf ;

                    // c'est plus bon �a -compare geometry
                    index_t new_c = find_corner( model_.vertex( tsurf_vertex_ptr[v_id] ) ) ;
                    if( new_c == NO_ID ) create_corner( tsurf_vertex_ptr[v_id] ) ;
                }
                /// 6. Read the Border information and store it
                else if( keyword == "BORDER" ) {
                    index_t id, p1, p2 ;
                    lis >> id >> p1 >> p2 ;
                    p1-- ;
                    p2-- ;

                    // Get the global corner id
                    index_t corner_id = find_corner( model_.vertex( tsurf_vertex_ptr[p1] ) ) ;
                    grgmesh_assert( corner_id != NO_ID ) ;

                    // Get the surface
                    index_t part_id = NO_ID ;
                    for( index_t i = 0; i < tface_vertex_start.size(); ++i ) {
                        if( p1 < tface_vertex_start[i] ) {
                            grgmesh_assert( p2 < tface_vertex_start[i] ) ;

                            // Get vertices ids in the surface
                            p1 += -tface_vertex_start[i - 1] ;
                            p2 += -tface_vertex_start[i - 1] ;

                            // i-1 is the id of the TFace in this TSurf
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    if( part_id == NO_ID ) {
                        // It is in the last built Tface
                        p1 += -tface_vertex_start[tface_vertex_start.size() - 1] ;
                        p2 += -tface_vertex_start[tface_vertex_start.size() - 1] ;
                        part_id = tface_vertex_start.size() - 1 ;
                    }
                    // The number of tfaces in previous tsurf is also to add
                    part_id += nb_tface_in_prev_tsurf ;

                    borders_to_build.push_back(
                        Border( part_id, corner_id, p1, p2 ) ) ;
                }
            }
        }
    
        make_vertices_unique() ;

        /// 7. Build the Lines
        build_lines( borders_to_build ) ;

        /// 8. Build the Contacts
        build_contacts() ;
    
        /// Eventual changes of the key facet orientation
        end_surfaces( change_key_facet ) ;
     
        // Finish up the model - CRASH if this failed
        grgmesh_assert( end_model() ) ;
        
        time( &end_load ) ;
#ifdef GRGMESH_DEBUG
        std::cout << "Info" << " Boundary model loading time"
            << difftime( end_load, start_load ) << " sec" << std::endl ;
#endif
    }
    

    /*!
     * @brief Find the facet which first 3 vertices are given 
     * @param[in] surface_id Index of the surface
     * @param[in] p0 First point coordinates
     * @param[in] p1 Second point coordinates
     * @param[in] p2 Third point coordinates
     * @param[out] same_sign Is true if the found facet has the same orientation than triangle p0p1p2
     * @return Index of the found facet, NO_ID if none found
     */
    index_t BoundaryModelBuilderGocad::find_key_facet( 
        index_t surface_id, const vec3& p0, const vec3& p1,
        const vec3& p2, bool& same_sign ) const 
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;
    
        for( index_t t = 0; t < surface.nb_cells(); ++t ){            
            const vec3& pp0 = model_.vertex( surface.model_vertex_id( t, 0 ) ) ;
            const vec3& pp1 = model_.vertex( surface.model_vertex_id( t, 1 ) ) ;
            const vec3& pp2 = model_.vertex( surface.model_vertex_id( t, 2 ) ) ;
            
            if( p0 == pp0 ) {
                if( p1 == pp1 && p2 == pp2 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp1 ) {
                    same_sign = false ;
                    return t ;
                }
            }
            if( p0 == pp1 ) {
                if( p1 == pp0 && p2 == pp2 ) {
                    same_sign = false ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp0 ) {
                    same_sign = true ;
                    return t ;
                }
            }
            if( p0 == pp2 ) {
                if( p1 == pp0 && p2 == pp1 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp1 && p2 == pp0 ) {
                    same_sign = false ;
                    return t ;
                }
            }
        }
        return NO_ID ;
    }


    /*!
     * @brief Verify that a surface key facet has an orientation consistent with the surface facets.
     * @param[in] surface_id Index of the surface
     * @return False if the key_facet orientation is not the same than the surface facets, else true.
     */
    bool BoundaryModelBuilderGocad::check_key_facet_orientation( index_t surface_id ) 
    {
        const Surface& S = model_.surface( surface_id ) ;
        const Surface::KeyFacet& key_facet = S.key_facet() ;

        if( key_facet.is_default() ) {
            set_surface_first_triangle_as_key( surface_id ) ;
            return true ;
        }
        else {
            const vec3& p0 = key_facet.p0_ ;
            const vec3& p1 = key_facet.p1_ ;
            const vec3& p2 = key_facet.p2_ ;
            bool same_sign = false ;

            index_t t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ;
            if( t == NO_ID ) {
                vec3 p00( p0 ); 
                p00.z *= -1 ;
                vec3 p10( p1 ) ;
                p10.z *= -1 ;
                vec3 p20( p2 ) ;
                p20.z *= -1 ;
                // It is because of the sign of Z that is not the same 
                t = find_key_facet( surface_id, p00, p10, p20, same_sign ) ;
            }
            grgmesh_assert( t != NO_ID ) ;
            return same_sign ;
        }
    }
       

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *   
     * WE ASSUME THAT THE STORAGE OF THE POINTS IS UNIQUE IN THE MODEL AND THAT 
     * SURFACES DO SHARE POINTS ON THEIR CONTACT LINES
     * make_vertices_unique() must have been called first
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_ids Indices of vertices on the Line (resized at 0 at the beginning)
     * @return Index of the Corner at which the Line ends
     */
    index_t BoundaryModelBuilderGocad::determine_line_vertices( 
        const Surface& S, 
        index_t id0, 
        index_t id1,
        std::vector< index_t >& border_vertex_model_ids 
        ) const 
    {
        grgmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_ids.resize( 0 ) ;
                 
        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        grgmesh_assert( f != Surface::NO_ID ) ;

        // Global ids at the model level 
        index_t p0 = S.model_vertex_id( id0 ) ;
        index_t p1 = S.model_vertex_id( id1 ) ;

        border_vertex_model_ids.push_back( p0 ) ;
        border_vertex_model_ids.push_back( p1 ) ;
            
        index_t p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        while( p1_corner == NO_ID ) {

            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id(f, id0), S.facet_vertex_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;

            grgmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;
            
            index_t next_id1 =  S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.model_vertex_id( next_id1 ) ;
            border_vertex_model_ids.push_back( p1 ) ;         
            p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        }
        return p1_corner ; 
    }
    
    /*!
     * @brief Creates all Lines for the model
     *
     * @param[in] borders Information on Surface boundaries gathered at .ml file reading
     */
    void BoundaryModelBuilderGocad::build_lines( const std::vector< Border >& borders )
    {      
        std::vector< index_t > global_ids ;

        for( index_t i = 0; i < borders.size(); ++i ) {
            const Border& b = borders[i] ;

            /// 1- Build the boundary : construct the vector
            /// of vertices on the border
            const Surface& S = model_.surface( b.part_id_) ;

            index_t end_corner_id = determine_line_vertices( 
                S, b.p0_, b.p1_, global_ids ) ;           

            /// 2 - Check if this border already exists
            index_t line_id = find_or_create_line( global_ids ) ;

            // Add the surface in which this line is
            add_element_in_boundary( BoundaryModelElement::LINE, line_id, b.part_id_ ) ;
        }
    }

    /*!
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces     
     */
    void BoundaryModelBuilderGocad::build_contacts()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< index_t > interfaces ;
            //std::vector< GEOL_FEATURE > types ;
            for( index_t j = 0; j < model_.line(i).nb_in_boundary(); ++j ) {
                index_t sp_id = model_.line(i).in_boundary_id( j ) ;
                const BoundaryModelElement& p = model_.surface(sp_id).parent() ;
                interfaces.insert( p.id() ) ;
                //types.push_back( p.geological_feature() ) ;
            }
            std::vector< index_t > toto( interfaces.begin(), interfaces.end() ) ;
            index_t contact_id = find_or_create_contact( toto ) ;
            add_child( BoundaryModelElement::CONTACT, contact_id, i ) ;
        }
    }

    /*!
     * @brief Finish up Surface 
     * @details Calls the end_surfaces() function and switch KeyFacet orientation of the required
     * surfaces.
     *
     * @param[in] change_orientation Indices of the surfaces in which KeyFacet orientation is not
     * consistent with the its facet orientation
     */
    void BoundaryModelBuilderGocad::end_surfaces(
        const std::vector< index_t >& change_orientation )
    {
        //end_surfaces() ;

        for( index_t i = 0; i < change_orientation.size(); i++ ) {
            index_t s_i = change_orientation[i] ; 
            
            // Change the key facet            
            set_surface_first_triangle_as_key( s_i ) ; 
            const Surface& S = model_.surface( s_i ) ;

            // Change the sign of this Surface in all the regions containing it
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {                
                BoundaryModelElement& R =
                    const_cast< BoundaryModelElement& >( S.in_boundary(j) ) ;              

                for( index_t b = 0; b < R.nb_boundaries(); ++b ){
                    if( R.boundary_id(b) == s_i ){
                        bool old_side = R.side( b ) ;
                        R.set_boundary( b, R.boundary_id(b), !old_side ) ;
                    }
                }
            }
        }
    }
    

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] key KeyFacet for this Surface
     */
    void BoundaryModelBuilderGocad::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const Surface::KeyFacet& key )
    {
        index_t parent = find_interface( interface_name ) ;
        if( interface_name != "" ) grgmesh_assert( parent != NO_ID ) ;

        index_t id = create_element( BoundaryModelElement::SURFACE ) ;
        set_parent( BoundaryModelElement::SURFACE, id, parent ) ;
        //set_element_geol_feature( BoundaryModelElement::SURFACE, id, t ) ;
        set_surface_key_facet( id, key ) ;
    }


    /*!
     * @brief Get the index of the Corner at given coordinates
     * @param[in] p Coordinates of the vertex
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilderGocad::find_corner( const vec3& p ) const
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).vertex() == p ) return i ;
        }
        return NO_ID ;
    }


    /*!
    * @brief Utility structure to build a BoundaryModel knowing only its surface
    * @details Store the vertices of a triangle that is on the boundary of a surface
    */
    struct BorderTriangle {
        /*!
         * @brief Constructor 
         * @param s Index of the surface 
         * @param f Index of the facet containing the 3 vertices
         * @param vi Indices in the BoundaryModel of the vertices defining the triangle
         *           the edge v0 - v1 is the one on the boundary
         */
        BorderTriangle( index_t s, index_t f, index_t v0, index_t v1, index_t v2 )                
            : s_( s ), v0_( v0 ), v1_( v1 ), v2_( v2 ) {};

        bool operator<( const BorderTriangle& rhs ) const {
            if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) 
                return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ ) ;
            if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) 
                return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ ) ;
            if( s_ != rhs.s_ ) return s_ < rhs.s_ ;
            if( f_ != rhs.f_ ) return f_ != rhs.f_ ;            
            return v2_ < rhs.v2_ ;
        }
        bool same_edge( const BorderTriangle& rhs ) const {
            return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ ) &&
                std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ ) ;
        }

        /// Indices of the points in the model. Triangle has the Surface orientation
        /// The edge v0v1 is, in surface s_, on the border.            
        index_t v0_ ; 
        index_t v1_ ;
        index_t v2_ ;
        // Index of the model surface
        index_t s_ ;
        // Index of the facet in the surface
        index_t f_ ;
    } ;

    /*! 
     * @brief Get the BorderTriangle corresponding to the next edge on border 
     * in the considered Surface
     */
    index_t get_next_border_triangle(
        const BoundaryModel& M, 
        const std::vector< BorderTriangle >& BT,
        index_t from, 
        bool backward = false 
    ){
        const BorderTriangle& in = BT[from] ;
        const Surface& S = M.surface( in.s_ ) ;
        index_t NO_ID (-1) ;

        // Get the next edge on border in the Surface
        index_t f = in.f_ ;
        index_t f_v0 = S.facet_id_from_model(f, in.v0_ ) ;
        index_t f_v1 = S.facet_id_from_model(f, in.v1_ ) ;
        grgmesh_assert( f_v0 != NO_ID && f_v1 != NO_ID ) ;

        index_t next_f = NO_ID ;
        index_t next_f_v0 = NO_ID ;
        index_t next_f_v1 = NO_ID ;

        if( !backward ) S.next_on_border( f, f_v0, f_v1, next_f, next_f_v0, next_f_v1 ) ;
        else            S.next_on_border( f, f_v1, f_v0, next_f, next_f_v0, next_f_v1 ) ;

        // Find the BorderTriangle that is correspond to this
        // It must exist and there is only one
        BorderTriangle bait( in.s_, next_f, S.model_vertex_id(next_f, next_f_v0), S.model_vertex_id(next_f, next_f_v1), NO_ID ) ;
        index_t result = std::lower_bound( BT.begin(), BT.end(), bait )- BT.begin() ;
       
        grgmesh_assert( result < BT.size() ) ;
        return result ;
    }

    void visit_border_triangle_on_same_edge( 
        const std::vector< BorderTriangle >& border_triangles,
        index_t i, 
        index_t line_id, 
        std::vector< index_t >& border_line_ids )
    {        
        index_t j = i+1 ;        
        while( j < border_triangles.size() && 
              border_triangles[i].same_edge( border_triangles[j] ) ) 
        {    
            border_line_ids[j] = line_id ;
            j++ ;            
        }       
        signed_index_t k = i -1 ;
        while( k > -1 && border_triangles[i].same_edge( border_triangles[k] ))
        {                             
            border_line_ids[k] = line_id ;
            k-- ;
        }
    }

    
    void get_adjacent_surfaces( 
        const std::vector< BorderTriangle >& border_triangles,
        index_t i, 
        std::vector< index_t >& adjacent_surfaces )
    {
        adjacent_surfaces.resize(0) ;
        adjacent_surfaces.push_back( border_triangles[i].s_ ) ;
        
        index_t j = i+1 ;        
        while( j < border_triangles.size() && 
              border_triangles[i].same_edge( border_triangles[j] ) ) 
        {
            adjacent_surfaces.push_back( border_triangles[j].s_ ) ;           
            j++ ;            
        }
        
        signed_index_t k = i -1 ;
        while( k > -1 && border_triangles[i].same_edge( border_triangles[k] )) 
        {                       
            adjacent_surfaces.push_back( border_triangles[k].s_ ) ;           
            k-- ;
        }
        // In rare cases - the same surface can appear twice around a contact
        // Make unique and sort the adjacent regions
        std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() ) ;                    
        adjacent_surfaces.resize(
            std::unique( adjacent_surfaces.begin(), adjacent_surfaces.end() )-adjacent_surfaces.begin() ) ;
    }


    /*! Computes and creates the corners, contacts and regions of the 
     *  using geometrical tests
     */
    void BoundaryModelBuilderSurface::build_model() {
     
        /// 1. Make the storage of the model vertices unique 
        /// So now we can make index comparison for edges :)
        make_vertices_unique() ;

        /// Store the edges on a boundary from the input mesh
        std::vector< BorderTriangle > border_triangles ;
       
        /// 1. Get for all Surface the triangles on the boundary      
        for( index_t i = 0; i < model_.nb_surfaces(); ++i ) {
            const Surface& S = model_.surface( i ) ;
            for( index_t j= 0; j < S.nb_cells(); ++j ) {
                for( index_t v = 0; v < S.nb_vertices_in_facet(j); ++v ){
                    if( S.is_on_border( j, v ) ) {
                        border_triangles.push_back( BorderTriangle( i, j,
                            S.model_vertex_id( j, v ), 
                            S.model_vertex_id( j, S.next_in_facet(j,v) ),
                            S.model_vertex_id( j, S.prev_in_facet(j,v) ) 
                        ) ) ;
                    }
                }
            }                      
        }
        // Sort them so that the ones sharing the same edges are the following
        std::sort( border_triangles.begin(), border_triangles.end() ) ;
   
        /// The index of the Line in which is each one of these BorderTriangle
        std::vector< index_t > border_line_ids( border_triangles.size(), NO_ID ) ;       

        
        /// 2. Build the contact parts 
        // For each contact part built keep 

        // Info to build the regions
        std::vector< ContactSort > regions_info;

        {
            // Model indices of the vertices defining the the contact 
            std::vector< index_t > vertices ;
            
            // The goal is to put all the edges that are on the border in one of Line
            index_t cur_line_id = 0 ;
            for( index_t i = 0; i < border_triangles.size(); ++i ) {
                if( border_line_ids[i] != NO_ID )
                {                       
                    border_line_ids[i] = cur_line_id ;

                    // Begin the gathering of information to create a new Line                    
                    vertices.clear() ;
                                     
                    // For all colocalized edges - so all the edges that follows
                    // this one in the vector and have the same edge v0v1
                    
                    // Get the indices of the adjacent surfaces
                    std::vector< index_t > adjacent ;

                    get_adjacent_surfaces ( border_triangles, i, adjacent ) ;
                    visit_border_triangle_on_same_edge( border_triangles, i, cur_line_id, border_line_ids ) ;
                    
                    // Get info to sort Surfaces around the contact
                    regions_info.push_back( ContactSort( cur_line_id ) ) ;
                    index_t j = i+1 ;
                    while( border_triangles[i].same_edge( border_triangles[j] )) {                       
                        regions_info.back().add_surface( 
                            border_triangles[j].s_ ,
                            model_.vertex( border_triangles[j].v0_ ) , 
                            model_.vertex( border_triangles[j].v1_ ) , 
                            model_.vertex( border_triangles[j].v2_ ) 
                            ) ;
                        j++ ;
                    }   

                    vertices.push_back( border_triangles[i].v0_ ) ;
                    vertices.push_back( border_triangles[i].v1_ ) ;
            
                    // Build the contact propating forward on the border of the Surface
                    // While the adjacent surfaces are the smae the vertices of the following edges are added
                    bool same_surfaces = true ;
                    index_t next_i = get_next_border_triangle( model_, border_triangles, i ) ;
                    do {                         
                        if( border_line_ids[ next_i ] == NO_ID ){
                            std::vector< index_t > adjacent_next ;
                            get_adjacent_surfaces( border_triangles, next_i, adjacent_next ) ;

                            if( adjacent.size() == adjacent_next.size() &&
                                std::equal( adjacent.begin(), adjacent.end(), adjacent_next.begin() )
                            ){
                                // The surfaces in contact are still the same
                                // Assign the current_line_id to the BorderTriangle sharing the edge
                                visit_border_triangle_on_same_edge( 
                                    border_triangles, next_i, cur_line_id, border_line_ids ) ;
                                
                                // Fill the Line vertices - Je suis pas dut tou sure que ce soit ça
                                if( border_triangles[next_i].v0_ == vertices.back() ) 
                                    vertices.push_back( border_triangles[next_i].v1_ ) ;
                                else {
                                    grgmesh_assert( border_triangles[next_i].v1_ == vertices.back() ) ;
                                    vertices.push_back( border_triangles[next_i].v0_ );
                                }
                            }
                            else same_surfaces = false ;                    
                        }
                        else same_surfaces = false ;            
                        next_i = get_next_border_triangle( model_, border_triangles, next_i ) ;

                    } while( same_surfaces && next_i != i ) ;

                    if( next_i != i ) {
                        // The boundary is not a loop
                        // Propagate backward to reach the other extremity 
                        same_surfaces = true ;            
                        // Going backwards
                        index_t prev_i = get_next_border_triangle( model_, border_triangles, i, true ) ;
                        do { 
                            grgmesh_assert( prev_i != i ) ;
                            if( border_line_ids[ prev_i ] == NO_ID ){
                                std::vector< index_t > adjacent_prev;
                                get_adjacent_surfaces( border_triangles, prev_i, adjacent_prev ) ;

                                if( adjacent.size() == adjacent_prev.size() &&
                                    std::equal( adjacent.begin(), adjacent.end(), adjacent_prev.begin() )
                                    ){
                                        // The surfaces in contact are still the same
                                        // Assign the current_line_id to the BorderTriangle sharing the edge
                                        visit_border_triangle_on_same_edge( 
                                            border_triangles, prev_i, cur_line_id, border_line_ids ) ;

                                        // Fill the Line vertices
                                if( border_triangles[prev_i].v0_ == vertices.front() ) 
                                    vertices.push_back( border_triangles[prev_i].v1_ ) ;
                                else {
                                    grgmesh_assert( border_triangles[prev_i].v1_ == vertices.front() ) ;
                                    vertices.push_back( border_triangles[prev_i].v0_ ) ;
                                }
                                }
                                else same_surfaces = false ;                    
                            }                        
                            else same_surfaces = false ;                   
                            prev_i = get_next_border_triangle( model_, border_triangles, prev_i, true ) ;

                        } while( same_surfaces ) ;
                    }

                    grgmesh_assert( vertices.size() > 1 )
               
                    // Now we have all the points on the contact 
                    // The first and the last are corners if the contact is not a loop
                    
                    // Find or create the corners
                    index_t corner0 = find_or_create_corner( vertices.front() ) ;
                    index_t corner1 = find_or_create_corner( vertices.back() ) ;
                        
                    // what happens with closed lines ???? 

                    // Create the current Line - corners must exist
                    index_t created = create_line( vertices ) ; 
                    grgmesh_assert( created == cur_line_id ) ;
                    
                   
                    for( index_t j = 0; j < adjacent.size(); ++j ) {
                        add_element_in_boundary( BME::LINE, created, adjacent[j] ) ;
                    }
                    cur_line_id++ ;
                }
            }
        }

   
        /// 3. Build the regions 

        // Sort surfaces around the contacts
        for( index_t i = 0 ; i < regions_info.size(); ++ i) {
            regions_info[i].sort() ;
        }        
        
        if( model_.nb_surfaces() == 1 ) {
            // Check that this surface is closed and define an interior
            // and exterior (universe) regions
            grgmesh_assert_not_reached ;            
        }
        else {
            index_t cur_region_id  = 0 ;
            // For each side of each Surface store the region in it is
            std::vector< index_t > surf_2_region ( 2*model_.nb_surfaces(), NO_ID ) ;

            // Start with the first Surface on its + side            
            std::stack< std::pair< index_t, bool > > S ;
            S.push( std::pair< index_t, bool > ( 0, true ) ) ;

            while( !S.empty() ){
                std::pair< index_t, bool > cur = S.top() ;
                S.pop() ;

                // Already visited
                if( surf_2_region[ cur.second == true ? 2*cur.first : 2*cur.first+1 ] != NO_ID ) continue ;
                                                
                // Create a new region
                create_region() ;

                std::stack< std::pair< index_t, bool > > SR ;
                SR.push( cur ) ;

                while( !SR.empty() ) {
                    std::pair< index_t, bool > s = SR.top() ;
                    SR.pop() ;

                    index_t s_id = s.second == true ? 2*s.first : 2*s.first+1 ;
                    // Already in a region
                    if( surf_2_region[ s_id ] != NO_ID ) continue ;
                    
                    // Add the surface to the current region
                    add_element_boundary( BME::REGION, cur_region_id, s.first, s.second ) ;
                    surf_2_region[ s_id ] = cur_region_id ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp = !s.second == true ? 2*s.first : 2*s.first+1 ;
                    if( surf_2_region[ s_id_opp ] == -1 ) 
                        S.push( std::pair< index_t, bool >(s.first, !s.second ) );                    

                    // For each contact, push the next oriented surface in the region
                    const BoundaryModelElement& surface_part = model_.surface( s.first ) ;
                    for( index_t i = 0 ; i < surface_part.nb_boundaries(); ++i ){
                        index_t cp_id = surface_part.boundary( i ).id() ;
                        grgmesh_assert( cp_id < regions_info.size() ) ;
                        const std::pair< index_t, bool >& n = regions_info[ cp_id ].next( s ) ;
                        index_t n_id =  n.second == true ? 2*n.first : 2*n.first+1 ;
                        if( surf_2_region[ n_id ] == -1 )
                            SR.push( n ) ;
                    }                
                }
                cur_region_id ++ ;
            }
            
            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            // TO DO implement the code to take into account included region building
            grgmesh_assert( std::count( surf_2_region.begin(), surf_2_region.end(), -1 ) == 0 ) ;        
        }

        // We need to remove form the regions_ the one corresponding
        // to the universe_ that is the one with the biggest volume
        double max_volume = -1. ;
        index_t universe_id = NO_ID ;
        for( index_t i = 0; i < model_.nb_regions(); ++i ){
            double cur_volume = BoundaryModelElementMeasure::size( &model_.region(i) ) ;
            if( cur_volume > max_volume ) {
                max_volume = cur_volume ;
                universe_id = i ;
            }
        }
        grgmesh_assert( universe_id != NO_ID ) ;
        
        const BoundaryModelElement& cur_region = model_.region(universe_id) ;
        std::vector< std::pair< index_t, bool > > univ_boundaries( cur_region.nb_boundaries() ) ;
        for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ){
            univ_boundaries[i].first = cur_region.boundary( i ).id() ;
            univ_boundaries[i].second = cur_region.side( i ) ;
        }
        set_universe( univ_boundaries ) ;

        // Decrease by one the ids of the regions that are after the
        // one converted to the universe
        // Remove the region converted to universe from the regions
        remove_universe_from_regions( universe_id ) ;                   

        end_model() ;

    } 

    

} //namespace
