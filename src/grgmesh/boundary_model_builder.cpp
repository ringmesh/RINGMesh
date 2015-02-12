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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr 
 *     Antoine.Mazuyer@univ-lorraine.fr 
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


/*! \author Jeanne Pellerin */


#include <grgmesh/boundary_model_builder.h>

#include <geogram/basic/line_stream.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>
#include <stack>

namespace GRGMesh {
    double read_double( GEO::LineInput& in, index_t field ) {
        double result ;
        GEO::String::from_string( in.field( field ), result ) ;
        return result ;
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
        model_.regions_.resize( from.nb_regions(), BoundaryModelElement( &model_, BME::REGION ) ) ;
        model_.layers_.resize( from.nb_layers(), BoundaryModelElement( &model_, BME::LAYER ) ) ;
        model_.contacts_.resize( from.nb_contacts(), BoundaryModelElement( &model_, BME::CONTACT ) ) ;
        model_.interfaces_.resize( from.nb_interfaces(), BoundaryModelElement( &model_, BME::INTERFACE ) ) ;
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
    index_t BoundaryModelBuilder::create_element( BME::TYPE type ) {
        index_t id = model_.nb_elements( type ) ;
        grgmesh_assert( id != NO_ID ) ;
        switch( type ) {
        case BME::CORNER:             
            model_.corners_.push_back( Corner( &model_, id ) ) ;
            break ;

        case BME::LINE:            
            model_.lines_.push_back( Line( &model_, id ) ) ;
            break ;

        case BME::SURFACE:            
            model_.surfaces_.push_back( Surface( &model_, id ) ) ;
            break ;

        case BME::REGION:            
            model_.regions_.push_back( BME( &model_, BME::REGION, id ) ) ;
            break ;

        case BME::CONTACT:            
            model_.contacts_.push_back( BME( &model_,BME::CONTACT, id ) ) ;
            break ;

        case BME::INTERFACE:            
            model_.interfaces_.push_back( BME( &model_,BME::INTERFACE, id ) ) ;
            break ;

        case BME::LAYER:            
            model_.layers_.push_back( BME( &model_, BME::LAYER, id ) ) ;
            break ;            
        default:   
            grgmesh_assert_not_reached ;
        }
        return id ;
    }

    /*!
     * @brief To use with extreme caution -  Erase one element of the BoundaryModel 
     * @details TO USE ONLY AFTER having removed all references to this element,
     * AND having updated the indices of the elements of the same type
     * AND having updated all references to these elements in their boundaries, 
     * in_boundaries, parent or children.
     *
     * \todo Implement a generic function remove all references to an element (its index)
     * update the indices of elements of the same type (after it in the vector) and update
     * references to these elements in all other elements.
     */
    void BoundaryModelBuilder::erase_element( BME::TYPE type, index_t id ) {
        switch( type ) {
        case BME::CORNER : 
            model_.corners_.erase( model_.corners_.begin()+id ) ;
            break ;     

        case BME::LINE:
            model_.lines_.erase( model_.lines_.begin()+id ) ;
            break ;

        case BME::SURFACE:
            model_.surfaces_.erase( model_.surfaces_.begin()+id ) ;
            break ;    

        case BME::REGION:
            model_.regions_.erase( model_.regions_.begin()+id ) ;
            break ;

        case BME::CONTACT:
            model_.contacts_.erase( model_.contacts_.begin()+id ) ;
            break ;

        case BME::INTERFACE:
            model_.interfaces_.erase( model_.interfaces_.begin()+id ) ;
            break ;

        case BME::LAYER:
            model_.layers_.erase( model_.layers_.begin()+id ) ;
            break ;

        default:   
            grgmesh_assert_not_reached ;
        }
    }


    void BoundaryModelBuilder::resize_elements( BME::TYPE type, index_t nb ) 
    {
         switch( type ) {
        case BME::CORNER : 
            model_.corners_.resize( nb, Corner( &model_) ) ;
            break ;     

        case BME::LINE:
            model_.lines_.resize( nb, Line( &model_) ) ;
            break ;

        case BME::SURFACE:
            model_.surfaces_.resize( nb, Surface( &model_) ) ;
            break ;    

        case BME::REGION:
            model_.regions_.resize( nb, BME( &model_, BME::REGION ) ) ;
            break ;

        case BME::CONTACT:
            model_.contacts_.resize( nb, BME( &model_, BME::CONTACT ) ) ;
            break ;

        case BME::INTERFACE:
            model_.interfaces_.resize( nb, BME( &model_, BME::INTERFACE ) ) ;
            break ;

        case BME::LAYER:
            model_.layers_.resize( nb, BME( &model_, BME::LAYER ) ) ;
            break ;
        default:   
            break ;
         }      
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
       index_t id = create_element( BME::CORNER ) ; 
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
        index_t id = create_element( BME::LINE ) ;
        set_line( id, points ) ;
               
        // Find the indices of the corner at both extremities
        index_t c0 = find_or_create_corner( points.front() ) ;
        index_t c1 = find_or_create_corner( points.back() ) ;       
        add_element_boundary( BME::LINE, id, c0 ) ;
        if( c1 != c0 ) add_element_boundary( BME::LINE, id, c1 ) ;         

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
        return create_element( BME::SURFACE ) ; 
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
        
        index_t id = create_element( BME::CONTACT ) ;
        set_element_name( BME::CONTACT, id, name ) ;
        
        /*for( index_t i = 0; i < interfaces.size(); ++i ) {
            add_element_in_boundary( BME::CONTACT, id, interfaces[i] ) ;
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
        BME::GEOL_FEATURE type )
    {
        index_t id = create_element( BME::INTERFACE ) ;
        set_element_geol_feature( BME::INTERFACE, id, type ) ;
        set_element_name( BME::INTERFACE, id, name ) ;
        return id ;
    }

     /*
    * @brief Adds an empty region to the model
    *
    *  Used in Geomodeling to convert a surface to a model
    */
    index_t BoundaryModelBuilder::create_region() {
        return create_element( BME::REGION ) ;
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
        index_t id = create_element( BME::REGION ) ;
        set_element_name( BME::REGION, id, name ) ;
        for( index_t i = 0; i < boundaries.size(); ++i ) {            
            add_element_boundary( BME::REGION, id, boundaries[i].first, boundaries[i].second ) ;
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
        index_t id = create_element( BME::LAYER ) ;
        set_element_name( BME::LAYER, id, name ) ;
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
        model_.universe_.set_element_type( BME::REGION ) ;
        model_.universe_.set_model( &model_ ) ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            grgmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary( boundaries[i].first,
                boundaries[i].second ) ;
            // If this surface have no type, set it at BME::VOI
            model_.surfaces_[boundaries[i].first].set_geological_feature( BME::VOI ) ;
        }
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
        if( E.element_type() == BME::NO_TYPE ) return false ;
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
        if( BME::boundary_allowed( T ) && T != BME::SURFACE ) {
            // A closed surface - bubble might have no boundary
            // The others Line - and Region must have one
            if( E.nb_boundaries() == 0 ){
                return false ;
            }            
        }
        // In_boundary
        if( BME::in_boundary_allowed( T ) ) {
            if( E.nb_in_boundary() == 0 ){
                return false ;
            }
        }
        // Parent - High level elements are not mandatory
        // But if the model has elements of the parent type, the element must have a parent
        if( BME::parent_allowed( T ) ) {
            if( E.parent_id() == NO_ID && 
                model_.nb_elements( BME::parent_type(T) ) > 0 ){
                    return false ;
            }
        }
        // Children
        if( BME::child_allowed( T ) ) {
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
                fill_elements_boundaries( BME::LINE ) ;
            }
            if( model_.line(0).nb_in_boundary() == 0 ){
                fill_elements_in_boundaries( BME::LINE ) ;    
            }
            if( model_.line(0).parent_id() == NO_ID && model_.nb_contacts() > 0 ){
                fill_elements_parent( BME::LINE ) ;
            }
        }
        // Corners
        if( model_.nb_corners() > 0 && model_.corner(0).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( BME::CORNER ) ;
        }    
        // Surfaces - There MUST be at least one
        if( model_.surface(0).nb_boundaries() == 0 ) {            
            fill_elements_boundaries( BME::SURFACE ) ;
        }        
        if( model_.surface(0).nb_in_boundary() == 0 ) {            
            fill_elements_in_boundaries( BME::SURFACE ) ;    
        }
        if( model_.surface(0).parent_id() == NO_ID ) {            
            fill_elements_parent( BME::SURFACE ) ;    
        }
        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region(0).nb_boundaries() == 0 ) {
                fill_elements_boundaries( BME::REGION ) ;
            }
            if( model_.region(0).parent_id() == NO_ID && model_.nb_layers() > 0 ){
                fill_elements_parent( BME::REGION ) ;
            }
        }
        // Contacts
        if( model_.nb_contacts() > 0 &&
            model_.contact(0).nb_children() == 0 ) {
            fill_elements_children( BME::CONTACT ) ;
        }
        // Interfaces
        if( model_.nb_interfaces() > 0 &&
            model_.one_interface(0).nb_children() == 0 ) {
            fill_elements_children( BME::INTERFACE ) ;
        }
        // Layers
        if( model_.nb_layers() > 0 &&
            model_.layer(0).nb_children() == 0 ) {
                fill_elements_children( BME::LAYER ) ;
        }
        return true ;
    }

    /*!
     * @brief Fill the boundaries of all elements of the given type
     *
     * @details If the boundary elements do not have any in_boundary
     * information, nothing is done, and model construction will eventually fail.
     */
    void BoundaryModelBuilder::fill_elements_boundaries( BME::TYPE type ) 
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
    void BoundaryModelBuilder::fill_elements_in_boundaries( BME::TYPE type ) 
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
    void BoundaryModelBuilder::fill_elements_parent( BME::TYPE type ) 
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
    void BoundaryModelBuilder::fill_elements_children( BME::TYPE type ) 
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
     * @brief Last function to call when building a model
     *
     * @details check that the model is correct and has all required information
     * Calls the complete_element_connectivity function 
     * Fills nb_elements_per_type_ vector
     * 
     * @return False if the model is not valid and cannot be fixed
     * otherwise returns true.
     *
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


    /**
    * \brief Structure used to build Line by BoundaryModelBuilderGocad
    */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1):
        part_id_(part), corner_id_(corner), p0_(p0), p1_(p1) {};

        // Id of the Surface owning this Border
        index_t part_id_ ;
        // Id of p0 in the BoundaryModel corner vector
        index_t corner_id_ ;

        // Ids of the starting corner and second vertex on the border in the Surface
        // to which this Border belong
        index_t p0_ ;
        index_t p1_ ;
    } ;   


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
    void BoundaryModelBuilderGocad::load_ml_file( const std::string& ml_file_name )
    {
        GEO::LineInput in( ml_file_name);
        if(!in.OK()) {
            return ;
        }        
        
        // Clear the model. Not sure this is actually useful
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

        while(!in.eof() && in.get_line()) {
            in.get_fields();
            if( in.nb_fields() > 0 ) 
            {
                if( read_model ) 
                {
                    if( strncmp( in.field(0), "name:", 5 )== 0 )
                    {
                        set_model_name( &in.field(0)[5] ) ;
                    }
                    else if( in.field_matches( 0, "TSURF" ) ) 
                    {
                        // 1. Create Interface its name
                        create_interface( in.field(1) ) ;
                        nb_tsurf++ ;
                    }
                    else if( in.field_matches( 0, "TFACE" ) ) 
                    {
                        // 2. Create the Surface from the name of its parent Interface
                        // its geological feature
                        index_t id = in.field_as_uint(1) ;
                        std::string geol = in.field(2) ;
                        std::string interface_name = in.field(3) ;
                        // And its key facet that give the orientation of the surface part
                        in.get_line() ; in.get_fields() ;              
                        vec3 p0 ( read_double( in, 0 ), read_double( in, 1 ), read_double( in, 2 ) ) ;
                        in.get_line() ; in.get_fields() ;              
                        vec3 p1 ( read_double( in, 0 ), read_double( in, 1 ), read_double( in, 2 ) ) ;
                        in.get_line() ; in.get_fields() ;              
                        vec3 p2 ( read_double( in, 0 ), read_double( in, 1 ), read_double( in, 2 ) ) ;
                    
                        create_surface( interface_name, geol, p0, p1, p2 ) ; 
                        nb_tface++ ;
                    }
                    else if( in.field_matches( 0, "REGION" ) ) 
                    {
                        // 3. Read Region information and create them from their name,
                        // the surfaces on their boundary                    
                        index_t id = in.field_as_uint(1) ;
                        std::string name = in.field(2) ;
                        
                        std::vector< std::pair< index_t, bool > > region_boundaries ;
                        bool end_region = false ;
                        while( !end_region ) {
                            in.get_line() ; in.get_fields() ;              
                            for( index_t i = 0; i < 5; ++i ) {
                                int s = in.field_as_int(i) ;                                
                                if( s == 0 ) {
                                    end_region = true ;
                                    break ;
                                }
                                bool side = s > 0 ;
                                if( s > 0 ) s -= 1 ;
                                else s = -s-1 ;
                                region_boundaries.push_back(
                                        std::pair< index_t, bool >( s, side ) ) ;                                
                            }
                        }
                        // The Universe is not a regular region                       
                        if( name == "Universe" ) set_universe( region_boundaries ) ;
                        else create_region( name, region_boundaries ) ;
                    }
                    else if( in.field_matches( 0, "LAYER" ) ) 
                    {
                        // 4. Build the volumetric layers from their name and 
                        // the ids of the regions they contain
                        index_t layer_id = create_layer( in.field(2) ) ;
                        bool end_layer = false ;
                        while( !end_layer ) {
                            in.get_line() ; in.get_fields() ;              
                            for( index_t i = 0; i < 5; ++i ) {
                                index_t region_id = in.field_as_uint(i) ;
                                if( region_id == 0 ) {
                                    end_layer = true ;
                                    break ;
                                } else {
                                    region_id -= nb_tface+1 ; // Remove Universe region
                                    // Correction because ids begin at 1 in the file
                                    add_child( BME::LAYER, layer_id, region_id-1 ) ;
                                }
                            }
                        }
                    } 
                    else if( in.field_matches( 0, "END" ) ) 
                    {
                        // End of the high level information on the model
                        // Switch to reading the geometry of the model surfaces
                        read_model = false ;
                        continue ;
                    }
                } 
                else {
                    if( in.field_matches( 0, "GOCAD") ) 
                    {
                        // This is the beginning of a new TSurf = Interface
                        tsurf_count++ ;
                    }
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) 
                    {
                        if( in.field_matches( 1, "Elevation" ) ) z_sign = 1 ;
                        else if( in.field_matches( 1, "Depth" ) ) z_sign = -1 ;
                        else grgmesh_assert_not_reached ;
                    } 
                    else if( in.field_matches( 0, "END" ) ) 
                    {
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
                    }
                    else if( in.field_matches( 0, "TFACE" ) ) 
                    {
                        // Beginning of a new TFace - Surface
                        if( tface_vertex_start.size() > 0 ) {
                            // End the previous TFace - Surface  (copy from line 1180)
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
                    // 4. Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    else if( in.field_matches( 0, "VRTX" ) || in.field_matches( 0, "PVRTX" ) ) 
                    {
                        const vec3 p( read_double( in, 2 ), read_double( in, 3 ), z_sign*read_double( in, 4 )) ;
                        tsurf_vertex_ptr.push_back( add_vertex( p ) ) ;
                    } 
                    else if( in.field_matches( 0, "PATOM" ) | in.field_matches( 0, "ATOM" ) ) 
                    {
                        tsurf_vertex_ptr.push_back( tsurf_vertex_ptr[in.field_as_uint(2) - 1] ) ;
                    } 
                    else if( in.field_matches( 0, "TRGL" ) ) 
                    {
                        // Read ids of the vertices of each triangle in the TSurf                        
                        // and switch to ids in the TFace
                        tface_facets.push_back( (index_t) in.field_as_uint(1)-tface_vertex_start.back()-1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(2)-tface_vertex_start.back()-1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(3)-tface_vertex_start.back()-1 ) ;
                        tface_facets_ptr.push_back( tface_facets.size() ) ;
                    }
                    // 5. Build the corners from their position and the surface parts
                    //    containing them
                    else if( in.field_matches( 0, "BSTONE" ) ) 
                    {
                        index_t v_id = in.field_as_uint(1)-1 ;                        
                        if( find_corner( model_.vertex( tsurf_vertex_ptr[v_id] ) ) == NO_ID ) 
                            create_corner( tsurf_vertex_ptr[v_id] ) ;
                    }
                    /// 6. Read the Border information and store it
                    else if( in.field_matches( 0, "BORDER" ) ) 
                    {
                        index_t p1 = in.field_as_uint(2)-1 ;
                        index_t p2 = in.field_as_uint(3)-1 ;

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
        }
    
        make_vertices_unique() ;

        /// 7. Build the Lines
        { 
            std::vector< index_t > global_ids ;
            for( index_t i = 0; i < borders_to_build.size(); ++i ) {
                const Border& b = borders_to_build[i] ;

                // 1- Build the boundary : construct the vector
                // of vertices on the border
                const Surface& S = model_.surface( b.part_id_) ;

                index_t end_corner_id = determine_line_vertices( 
                    S, b.p0_, b.p1_, global_ids ) ;           

                // 2 - Check if this border already exists
                index_t line_id = find_or_create_line( global_ids ) ;

                // Add the surface in which this line is
                add_element_in_boundary( BME::LINE, line_id, b.part_id_ ) ;
            }
        }
        
        /// 8. Build the Contacts
        build_contacts() ;
            
        // Modify in the Region the side of the Surface for which the key facet 
        // orientation was not the same than their facet orientations        
        for( index_t i = 0; i < change_key_facet.size(); i++ ) {   
            const Surface& S = model_.surface( change_key_facet[i] ) ;
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {                
                BoundaryModelElement& R = element( BME::REGION, S.in_boundary_id(j) ) ;              
                for( index_t b = 0; b < R.nb_boundaries(); ++b ){
                    if( R.boundary_id(b) == change_key_facet[i] ){
                        bool old_side = R.side(b) ;
                        R.set_boundary( b, R.boundary_id(b), !old_side ) ;
                    }
                }
            }
        }        
     
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
        const KeyFacet& key_facet = key_facets_[surface_id] ;
      
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
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces     
     */
    void BoundaryModelBuilderGocad::build_contacts()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< index_t > interfaces ;
            for( index_t j = 0; j < model_.line(i).nb_in_boundary(); ++j ) {
                index_t sp_id = model_.line(i).in_boundary_id( j ) ;
                const BoundaryModelElement& p = model_.surface(sp_id).parent() ;
                interfaces.insert( p.id() ) ;
            }
            std::vector< index_t > toto( interfaces.begin(), interfaces.end() ) ;
            index_t contact_id = find_or_create_contact( toto ) ;
            add_child( BME::CONTACT, contact_id, i ) ;
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
        const vec3& p0, 
        const vec3& p1,
        const vec3& p2 )
    {
        index_t parent = find_interface( interface_name ) ;
        if( interface_name != "" ) grgmesh_assert( parent != NO_ID ) ;

        index_t id = create_element( BME::SURFACE ) ;
        set_parent( BME::SURFACE, id, parent ) ;
        key_facets_.push_back( KeyFacet( p0, p1, p2 ) ) ;
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


    bool BoundaryModelBuilderBM::load_file( const std::string& bm_file_name ) {
        GEO::LineInput in( bm_file_name);
        if(!in.OK()) {
            return false ;
        }        
        while(!in.eof() && in.get_line()) {
            in.get_fields();
            if( in.nb_fields() > 0 ) {
                // Name of the model
                if( in.field_matches( 0, "NAME" ) ) {
                    if( in.nb_fields() > 1 ) set_model_name( in.field(1) ) ;
                }
                // Number of elements of a given type
                else if( match_nb_elements( in.field( 0 ) ) != BME::NO_TYPE ) 
                {
                    // Allocate the space
                    if( in.nb_fields() > 1 ) 
                        resize_elements( match_nb_elements(in.field(0)), in.field_as_uint(1) ) ;
                }
                // High-level elements
                else if( match_high_level_type( in.field(0) ) )
                {
                    // Read this element 
                    // First line id - name - geol_feature
                    if( in.nb_fields() < 4 ) {
                        std::cout << "I/O Error File line " << in.line_number() << std::endl;
                        return false;
                    }
                    BME::TYPE t = match_type( in.field(0) ) ;
                    index_t id = in.field_as_uint(1) ;
                    set_element_index( t, id ) ;
                    set_element_name( t, id, in.field(2) ) ;
                    set_element_geol_feature( t, id, BME::determine_geological_type( in.field(3) ) ) ;
                    
                    // Second line - indices of its children
                    in.get_line() ; in.get_fields() ;
                    for(index_t c = 0; c < in.nb_fields(); c++) {
                        add_child( t, id, in.field_as_uint(c) ) ;
                    }                    
                }
                // Regions
                else if( match_type( in.field(0) ) == BME::REGION )
                {
                    // First line id - name
                    if( in.nb_fields() < 3 ) {
                        std::cout << "I/O Error File line " << in.line_number() << std::endl;
                        return false;
                    }
                    index_t id = in.field_as_uint(1) ;
                    set_element_index( BME::REGION, id ) ;
                    set_element_name( BME::REGION, id, in.field(2) ) ;

                    // Second line - signed indices of boundaries 
                    in.get_line() ;  in.get_fields() ;
                    for(index_t c = 0; c < in.nb_fields(); c++) {
                        bool side = false ;
                        if( strncmp( in.field(c), "+", 1 ) == 0 ) side = true ;
                        else {
                            // to remove after debug
                            grgmesh_assert( strncmp( in.field(c), "-", 1 ) == 0 ) ;
                        }
                       index_t s ;
                       GEO::String::from_string( &in.field(c)[1], s ) ;
                        
                        
                        add_element_boundary( BME::REGION, id, s, side ) ;
                    }                    
                }
                // Universe
                else if( in.field_matches(0, "UNIVERSE") )
                {
                    std::vector< std::pair< index_t, bool > > b_universe ;
                    // Second line - signed indices of boundaries 
                    in.get_line() ;  in.get_fields() ;
                    for(index_t c = 0; c < in.nb_fields(); c++) {
                        bool side = false ;
                        if( strncmp( in.field(c), "+", 1 ) == 0 ) side = true ;  
                        index_t s ;
                        GEO::String::from_string( &in.field(c)[1], s ) ;
                        
                        b_universe.push_back( std::pair< index_t, bool >( s, side ) ) ;
                    }
                    set_universe( b_universe ) ;
                }

                // Model vertices
                else if( in.field_matches( 0, "MODEL_VERTICES" ) ) 
                {
                    index_t nb_vertices = in.field_as_uint(1) ;
                    reserve_vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; ++i ){
                        in.get_line() ; in.get_fields() ;              
                        add_vertex( vec3( 
                            read_double( in, 0 ), read_double( in, 1 ), read_double( in, 2 ) ) ) ;
                    }
                }
                // Corners
                else if( match_type( in.field(0) ) == BME::CORNER )
                {
                    // One line id - vertex id
                    if( in.nb_fields() < 3 ) {
                        std::cout << "I/O Error File line " << in.line_number() << std::endl;
                        return false;
                    }
                    index_t id = in.field_as_uint(1) ;                    
                    set_element_index( BME::CORNER, id ) ;
                    set_element_vertex( BME::CORNER, id, 0, in.field_as_uint(2) ) ;
                }
                // Lines
                else if( match_type( in.field(0) ) == BME::LINE )
                {
                    index_t id = in.field_as_uint(1) ;
                    Line& L = dynamic_cast< Line& >( element( BME::LINE, id ) ) ;
                    L.set_id( id ) ;
                    
                    // Following information - vertices of the lines
                    in.get_line() ;  in.get_fields() ;
                    grgmesh_assert( in.field_matches( 0,"LINE_VERTICES" ) ) ;
                    index_t nb_vertices = in.field_as_uint(1) ;

                    // Read the vertices indices
                    in.get_line() ;  in.get_fields() ;                    
                    std::vector< index_t > vertices( nb_vertices ) ;
                    int count = 0 ;
                    for( index_t i = 0; i < nb_vertices; i++ ){
                        vertices[i] = in.field_as_uint( count ) ;
                        count++ ;
                        if( count == 20 && i+1 < nb_vertices ) {
                            count = 0 ;
                            in.get_line() ;  in.get_fields() ;
                        }
                    }
                    L.set_vertices( vertices ) ;
                    // Set the corners
                    index_t c0 = find_corner( vertices.front() ) ;
                    index_t c1 = find_corner( vertices.back() ) ;       
                    add_element_boundary( BME::LINE, id, c0 ) ;
                    if( c1 != c0 ) add_element_boundary( BME::LINE, id, c1 ) ; 
                    
                    // Finally we have the in_boundary information
                    in.get_line() ; in.get_fields() ;
                    grgmesh_assert( in.field_matches( 0,"IN_BOUNDARY" ) ) ;
                    for( index_t b = 1 ; b < in.nb_fields(); b++ ) {
                        L.add_in_boundary( in.field_as_uint( b ) ) ; 
                    }   
                }
                // Surfaces
                else if( match_type( in.field(0) ) == BME::SURFACE )
                {
                    index_t id = in.field_as_uint(1) ;
                    Surface& S = dynamic_cast< Surface& >( element( BME::SURFACE, id ) ) ;
                    S.set_id( id ) ;

                    // Read the surface facets
                    in.get_line() ; in.get_fields() ;
                    grgmesh_assert( in.field_matches( 0,"SURFACE_CORNERS" ) ) ;
                    index_t nb_corners = in.field_as_uint(1) ;

                    in.get_line() ; in.get_fields() ;
                    grgmesh_assert( in.field_matches( 0,"SURFACE_FACETS" ) ) ;
                    index_t nb_facets = in.field_as_uint(1) ;

                    std::vector< index_t > corners( nb_corners ) ;
                    std::vector< index_t > facet_ptr( nb_facets+1, 0 ) ;
                    index_t count_facets = 0 ;
                    for( index_t f = 0 ; f < nb_facets ; f++ ) {
                        in.get_line() ; in.get_fields() ;
                        
                        index_t nb_v = in.field_as_uint(0) ;
                        
                        for( index_t v = 0; v < nb_v; ++v ) {
                            corners[ count_facets+v ] = in.field_as_uint(v+1) ; 
                        }
                        count_facets += nb_v ;
                        facet_ptr[ f+1 ] = count_facets ;
                    }   

                    S.set_geometry( corners, facet_ptr ) ;
                    set_surface_adjacencies( id ) ;
                }
            }
        }
        grgmesh_assert( end_model() ) ;       
        return true ;
    }  

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_nb_elements( const char* s )
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3) != 0 ) return BME::NO_TYPE ;
        else {
            for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
                BME::TYPE type = (BME::TYPE) i ; 
                if( strstr( s, BME::type_name( type ).data() ) != NULL ) {
                    return type ;
                }
            }
            return BME::NO_TYPE ;
        }
    }

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_type( const char* s ) {
       for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
            BME::TYPE type = (BME::TYPE) i ; 
            if( strcmp( s, BME::type_name( type ).data() ) == 0 ) {
                return type ;
            }
        }
        return BME::NO_TYPE ;
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
            : s_( s ), f_( f ), v0_( v0 ), v1_( v1 ), v2_( v2 ) {};

        bool operator<( const BorderTriangle& rhs ) const {
            if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) 
                return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ ) ;
            if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) 
                return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ ) ;
            if( s_ != rhs.s_ ) return s_ < rhs.s_ ;
            if( f_ != rhs.f_ ) return f_ < rhs.f_ ;            
            return rhs.v2_ == index_t(-1) ? false : v2_ < rhs.v2_ ;
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
        BorderTriangle bait( in.s_, next_f, S.model_vertex_id(next_f, next_f_v0),
            S.model_vertex_id(next_f, next_f_v1), NO_ID ) ;
        // lower_bound returns an iterator pointing to the first element in the range [first,last) 
        // which does not compare less than the given val.
        // See operator< on BorderTriangle
        index_t result = std::lower_bound( BT.begin(), BT.end(), bait )-BT.begin() ;
       
        grgmesh_assert( result < BT.size() ) ;
        return result ;
    }

     /*! 
     * @brief Mark as visited all BorderTriangle which first edge is the same than
     * the first edge of i.
     *
     * @param[in] border_triangles Information on triangles MUST be sorted so that 
     *            BorderTriangle having the same boundary edge are adjacent
     *
     */
    void visit_border_triangle_on_same_edge( 
        const std::vector< BorderTriangle >& border_triangles,
        index_t i, 
        std::vector< bool >& visited )
    {        
        index_t j = i ;        
        while( j < border_triangles.size() && 
              border_triangles[i].same_edge( border_triangles[j] ) ) 
        {    
            visited[j] = true ;
            j++ ;            
        }       
        signed_index_t k = i -1 ;
        while( k > -1 && border_triangles[i].same_edge( border_triangles[k] ))
        {                             
            visited[k] = true ;
            k-- ;
        }
    }

    
    /*! 
     * @brief Get the indices of the Surface adjacent to the first edge of a BorderTriangle
     * 
     * @param[in] border_triangles Information on triangles MUST be sorted so that 
     *          BorderTriangle having the same boundary edge are adjacent
     * @param[in] i Index of the BorderTriangle
     * @param[out] adjacent_surfaces Indices of the Surface stored by the BorderTriangle sharing 
     *             the first edge of i
     */
    void get_adjacent_surfaces( 
        const std::vector< BorderTriangle >& border_triangles,
        index_t i, 
        std::vector< index_t >& adjacent_surfaces )
    {
        adjacent_surfaces.resize(0) ;
        
        index_t j = i ;        
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


    /*! 
     * @brief From a BoundaryModel in which only Surface are defined, create 
     * corners, contacts and regions.
     *
     */
    void BoundaryModelBuilderSurface::build_model() {

        grgmesh_assert( model_.nb_surfaces() > 0 ) ;
     
        /// 1. Make the storage of the model vertices unique 
        /// So now we can make index comparison to find colocated edges
        make_vertices_unique() ;
       
        /// 2.1 Get for all Surface, the triangles that have an edge
        /// on the boundary.
        std::vector< BorderTriangle > border_triangles ;
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
        /// 2.2 Sort these triangles so that triangles sharing the same edge follow one another
        std::sort( border_triangles.begin(), border_triangles.end() ) ;
           
        /// 3. Build the Lines and gather information to build the regions
        std::vector< SortTriangleAroundEdge > regions_info;
        // The goal is to visit all BorderTriangle and propagate to get each Line vertices
        std::vector< bool > visited ( border_triangles.size(), false ) ;                   
        for( index_t i = 0; i < border_triangles.size(); ++i ) {
            if( !visited[i] ) {                       
                // This is a new Line
                std::vector< index_t > vertices ;
                // Get the indices of the Surfaces around this Line     
                std::vector< index_t > adjacent ;
                get_adjacent_surfaces ( border_triangles, i, adjacent ) ;
                // Mark as visited the BorderTriangle around the same first edge
                visit_border_triangle_on_same_edge( border_triangles, i, visited ) ;
                    
                // Gather information to sort triangles around the contact
                regions_info.push_back( SortTriangleAroundEdge() ) ;
                index_t j = i ;
                while( j < border_triangles.size() &&
                        border_triangles[i].same_edge( border_triangles[j] )) 
                {                       
                    regions_info.back().add_triangle( 
                        border_triangles[j].s_ ,
                        model_.vertex( border_triangles[j].v0_ ) , 
                        model_.vertex( border_triangles[j].v1_ ) , 
                        model_.vertex( border_triangles[j].v2_ ) 
                        ) ;
                    j++ ;
                }   

                // Add vertices to the Line
                vertices.push_back( border_triangles[i].v0_ ) ;
                vertices.push_back( border_triangles[i].v1_ ) ;
            
                // Build the contact propating forward on the border of the Surface
                // While the adjacent surfaces are the same the vertices the next edge on the 
                // boundary of the Surface are added
                bool same_surfaces = true ;
                index_t next_i = get_next_border_triangle( model_, border_triangles, i ) ;
                do {                         
                    grgmesh_assert( next_i != NO_ID ) ;
                    if( !visited[ next_i ] ){
                        std::vector< index_t > adjacent_next ;
                        get_adjacent_surfaces( border_triangles, next_i, adjacent_next ) ;

                        if( adjacent.size() == adjacent_next.size() &&
                            std::equal( adjacent.begin(), adjacent.end(), adjacent_next.begin() )
                        ){                            
                            visit_border_triangle_on_same_edge( 
                                border_triangles, next_i, visited ) ;
                            // Add the next vertex    
                            if( border_triangles[next_i].v0_ == vertices.back() ) 
                                vertices.push_back( border_triangles[next_i].v1_ ) ;
                            else {
                                grgmesh_assert( border_triangles[next_i].v1_ == vertices.back() ) ;
                                vertices.push_back( border_triangles[next_i].v0_ );
                            }
                        }  else same_surfaces = false ;                    
                    } else same_surfaces = false ;            
                    next_i = get_next_border_triangle( model_, border_triangles, next_i ) ;
                } while( same_surfaces && next_i != i ) ;

                if( next_i != i ) {                    
                    // Propagate backward to reach the other extremity 
                    same_surfaces = true ;
                    index_t prev_i = get_next_border_triangle( model_, border_triangles, i, true ) ;
                    do { 
                        grgmesh_assert( prev_i != NO_ID && prev_i != i ) ;
                        if( !visited[ prev_i ] ){
                            std::vector< index_t > adjacent_prev ;
                            get_adjacent_surfaces( border_triangles, prev_i, adjacent_prev ) ;

                            if( adjacent.size() == adjacent_prev.size() &&
                                std::equal( adjacent.begin(), adjacent.end(), adjacent_prev.begin() )
                                )
                            {
                                visit_border_triangle_on_same_edge( 
                                    border_triangles, prev_i, visited ) ;
                                // Fill the Line vertices
                                if( border_triangles[prev_i].v0_ == vertices.front() ) 
                                    vertices.insert( vertices.begin(), border_triangles[prev_i].v1_ ) ;
                                else {
                                    grgmesh_assert( border_triangles[prev_i].v1_ == vertices.front() ) ;
                                    vertices.insert( vertices.begin(), border_triangles[prev_i].v0_ ) ;
                                }
                            }  else same_surfaces = false ;                    
                        } else same_surfaces = false ;                   
                        prev_i = get_next_border_triangle( model_, border_triangles, prev_i, true ) ;
                    } while( same_surfaces ) ;
                }

                grgmesh_assert( vertices.size() > 1 )
               
                // At last create the Line
                index_t created = create_line( vertices ) ; 
                for( index_t j = 0; j < adjacent.size(); ++j ) {
                    add_element_in_boundary( BME::LINE, created, adjacent[j] ) ;
                }
            }
        }     
   
        /// 4. Build the regions 
        
        // Complete boundary information for surfaces 
        // We need it to compute volumetric regions
        fill_elements_boundaries( BME::SURFACE ) ;
        
        /// 4.1 Sort surfaces around the contacts
        for( index_t i = 0 ; i < regions_info.size(); ++ i) {
            regions_info[i].sort() ;
        }        
        
        if( model_.nb_surfaces() == 1 ) {
            /// \todo Build a Region when a BoundaryModel has only one Surface 
            // Check that this surface is closed and define an interior
            // and exterior (universe) regions
            grgmesh_assert_not_reached ;            
        }
        else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2*model_.nb_surfaces(), NO_ID ) ;

            // Start with the first Surface on its + side
            std::stack< std::pair< index_t, bool > > S ;
            S.push( std::pair< index_t, bool > ( 0, true ) ) ;

            while( !S.empty() ){
                std::pair< index_t, bool > cur = S.top() ;
                S.pop() ;

                // This side is already assigned
                if( surf_2_region[ cur.second == true ? 2*cur.first : 2*cur.first+1 ] != NO_ID ) continue ;
                                                
                // Create a new region
                index_t cur_region_id = create_region() ;

                std::stack< std::pair< index_t, bool > > SR ;
                SR.push( cur ) ;
                while( !SR.empty() ) {
                    std::pair< index_t, bool > s = SR.top() ;
                    SR.pop() ;
                    index_t s_id = s.second == true ? 2*s.first : 2*s.first+1 ;
                    
                    // This side is already assigned
                    if( surf_2_region[ s_id ] != NO_ID ) continue ;
                    
                    // Add the surface to the current region
                    add_element_boundary( BME::REGION, cur_region_id, s.first, s.second ) ;
                    surf_2_region[ s_id ] = cur_region_id ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp = !s.second == true ? 2*s.first : 2*s.first+1 ;
                    if( surf_2_region[ s_id_opp ] == NO_ID ) 
                        S.push( std::pair< index_t, bool >(s.first, !s.second) );                    

                    // For each contact, push the next oriented surface that is in the same region
                    const BoundaryModelElement& surface = model_.surface( s.first ) ;
                    for( index_t i = 0 ; i < surface.nb_boundaries(); ++i ){
                        const std::pair< index_t, bool >& n =
                            regions_info[ surface.boundary_id(i) ].next( s ) ;
                        index_t n_id =  n.second == true ? 2*n.first : 2*n.first+1 ;
                        
                        if( surf_2_region[ n_id ] == NO_ID ) SR.push( n ) ;
                    }                
                }
            }
            
            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            /// \todo Implement the code to take into regions included in others (bubbles)
            grgmesh_assert( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID ) == 0 ) ;        
        }

        // We need to remove from the regions_ the one corresponding
        // to the universe_, the one with the biggest volume
        double max_volume = -1. ;
        index_t universe_id = NO_ID ;
        for( index_t i = 0; i < model_.nb_regions(); ++i ){
            double cur_volume = BoundaryModelElementMeasure::size( &model_.region(i) ) ;
            if( cur_volume > max_volume ) {
                max_volume = cur_volume ;
                universe_id = i ;
            }
        }
        
        const BoundaryModelElement& cur_region = model_.region(universe_id) ;
        std::vector< std::pair< index_t, bool > > univ_boundaries( cur_region.nb_boundaries() ) ;
        for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ){
            univ_boundaries[i].first = cur_region.boundary( i ).id() ;
            univ_boundaries[i].second = cur_region.side( i ) ;
        }
        set_universe( univ_boundaries ) ;

        // Decrease by one the ids of the regions that are after the
        // one converted to the universe
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            index_t cur_id = model_.region(i).id() ;            
            if( i > universe_id ) element( BME::REGION, i ).set_id( cur_id-1 ) ;
        }
        // Remove the region converted to universe from the regions
        erase_element( BME::REGION, universe_id ) ;
        // We are not in trouble since the boundaries of surface are not yet set
        // And we have no layer in the model

        grgmesh_assert( end_model() ) ;
    } 

    

} //namespace

