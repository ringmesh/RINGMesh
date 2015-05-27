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

#include <ringmesh/boundary_model_builder.h>
#include <ringmesh/utils.h>

#include <geogram/basic/line_stream.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_preprocessing.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>
#include <stack>

namespace {
    using namespace RINGMesh;

    double read_double( GEO::LineInput& in, index_t field )
    {
        double result ;
        std::istringstream iss( in.field( field ) ) ;
        iss >> result >> std::ws ;
        return result ;
    }

    /*!
     * @brief Fill the geological info from parent or children
     *
     * @details Set the geological feature of a BME to the one of its parent
     * if it has a parent that has a geol. feature,
     * or to the one of its first child if it has one with a geol. feature.
     *
     * @note The geol feature of \b E whatever its initial value
     */
    void fill_element_geological_feature( BoundaryModelElement& E )
    {
        if( E.has_parent() && E.parent().has_geological_feature() ) {
            E.set_geological_feature( E.parent().geological_feature() ) ;
        } else if( E.nb_children() > 0 && E.child( 0 ).has_geological_feature() ) {
            E.set_geological_feature( E.child( 0 ).geological_feature() ) ;
        }

        // Paranoia : we should perhaps verify that all children
        // have the same geological features
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
    *
    * \brief Tests whether a facet is degenerate.
    * \param[in] M the mesh that the facet belongs to
    * \param[in] f the index of the facet in \p M
    * \return true if facet \p f has duplicated vertices,
    *  false otherwise
    */
    bool facet_is_degenerate( 
        const GEO::Mesh& M,
        index_t f,
        GEO::vector<index_t>& colocated_vertices )
    {
        index_t nb_vertices = M.facets.nb_vertices( f );
        if( nb_vertices != 3 ) {
            index_t* vertices = (index_t*)alloca( nb_vertices*sizeof( index_t ) );
            for( index_t lv = 0; lv<nb_vertices; ++lv ) {
                vertices[ lv ] = colocated_vertices[ M.facets.vertex( f, lv ) ];
            }
            std::sort( vertices, vertices + nb_vertices );
            return std::unique(
                vertices, vertices + nb_vertices
                ) != vertices + nb_vertices;
        }
        index_t c1 = M.facets.corners_begin( f );
        index_t c2 = c1 + 1;
        index_t c3 = c2 + 1;
        index_t v1 = colocated_vertices[ M.facet_corners.vertex( c1 ) ];
        index_t v2 = colocated_vertices[ M.facet_corners.vertex( c2 ) ];
        index_t v3 = colocated_vertices[ M.facet_corners.vertex( c3 ) ];
        return v1 == v2 || v2 == v3 || v3 == v1;
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
     */
    void mesh_detect_degenerate_facets(
        const GEO::Mesh& M, 
        GEO::vector<index_t>& f_is_degenerate,
        GEO::vector<index_t>& colocated_vertices
        )
    {
        f_is_degenerate.resize( M.facets.nb() );
        for( index_t f = 0; f<M.facets.nb(); ++f ) {
            f_is_degenerate[ f ] = facet_is_degenerate( M, f, colocated_vertices );
        }
    }

    /*! 
     * @brief Detect and remove degenerated facets in a Mesh
     */
    index_t detect_degenerate_facets( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_facets( M, degenerate, colocated ) ;
        return std::count( degenerate.begin(), degenerate.end(), 1 ) ;
    }

    bool edge_is_degenerate(
        const GEO::Mesh& M,
        index_t e,
        GEO::vector<index_t>& colocated_vertices )
    {
        index_t v1 = colocated_vertices[ M.edges.vertex( e, 0 ) ];
        index_t v2 = colocated_vertices[ M.edges.vertex( e, 1 ) ];
        return v1 == v2 ;
    }

    void mesh_detect_degenerate_edges(
        const GEO::Mesh& M,
        GEO::vector<index_t>& e_is_degenerate,
        GEO::vector<index_t>& colocated_vertices
        )
    {
        e_is_degenerate.resize( M.edges.nb() );
        for( index_t e = 0; e<M.edges.nb(); ++e ) {
            e_is_degenerate[ e ] = edge_is_degenerate( M, e, colocated_vertices );
        }
    }

    /*!
    * @brief Detect and remove degenerated edges in a Mesh
    */
    index_t repair_line_mesh( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_edges( M, degenerate, colocated ) ;
        index_t nb = std::count( degenerate.begin(), degenerate.end(), 1 ) ;
        M.edges.delete_elements( degenerate ) ;
        return nb ;
    }
  


}

namespace RINGMesh {
  

    /*!
     * @brief Update the indices stored by each element of the model \
     * according to its actual position in the corresponding vector in the model
     */
    void BoundaryModelBuilder::update_all_ids()
    {
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            model_.corners_[co]->set_id( co ) ;
        }
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            model_.lines_[cp]->set_id( cp ) ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp]->set_id( sp ) ;
        }
        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            model_.contacts_[c]->set_id( c ) ;
        }
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            model_.interfaces_[s]->set_id( s ) ;
        }
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            model_.regions_[r]->set_id( r ) ;
        }
        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            model_.layers_[l]->set_id( l ) ;
        }
    }

    /*!
     * @brief Creates a element of the given type and add it to the correct vector
     * The BoundaryModelElement is created from its type and its index
     *
     * @param[in] type Type of the element to create
     * @return The index of the created element
     */
    BME::bme_t BoundaryModelBuilder::create_element( BME::TYPE type )
    {
        index_t id = model_.nb_elements( type ) ;
        ringmesh_assert( id != NO_ID ) ;
        switch( type ) {
            case BME::CORNER:
                model_.corners_.push_back( new Corner( &model_, id ) ) ;
                break ;

            case BME::LINE:
                model_.lines_.push_back( new Line( &model_, id ) ) ;
                break ;

            case BME::SURFACE:
                model_.surfaces_.push_back( new Surface( &model_, id ) ) ;
                break ;

            case BME::REGION:
                model_.regions_.push_back( new BME( &model_, BME::REGION, id ) ) ;
                break ;

            case BME::CONTACT:
                model_.contacts_.push_back( new BME( &model_, BME::CONTACT, id ) ) ;
                break ;

            case BME::INTERFACE:
                model_.interfaces_.push_back(
                    new BME( &model_, BME::INTERFACE, id ) ) ;
                break ;

            case BME::LAYER:
                model_.layers_.push_back( new BME( &model_, BME::LAYER, id ) ) ;
                break ;
            default:
                ringmesh_assert_not_reached;
                break ;
            }
        return BME::bme_t( type, id ) ;
    }

    /*!
     * @brief Use with EXTREME caution -  Erase one element of the BoundaryModel
     * @details TO USE ONLY AFTER having removed all references to this element,
     * AND having updated the indices of the elements of the same type
     * AND having updated all references to these elements in their boundaries,
     * in_boundaries, parent or children.
     *
     * \todo Implement a generic function remove all references to an element (its index)
     * update the indices of elements of the same type (after it in the vector) and update
     * references to these elements in all other elements.
     */
    void BoundaryModelBuilder::erase_element( const BME::bme_t& t )
    {
        switch( t.type ) {
            case BME::CORNER:
                delete model_.corners_[t.index] ;
                model_.corners_.erase( model_.corners_.begin() + t.index ) ;
                break ;

            case BME::LINE:
                delete model_.lines_[t.index] ;
                model_.lines_.erase( model_.lines_.begin() + t.index ) ;
                break ;

            case BME::SURFACE:
                delete model_.surfaces_[t.index] ;
                model_.surfaces_.erase( model_.surfaces_.begin() + t.index ) ;
                break ;

            case BME::REGION:
                delete model_.regions_[t.index] ;
                model_.regions_.erase( model_.regions_.begin() + t.index ) ;
                break ;

            case BME::CONTACT:
                delete model_.contacts_[t.index] ;
                model_.contacts_.erase( model_.contacts_.begin() + t.index ) ;
                break ;

            case BME::INTERFACE:
                delete model_.interfaces_[t.index] ;
                model_.interfaces_.erase( model_.interfaces_.begin() + t.index ) ;
                break ;

            case BME::LAYER:
                delete model_.layers_[t.index] ;
                model_.layers_.erase( model_.layers_.begin() + t.index ) ;
                break ;

            default:
                ringmesh_assert_not_reached;
                break ;
            }
        }

    void BoundaryModelBuilder::resize_elements( BME::TYPE type, index_t nb )
    {
        switch( type ) {
            case BME::CORNER:
                model_.corners_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.corners_[i] ) {
                        model_.corners_[i] = new Corner( &model_ ) ;
                    }
                }
                break ;

            case BME::LINE:
                model_.lines_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.lines_[i] ) {
                        model_.lines_[i] = new Line( &model_ ) ;
                    }
                }
                break ;

            case BME::SURFACE:
                model_.surfaces_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.surfaces_[i] ) {
                        model_.surfaces_[i] = new Surface( &model_ ) ;
                    }
                }
                break ;

            case BME::REGION:
                model_.regions_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.regions_[i] ) {
                        model_.regions_[i] = new BME( &model_, BME::REGION ) ;
                    }
                }
                break ;

            case BME::CONTACT:
                model_.contacts_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.contacts_[i] ) {
                        model_.contacts_[i] = new BME( &model_, BME::CONTACT ) ;
                    }
                }
                break ;

            case BME::INTERFACE:
                model_.interfaces_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.interfaces_[i] ) {
                        model_.interfaces_[i] = new BME( &model_, BME::INTERFACE ) ;
                    }
                }
                break ;

            case BME::LAYER:
                model_.layers_.resize( nb, nil ) ;
                for( index_t i = 0; i < nb; i++ ) {
                    if( !model_.layers_[i] ) {
                        model_.layers_[i] = new BME( &model_, BME::LAYER ) ;
                    }
                }
                break ;
            default:
                break ;
        }
    }

    /*!
     * @brief Get the index of the Corner for a given point
     * @param[in] point Geometric location to look for 
     * @return NO_ID or the index of the Corner
     */
    BME::bme_t BoundaryModelBuilder::find_corner( const vec3& point ) const
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner( i ).vertex() == point ) {
                return BME::bme_t( BME::CORNER, i ) ;
            }
        }
        return BME::bme_t() ;
    }

    /*!
     * @brief Get the index of the Corner at a given model point
     * @param[in] model_point_id Index of the point in the BoudaryModel
     * @return NO_ID or the index of the Corner
     */
    BME::bme_t BoundaryModelBuilder::find_corner( index_t model_point_id ) const
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner( i ).model_vertex_id() == model_point_id ) {
                return BME::bme_t( BME::CORNER, i ) ;
            }
        }
        return BME::bme_t() ;
    }

    /*!
     * @brief Create a corner at given coordinates.
     *
     * @param[in] point Geometric location of the new Corner
     * @return Index of the Corner
     */
    BME::bme_t BoundaryModelBuilder::create_corner( const vec3& point )
    {
        BME::bme_t id = create_element( BME::CORNER ) ;
        set_corner( id, point ) ;
        return id ;
    }

    /*!
     * @brief Find or create a corner at given coordinates.
     *
     * @param[in] point Geometric location of the Corner
     * @return Index of the Corner
     */
    BME::bme_t BoundaryModelBuilder::find_or_create_corner( const vec3& point )
    {
        BME::bme_t result = find_corner( point ) ;
        if( result.is_defined() ) {
            return result ;
        } else {
            return create_corner( point ) ;
        }
    }

    /*!
     * @brief Looks for a line in the model
     *
     * @param[in] vertices Coordinates of the vertices of the line
     * @return NO_ID or the index of the Line
     */
    BME::bme_t BoundaryModelBuilder::find_line(
        const std::vector< vec3 >& vertices ) const
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            if( model_.line( i ).equal( vertices ) ) {
                return BME::bme_t( BME::LINE, i ) ;
            }
        }
        return BME::bme_t() ;
    }

    /*!
     * @brief Add a Line knowing from teh coordinates of its vertices.
     * The corners are created if they do not exist.
     */
    BME::bme_t BoundaryModelBuilder::create_line( const std::vector< vec3 >& points )
    {
        BME::bme_t id = create_element( BME::LINE ) ;
        set_line( id, points ) ;

        // Find the indices of the corner at both extremities
        // Both must be defined to have a valid LINE
        add_element_boundary( id, find_or_create_corner( points.front() ) ) ;
        add_element_boundary( id, find_or_create_corner( points.back() ) ) ;

        return id ;
    }

    /*!
     * @brief Find or create a line
     *
     * @param[in] vertices Coordinates of the vertices of the line
     * @return Index of the Line
     */
    BME::bme_t BoundaryModelBuilder::find_or_create_line(
        const std::vector< vec3 >& vertices )
    {
        BME::bme_t result = find_line( vertices ) ;
        if( result.is_defined() ) {
            return result ;
        } else {
            return create_line( vertices ) ;
        }
    }

    /*!
     * @brief Create a surface
     *
     * @return Index of the Surface in the surfaces_ vector
     */
    BME::bme_t BoundaryModelBuilder::create_surface()
    {
        return create_element( BME::SURFACE ) ;
    }

    /*!
     * @brief Find a Contact
     * @param[in] interfaces Indices of the Interfaces determining the contact
     * @return NO_ID or index of the contact
     */
    BME::bme_t BoundaryModelBuilder::find_contact(
        const std::vector< index_t >& interfaces ) const
    {
        std::vector< const BoundaryModelElement* > comp( interfaces.size() ) ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            comp[i] = &model_.one_interface( interfaces[i] ) ;
        }

        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            if( comp.size() == model_.contact( i ).nb_in_boundary() ) {
                bool equal = true ;
                for( index_t j = 0; j < model_.contact( i ).nb_in_boundary(); j++ ) {
                    if( comp[j] != &model_.contact( i ).in_boundary( j ) ) {
                        equal = false ;
                        break ;
                    }
                }
                if( equal ) {
                    return BME::bme_t( BME::CONTACT, i ) ;
                }
            }
        }
        return BME::bme_t() ;
    }

    /*!
     * @brief Create a contact between the given Interfaces
     * The name of the contact is determined from the names of the interfaces.
     *
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    BME::bme_t BoundaryModelBuilder::create_contact(
        const std::vector< index_t >& interfaces )
    {
        // Create a name for this contact
        std::string name = "contact_" ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            name += model_.one_interface( interfaces[i] ).name() ;
            name += "_" ;
        }

        BME::bme_t id = create_element( BME::CONTACT ) ;
        set_element_name( id, name ) ;

        return id ;
    }

    /*!
     * @brief Find or create a contact between given Interfaces
     *
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    BME::bme_t BoundaryModelBuilder::find_or_create_contact(
        const std::vector< index_t >& interfaces )
    {
        BME::bme_t result = find_contact( interfaces ) ;
        if( result.is_defined() ) {
            return result ;
        } else {
            return create_contact( interfaces ) ;
        }
    }

    /*!
     * @brief Get the index of an Interface from its name
     *
     * @param[in] name Name of the Interface
     * @return Index of the interface in the model, NO_ID if not found.
     */
    BME::bme_t BoundaryModelBuilder::find_interface( const std::string& name ) const
    {
        for( index_t i = 0; i < model_.nb_interfaces(); ++i ) {
            if( model_.one_interface( i ).name() == name ) {
                return BME::bme_t( BME::INTERFACE, i ) ;
            }
        }
        return BME::bme_t() ;
    }

    /*!
     * @brief Create a new Interface
     *
     * @param[in] name Name of the interface
     * @param[in] type Type of the interface
     * @return The Interface index.
     */
    BME::bme_t BoundaryModelBuilder::create_interface(
        const std::string& name,
        BME::GEOL_FEATURE type )
    {
        BME::bme_t id = create_element( BME::INTERFACE ) ;
        set_element_geol_feature( id, type ) ;
        set_element_name( id, name ) ;
        return id ;
    }

    /*
     * @brief Adds an empty region to the model
     *
     *  Used in Geomodeling to convert a surface to a model
     */
    BME::bme_t BoundaryModelBuilder::create_region()
    {
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
    BME::bme_t BoundaryModelBuilder::create_region(
        const std::string& name,
        const std::vector< std::pair< index_t, bool > >& boundaries )
    {
        BME::bme_t id = create_element( BME::REGION ) ;
        set_element_name( id, name ) ;
        for( index_t i = 0; i < boundaries.size(); ++i ) {
            add_element_boundary( id,
                BME::bme_t( BME::SURFACE, boundaries[i].first ),
                boundaries[i].second ) ;
        }
        return id ;
    }

    /*!
     * @brief Creates a new empty Layer with the given name
     *
     * @param[in] name Name of the layer
     * @return The layer index
     */
    BME::bme_t BoundaryModelBuilder::create_layer( const std::string& name )
    {
        BME::bme_t id = create_element( BME::LAYER ) ;
        set_element_name( id, name ) ;
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
            ringmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary(
                BME::bme_t( BME::SURFACE, boundaries[i].first ),
                boundaries[i].second ) ;
        }
    }

    /*!
     * @brief Set the geometric location of a Corner
     *
     * @param[in] corner_id Index of the corner
     * @param[in] point Coordinates of the vertex
     */
    void BoundaryModelBuilder::set_corner(
        const BME::bme_t& corner_id,
        const vec3& point )
    {
        ringmesh_assert( corner_id.index < model_.nb_corners() ) ;
        model_.corners_[corner_id.index]->set_vertex( point, false ) ;
    }

    /*!
     * @brief Set one Line points
     *
     * @param[in] id Line index
     * @param[in] vertices Coordinates of the vertices on the line
     */
    void BoundaryModelBuilder::set_line(
        const BME::bme_t& id,
        const std::vector< vec3 >& vertices )
    {
        ringmesh_assert( id.index < model_.nb_lines() ) ;
        model_.lines_[id.index]->set_vertices( vertices ) ;
    }

    /*!
     * @brief Set the points and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] points Coordinates of the vertices
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const BME::bme_t& surface_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }

        model_.surfaces_[surface_id.index]->
            set_geometry( points, facets, facet_ptr ) ;

        // Test - remove degenerated facets 
      /*  GEO::Mesh& M = model_.surface( surface_id.index ).mesh() ;
        index_t nb = repair_surface_mesh( M ) ;
        if( nb > 0 ) {
            GEO::Logger::out( "BoundaryModel" )
                << " Removed " << nb << " degenerated facets in Surface "
                << surface_id.index << std::endl ;
        }
      */
        set_surface_adjacencies( surface_id ) ;
    }

    /*!
     * @brief Add a point to the BoundaryModel and not to one of its elements
     * @details To use when adding the points to the model before building its elements
     */
    index_t BoundaryModelBuilder::add_unique_vertex( const vec3& p )
    {
        return model_.vertices.add_unique_vertex( p ) ;
    }

    /*!
     * @brief Set the vertex for a Corner. Store the info in the BM vertices
     *
     * @param[in] corner_id Index of the corner
     * @param[in] unique_vertex Index of the vertex in the model
     */
    void BoundaryModelBuilder::set_corner(
        const BME::bme_t& corner_id,
        index_t unique_vertex )
    {
        ringmesh_assert( corner_id.index < model_.nb_corners() ) ;
        model_.corners_[corner_id.index]->set_vertex( unique_vertex ) ;
    }

    /*!
     * @brief Set one Line vertices. Store the info in the BM vertices
     *
     * @param[in] id Line index
     * @param[in] unique_vertices Indices in the model of the unique vertices with which to build the Line
     */
    void BoundaryModelBuilder::set_line(
        const BME::bme_t& id,
        const std::vector< index_t >& unique_vertices )
    {
        ringmesh_assert( id.index < model_.nb_lines() ) ;
        model_.lines_[id.index]->set_vertices( unique_vertices ) ;
    }

    /*!
     * @brief Set the vertices and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] model_vertex_ids Indices of unique vertices in the BoundaryModel
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const BME::bme_t& surface_id,
        const std::vector< index_t >& model_vertex_ids,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }

        model_.surfaces_[ surface_id.index ]->
            set_geometry( model_vertex_ids, facets, facet_ptr ) ;

        set_surface_adjacencies( surface_id ) ;
    }

    /*!
     * @brief Set the facets of a surface
     *
     * @param[in] surface_id Index of the surface
     * @param[in] facets Indices of the model vertices defining the facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets     
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const BME::bme_t& surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }

        // Compute the vertices from the corners
        // This is quite stupid !! The real solution would be to remove
        // the vertices vector from the Surface
        std::map< index_t, index_t > old_2_new ;

        std::vector< index_t > vertices ;
        std::vector< index_t > facets_local( facets.size() ) ;
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t c = facets[i] ;
            std::map< index_t, index_t >::iterator it = old_2_new.find( c ) ;
            index_t new_corner_id = NO_ID ;

            if( it == old_2_new.end() ) {
                new_corner_id = vertices.size() ;
                old_2_new[c] = new_corner_id ;

                // Not so great to push back, but whatever
                vertices.push_back( c ) ;
            } else {
                new_corner_id = old_2_new[c] ;
            }
            facets_local[i] = new_corner_id ;
        }

        set_surface_geometry( surface_id, vertices, facets_local, facet_ptr ) ;
    }

    /*!
     * @brief Compute and set the adjacencies between the facets
     * @details The adjacent facet is given for each vertex of each facet for the edge
     * starting at this vertex.
     * If there is no neighbor inside the same Surface adjacent is set to NO_ADJACENT
     *
     * @param[in] surface_id Index of the surface
     */
    void BoundaryModelBuilder::set_surface_adjacencies(
        const BME::bme_t& surface_id )
    {
        Surface& S = *model_.surfaces_[surface_id.index] ;
        ringmesh_assert( S.nb_cells() > 0 ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( S.nb_cells() - 1 ), Surface::NO_ADJACENT ) ;

        index_t nb_facets = S.nb_cells() ;
        index_t nb_vertices = S.nb_vertices() ;

        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > toto ;
        toto.reserve( 10 ) ;
        std::vector< std::vector< index_t > > vertex_to_facets( nb_vertices, toto ) ;

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                vertex_to_facets[S.surf_vertex_id( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ) {
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                index_t cur = S.surf_vertex_id( f, v ) ;
                index_t prev = S.surf_vertex_id( f, S.prev_in_facet( f, v ) ) ;

                const std::vector< index_t >& f_prev = vertex_to_facets[prev] ;
                const std::vector< index_t >& f_cur = vertex_to_facets[cur] ;

                std::vector< index_t > inter(
                    std::min( f_prev.size(), f_cur.size() ) ) ;
                index_t end = narrow_cast< index_t >(
                    std::set_intersection( f_prev.begin(), f_prev.end(),
                        f_cur.begin(), f_cur.end(), inter.begin() )
                        - inter.begin() ) ;

                if( end == 2 ) {
                    // There is one neighbor
                    index_t f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[S.facet_begin( f ) + S.prev_in_facet( f, v )] = f2 ;
                } else {
                    ringmesh_debug_assert( end == 1 ) ;
                }
            }
        }
        S.set_adjacent( adjacent ) ;
    }

    /*!
     * @brief Complete missing information in BoundaryModelElements
     * boundaries - in_boundary - parent - children
     *
     * @details For all 7 types of elements, check what information is available
     * for the first one and fill the elements of the same type accordingly
     * THIS MEANS that the all the elements of the same type have been initialized with
     * the same information.
     */
    bool BoundaryModelBuilder::complete_element_connectivity()
    {
        // Lines
        if( model_.nb_lines() > 0 ) {
            if( model_.line( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( BME::LINE ) ;
            }
            if( model_.line( 0 ).nb_in_boundary() == 0 ) {
                fill_elements_in_boundaries( BME::LINE ) ;
            }
            if( !model_.line( 0 ).parent_id().is_defined()
                && model_.nb_contacts() > 0 ) {
                fill_elements_parent( BME::LINE ) ;
            }
        }

        // Corners
        if( model_.nb_corners() > 0 && model_.corner( 0 ).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( BME::CORNER ) ;
        }

        // Surfaces - There MUST be at least one
        if( model_.surface( 0 ).nb_boundaries() == 0 ) {
            fill_elements_boundaries( BME::SURFACE ) ;
        }
        if( model_.surface( 0 ).nb_in_boundary() == 0 ) {
            fill_elements_in_boundaries( BME::SURFACE ) ;
        }
        if( !model_.surface( 0 ).parent_id().is_defined() ) {
            fill_elements_parent( BME::SURFACE ) ;
        }

        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( BME::REGION ) ;
            }
            if( !model_.region( 0 ).parent_id().is_defined()
                && model_.nb_layers() > 0 ) {
                fill_elements_parent( BME::REGION ) ;
            }
        }

        // Contacts
        if( model_.nb_contacts() > 0 && model_.contact( 0 ).nb_children() == 0 ) {
            fill_elements_children( BME::CONTACT ) ;
        }

        // Interfaces
        if( model_.nb_interfaces() > 0
            && model_.one_interface( 0 ).nb_children() == 0 ) {
            fill_elements_children( BME::INTERFACE ) ;
        }

        // Layers
        if( model_.nb_layers() > 0 && model_.layer( 0 ).nb_children() == 0 ) {
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
        ringmesh_debug_assert( type != BME::REGION ) ;

        BME::TYPE b_type = BME::boundary_type( type ) ;
        if( b_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < model_.nb_elements( b_type ); ++i ) {
                const BME& b = model_.element( BME::bme_t( b_type, i ) ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_element_boundary( b.in_boundary_id( j ),
                        BME::bme_t( b_type, i ) ) ;
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
                const BME& in_b = model_.element( BME::bme_t( in_b_type, i ) ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_element_in_boundary( in_b.boundary_id( j ),
                        BME::bme_t( in_b_type, i ) ) ;
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
                const BME& p = model_.element( BME::bme_t( p_type, i ) ) ;
                for( index_t j = 0; j < p.nb_children(); ++j ) {
                    set_parent( p.child_id( j ), BME::bme_t( p_type, i ) ) ;
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
                BME::bme_t cur_child = BME::bme_t( c_type, i ) ;
                const BME::bme_t& parent = model_.element( cur_child ).parent_id() ;
                if( parent.is_defined() ) {
                    add_child( parent, cur_child ) ;
                }
            }
        }
    }

    void BoundaryModelBuilder::remove_degenerate_facet_and_edges()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            index_t nb = repair_line_mesh( model_.line( i ).mesh() ) ;
            if( nb > 0 ) {
                GEO::Logger::out( "BoundaryModel" )
                    << nb << " degenerated edges removed in LINE "
                    << i << std::endl ;

                /// @todo The line may be empty now - remove it from the model
            }

        }

        for( index_t i = 0; i < model_.nb_surfaces(); ++i ) {
            GEO::Mesh& M = model_.surface( i ).mesh() ;
            index_t nb = detect_degenerate_facets( M ) ;
            if( nb > 0 ) {
                // If there are some degenerated facets 
                // We need to repair the model 
                GEO::Logger::out( "BoundaryModel" )
                    << nb << " degenerated facets in SURFACE "
                    << i << std::endl ;
                // Using repair function of geogram
                // Warning - This triangulates the mesh
                GEO::mesh_repair( M ) ;
                                
                // This might create some small components - remove them
                // How to choose the epsilon ? and the maximum number of facets ?
                GEO::remove_small_connected_components( M, epsilon_sq, 3 ) ;

                // Ok this is a bit of an overkill
                // The alternative is to copy mesh_repair and change it
                GEO::mesh_repair( M ) ;

                // If the Surface has internal boundaries, we need to 
                // re-cut the Surface along these lines
                Surface& S = *model_.surfaces_[ i ] ;
                for( index_t l = 0; l < S.nb_boundaries(); ++l ) {
                    const Line& L = model_.line( S.boundary_id( l ).index ) ;
                    if( L.is_inside_border( S ) ) {
                        S.cut_by_line( L ) ;
                    }
                }
            }
        }
        
    }

    /*!
     * @brief Fills the model nb_elements_per_type_ vector
     * @details See global element access with BoundaryModel::element( BME::TYPE, index_t )
     */
    void BoundaryModelBuilder::init_global_model_element_access()
    {
        index_t count = 0 ;
        model_.nb_elements_per_type_.push_back( count ) ;
        for( index_t type = BME::CORNER; type < BME::NO_TYPE; type++ ) {
            count += model_.nb_elements( (BME::TYPE) type ) ;
            model_.nb_elements_per_type_.push_back( count ) ;
        }
    }

    /*!
     * @brief This function MUST be the last function called when building a BoundaryModel
     *
     * @details Check that the model is correct and has all required information
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
        if( model_.name() == "" ) {
            set_model_name( "model_default_name" ) ;
        }

        init_global_model_element_access() ;

        complete_element_connectivity() ;

        // Fill geological feature if missing
        for( index_t i = 0; i < model_.nb_elements( BME::ALL_TYPES ); ++i ) {
            BME& E = element( BME::bme_t( BME::ALL_TYPES, i ) ) ;
            if( !E.has_geological_feature() ) {
                fill_element_geological_feature( E ) ;
            }
        }

        // Mesh repair for surfaces and lines
        // Not activated now - because it might create empty elements
        // That must be removed of the BoundaryModel
    
        remove_degenerate_facet_and_edges() ;

        GEO::Logger::out( "BoundaryModel" ) << "Model " << model_.name() << " has "
            << std::endl << std::setw( 10 ) << std::left << model_.nb_vertices()
            << " vertices " << std::endl << std::setw( 10 ) << std::left
            << model_.nb_facets() << " facets " << std::endl << std::endl
            << std::setw( 10 ) << std::left << model_.nb_regions() << " regions "
            << std::endl << std::setw( 10 ) << std::left << model_.nb_surfaces()
            << " surfaces " << std::endl << std::setw( 10 ) << std::left
            << model_.nb_lines() << " lines " << std::endl << std::setw( 10 )
            << std::left << model_.nb_corners() << " corners " << std::endl
            << std::endl << std::setw( 10 ) << std::left << model_.nb_contacts()
            << " contacts " << std::endl << std::setw( 10 ) << std::left
            << model_.nb_interfaces() << " interfaces " << std::endl
            << std::setw( 10 ) << std::left << model_.nb_layers() << " layers "
            << std::endl << std::endl ;

        // What do we do if the model is not valid ?
        if( model_.check_model_validity() ) {
            GEO::Logger::out( "BoundaryModel" ) << "Model " << model_.name()
                << " is valid " << std::endl ;
        } else {
            GEO::Logger::out( "BoundaryModel" ) << "Model " << model_.name()
                << " is invalid " << std::endl ;
        }

        return true ;
    }

    /*!
     * @brief Structure used to build Line by BoundaryModelBuilderGocad
     */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1 )
            : part_id_( part ), corner_id_( corner ), p0_( p0 ), p1_( p1 )
        {
        }

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
     *  - Gocad TSurf  <-> BoundaryModel Interface
     *  - Gocad TFace  <-> BoundaryModel Surface
     *  - Gocad Region <-> BoundaryModel Region
     *  - Gocad Layer  <-> BoundaryModel Layer
     *
     * @param[in] ml_file_name Input .ml file stream
     */
    void BoundaryModelBuilderGocad::load_ml_file( const std::string& ml_file_name )
    {
        GEO::LineInput in( ml_file_name ) ;
        if( !in.OK() ) {
            return ;
        }

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

        /// The file contains 2 parts and is read in 2 steps
        /// 1. Read global information on model elements
        /// 2. Read surface geometries and info to build corners and contacts
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and betweeen TSurfs

        // Indices of the vertices of the currently built TSurf in the model
        std::vector< vec3 > tsurf_vertices ;

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

        while( !in.eof() && in.get_line() ) {
            in.get_fields() ;
            if( in.nb_fields() > 0 ) {
                if( read_model ) {
                    if( strncmp( in.field( 0 ), "name:", 5 ) == 0 ) {
                        set_model_name( &in.field( 0 )[5] ) ;
                    } else if( in.field_matches( 0, "TSURF" ) ) {
                        /// 1.1 Create Interface from its name
                        index_t f = 1 ;
                        std::ostringstream oss ;
                        do {
                            oss << in.field( f++ ) ;
                        } while( f < in.nb_fields() ) ;
                        create_interface( oss.str() ) ;
                        nb_tsurf++ ;
                    } else if( in.field_matches( 0, "TFACE" ) ) {
                        /// 1.2 Create Surface from the name of its parent Interface
                        /// and its geological feature
                        index_t id = in.field_as_uint( 1 ) ;
                        std::string geol = in.field( 2 ) ;
                        index_t f = 3 ;
                        std::ostringstream oss ;
                        do {
                            oss << in.field( f++ ) ;
                        } while( f < in.nb_fields() ) ;
                        std::string interface_name = oss.str() ;

                        // And its key facet that give the orientation of the surface part
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p0( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p1( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p2( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;

                        create_surface( interface_name, geol, p0, p1, p2 ) ;
                        nb_tface++ ;
                    } else if( in.field_matches( 0, "REGION" ) ) {
                        /// 1.3 Read Region information and create them from their name,
                        /// and the surfaces on their boundary
                        index_t id = in.field_as_uint( 1 ) ;
                        std::string name = in.field( 2 ) ;

                        std::vector< std::pair< index_t, bool > > region_boundaries ;
                        bool end_region = false ;
                        while( !end_region ) {
                            in.get_line() ;
                            in.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                int s = in.field_as_int( i ) ;
                                if( s == 0 ) {
                                    end_region = true ;
                                    break ;
                                }
                                bool side = s > 0 ;
                                if( s > 0 ) {
                                    s -= 1 ;
                                } else {
                                    s = -s - 1 ;
                                }
                                region_boundaries.push_back(
                                    std::pair< index_t, bool >( s, side ) ) ;
                            }
                        }

                        // The Universe is not a regular region
                        if( name == "Universe" ) {
                            set_universe( region_boundaries ) ;
                        } else {
                            create_region( name, region_boundaries ) ;
                        }
                    } else if( in.field_matches( 0, "LAYER" ) ) {
                        /// 1.4 Build the volumetric layers from their name and
                        /// the ids of the regions they contain
                        BME::bme_t layer_id = create_layer( in.field( 1 ) ) ;
                        bool end_layer = false ;
                        while( !end_layer ) {
                            in.get_line() ;
                            in.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                index_t region_id = in.field_as_uint( i ) ;
                                if( region_id == 0 ) {
                                    end_layer = true ;
                                    break ;
                                } else {
                                    region_id -= nb_tface + 1 ; // Remove Universe region
                                    // Correction because ids begin at 1 in the file
                                    add_child( layer_id,
                                        BME::bme_t( BME::REGION, region_id - 1 ) ) ;
                                }
                            }
                        }
                    } else if( in.field_matches( 0, "END" ) ) {
                        // End of the high level information on the model
                        // Switch to reading the geometry of the model surfaces
                        read_model = false ;
                        continue ;
                    }
                } else {
                    if( in.field_matches( 0, "GOCAD" ) ) {
                        // This is the beginning of a new TSurf = Interface
                        tsurf_count++ ;
                    }
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( in.field_matches( 1, "Elevation" ) ) {
                            z_sign = 1 ;
                        } else if( in.field_matches( 1, "Depth" ) ) {
                            z_sign = -1 ;
                        } else {
                            ringmesh_assert_not_reached;}
                    } else if( in.field_matches( 0, "END" ) ) {
                        // This the END of a TSurf
                        if( tsurf_count > 0 ) {
                            // End the last TFace - Surface of this TSurf
                            set_surface_geometry(
                                BME::bme_t( BME::SURFACE, tface_count - 1 ),
                                std::vector< vec3 >(
                                    tsurf_vertices.begin() +
                                    tface_vertex_start.back(),
                                    tsurf_vertices.end() ),
                                tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;

                            // End this TSurf - Interface
                            nb_tface_in_prev_tsurf += tface_vertex_start.size() ;
                            tsurf_vertices.clear() ;
                            tface_vertex_start.clear() ;
                        }
                    } else if( in.field_matches( 0, "TFACE" ) ) {
                        // Beginning of a new TFace - Surface
                        if( tface_vertex_start.size() > 0 ) {
                            // End the previous TFace - Surface  (copy from line 1180)
                            set_surface_geometry(
                                BME::bme_t( BME::SURFACE, tface_count - 1),
                                std::vector< vec3 >(
                                    tsurf_vertices.begin() +
                                    tface_vertex_start.back(),
                                    tsurf_vertices.end() ),
                                tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;
                        }

                        // Register where begin the new TFace vertices
                        tface_vertex_start.push_back( tsurf_vertices.size() ) ;

                        tface_count++ ;
                    }

                    /// 2.1 Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    else if( in.field_matches( 0,
                            "VRTX" ) || in.field_matches( 0, "PVRTX" ) )
                    {
                        const vec3 p( read_double( in, 2 ), read_double( in,
                                3 ), z_sign * read_double( in, 4 ) ) ;
                        tsurf_vertices.push_back( p ) ;
                    } else if( in.field_matches( 0,
                            "PATOM" ) | in.field_matches( 0, "ATOM" ) )
                    {
                        tsurf_vertices.push_back( tsurf_vertices[ in.
                            field_as_uint(
                                2 ) - 1 ] ) ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        // Read ids of the vertices of each triangle in the TSurf
                        // and switch to ids in the TFace
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                1 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                2 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                3 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets_ptr.push_back( tface_facets.size() ) ;
                    }

                    // 2.2 Build the corners from their position and the surface parts
                    //    containing them
                    else if( in.field_matches( 0, "BSTONE" ) ) {
                        index_t v_id = in.field_as_uint( 1 ) - 1 ;
                        if( !find_corner(tsurf_vertices[v_id]).is_defined() )
                        {
                            create_corner( tsurf_vertices[ v_id ] ) ;
                        }
                    }

                    /// 2.3 Read the Border information and store it
                    else if( in.field_matches( 0, "BORDER" ) ) {
                        index_t p1 = in.field_as_uint( 2 ) - 1 ;
                        index_t p2 = in.field_as_uint( 3 ) - 1 ;

                        // Get the global corner id
                        BME::bme_t corner_id =
                        find_corner( tsurf_vertices[ p1 ] ) ;
                        ringmesh_assert( corner_id.is_defined() ) ;

                        // Get the surface
                        index_t part_id = NO_ID ;
                        for( index_t i = 0 ; i < tface_vertex_start.size() ; ++i ) {
                            if( p1 < tface_vertex_start[ i ] ) {
                                ringmesh_assert( p2 < tface_vertex_start[ i ] ) ;

                                // Get vertices ids in the surface
                                p1 = p1 - tface_vertex_start[ i - 1 ] ;
                                p2 = p2 - tface_vertex_start[ i - 1 ] ;

                                // i-1 is the id of the TFace in this TSurf
                                part_id = i - 1 ;
                                break ;
                            }
                        }
                        if( part_id == NO_ID ) {
                            // It is in the last built Tface
                            p1 = p1 - tface_vertex_start[ tface_vertex_start.size() - 1 ] ;
                            p2 = p2 - tface_vertex_start[ tface_vertex_start.size() - 1 ] ;

                            part_id = tface_vertex_start.size() - 1 ;
                        }

                        // The number of tfaces in previous tsurf is also to add
                        part_id += nb_tface_in_prev_tsurf ;

                        borders_to_build.push_back(
                            Border( part_id, corner_id.index, p1, p2 ) ) ;
                    }
                }
            }
        }

        // I agree that we do not need to compute the BoundaryModelVertices here
        // But perhaps the computation of Lines would be faster and safer - Jeanne

        /// 3. Build the Lines
        {
            std::vector< vec3 > line_vertices ;
            for( index_t i = 0; i < borders_to_build.size(); ++i ) {
                const Border& b = borders_to_build[i] ;

                // 1- Build the boundary : construct the vector
                // of vertices on the border
                const Surface& S = model_.surface( b.part_id_ ) ;

                BME::bme_t end_corner_id = determine_line_vertices( S, b.p0_, b.p1_,
                    line_vertices ) ;

                // 2 - Check if this border already exists
                BME::bme_t line_id = find_or_create_line( line_vertices ) ;

                // Add the surface in which this line is
                add_element_in_boundary( line_id,
                    BME::bme_t( BME::SURFACE, b.part_id_ ) ) ;
            }
        }

        /// 4. Build the Contacts
        build_contacts() ;

        // Modify in the Region the side of the Surface for which the key facet
        // orientation was not the same than their facet orientations
        for( index_t i = 0; i < change_key_facet.size(); i++ ) {
            const Surface& S = model_.surface( change_key_facet[i] ) ;
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {
                BoundaryModelElement& R = element( S.in_boundary_id( j ) ) ;
                for( index_t b = 0; b < R.nb_boundaries(); ++b ) {
                    if( R.boundary_id( b ).index == change_key_facet[i] ) {
                        bool old_side = R.side( b ) ;
                        R.set_boundary( b, R.boundary_id( b ), !old_side ) ;
                    }
                }
            }
        }

        /// 5. Fill missing information and check model validity
        if( !end_model() ) {
            GEO::Logger::err("BoundaryModel") << " Model " << model_.name() 
                << " is not a valid boundary representation. " << std::endl ;
        }

        time( &end_load ) ;
        // Output of loading time only in debug mode has no meaning (JP)
        GEO::Logger::out("I/O") << " Model loading time "
            << difftime( end_load, start_load ) << " sec" << std::endl ;
    }

    /*!
     * @brief Find the facet which first 3 vertices are given
     * 
     * @param[in] surface_id Index of the surface
     * @param[in] p0 First point coordinates
     * @param[in] p1 Second point coordinates
     * @param[in] p2 Third point coordinates
     * @param[out] same_sign Is true if the found facet has the same orientation than triangle p0p1p2
     * @return Index of the found facet, NO_ID if none found
     */
    index_t BoundaryModelBuilderGocad::find_key_facet(
        index_t surface_id,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool& same_sign ) const
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;

        for( index_t t = 0; t < surface.nb_cells(); ++t ) {
            const vec3& pp0 = surface.vertex( t, 0 ) ;
            const vec3& pp1 = surface.vertex( t, 1 ) ;
            const vec3& pp2 = surface.vertex( t, 2 ) ;

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
     * 
     * @param[in] surface_id Index of the surface
     * @return False if the key_facet orientation is not the same than the surface facets, else true.
     */
    bool BoundaryModelBuilderGocad::check_key_facet_orientation(
        index_t surface_id ) const
    {
        const KeyFacet& key_facet = key_facets_[surface_id] ;

        const vec3& p0 = key_facet.p0_ ;
        const vec3& p1 = key_facet.p1_ ;
        const vec3& p2 = key_facet.p2_ ;
        bool same_sign = false ;

        index_t t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ;
        if( t == NO_ID ) {
            vec3 p00( p0 ) ;
            p00.z *= -1 ;
            vec3 p10( p1 ) ;
            p10.z *= -1 ;
            vec3 p20( p2 ) ;
            p20.z *= -1 ;

            // It is because of the sign of Z that is not the same
            t = find_key_facet( surface_id, p00, p10, p20, same_sign ) ;
        }
        ringmesh_assert( t != NO_ID ) ;
        return same_sign ;
    }

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_vertices Coordinates of the vertices on the Line (emptied and filled again)
     * @return Index of the Corner at which the Line ends
     */
    BME::bme_t BoundaryModelBuilderGocad::determine_line_vertices(
        const Surface& S,
        index_t id0,
        index_t id1,
        std::vector< vec3 >& border_vertex_model_vertices ) const
    {
        ringmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_vertices.resize( 0 ) ;

        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        ringmesh_assert( f != Surface::NO_ID ) ;

        vec3 p0 = S.vertex( id0 ) ;
        vec3 p1 = S.vertex( id1 ) ;

        border_vertex_model_vertices.push_back( p0 ) ;
        border_vertex_model_vertices.push_back( p1 ) ;

        BME::bme_t p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        while( !p1_corner.is_defined() ) {
            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;

            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.vertex( next_id1 ) ;
            border_vertex_model_vertices.push_back( p1 ) ;
            p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        }
        return p1_corner ;
    }

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *
     * WE ASSUME THAT THE MODEL VERTICES ARE AVAILABLE AND CORRECT
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_ids Indices of vertices on the Line (resized at 0 at the beginning)
     * @return Index of the Corner at which the Line ends
     */
    BME::bme_t BoundaryModelBuilderGocad::determine_line_vertices(
        const Surface& S,
        index_t id0,
        index_t id1,
        std::vector< index_t >& border_vertex_model_ids ) const
    {
        ringmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_ids.resize( 0 ) ;

        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        ringmesh_assert( f != Surface::NO_ID ) ;

        // Global ids at the model level
        index_t p0 = S.model_vertex_id( id0 ) ;
        index_t p1 = S.model_vertex_id( id1 ) ;

        border_vertex_model_ids.push_back( p0 ) ;
        border_vertex_model_ids.push_back( p1 ) ;

        BME::bme_t p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        while( !p1_corner.is_defined() ) {
            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;

            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

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
            for( index_t j = 0; j < model_.line( i ).nb_in_boundary(); ++j ) {
                index_t sp_id = model_.line( i ).in_boundary_id( j ).index ;
                const BoundaryModelElement& p = model_.surface( sp_id ).parent() ;
                interfaces.insert( p.bme_id().index ) ;
            }
            std::vector< index_t > toto( interfaces.begin(), interfaces.end() ) ;
            BME::bme_t contact_id = find_or_create_contact( toto ) ;
            add_child( contact_id, BME::bme_t( BME::LINE, i ) ) ;
        }
    }

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] p0 Coordinates of the 1 point of the TFace key facet 
     * @param[in] p1 Coordinates of the 2 point of the TFace key facet 
     * @param[in] p2 Coordinates of the 3 point of the TFace key facet 
     */
    void BoundaryModelBuilderGocad::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        BME::bme_t parent = find_interface( interface_name ) ;
        if( interface_name != "" ) {
            ringmesh_assert( parent.is_defined() ) ;
        }

        BME::bme_t id = create_element( BME::SURFACE ) ;
        set_parent( id, parent ) ;
        set_element_geol_feature( parent, BME::determine_geological_type( type ) ) ;
        key_facets_.push_back( KeyFacet( p0, p1, p2 ) ) ;
    }

    bool BoundaryModelBuilderBM::load_file( const std::string& bm_file_name )
    {
        GEO::LineInput in( bm_file_name ) ;
        if( !in.OK() ) {
            return false ;
        }
        while( !in.eof() && in.get_line() ) {
            in.get_fields() ;
            if( in.nb_fields() > 0 ) {
                // Name of the model
                if( in.field_matches( 0, "NAME" ) ) {
                    if( in.nb_fields() > 1 ) {
                        set_model_name( in.field( 1 ) ) ;
                    }
                }

                // Number of elements of a given type
                else if( match_nb_elements( in.field( 0 ) ) != BME::NO_TYPE ) {
                    // Allocate the space
                    if( in.nb_fields() > 1 ) {
                        resize_elements( match_nb_elements( in.field( 0 ) ),
                            in.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level elements
                else if( match_high_level_type( in.field( 0 ) ) ) {
                    // Read this element
                    // First line id - name - geol_feature
                    if( in.nb_fields() < 4 ) {
                        std::cout << "I/O Error File line " << in.line_number()
                            << std::endl ;
                        return false ;
                    }
                    BME::TYPE t = match_type( in.field( 0 ) ) ;
                    index_t id = in.field_as_uint( 1 ) ;
                    BME::bme_t element( t, id ) ;
                    set_element_index( element ) ;
                    set_element_name( element, in.field( 2 ) ) ;
                    set_element_geol_feature( element,
                        BME::determine_geological_type( in.field( 3 ) ) ) ;

                    // Second line - indices of its children
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        add_child( element,
                            BME::bme_t( BME::child_type( t ),
                                in.field_as_uint( c ) ) ) ;
                    }
                }

                // Regions
                else if( match_type( in.field( 0 ) ) == BME::REGION ) {
                    // First line id - name
                    if( in.nb_fields() < 3 ) {
                        std::cout << "I/O Error File line " << in.line_number()
                            << std::endl ;
                        return false ;
                    }
                    index_t id = in.field_as_uint( 1 ) ;
                    BME::bme_t element( BME::REGION, id ) ;
                    set_element_index( element ) ;
                    set_element_name( element, in.field( 2 ) ) ;

                    // Second line - signed indices of boundaries
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( in.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &in.field( c )[1], s ) ;

                        add_element_boundary( element, BME::bme_t( BME::SURFACE, s ),
                            side ) ;
                    }
                }

                // Universe
                else if( in.field_matches( 0, "UNIVERSE" ) ) {
                    std::vector< std::pair< index_t, bool > > b_universe ;

                    // Second line - signed indices of boundaries
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( in.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &in.field( c )[1], s ) ;

                        b_universe.push_back(
                            std::pair< index_t, bool >( s, side ) ) ;
                    }
                    set_universe( b_universe ) ;
                }

                // Model vertices
//                else if( in.field_matches( 0, "MODEL_VERTICES" ) ) {
//                    index_t nb_vertices = in.field_as_uint( 1 ) ;
//
//                    // Attributes
//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "MODEL_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute< BoundaryModel::VERTEX > >
//                    vertex_attribs( nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            model_.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }
//                    for( index_t i = 0; i < nb_vertices; ++i ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        add_vertex( vec3(
//                                read_double( in,
//                                    0 ), read_double( in, 1 ), read_double( in, 2 ) ) ) ;
//                        serialize_read_attributes( in, 3, i, vertex_attribs ) ;
//                    }
//                }

                // Corners
                else if( match_type( in.field( 0 ) ) == BME::CORNER ) {
                    // One line id - vertex id
                    if( in.nb_fields() < 3 ) {
                        std::cout << "I/O Error File line " << in.line_number()
                            << std::endl ;
                        return false ;
                    }
                    index_t id = in.field_as_uint( 1 ) ;
                    BME::bme_t element( BME::CORNER, id ) ;
                    set_element_index( element ) ;
                    vec3 point( read_double( in, 2 ), read_double( in, 3 ),
                        read_double( in, 4 ) ) ;
                    set_element_vertex( element, 0, point ) ;
                }

                // Lines
                else if( match_type( in.field( 0 ) ) == BME::LINE ) {
                    index_t id = in.field_as_uint( 1 ) ;
                    BME::bme_t cur_element( BME::LINE, id ) ;
                    Line& L = dynamic_cast< Line& >( element( cur_element ) ) ;
                    L.set_id( id ) ;

                    // Following information - vertices of the lines
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "LINE_VERTICES" ) ) ;
                    index_t nb_vertices = in.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 point( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        vertices[i] = point ;
                    }

                    // Set the line points
                    L.set_vertices( vertices ) ;

                    // Attributes on line vertices
//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "LINE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute > vertex_attribs(
//                        nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            L.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices indices and attributes on vertices
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        serialize_read_attributes( in, 1, i, vertex_attribs ) ;
//                    }

                    // Read attributes on line segments
//                    in.get_line() ;
//                    in.get_fields() ;
//                   ringmesh_assert( in.field_matches( 0, "LINE_SEGMENT_ATTRIBUTES" ) ) ;
//                    index_t nb_segment_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    if( nb_segment_attribs > 0 ) {
//                        std::vector< SerializedAttribute< BME::FACET > >
//                        segment_attribs( nb_segment_attribs ) ;
//                        for( index_t i = 0; i < nb_segment_attribs; i++ ) {
//                            segment_attribs[ i ].bind(
//                                L.facet_attribute_manager(), in.field(
//                                    1 + 2 * i ), in.field( 2 + 2 * i ), L.nb_cells() ) ;
//                        }
//                        for( index_t i = 0; i < L.nb_cells(); i++ ) {
//                            in.get_line() ;
//                            in.get_fields() ;
//                            serialize_read_attributes( in, 1, in.field_as_uint(
//                                    0 ), segment_attribs ) ;
//                        }
//                    }

                    // Set the corners - they can be the same
                    add_element_boundary( cur_element,
                        find_corner( vertices.front() ) ) ;
                    add_element_boundary( cur_element,
                        find_corner( vertices.back() ) ) ;

                    // Finally we have the in_boundary information
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "IN_BOUNDARY" ) ) ;
                    for( index_t b = 1; b < in.nb_fields(); b++ ) {
                        L.add_in_boundary(
                            BME::bme_t( BME::SURFACE, in.field_as_uint( b ) ) ) ;
                    }
                }

                // Surfaces
                else if( match_type( in.field( 0 ) ) == BME::SURFACE ) {
                    index_t id = in.field_as_uint( 1 ) ;
                    BME::bme_t cur_element( BME::SURFACE, id ) ;
                    Surface& S = dynamic_cast< Surface& >( element( cur_element ) ) ;
                    S.set_id( id ) ;

                    // Read the surface vertices and their attributes
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_VERTICES" ) ) ;
                    index_t nb_vertices = in.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 point( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        vertices[i] = point ;
                    }

//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0,
//                            "SURFACE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_vertex_attribs = ( in.nb_fields() - 1 ) / 2 ;
//
                    // Bind the vertex attributes
//                    std::vector< SerializedAttribute< BME::VERTEX > > vertex_attribs(
//                        nb_vertex_attribs ) ;
//                    for( index_t i = 0; i < nb_vertex_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            S.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices global ids and attributes
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        serialize_read_attributes( in, 1, i, vertex_attribs ) ;
//                    }

                    // Read the surface facets
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_CORNERS" ) ) ;
                    index_t nb_corners = in.field_as_uint( 1 ) ;

                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_FACETS" ) ) ;
                    index_t nb_facets = in.field_as_uint( 1 ) ;

//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "SURFACE_FACET_ATTRIBUTES" ) ) ;
//                    index_t nb_facet_attribs = ( in.nb_fields() - 1 ) / 2 ;

                    // Bind the facet attributes
//                    std::vector< SerializedAttribute< BME::FACET > > facet_attribs(
//                        nb_facet_attribs ) ;
//                    for( index_t i = 0; i < nb_facet_attribs; i++ ) {
//                        facet_attribs[ i ].bind(
//                            S.facet_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_facets ) ;
//                    }

                    std::vector< index_t > corners( nb_corners ) ;
                    std::vector< index_t > facet_ptr( nb_facets + 1, 0 ) ;
                    index_t count_facets = 0 ;
                    for( index_t f = 0; f < nb_facets; f++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        index_t nb_v = in.field_as_uint( 0 ) ;
                        for( index_t v = 0; v < nb_v; ++v ) {
                            corners[count_facets + v] = in.field_as_uint( v + 1 ) ;
                        }
                        count_facets += nb_v ;
                        facet_ptr[f + 1] = count_facets ;
//                        serialize_read_attributes( in, nb_v + 1, f, facet_attribs ) ;
                    }

                    S.set_geometry( vertices, corners, facet_ptr ) ;
                    set_surface_adjacencies( cur_element ) ;
                }
            }
        }
        if( !end_model() ) {
            std::cout << "Invalid BoundaryModel loaded" << std::endl ;
        }
        return true ;
    }

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_nb_elements(
        const char* s )
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3 ) != 0 ) {
            return BME::NO_TYPE ;
        } else {
            for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
                BME::TYPE type = (BME::TYPE) i ;
                if( strstr( s, BME::type_name( type ).data() ) != NULL ) {
                    return type ;
                }
            }
            return BME::NO_TYPE ;
        }
    }

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_type( const char* s )
    {
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
         * @param v0, v1, v2 Indices in the BoundaryModel of the vertices defining the triangle
         *           the edge v0 - v1 is the one on the boundary
         */
        BorderTriangle( index_t s, index_t f, index_t v0, index_t v1, index_t v2 )
            : v0_( v0 ), v1_( v1 ), v2_( v2 ), s_( s ), f_( f )
        {
        }

        bool operator<( const BorderTriangle& rhs ) const
        {
            if( s_ != rhs.s_ ) {
                return s_ < rhs.s_ ;
            }
            if( f_ != rhs.f_ ) {
                return f_ < rhs.f_ ;
            }
            if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) {
                return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ ) ;
            }
            if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) {
                return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ ) ;
            }
            return rhs.v2_ == index_t( -1 ) ? false : v2_ < rhs.v2_ ;
        }

        bool same_edge( const BorderTriangle& rhs ) const
        {
            return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ )
                && std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ ) ;
        }

        /// Indices of the points in the surface. Triangle has the Surface orientation
        /// The edge v0v1 is, in surface s_, on the border.
        index_t v0_ ;
        index_t v1_ ;
        index_t v2_ ;

        /// Index of the model surface
        index_t s_ ;

        /// Index of the facet in the surface
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
        bool backward = false )
    {
        const BorderTriangle& in = BT[from] ;
        const Surface& S = M.surface( in.s_ ) ;
        index_t NO_ID( -1 ) ;

        // Get the next edge on border in the Surface
        index_t f = in.f_ ;
        index_t f_v0 = in.v0_ ;
        index_t f_v1 = in.v1_ ;
        ringmesh_assert( f_v0 != NO_ID && f_v1 != NO_ID ) ;

        index_t next_f = NO_ID ;
        index_t next_f_v0 = NO_ID ;
        index_t next_f_v1 = NO_ID ;

        if( !backward ) {
            S.next_on_border( f, f_v0, f_v1, next_f, next_f_v0, next_f_v1 ) ;
        } else {
            S.next_on_border( f, f_v1, f_v0, next_f, next_f_v0, next_f_v1 ) ;
        }

        // Find the BorderTriangle that is correspond to this
        // It must exist and there is only one
        BorderTriangle bait( in.s_, next_f, S.surf_vertex_id( next_f, next_f_v0 ),
            S.surf_vertex_id( next_f, next_f_v1 ), NO_ID ) ;

        // lower_bound returns an iterator pointing to the first element in the range [first,last)
        // which does not compare less than the given val.
        // See operator< on BorderTriangle
        index_t result = narrow_cast< index_t >(
            std::lower_bound( BT.begin(), BT.end(), bait ) - BT.begin() ) ;

        ringmesh_assert( result < BT.size() ) ;
        return result ;
    }

    /*!
     * @brief Mark as visited all BorderTriangle which first edge is the same than
     * the first edge of i.
     *
     * @param[in] border_triangles Information on triangles MUST be sorted so that
     *            BorderTriangle having the same boundary edge are adjacent
     * @param[in] i Index of reference BorderTriangle in border_triangles 
     * @param[out] visited Stores which of the border_triangles have a matching edge
     */
    void visit_border_triangle_on_same_edge(
        const std::vector< BorderTriangle >& border_triangles,
        index_t i,
        std::vector< bool >& visited )
    {
        index_t j = i ;
        while( j < border_triangles.size()
            && border_triangles[i].same_edge( border_triangles[j] ) ) {
            visited[j] = true ;
            j++ ;
        }
        signed_index_t k = i - 1 ;
        while( k > -1 && border_triangles[i].same_edge( border_triangles[k] ) ) {
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
        adjacent_surfaces.resize( 0 ) ;

        index_t j = i ;
        while( j < border_triangles.size()
            && border_triangles[i].same_edge( border_triangles[j] ) ) {
            adjacent_surfaces.push_back( border_triangles[j].s_ ) ;
            j++ ;
        }

        signed_index_t k = i - 1 ;
        while( k > -1 && border_triangles[i].same_edge( border_triangles[k] ) ) {
            adjacent_surfaces.push_back( border_triangles[k].s_ ) ;
            k-- ;
        }

        // In rare cases - the same surface can appear twice around a contact
        // Make unique and sort the adjacent regions
        std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() ) ;
        adjacent_surfaces.resize(
            std::unique( adjacent_surfaces.begin(), adjacent_surfaces.end() )
                - adjacent_surfaces.begin() ) ;
    }

    /*!
     * @brief From a BoundaryModel in which only Surface are defined, create
     * corners, contacts and regions.
     *
     */
    void BoundaryModelBuilderSurface::build_model()
    {
        ringmesh_assert( model_.nb_surfaces() > 0 ) ;

        /// 1. Make the storage of the model vertices unique
        /// So now we can make index comparison to find colocated edges

        // Force the computation of the model vertices to avoid troubles 
        model_.vertices.nb_unique_vertices() ;

        /// 2.1 Get for all Surface, the triangles that have an edge
        /// on the boundary.
        std::vector< BorderTriangle > border_triangles ;
        for( index_t s = 0; s < model_.nb_surfaces(); ++s ) {
            const Surface& S = model_.surface( s ) ;
            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                for( index_t v = 0; v < S.nb_vertices_in_facet( f ); ++v ) {
                    if( S.is_on_border( f, v ) ) {
                        border_triangles.push_back(
                            BorderTriangle( s, f, S.model_vertex_id( f, v ),
                                S.model_vertex_id( f, S.next_in_facet( f, v ) ),
                                S.model_vertex_id( f, S.prev_in_facet( f, v ) ) ) ) ;
                    }
                }
            }
        }

        /// 2.2 Sort these triangles so that triangles sharing the same edge follow one another
        std::sort( border_triangles.begin(), border_triangles.end() ) ;

        /// 3. Build the Lines and gather information to build the regions
        std::vector< SortTriangleAroundEdge > regions_info ;

        // The goal is to visit all BorderTriangle and propagate to get each Line vertices
        std::vector< bool > visited( border_triangles.size(), false ) ;
        for( index_t i = 0; i < border_triangles.size(); ++i ) {
            if( !visited[i] ) {
                // This is a new Line
                std::vector< vec3 > vertices ;

                // Get the indices of the Surfaces around this Line
                std::vector< index_t > adjacent ;
                get_adjacent_surfaces( border_triangles, i, adjacent ) ;

                // Mark as visited the BorderTriangle around the same first edge
                visit_border_triangle_on_same_edge( border_triangles, i, visited ) ;

                // Gather information to sort triangles around the contact
                regions_info.push_back( SortTriangleAroundEdge() ) ;
                index_t j = i ;
                while( j < border_triangles.size()
                    && border_triangles[i].same_edge( border_triangles[j] ) ) {
                    index_t cur_surface_id = border_triangles[j].s_ ;
                    const Surface& cur_surface = model_.surface( cur_surface_id ) ;
                    regions_info.back().add_triangle( cur_surface_id,
                        cur_surface.vertex( border_triangles[j].v0_ ),
                        cur_surface.vertex( border_triangles[j].v1_ ),
                        cur_surface.vertex( border_triangles[j].v2_ ) ) ;
                    j++ ;
                }

                // Add vertices to the Lin
                index_t surface_id = border_triangles[i].s_ ;
                const Surface& surface = model_.surface( surface_id ) ;
                vertices.push_back( surface.vertex( border_triangles[i].v0_ ) ) ;
                vertices.push_back( surface.vertex( border_triangles[i].v1_ ) ) ;

                // Build the contact propating forward on the border of the Surface
                // While the adjacent surfaces are the same the vertices the next edge on the
                // boundary of the Surface are added
                bool same_surfaces = true ;
                index_t next_i = get_next_border_triangle( model_, border_triangles,
                    i ) ;
                do {
                    ringmesh_assert( next_i != NO_ID ) ;
                    if( !visited[next_i] ) {
                        index_t cur_surface_id = border_triangles[next_i].s_ ;
                        const Surface& cur_surface = model_.surface(
                            cur_surface_id ) ;
                        std::vector< index_t > adjacent_next ;
                        get_adjacent_surfaces( border_triangles, next_i,
                            adjacent_next ) ;

                        if( adjacent.size() == adjacent_next.size()
                            && std::equal( adjacent.begin(), adjacent.end(),
                                adjacent_next.begin() ) ) {
                            visit_border_triangle_on_same_edge( border_triangles,
                                next_i, visited ) ;

                            // Add the next vertex
                            if( cur_surface.vertex( border_triangles[next_i].v0_ )
                                == vertices.back() ) {
                                vertices.push_back(
                                    cur_surface.vertex(
                                        border_triangles[next_i].v1_ ) ) ;
                            } else {
                                ringmesh_assert(
                                    cur_surface.vertex(
                                        border_triangles[next_i].v1_ )
                                        == vertices.back() ) ;
                                vertices.push_back(
                                    cur_surface.vertex(
                                        border_triangles[next_i].v0_ ) ) ;
                            }
                        } else {
                            same_surfaces = false ;
                        }
                    } else {
                        same_surfaces = false ;
                    }
                    next_i = get_next_border_triangle( model_, border_triangles,
                        next_i ) ;
                } while( same_surfaces && next_i != i ) ;

                if( next_i != i ) {
                    // Propagate backward to reach the other extremity
                    same_surfaces = true ;
                    index_t prev_i = get_next_border_triangle( model_,
                        border_triangles, i, true ) ;
                    do {
                        ringmesh_assert( prev_i != NO_ID && prev_i != i ) ;
                        if( !visited[prev_i] ) {
                            index_t cur_surface_id = border_triangles[prev_i].s_ ;
                            const Surface& cur_surface = model_.surface(
                                cur_surface_id ) ;
                            std::vector< index_t > adjacent_prev ;
                            get_adjacent_surfaces( border_triangles, prev_i,
                                adjacent_prev ) ;

                            if( adjacent.size() == adjacent_prev.size()
                                && std::equal( adjacent.begin(), adjacent.end(),
                                    adjacent_prev.begin() ) ) {
                                visit_border_triangle_on_same_edge( border_triangles,
                                    prev_i, visited ) ;

                                // Fill the Line vertices
                                if( cur_surface.vertex(
                                    border_triangles[prev_i].v0_ )
                                    == vertices.front() ) {
                                    vertices.insert( vertices.begin(),
                                        cur_surface.vertex(
                                            border_triangles[prev_i].v1_ ) ) ;
                                } else {
                                    ringmesh_assert(
                                        cur_surface.vertex(
                                            border_triangles[prev_i].v1_ )
                                            == vertices.front() ) ;
                                    vertices.insert( vertices.begin(),
                                        cur_surface.vertex(
                                            border_triangles[prev_i].v0_ ) ) ;
                                }
                            } else {
                                same_surfaces = false ;
                            }
                        } else {
                            same_surfaces = false ;
                        }
                        prev_i = get_next_border_triangle( model_, border_triangles,
                            prev_i, true ) ;
                    } while( same_surfaces ) ;
                }

                ringmesh_assert( vertices.size() > 1 )

                // At last create the Line
                BME::bme_t created = create_line( vertices ) ;
                for( index_t j = 0; j < adjacent.size(); ++j ) {
                    add_element_in_boundary( created,
                        BME::bme_t( BME::SURFACE, adjacent[j] ) ) ;
                }
            }
        }

        /// 4. Build the regions

        // Complete boundary information for surfaces
        // We need it to compute volumetric regions
        fill_elements_boundaries( BME::SURFACE ) ;

        /// 4.1 Sort surfaces around the contacts
        for( index_t i = 0; i < regions_info.size(); ++i ) {
            regions_info[i].sort() ;
        }

        if( model_.nb_surfaces() == 1 ) {
            /// \todo Build a Region when a BoundaryModel has only one Surface
            // Check that this surface is closed and define an interior
            // and exterior (universe) regions
            ringmesh_assert_not_reached;
        } else {
            // Each side of each Surface is in one Region( +side is first )
            std::vector< index_t > surf_2_region( 2 * model_.nb_surfaces(), NO_ID ) ;

            // Start with the first Surface on its + side
            std::stack< std::pair< index_t, bool > > S ;
            S.push( std::pair< index_t, bool > ( 0, true ) ) ;

            while( !S.empty() ) {
                std::pair< index_t, bool > cur = S.top() ;
                S.pop() ;

                // This side is already assigned
                if( surf_2_region[ cur.second == true ? 2 * cur.first : 2 *
                    cur.first + 1 ] != NO_ID ) {continue ;}

                // Create a new region
                BME::bme_t cur_region_id = create_region() ;

                std::stack< std::pair< index_t, bool > > SR ;
                SR.push( cur ) ;
                while( !SR.empty() ) {
                    std::pair< index_t, bool > s = SR.top() ;
                    SR.pop() ;
                    index_t s_id = s.second == true ? 2 * s.first : 2 * s.first + 1 ;

                    // This side is already assigned
                    if( surf_2_region[ s_id ] != NO_ID ) {continue ;}

                    // Add the surface to the current region
                    add_element_boundary( cur_region_id, BME::bme_t( BME::SURFACE, s.first ),
                        s.second ) ;
                    surf_2_region[ s_id ] = cur_region_id.index ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp = !s.second == true ? 2 * s.first : 2 *
                    s.first + 1 ;
                    if( surf_2_region[ s_id_opp ] == NO_ID ) {
                        S.push( std::pair< index_t, bool >( s.first, !s.second ) ) ;
                    }

                    // For each contact, push the next oriented surface that is in the same region
                    const BoundaryModelElement& surface = model_.surface( s.first ) ;
                    for( index_t i = 0 ; i < surface.nb_boundaries() ; ++i ) {
                        const std::pair< index_t, bool >& n =
                        regions_info[ surface.boundary_id( i ).index ].next( s ) ;
                        index_t n_id = n.second == true ? 2 * n.first : 2 *
                        n.first + 1 ;

                        if( surf_2_region[ n_id ] == NO_ID ) {SR.push( n ) ;}
                    }
                }
            }

            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            /// \todo Implement the code to take into regions included in others (bubbles)
            ringmesh_assert( std::count( surf_2_region.begin(), surf_2_region.end(),
                    NO_ID ) == 0 ) ;
        }

        // We need to remove from the regions_ the one corresponding
        // to the universe_, the one with the biggest volume
        double max_volume = -1. ;
        index_t universe_id = NO_ID ;
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            double cur_volume = BoundaryModelElementMeasure::size(
                &model_.region( i ) ) ;
            if( cur_volume > max_volume ) {
                max_volume = cur_volume ;
                universe_id = i ;
            }
        }

        const BoundaryModelElement& cur_region = model_.region( universe_id ) ;
        std::vector< std::pair< index_t, bool > > univ_boundaries(
            cur_region.nb_boundaries() ) ;
        for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ) {
            univ_boundaries[i].first = cur_region.boundary( i ).bme_id().index ;
            univ_boundaries[i].second = cur_region.side( i ) ;
        }
        set_universe( univ_boundaries ) ;

        // Decrease by one the ids of the regions that are after the
        // one converted to the universe
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            index_t cur_id = model_.region( i ).bme_id().index ;
            if( i > universe_id ) {
                element( BME::bme_t( BME::REGION, i ) ).set_id( cur_id - 1 ) ;
            }
        }

        // Remove the region converted to universe from the regions
        erase_element( BME::bme_t( BME::REGION, universe_id ) ) ;

        // We are not in trouble since the boundaries of surface are not yet set
        // And we have no layer in the model

        ringmesh_assert( end_model() ) ;
    }
} // namespace
