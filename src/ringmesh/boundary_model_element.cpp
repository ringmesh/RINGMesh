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

/*! \author Jeanne Pellerin and Arnaud Botella */

#include <ringmesh/boundary_model_element.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/utils.h>

#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_AABB.h>

#include <set>
#include <stack>
#include <fstream>

namespace RINGMesh {
    /*!
     * @brief Map the name of a geological type with a value of GEOL_FEATURE
     *
     * @param[in] in Name of the feature
     * @return The geological feature index
     *
     * \todo Keep all the information ( add new GEOL_FEATURE) instead of simplfying it.
     */
    BoundaryModelElement::GEOL_FEATURE BoundaryModelElement::
    determine_geological_type( const std::string& in )
    {
        if( in == "" ) {return NO_GEOL ;}
        if( in == "reverse_fault" ) {return FAULT ;}
        if( in == "normal_fault" ) {return FAULT ;}
        if( in == "fault" ) {return FAULT ;}
        if( in == "top" ) {return STRATI ;}
        if( in == "none" ) {return STRATI ;}
        if( in == "unconformity" ) {return STRATI ;}
        if( in == "boundary" ) {return VOI ;}

        std::cout << "ERROR" << "Unexpected type in the model file " << in
                  << std::endl ;
        return NO_GEOL ;
    }


    /*!
     * @brief Compute an intersection type
     *
     * @param[in] types Type that intersect
     * @return Intersection type
     */
    BoundaryModelElement::GEOL_FEATURE BoundaryModelElement::determine_type(
        const std::vector< GEOL_FEATURE >& types )
    {
        if( types.size() == 0 ) {return NO_GEOL ;}

        // Sort and remove duplicates form the in types
        std::vector< GEOL_FEATURE > in = types ;
        std::sort( in.begin(), in.end() ) ;
        index_t new_size = narrow_cast< index_t >(std::unique( in.begin(), in.end() ) - in.begin()) ;
        in.resize( new_size ) ;

        if( in.size() == 1 ) {return in[ 0 ] ;}

        if( in.size() == 2 ) {
            if( in[ 0 ] == NO_GEOL ) {return NO_GEOL ;}
            if( in[ 0 ] == STRATI ) {
                if( in[ 1 ] == FAULT ) {return STRATI_FAULT ;}
                if( in[ 1 ] == VOI ) {return STRATI_VOI ;}
            } else if( in[ 0 ] == FAULT ) {
                if( in[ 1 ] == VOI ) {return FAULT_VOI ;}
            }

            // Other cases ? for corners ? what is the vertex ?
            return NO_GEOL ;
        }
        return NO_GEOL ;
    }


    std::string BoundaryModelElement::type_name( BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case CORNER    : return "CORNER" ;
             case LINE      : return "LINE" ;
             case SURFACE   : return "SURFACE" ;
             case REGION    : return "REGION" ;
             case CONTACT   : return "CONTACT" ;
             case INTERFACE : return "INTERFACE" ;
             case LAYER     : return "LAYER" ;
             default        : return "NO_TYPE_NAME" ;
        }
    }


    std::string BoundaryModelElement::geol_name(
        BoundaryModelElement::GEOL_FEATURE t )
    {
        switch( t ) {
             case STRATI  : return "top" ;
             case FAULT   : return "fault" ;
             case VOI     : return "boundary" ;
             case NO_GEOL : return "none" ;
             default      : return "none" ;
                 break ;
        }
    }


    /*!
     * @brief Define the type of the parent of an element of type @param t
     *        If no parent is allowed returns NO_TYPE
     * @details The elements that can have a parent are LINE, SURFACE, and REGION
     */
    BoundaryModelElement::TYPE BoundaryModelElement::parent_type(
        BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case LINE      : return CONTACT ;
             case SURFACE   : return INTERFACE ;
             case REGION    : return LAYER ;
             default :

                 // The others have no parent
                 return NO_TYPE ;
        }
    }


    /*!
     * @brief Define the type of a child of an element of type @param t
     *        If no child is allowed returns NO_TYPE
     * @details The elements that can have a parent are CONTACT, INTERFACE, and LAYER
     */
    BoundaryModelElement::TYPE BoundaryModelElement::child_type(
        BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case CONTACT   : return LINE  ;
             case INTERFACE : return SURFACE ;
             case LAYER     : return REGION ;
             default :
                 return NO_TYPE ;
        }
    }


    /*!
     * @brief Define the type of an element on the boundary of an element of type @param t
     *        If no boundary is allowed returns NO_TYPE
     * @details The elements that can have a boundary are LINE, SURFACE, and REGION
     */
    BoundaryModelElement::TYPE BoundaryModelElement::boundary_type(
        BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case LINE      : return CORNER ;
             case SURFACE   : return LINE ;
             case REGION    : return SURFACE ;
             default :
                 return NO_TYPE ;
        }
    }


    /*!
     * @brief Define the type of an element into which boundary an element of type @param t can be
     *        If no in_boundary is allowed returns NO_TYPE
     * @details The elements that can be in the boudanry of another are CORNER, LINE, and SURFACE
     */
    BoundaryModelElement::TYPE BoundaryModelElement::in_boundary_type(
        BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case CORNER    : return LINE ;
             case LINE      : return SURFACE ;
             case SURFACE   : return REGION ;
             default :
                 return NO_TYPE ;
        }
    }


    /*!
     * @brief Dimension 0, 1, 2, or 3 of anelement of type @param t
     */
    index_t BoundaryModelElement::dimension( BoundaryModelElement::TYPE t )
    {
        switch( t ) {
             case CORNER    : return 0 ;
             case LINE      : return 1 ;
             case SURFACE   : return 2 ;
             case REGION    : return 3 ;
             case CONTACT   : return 1 ;
             case INTERFACE : return 2 ;
             case LAYER     : return 3 ;
             default        : return NO_ID ;
        }
    }


    bool BoundaryModelElement::operator==( const BoundaryModelElement& rhs ) const
    {
        if( model_ != rhs.model_ ) {return false ;}
        if( type_ != rhs.type_ ) {return false ;}
        if( id_ != rhs.id_ ) {return false ;}
        if( name_ != rhs.name_ ) {return false ;}
        if( geol_feature_ != rhs.geol_feature_ ) {return false ;}
        if( nb_boundaries() != rhs.nb_boundaries() ) {return false ;}
        if( !std::equal( boundaries_.begin(), boundaries_.end(),
                rhs.boundaries_.begin() ) ) { return false ;}
        if( !std::equal( sides_.begin(), sides_.end(),
                rhs.sides_.begin() ) ) { return false ;}
        if( nb_in_boundary() != rhs.nb_in_boundary() ) {return false ;}
        if( !std::equal( in_boundary_.begin(), in_boundary_.end(),
                rhs.in_boundary_.begin() ) ) { return false ;}
        if( parent_ != rhs.parent_ ) {return false ;}
        if( nb_children() != rhs.nb_children() ) {return false ;}
        if( !std::equal( children_.begin(), children_.end(),
                rhs.children_.begin() ) ) { return false ;}

        return true ;
    }


    /*!
     * @return Assert that the parent exists and returns it.
     */
    const BoundaryModelElement& BoundaryModelElement::parent() const
    {
        ringmesh_assert( parent_id() != NO_ID ) ;
        return model_->element( parent_type( type_ ), parent_id() ) ;
    }


    /*!
     *
     * @param[in] x Index of the boundary element
     * @return Assert that is exits and return the element on the boundary
     */
    const BoundaryModelElement& BoundaryModelElement::boundary( index_t x ) const
    {
        ringmesh_assert( x < nb_boundaries() ) ;
        return model_->element( boundary_type( type_ ), boundary_id( x ) ) ;
    }


    /*!
     *
     * @param[in] x Index of the in_boundary element
     * @return Assert that it exist and return the element in in_boundary.
     */
    const BoundaryModelElement& BoundaryModelElement::in_boundary( index_t x ) const
    {
        ringmesh_assert( x < nb_in_boundary() ) ;
        return model_->element( in_boundary_type( type_ ), in_boundary_id( x ) ) ;
    }


    /*!
     *
     * @param[in] x Index of the child
     * @return Assert that the child exists and returns it.
     */
    const BoundaryModelElement& BoundaryModelElement::child( index_t x ) const
    {
        ringmesh_assert( x < nb_children() ) ;
        return model_->element( child_type( type_ ), child_id( x ) ) ;
    }


    /*!
     * @brief Copy all attribute except model_ from @param rhs to this element
     * @param[in] rhs To copy from
     * @param[in] model Model to associate to this element
     */
    void BoundaryModelElement::copy_macro_topology(
        const BoundaryModelElement& rhs,
        BoundaryModel& model )
    {
        model_        = &model ;
        type_         = rhs.type_ ;
        id_           = rhs.id_ ;
        name_         = rhs.name_ ;
        geol_feature_ = rhs.geol_feature_ ;
        boundaries_   = rhs.boundaries_ ;
        sides_        = rhs.sides_ ;
        in_boundary_  = rhs.in_boundary_ ;
        parent_       = rhs.parent_ ;
        children_     = rhs.children_ ;
    }


    /*!
     * @brief Checks if this element or one of the element containing it
     * determines the model Volume Of Interest
     * @details This is known with the type of an element
     */
    bool BoundaryModelElement::is_on_voi() const
    {
        if( geol_feature_ == NO_GEOL ) {
            for( index_t j = 0; j < nb_in_boundary(); ++j ) {
                GEOL_FEATURE t = in_boundary( j ).geological_feature() ;
                if( t == VOI || t == STRATI_VOI || t == FAULT_VOI ) {
                    return true ;
                }
            }
        } else if( geol_feature_ == VOI        ||
                   geol_feature_ == STRATI_VOI ||
                   geol_feature_ == FAULT_VOI )
        {
            return true ;
        }
        return false ;
    }

    index_t Corner::model_vertex_id( index_t p ) const
    {
       // return model_->vertices.unique_vertex_id( CORNER, id(), p ) ;
       ringmesh_debug_assert( model_vertex_id_[0] != NO_ID ) ;

       return model_vertex_id_[0] ;
    }

    
    void Corner::set_vertex( index_t model_point_id ) 
    {
        mesh_.vertices.point( 0 ) = model_->vertex( model_point_id ) ;
        model_vertex_id_[0] = model_point_id ;
        model_->vertices.add_unique_to_bme( model_point_id, element_type(), id(), 0 ) ;
    }

    void Corner::set_vertex( index_t index, const vec3& point, bool update )
    {
        if( update )
            model_->vertices.update_point(
                model_->vertices.unique_vertex_id( CORNER, id(), index ), point ) ;
        else
            mesh_.vertices.point( 0 ) = point ;
    }

    /*!
     * @brief Construct a Line
     *
     * @param[in] model The parent model
     * @param[in] id The index of the line in the lines_ vector of the parent model
     */
    Line::Line(
        BoundaryModel* model,
        index_t id ) :
        BoundaryModelElement( model, LINE, id )
    {
        model_vertex_id_.bind( mesh_.vertices.attributes(), "model_vertex_id" ) ; 
    }


    /*!
     * @brief Construct a Line knowing its vertices
     *
     * @param[in] model  The parent model
     * @param[in] id The index of the line in the lines_ vector of the parent model
     * @param[in] vertices Indices (in the model) of the vertices defining this Line
     */
    Line::Line(
        BoundaryModel* model,
        index_t id,
        const std::vector< vec3 >& vertices )
          : BoundaryModelElement( model, LINE, id )
    {
        mesh_.vertices.create_vertices( vertices.size() ) ;
        for( index_t v = 0; v < vertices.size(); v++ ) {
            mesh_.vertices.point( v ) = vertices[v] ;
        }
        model_vertex_id_.bind( mesh_.vertices.attributes(), "model_vertex_id" ) ; 
    }


    /*!
     * @brief Construct a Line knowing its vertices
     *
     * @param[in] model  The parent model
     * @param[in] id The index of the line in the lines_ vector of the parent model
     * @param[in] vertices Indices (in the model) of the vertices defining this Line
     * @param[in] corner0 Index of the starting corner
     * @param[in] corner1 Index of the ending corner
     */
    Line::Line(
        BoundaryModel* model,
        index_t id,
        index_t corner0,
        index_t corner1,
        const std::vector< vec3 >& vertices )
        : BoundaryModelElement( model, LINE, id )
    {
        mesh_.vertices.create_vertices( vertices.size() ) ;
        for( index_t v = 0; v < vertices.size(); v++ ) {
            mesh_.vertices.point( v ) = vertices[v] ;
        }
        model_vertex_id_.bind( mesh_.vertices.attributes(), "model_vertex_id" ) ; 

        boundaries_.push_back( corner0 ) ;
        boundaries_.push_back( corner1 ) ;
    }

    index_t Line::model_vertex_id( index_t p ) const
    {
        ringmesh_debug_assert( model_vertex_id_[p] != NO_ID ) ;
        return model_vertex_id_[p] ;
    }

    void Line::set_vertex( index_t index, const vec3& point, bool update )
    {
        if( update )
            model_->vertices.update_point(
                model_->vertices.unique_vertex_id( CORNER, id(), index ),
                point ) ;
        else
            mesh_.vertices.point( index ) = point ;
    }

    void Line::set_vertices( const std::vector< index_t >& model_vertex_ids )
    {
        mesh_.clear( true, true ) ;
        mesh_.vertices.create_vertices( model_vertex_ids.size() ) ;
        for( index_t v = 0; v < model_vertex_ids.size(); v++ ) {
            index_t cur = model_vertex_ids[v] ;
            mesh_.vertices.point( v ) = model_->vertex( cur ) ;
            model_vertex_id_[v] = cur ;
            model_->vertices.add_unique_to_bme( cur, element_type(), id(), v ) ;
        }
    }

    void Line::set_model_vertex_id( index_t line_id, index_t model_id ) {
        model_vertex_id_[line_id] = model_id ;
    }

    /*!
     * @brief Check if the Line is twice on the boundary of a surface
     *
     * @param[in] surface The surface to test
     */
    bool Line::is_inside_border( const BoundaryModelElement& surface ) const
    {
        // Find out if this surface is twice in the in_boundary vector
        return std::count( in_boundary_.begin(), in_boundary_.end(),
            surface.id() ) > 1 ;
    }


    /*!
     *
     * @param[in] s Segment index
     * @return The coordinates of the barycenter of the segment
     */
    vec3 Line::segment_barycenter( index_t s ) const
    {
        return 0.5 * ( vertex( s ) + vertex( s + 1 ) ) ;
    }


    /*!
     *
     * @param[in] s Segment index
     * @return The length of the segment
     */
    double Line::segment_length( index_t s ) const
    {
        return length( vertex( s + 1 ) - vertex( s ) ) ;
    }


    /*!
     *
     * @return The length of the Line
     */
    double Line::total_length() const
    {
        double result = 0 ;
        for( index_t s = 0; s < nb_cells(); s++ ) {
            result += segment_length( s ) ;
        }
        return result ;
    }


    /*!
     * @brief Returns true if the Line has exaclty the given vertices
     *
     * @param[in] rhs_vertices Vertices to compare to
     */
    bool Line::equal( const std::vector< vec3 >& rhs_vertices ) const
    {
        if( nb_vertices() != rhs_vertices.size() ) {
            return false ;
        }

        bool equal = true ;
        for( index_t i = 0; i < nb_vertices(); i++ ) {
            if( rhs_vertices[i] != mesh_.vertices.point( i ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) return true ;

        equal = true ;
        for( index_t i = 0; i < nb_vertices(); i++ ) {
            if( rhs_vertices[i] != mesh_.vertices.point( nb_vertices()-i-1 ) ) {
                equal = false ;
                break ;
            }
        }
        if( equal ) return true ;

        return false ;
    }

    Surface::~Surface()
    {
    }

    /*!
     * @param[in] f Facet index
     * @param[in] v Vertex index in the facet
     * @return The coordinates of the vertex
     */
    const vec3& Surface::vertex(
        index_t f,
        index_t v ) const
    {
        ringmesh_debug_assert( v < nb_vertices_in_facet( f ) ) ;
        return vertex( surf_vertex_id( f, v ) ) ;
    }

    void Surface::set_vertex(
        index_t index,
        const vec3& point,
        bool update )
    {
        ringmesh_debug_assert( index < nb_vertices() ) ;
        if( update )
            model_->vertices.update_point(
                 model_->vertices.unique_vertex_id( SURFACE, id(), index ),
                point ) ;
        else
            mesh_.vertices.point( index ) = point ;
    }

    index_t Surface::model_vertex_id( index_t p ) const {
        ringmesh_debug_assert( model_vertex_id_[p] != NO_ID ) ;
        return model_vertex_id_[p] ;
    }

    /*!
     * @param[in] surf_vertex_id Index of the vertex in the surface
     * @return The coordinates of the vertex
     */
    const vec3& Surface::vertex( index_t surf_vertex_id ) const
    {
        return mesh_.vertices.point( surf_vertex_id ) ;
    }

    index_t Surface::surf_vertex_id( index_t model_vertex_id ) const
    {
        const std::vector< BoundaryModelVertices::VertexInBME >& reverse_db =
            model_->vertices.bme_vertices( model_vertex_id ) ;
        for( index_t i = 0; i < reverse_db.size(); i++ ) {
            const BoundaryModelVertices::VertexInBME& info = reverse_db[i] ;
            if( info.bme_type == SURFACE && info.bme_id == id() ) {
                return info.v_id ;
            }
        }
        return NO_ID ;
    }


    void Surface::set_model_vertex_id( index_t surf_id, index_t model_id ) {
        model_vertex_id_[surf_id] = model_id ;
    }

    void Surface::set_geometry(
        const std::vector< vec3 >& vertices,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        mesh_.clear( true, true ) ;
        mesh_.vertices.create_vertices( vertices.size() ) ;
        for( index_t v = 0; v < vertices.size() ; v++ ) {
            mesh_.vertices.point( v ) = vertices[v] ;
        }
        set_geometry( facets, facet_ptr ) ;
    }

    void Surface::set_geometry(
        const std::vector< index_t >& model_vertex_ids,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        mesh_.clear( true, true ) ;
        mesh_.vertices.create_vertices( model_vertex_ids.size() ) ;
        for( index_t v = 0; v < model_vertex_ids.size() ; v++ ) {
            index_t cur = model_vertex_ids[v] ;
            mesh_.vertices.point( v ) = model_->vertex( cur ) ;
            model_vertex_id_[v] = cur ;
            model_->vertices.add_unique_to_bme( cur, element_type(), id(), v ) ;
        }
        set_geometry( facets, facet_ptr ) ;
    }

    void Surface::set_geometry(
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        for( index_t f = 0; f < facet_ptr.size()-1; f++ ) {
            index_t size = facet_ptr[f+1] - facet_ptr[f] ;
            GEO::vector< index_t > facet_vertices( size ) ;
            index_t start = facet_ptr[f] ;
            for( index_t lv = 0; lv < size; lv++ ) {
                facet_vertices[lv] = facets[start++] ;
            }
            mesh_.facets.create_polygon( facet_vertices ) ;
        }
    }
    /*!
     * @brief Traversal of a surface border
     * @details From the input facet f, get the facet that share vertex v and
     * get the indices of vertex v and of the following vertex in this new facet.
     * The next facet next_f may be the same, and @param is required to avoid going back.
     *
     * @param[in] f Index of the facet
     * @param[in] from Index in the facet of the previous point on the border - gives the direction
     * @param[in] v Index in the facet of the point for which we want the next point on border
     * @param[out] next_f Index of the facet containing the next point on border
     * @param[out] v_in_next Index of vertex v in facet next_f
     * @param[out] next_in_next Index of the next vertex on border in facet v_in_next
     */
    void Surface::next_on_border(
        index_t f,
        index_t from,
        index_t v,
        index_t& next_f,
        index_t& v_in_next,
        index_t& next_in_next ) const
    {
        ringmesh_assert( v < nb_vertices_in_facet( f ) ) ;
        ringmesh_assert( is_on_border( f, v ) || is_on_border( f, from ) ) ;

        index_t V = surf_vertex_id( f, v ) ;

        // We want the next triangle that is on the boundary and share V
        // If there is no such triangle, the next vertex on the boundary
        // is the vertex of F neighbor of V that is not from

        // Get the facets around the shared vertex that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets ;
        index_t nb_around = facets_around_vertex( V, facets, true, f ) ;
        ringmesh_assert( nb_around < 3 && nb_around > 0 ) ;

        next_f = facets[ 0 ] ;

        if( nb_around == 2 ) {
            if( next_f == f ) {next_f = facets[ 1 ] ;}
            ringmesh_debug_assert( next_f != NO_ID ) ;

            // Now get the other vertex that is on the boundary opposite to p1
            v_in_next = facet_vertex_id( next_f, V ) ;
            ringmesh_assert( v_in_next != NO_ID ) ;

            // The edges containing V in next_f are
            // the edge starting at v_in_next and the one ending there
            index_t prev_v_in_next = prev_in_facet( next_f, v_in_next )  ;

            bool e0_on_boundary = is_on_border( next_f, v_in_next ) ;
            bool e1_on_boundary = is_on_border( next_f, prev_v_in_next ) ;

            // Only one must be on the boundary otherwise there is a corner missing
            ringmesh_assert( e0_on_boundary != e1_on_boundary ) ;

            // From the edge that is on boundary get the next vertex on this boundary
            // If the edge starting at p_in_next is on boundary, new_vertex is its next
            // If the edge ending at p_in_next is on boundary, new vertex is its prev
            next_in_next =
                e0_on_boundary ? next_in_facet( next_f, v_in_next ) : prev_v_in_next ;
        } else if( nb_around == 1 ) {
            // V must be in two border edges of facet f
            // Get the id in the facet of the vertex neighbor of v1 that is not v0
            v_in_next = v ;
            if( prev_in_facet( f, v ) == from  ) {
                ringmesh_debug_assert( is_on_border( f, v ) ) ;
                next_in_next = next_in_facet( f, v ) ;
            } else {
                ringmesh_debug_assert( is_on_border( f, prev_in_facet( f, v ) ) ) ;
                next_in_next = prev_in_facet( f, v ) ;
            }
        }
    }


    /*!
     * @brief Get the next edge on the border
     * @param[in] f Input facet index
     * @param[in] e Edge index in the facet
     * @param[out] next_f Next facet index
     * @param[out] next_e Next edge index in the facet
     */
    void Surface::next_on_border(
        index_t f,
        index_t e,
        index_t& next_f,
        index_t& next_e ) const
    {
        index_t v = next_in_facet( f, e ) ;
        return next_on_border( f, e, v, next_f, next_e ) ;
    }

    /*!
     * @brief Get the first facet of the surface that has an edge linking the two vertices (ids in the surface)
     *
     * @param[in] in0 Index of the first vertex in the surface
     * @param[in] in1 Index of the second vertex in the surface
     * @return NO_ID or the index of the facet
     */
    index_t Surface::facet_from_surface_vertex_ids(
        index_t in0,
        index_t in1 ) const
    {
        ringmesh_debug_assert( in0 < mesh_.vertices.nb() && in1 < mesh_.vertices.nb() ) ;

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for( index_t f = 0; f < nb_cells(); ++f ) {
            bool found = false ;
            index_t prev = surf_vertex_id( f, nb_vertices_in_facet( f ) - 1 ) ;
            for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
                index_t p = surf_vertex_id( f, v ) ;
                if( ( prev == in0 && p == in1 ) ||
                    ( prev == in1 && p == in0 ) )
                {
                    found = true ;
                    break ;
                }
                prev = p ;
            }
            if( found ) {
                return f ;
            }
        }
        return NO_ID ;
    }


    /*!
     * @brief Get the first facet of the surface that has an edge linking the
     * two vertices (ids in the model)
     *
     * @param[in] in0 Index of the first vertex in the model
     * @param[in] in1 Index of the second vertex in the model
     * @return NO_ID or the index of the facet
     */
    index_t Surface::facet_from_model_vertex_ids(
        index_t i0,
        index_t i1 ) const
    {
        index_t facet = NO_ID ;
        index_t edge = NO_ID ;
        edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
        return facet ;
    }


    /*!
     * @brief Determine the facet and the edge in this facet linking the 2 vertices
     * @details There might be two pairs facet-edge. This only gets the first.
     *
     * @param[in] i0 First vertex index in the model
     * @param[in] i1 Second vertex index in the model
     * @param[out] facet NO_ID or facet index in the surface
     * @param[out] edge NO_ID or edge index in the facet
     */
    void Surface::edge_from_model_vertex_ids(
        index_t i0,
        index_t i1,
        index_t& facet,
        index_t& edge ) const
    {
        edge = NO_ID ;

        // If a facet is given, look for the edge in this facet only
        if( facet != NO_ID ) {
            for( index_t v = 0; v < nb_vertices_in_facet( facet ); ++v ) {
                index_t prev = model_vertex_id( facet, prev_in_facet( facet, v )  ) ;
                index_t p = model_vertex_id( facet, v ) ;
                if( ( prev == i0 && p == i1 ) ||
                    ( prev == i1 && p == i0 ) )
                {
                    edge = prev_in_facet( facet, v ) ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) {return ;}
            }
        }

        // If we get here, no facet was found get out
        facet = NO_ID ;
        edge = NO_ID ;
    }


    /*!
     * @brief Determine the facet and the edge linking the 2 vertices with the same orientation
     * @details There might be two pairs facet-edge. This only gets the first.
     *
     * @param[in] i0 First vertex index in the model
     * @param[in] i1 Second vertex index in the model
     * @param[out] facet NO_ID or facet index in the surface
     * @param[out] edge NO_ID or edge index in the facet
     */
    void Surface::oriented_edge_from_model_vertex_ids(
        index_t i0,
        index_t i1,
        index_t& facet,
        index_t& edge ) const
    {
        // Copy from above .. tant pis
        edge = NO_ID ;

        // If a facet is given, look for the oriented edge in this facet only
        if( facet != NO_ID ) {
            for( index_t v = 0; v < nb_vertices_in_facet( facet ); ++v ) {
                index_t p = model_vertex_id( facet, v ) ;
                index_t next = model_vertex_id( facet, next_in_facet( facet, v )  ) ;

                if( p == i0 && next == i1 ) {
                    edge = v ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                oriented_edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) {return ;}
            }
        }
        facet = NO_ID ;
    }


    /*!
     * @brief Convert vertex surface index to an index in a facet
     * @param[in] f Index of the facet
     * @param[in] surf_vertex_id_in Index of the vertex in the surface
     * @return NO_ID or index of the vertex in the facet
     */
    index_t Surface::facet_vertex_id(
        index_t f,
        index_t surf_vertex_id_in ) const
    {
        for( index_t v = 0; v < nb_vertices_in_facet( f ); v++ ) {
            if( surf_vertex_id( f, v ) == surf_vertex_id_in ) {return v ;}
        }
        return NO_ID ;
    }


    /*!
     * @brief Convert model vertex index to an index in a facet
     * @param[in] f Index of the facet
     * @param[in] model_v_id Index of the vertex in the BoundaryModel
     * @return NO_ID or index of the vertex in the facet
     */
    index_t Surface::facet_id_from_model(
        index_t f,
        index_t model_v_id ) const
    {
        for( index_t v = 0; v < nb_vertices_in_facet( f ); v++ ) {
            if( model_vertex_id( f, v ) == model_v_id ) {return v ;}
        }
        return NO_ID ;
    }


    /*!
     * @brief Comparator of two vec3
     */
    struct comp_vec3bis {
        bool operator()(
            const vec3& l,
            const vec3& r ) const
        {
            if( l.x != r.x ) {return l.x < r.x ;}
            if( l.y != r.y ) {return l.y < r.y ;}
            return l.z < r.z ;
        }
    } ;

    /*!
     * @brief Determines the facets around a vertex
     *
     * @param[in] shared_vertex Index ot the vertex in the surface
     * @param[in] result Indices of the facets containing @param shared_vertex
     * @param[in] border_only If true only facets on the border are considered
     * @return The number of facet found
     *
     * \todo Evaluate if this is fast enough !!
     */
    index_t Surface::facets_around_vertex(
        index_t shared_vertex,
        std::vector< index_t >& result,
        bool border_only ) const
    {
        result.resize( 0 ) ;
        for( index_t t = 0; t < nb_cells(); ++t ) {
            for( index_t v = 0; v < nb_vertices_in_facet( t ); v++ ) {
                if( surf_vertex_id( t, v ) == shared_vertex ) {
                    return facets_around_vertex( shared_vertex, result,
                        border_only, t ) ;
                }
            }
        }
        ringmesh_assert_not_reached ;
        return dummy_index_t ;
    }


    /*!
     * @brief Determines the facets around a vertex
     *
     * @param[in] P Index ot the vertex in the surface
     * @param[in] result Indices of the facets containing @param P
     * @param[in] border_only If true only facets on the border are considered
     * @param[in] f0 Index of one facet containing the vertex @param P
     * @return The number of facet found
     *
     * \todo Evaluate if this is fast enough !!
     */
    index_t Surface::facets_around_vertex(
        index_t P,
        std::vector< index_t >& result,
        bool border_only,
        index_t f0 ) const
    {
        result.resize( 0 ) ;

        // Flag the visited facets
        std::vector< index_t > visited ;
        visited.reserve( 10 ) ;

        // Stack of the adjacent facets
        std::stack< index_t > S ;
        S.push( f0 ) ;
        visited.push_back( f0 ) ;

        do {
            index_t t = S.top() ;
            S.pop() ;

            for( index_t v = 0; v < nb_vertices_in_facet( t ); ++v ) {
                if( surf_vertex_id( t, v ) == P ) {
                    index_t adj_P = adjacent( t, v ) ;
                    index_t prev = prev_in_facet( t, v ) ;
                    index_t adj_prev = adjacent( t, prev ) ;

                    if( adj_P != NO_ADJACENT ) {
                        // The edge starting at P is not on the boundary
                        if( !Utils::contains( visited, adj_P ) ) {
                            S.push( adj_P ) ;
                            visited.push_back( adj_P ) ;
                        }
                    }
                    if( adj_prev != NO_ADJACENT ) {
                        // The edge ending at P is not on the boundary
                        if( !Utils::contains( visited, adj_prev ) ) {
                            S.push( adj_prev ) ;
                            visited.push_back( adj_prev  ) ;
                        }
                    }

                    if( border_only ) {
                        if( adj_P == NO_ADJACENT || adj_prev == NO_ADJACENT ) {
                            result.push_back( t ) ;
                        }
                    } else {result.push_back( t ) ;}

                    // We are done with this facet
                    break ;
                }
            }
        } while( !S.empty() ) ;

        return result.size() ;
    }


    /*!
     * @brief Compute the barycenter of a facet
     * @param[in] f Facet index in the surface
     * @return The coordinates of the facet barycenter
     */
    vec3 Surface::facet_barycenter( index_t f ) const
    {
        vec3 barycenter( 0., 0., 0. ) ;
        for( index_t i = 0; i < nb_vertices_in_facet( f ); i++ ) {
            barycenter += vertex( f, i ) ;
        }
        return barycenter / nb_vertices_in_facet( f ) ;
    }


    /*!
     * @brief Compute the area of a facet
     * @param[in] f Facet index in the surface
     * @return The area of the facet
     */
    double Surface::facet_area( index_t f ) const
    {
        double result = 0 ;
        for( index_t i = 1; i + 1 < nb_vertices_in_facet( f ); i++ ) {
            result += Utils::triangle_area(
                vertex( f, 0 ), vertex( f, i ), vertex( f, i + 1 ) ) ;
        }
        return result ;
    }


    /*!
     *
     * @param[in] f Facet index
     * @return Normal to the triangle made by the first 3 vertices
     * of the facet
     *
     * WARNING : if the facet is not planar calling this has no meaning
     */
    vec3 Surface::facet_normal( index_t f ) const
    {
        const vec3& p0 = vertex( f, 0 )  ;
        const vec3& p1 = vertex( f, 1 )  ;
        const vec3& p2 = vertex( f, 2 )  ;
        vec3 c0 = cross( p0 - p2, p1 - p2 ) ;
        return normalize( c0 ) ;
    }


    /*!
     * @brief Compute the normal to the surface vertices
     * @details The normal at a point is computed as the mean of the normal
     * to its adjacent facets.
     *
     * @param[out] normals Coordinates of the normal vectors to the vertices
     */
    void Surface::vertex_normals( std::vector< vec3 >& normals ) const
    {
        normals.resize( nb_vertices() ) ;
        for( index_t f = 0; f < nb_cells(); f++ ) {
            vec3 normal = facet_normal( f ) ;
            for( index_t p = 0; p < nb_vertices_in_facet( f ); p++ ) {
                index_t id = surf_vertex_id( f, p ) ;
                normals[ id ] += normal ;
            }
        }
        for( index_t p = 0; p < nb_vertices(); p++ ) {
            normals[ p ] = normalize( normals[ p ] ) ;
        }
    }


    /*!
     * @brief Compute closest vertex in a facet to a point
     * @param[in] f Facet index
     * @param[in] v Coordinates of the point to which distance is measured
     * @return Index of the vertex of @param f closest to @param v
     */
    index_t Surface::closest_vertex_in_facet(
        index_t f,
        const vec3& v ) const
    {
        index_t result = 0 ;
        double dist = DBL_MAX ;
        for( index_t p = 0; p < nb_vertices_in_facet( f ); p++ ) {
            double distance = length2( v - vertex( f, p ) ) ;
            if( dist > distance ) {
                dist = distance ;
                result = p ;
            }
        }
        return result ;
    }


    vec3 Surface::edge_barycenter( index_t c ) const
    {
        vec3 result( 0, 0, 0 ) ;

        // Get the facet index
        index_t f = NO_ID ;
        for( ; f < mesh_.facets.nb(); f++ ) {
            if( mesh_.facets.corners_begin(f ) < c  ) {
                break ;
            }
        }
        ringmesh_debug_assert( f != NO_ID ) ;
//        index_t f = narrow_cast< index_t >( std::lower_bound(
//            facet_ptr_.begin(), facet_ptr_.end(), c ) - facet_ptr_.begin() ) ;
        index_t v = c - facet_begin( f ) ;
        result += vertex( f, v ) ;
        result += vertex( f, next_in_facet( f, v ) ) ;
        return .5 * result ;
    }

    SurfaceTools::SurfaceTools( const Surface& surface )
        : surface_( surface ), aabb_( nil ), ann_( nil )
    {
    }
    SurfaceTools::~SurfaceTools()
    {
        if( aabb_ ) delete aabb_ ;
        if( ann_ ) delete ann_ ;
    }

    const GEO::MeshFacetsAABB& SurfaceTools::aabb() const
    {
        if( !aabb_ ) {
            SurfaceTools* this_not_const = const_cast< SurfaceTools* >( this ) ;
            this_not_const->aabb_ = new GEO::MeshFacetsAABB(
                const_cast< GEO::Mesh& >( surface_.mesh() ) ) ;
            surface_.model_->vertices.clear() ;
            if( ann_ ) {
                delete ann_ ;
                this_not_const->ann_ = nil ;
            }
        }
        return *aabb_ ;
    }

    const ColocaterANN& SurfaceTools::ann() const
    {
        if( !ann_ ) {
            const_cast< SurfaceTools* >( this )->ann_ = new ColocaterANN(
                surface_.mesh(), ColocaterANN::VERTICES ) ;
        }
        return *ann_ ;
    }



    /*!
     * @brief Cut the Surface along the line
     * @details First modify to NO_ADJACENT the neighbors the edges that are along the line
     * and then duplicate the points along this new boundary
     * Corners are not duplicated - maybe they should be in some cases but not in general..
     *
     * @param[in] L The Line
     */
    void SurfaceMutator::cut_by_line( const Line& L )
    {
        for( index_t i = 0; i + 1 < L.nb_vertices(); ++i ) {
            index_t p0 = L.model_vertex_id( i ) ;
            index_t p1 =
                ( i == L.nb_vertices() -
                  1 ) ? L.model_vertex_id( 0 ) : L.model_vertex_id(
                    i + 1 ) ;

            index_t f = Surface::NO_ID ;
            index_t v = Surface::NO_ID ;
            S_.edge_from_model_vertex_ids( p0, p1, f, v ) ;
            ringmesh_debug_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

            index_t f2 = S_.adjacent( f, v ) ;
            index_t v2 = Surface::NO_ID ;
            ringmesh_debug_assert( f2 != Surface::NO_ADJACENT ) ;
            S_.edge_from_model_vertex_ids( p0, p1, f2, v2 ) ;
            ringmesh_debug_assert( v2 != Surface::NO_ID ) ;

            // Virtual cut - set adjacencies to NO_ADJACENT
            S_.set_adjacent( f, v, Surface::NO_ADJACENT ) ;
            S_.set_adjacent( f2, v2, Surface::NO_ADJACENT ) ;
        }

        // Now travel on one side of the "faked" boundary and actually duplicate
        // the vertices in the surface
        // Get started in the surface - find (again) one of the edge that contains
        // the first two vertices of the line
        index_t f = Surface::NO_ID ;
        index_t v = Surface::NO_ID ;
        S_.oriented_edge_from_model_vertex_ids( L.model_vertex_id(
                0 ), L.model_vertex_id( 1 ), f, v ) ;
        ringmesh_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

        index_t id0 = S_.surf_vertex_id( f, v ) ;
        index_t id1 = S_.surf_vertex_id( f, S_.next_in_facet( f, v ) ) ;

        // Stopping criterion
        index_t last_vertex = L.model_vertex_id( L.nb_vertices() - 1 ) ;

        // Hopefully we have all the vertices on the Line..
        /// \todo Check that all vertices on the line are recovered
        while( S_.model_vertex_id( id1 ) != last_vertex ) {
            // Get the next vertex on the border
            // Same algorithm than in determine_line_vertices function
            index_t next_f = Surface::NO_ID ;
            index_t id1_in_next = Surface::NO_ID ;
            index_t next_id1_in_next = Surface::NO_ID ;

            // Get the next facet and next triangle on this boundary
            S_.next_on_border( f,
                S_.facet_vertex_id( f, id0 ), S_.facet_vertex_id( f, id1 ),
                next_f, id1_in_next, next_id1_in_next ) ;
            ringmesh_assert(
                next_f != Surface::NO_ID && id1_in_next != Surface::NO_ID
                && next_id1_in_next != Surface::NO_ID ) ;

            index_t next_id1 = S_.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Duplicate the vertex at id1
            // After having determined the next 1 we can probably get both at the same time
            // but I am lazy, and we must be careful not to break next_on_border function (Jeanne)
            std::vector< index_t > facets_around_id1 ;
            S_.facets_around_vertex( id1, facets_around_id1, false, f ) ;

            S_.mesh_.vertices.create_vertex( S_.vertex( id1 ).data() ) ;
            ringmesh_debug_assert( S_.nb_vertices() > 0 ) ;
            index_t new_id1 = S_.nb_vertices() - 1 ;

            for( index_t i = 0; i < facets_around_id1.size(); ++i ) {
                index_t cur_f = facets_around_id1[ i ] ;
                for( index_t cur_v = 0;
                     cur_v < S_.nb_vertices_in_facet( cur_f );
                     cur_v++ )
                {
                    if( S_.surf_vertex_id( cur_f, cur_v ) == id1 ) {
                        S_.mesh_.facets.set_vertex( cur_f, cur_v, new_id1 ) ;
                        break ;
                    }
                }
            }

            // Update
            f = next_f ;
            id0 = new_id1 ;
            id1 = next_id1 ;
        }

        /// \todo Check qu'on ne coupe pas compl�tement la surface, si on a 2 surfaces � la fin c'est la merde
    }


    /*!
     * @brief Compute the size (volume, area, length) of an Element
     *
     * @param[in] E Element to evaluate
     */
    double BoundaryModelElementMeasure::size( const BoundaryModelElement* E )
    {
        double result = 0. ;

        /// If this element has children sum up their sizes
        for( index_t i = 0; i < E->nb_children(); ++i ) {
            result += BoundaryModelElementMeasure::size( &E->child( i ) )  ;
        }
        if( result != 0 ) {return result ;}

        /// Else it is a base element and its size is computed

        // If this is a region
        if( E->element_type() == BoundaryModelElement::REGION ) {
            // Compute the volume if this is a region
            for( index_t i = 0; i < E->nb_boundaries(); i++ ) {
                const Surface& surface =
                    dynamic_cast< const Surface& >( E->boundary( i ) ) ;

                for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                    const vec3& p0 = surface.vertex( t, 0 ) ;
                    for( index_t v = 1;
                         v + 1 < surface.nb_vertices_in_facet( t );
                         ++v )
                    {
                        double cur_volume = ( dot( p0,
                                                  cross( surface.vertex( t,
                                                          v ),
                                                      surface.vertex( t, v + 1 ) ) ) )
                                            / static_cast< double >( 6 ) ;
                        E->side( i ) ? result -= cur_volume : result += cur_volume ;
                    }
                }
            }
            return fabs( result ) ;
        } else if( E->element_type() == BoundaryModelElement::CORNER ) {
            return 0 ;
        } else if( E->element_type() == BoundaryModelElement::LINE ) {
            const Line* L = dynamic_cast< const Line* >( E ) ;
            ringmesh_assert( L != nil ) ;
            for( index_t i = 1; i < E->nb_vertices(); ++i ) {
                result += GEO::Geom::distance( E->vertex( i ), E->vertex( i - 1 ) ) ;
            }
            return result ;
        } else if( E->element_type() == BoundaryModelElement::SURFACE ) {
            const Surface* S = dynamic_cast< const Surface* >( E ) ;
            ringmesh_assert( S != nil ) ;

            for( index_t i = 0; i < S->nb_cells(); i++ ) {
                result += S->facet_area( i ) ;
            }
            return result ;
        }
        ringmesh_assert_not_reached ;
        return result ;
    }


    /*!
     * @brief Compute the barycenter of a part of a BoundaryModelElement
     * Only implemented for Surface and Line
     *
     * @param[in] E Pointer to the element
     * @param[in] cells Indices of the segments/facets to consider
     * @return The coordinates of the barycenter of the @param cells
     */
    vec3 BoundaryModelElementMeasure::barycenter(
        const BoundaryModelElement* E,
        const std::vector< index_t >& cells )
    {
        vec3 result( 0., 0., 0. ) ;
        double size = 0 ;

        const Line* L = dynamic_cast< const Line* >( E ) ;
        if( L != nil ) {
            for( index_t i = 0; i < cells.size(); ++i ) {
                result += L->segment_length( cells[ i ] ) * L->segment_barycenter(
                    cells[ i ] ) ;
                size   += L->segment_length( cells[ i ] ) ;
            }
            return size > epsilon ? result / size : result ;
        }
        const Surface* S = dynamic_cast< const Surface* >( E ) ;
        if( S != nil ) {
            for( index_t i = 0; i < cells.size(); ++i ) {
                result += S->facet_area( cells[ i ] ) *
                          S->facet_barycenter( cells[ i ] ) ;
                size   += S->facet_area( cells[ i ] ) ;
            }
            return size > epsilon ? result / size : result ;
        }
        ringmesh_assert_not_reached ;
        return result ;
    }


    /*!
     * @brief Measures the minimal distance between an element and a point
     * Implement only for Surface, Line and Corner
     *
     * @param[in] E Pointer to the element
     * @param[in] p Coordinates of the point to which distance is measured
     */
    double BoundaryModelElementMeasure::distance(
        const BoundaryModelElement* E,
        const vec3& p )
    {
        double result = FLT_MAX ;
        const Line* L = dynamic_cast< const Line* >( E ) ;
        if( L != nil ) {
            for( index_t i = 1; i < L->nb_vertices(); ++i ) {
                // Distance between a vertex and a segment
                const vec3& p0 = L->vertex( i - 1 ) ;
                const vec3& p1 = L->vertex( i ) ;

                double distance_pt_2_segment  = FLT_MAX ;
                vec3 c = ( p1 - p0 ) / 2 ;
                double half = GEO::distance( p1, c ) ;
                double cp_dot_p0p1 = dot( p - c, p1 - p0 ) ;

                if( cp_dot_p0p1 < - half ) {
                    distance_pt_2_segment =  GEO::distance(
                        p0,
                        p ) ;
                } else if( cp_dot_p0p1 >
                           half )
                {
                    distance_pt_2_segment = GEO::distance( p1, p ) ;
                } else {
                    vec3 projection = c + cp_dot_p0p1 * ( p1 - p0 ) ;
                    distance_pt_2_segment = GEO::distance( projection, p ) ;
                }
                result = distance_pt_2_segment <
                         result ? distance_pt_2_segment : result ;
            }
            return result ;
        }
        const Surface* S = dynamic_cast< const Surface* >( E ) ;
        if( S != nil ) {
            for( index_t i = 0; i < S->nb_cells(); i++ ) {
                for( index_t j = 1; j + 1 < S->nb_vertices_in_facet( i ); ++j ) {
                    double cur = GEO::Geom::point_triangle_squared_distance(
                        p, S->vertex( i, 0 ), S->vertex( i, j ), S->vertex( i, j + 1 ) ) ;
                    if( cur < result ) {result = cur ;}
                }
            }
            if( result != FLT_MAX ) {result = sqrt( result ) ;}
            return result ;
        }

        const Corner* C = dynamic_cast < const Corner* >( E ) ;
        if( C != nil ) {
            return GEO::distance( C->vertex(), p ) ;
        }

        // If it is not one of the basic types - compute it for the children
        // if any
        if( E->nb_children() == 0 ) {
            ringmesh_assert_not_reached ;
            return result ;
        } else {
            for( index_t i = 0; i < E->nb_children(); ++i ) {
                double dist = distance( &E->child( i ), p ) ;
                result = ( dist < result ) ? dist : result ;
            }
            return result ;
        }
    }
}
