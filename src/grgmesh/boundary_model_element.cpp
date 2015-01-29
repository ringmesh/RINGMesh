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

/*! \author Jeanne Pellerin and Arnaud Botella */

#include <grgmesh/boundary_model_element.h>
#include <grgmesh/boundary_model.h>
#include <grgmesh/utils.h>

#include <geogram/basic/geometry_nd.h>

#include <set>
#include <stack>
#include <fstream>


namespace GRGMesh {

    /*!
     *
     * @return Assert that the parent exist and returns it.
     */
    const BoundaryModelElement& BoundaryModelElement::parent() const
    {
        grgmesh_debug_assert( has_parent() ) ;
        switch( element_type() ) {
            case BM_LINE :
                return model_->contact( parent_ ) ;
            case BM_SURFACE :
                return model_->one_interface( parent_ ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    /*!
     *
     * @param[in] x Index of the boudnary element
     * @return Assert that is exits and return the element on the boundary
     */
    const BoundaryModelElement& BoundaryModelElement::boundary( index_t x ) const
    {
        grgmesh_debug_assert( x < nb_boundaries() ) ;
        index_t id = boundaries_.at(x) ;
        switch( element_type() ) {
            case BM_LINE     : return model_->corner( id ) ;
            case BM_SURFACE  : return model_->line( id ) ;
            case BM_INTERFACE: return model_->contact( id ) ;
            case BM_REGION   : return model_->surface( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    /*!
     *
     * @param[in] x Index of the in_boundary element
     * @return Assert that it exist and return the element in in_boundary.
     */
    const BoundaryModelElement& BoundaryModelElement::in_boundary( index_t x ) const
    {
        grgmesh_debug_assert( x < nb_in_boundary() ) ;
        index_t id = in_boundary_.at(x) ;
        switch( element_type() ) {
            case BM_CORNER      : return model_->line( id ) ;
            case BM_LINE        : return model_->surface( id ) ;
            case BM_CONTACT     : return model_->one_interface( id ) ;
            case BM_SURFACE     : return model_->region( id ) ;
            case BM_INTERFACE   : return model_->layer( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    /*!
     *
     * @param[in] x Index of the child
     * @return Assert that the chil exist and return it.
     */
    const BoundaryModelElement& BoundaryModelElement::child( index_t x ) const
    {
        grgmesh_debug_assert( x < nb_children() ) ;
        index_t id = children_.at(x) ;
        switch( element_type() ) {
            case BM_CONTACT     : return model_->line( id ) ;
            case BM_INTERFACE   : return model_->surface( id ) ;
            case BM_LAYER       : return model_->region( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    /*!
     * @brief Copy all attributes (except the model_) from \param rhs to this element
     * @param[in] rhs To copy from
     * @param[in] model Model to associate to this element
     */
    void BoundaryModelElement::copy_macro_topology(
        const BoundaryModelElement& rhs,
        BoundaryModel& model )
    {
        model_        = &model ;
        name_         = rhs.name_ ;
        id_           = rhs.id_ ;
        type_         = rhs.type_;
        geol_feature_ = rhs.geol_feature_;
        parent_       = rhs.parent_ ;
        boundaries_   = rhs.boundaries_ ;
        sides_        = rhs.sides_ ;
        in_boundary_  = rhs.in_boundary_ ;
        children_     = rhs.children_ ;
    }
  
    /*!
     * @brief Checks if this element or one of the element containing it 
     * determines the model Volume Of Interest
     * @details This is known with the type of an element
     */
    bool BoundaryModelElement::is_on_voi() const
    {
        if( geol_feature_ == ALL ) 
        {
            for( index_t j = 0; j < nb_in_boundary(); ++j ) {
                GEOL_FEATURE t = in_boundary( j ).geological_feature() ;
                if( t == VOI || t == STRATI_VOI || t == FAULT_VOI ){
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
    
    /*!
     *
     * @return The coordinates of the point at this corner.
     */
    const vec3& Corner::vertex( index_t p ) const
    {
        return model_->vertex( vertex_ ) ;
    }


    /*!
     * @brief Construct a Line
     * 
     * @param[in] model The parent model
     * @param[in] id The index of the line in the lines_ vector of the parent model
     */
    Line::Line( BoundaryModel* model, index_t id ):
        BoundaryModelElement( model, BM_LINE, id )
    { 
        boundaries_.resize( 2, nil) ; 
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
        const std::vector< index_t >& vertices )
        : BoundaryModelElement( model, BM_LINE, id ), vertices_( vertices )
    {
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
        const std::vector< index_t >& vertices
    ):  BoundaryModelElement( model, BM_LINE, id ),
        vertices_( vertices )
    {
        boundaries_.push_back( corner0 ) ;
        boundaries_.push_back( corner1 ) ;
    } ;
    
    /*!
     * @brief Check is the Line is twice on the boundary of a surface
     *
     * @param[in] surface The surface to test
     */
    bool Line::is_inside_border(
        const BoundaryModelElement& surface ) const
    {
        // Find out if this surface is twice in the in_boundary vector
        return std::count( in_boundary_.begin(), in_boundary_.end(), surface.id() ) > 1 ;
    }

    /*!
     * @return The coordinates of the \param line_vertex_id th vertex on the Line
     */
    const vec3& Line::vertex( index_t line_vertex_id ) const {
        return model_->vertex( vertices_.at( line_vertex_id ) ) ;
    }

    /*!
     * 
     * @param[in] s Segment index
     * @return The coordinates of the barycenter of the segment
     */
    vec3 Line::segment_barycenter( index_t s ) const {
        return 0.5*( vertex(s) + vertex(s+1) ) ;
    }
    
    /*!
     *
     * @param[in] s Segment index
     * @return The length of the segment
     */
    double Line::segment_length( index_t s ) const {
        return length( vertex(s+1)-vertex(s) ) ;
    }

    /*!
     *
     * @return The length of the Line
     */
    double Line::total_length() const {
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
    bool Line::equal( const std::vector< index_t >& rhs_vertices ) const {        
        if( nb_vertices() != rhs_vertices.size() ) return false ; 
        
        if( std::equal( rhs_vertices.begin(), rhs_vertices.end(), vertices_.begin() ) )
            return true ;
            
        if( std::equal( rhs_vertices.begin(), rhs_vertices.end(), vertices_.rbegin() ) )
            return true ;
            
        return false ;
    }
   
    /*!
     * \todo Check LineMutato function implementation and comment them
     * @param id
     * @param p
     *
    void LineMutator::set_vertex( index_t id, const vec3& p ) {
        M_.model_->vertices_[M_.vertices_[id]] = p ;
    }

    /*!
     *
     * @param p
     * @return
     *
    vec3& LineMutator::vertex( index_t p ) const
    {
        return M_.model_->vertices_[ M_.vertices_[p] ] ;
    }


    /*!
     * @param[in] f Facet index
     * @param[in] v Vertex index in the facet
     * @return The coordinates of the vertex
     */
    const vec3& Surface::vertex( index_t f, index_t v ) const
    {
        grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
        return vertex( surf_vertex_id( f, v) ) ;
    }

    /*!
     * @param[in] surf_vertex_id Index of the vertex in the surface
     * @return The coordinates of the vertex
     */
    const vec3& Surface::vertex( index_t surf_vertex_id ) const {
        return model_->vertex( vertices_.at(surf_vertex_id) ) ;
    }

    /*!
     * @param[in] f Facet index in the surface
     * @return Facet index in the parent model
     */
    index_t Surface::model_facet_id( index_t f ) const {
        return model_->model_facet( id_, f ) ;
    }

    /*!
     * @brief Initialize the KeyFacet to be the first 3 vertices of the surface
     */
    void Surface::set_first_triangle_as_key()
    {
        // I guess it should'nt be a problem if the first facet is not a triangle
        // that's only a guess (Jeanne)        
        key_facet_ = KeyFacet( model_->vertex( vertices_[facets_[0]] ),
            model_->vertex( vertices_[facets_[1]] ),
            model_->vertex( vertices_[facets_[2]] ) ) ;
    }

         
    /*!
     * @brief Traversal of a surface border
     * @details From the input facet @param f, get the facet that share vertex @param v and 
     * get the indices of vertex @param v and of the following vertex in this new facet.
     * The next facet @param next_f may be the same, and @param is required to avoid going back.
     *
     * @param[in] f Index of the facet
     * @param[in] from Index in the facet of the previous point on the border - gives the direction
     * @param[in] v Index in the facet of the point for which we want the next point on border
     * @param[out] next_f Index of the facet containing the next point on border
     * @param[out] v_in_next Index of vertex @param v in facet @param next_f
     * @param[out] next_in_next Index of the next vertex on border in facet @param v_in_next
     */
    void Surface::next_on_border( 
        index_t f, index_t from, index_t v, 
        index_t& next_f, index_t& v_in_next, index_t& next_in_next ) const
    {
        grgmesh_assert( v < nb_vertices_in_facet( f ) ) ;
        grgmesh_assert( is_on_border( f, v ) || is_on_border( f, from ) ) ;

        index_t V = surf_vertex_id( f, v ) ;

        // We want the next triangle that is on the boundary and share V
        // If there is no such triangle, the next vertex on the boundary 
        // is the vertex of F neighbor of V that is not from         

        // Get the facets around the shared vertex that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)        
        std::vector< index_t > facets ;  
        index_t nb_around = facets_around_vertex( V, facets, true, f ) ;
        grgmesh_assert( nb_around < 3 && nb_around > 0 ) ;

        next_f = facets[0] ;
        
        if( nb_around == 2 ) {
            if( next_f == f ) next_f = facets[1] ;
            grgmesh_debug_assert( next_f != NO_ID ) ;

            // Now get the other vertex that is on the boundary opposite to p1
            v_in_next = facet_vertex_id( next_f, V ) ;
            grgmesh_assert( v_in_next != NO_ID ) ;

            // The edges containing V in next_f are
            // the edge starting at v_in_next and the one ending there
            index_t prev_v_in_next = prev_in_facet( next_f, v_in_next )  ;

            bool e0_on_boundary = is_on_border( next_f, v_in_next ) ;
            bool e1_on_boundary = is_on_border( next_f, prev_v_in_next ) ;

            // Only one must be on the boundary otherwise there is a corner missing
            grgmesh_assert( e0_on_boundary != e1_on_boundary ) ;

            // From the edge that is on boundary get the next vertex on this boundary
            // If the edge starting at p_in_next is on boundary, new_vertex is its next
            // If the edge ending at p_in_next is on boundary, new vertex is its prev
            next_in_next = e0_on_boundary ? next_in_facet( next_f, v_in_next ) : prev_v_in_next ;

        } else if( nb_around == 1 ) {
            // V must be in two border edges of facet f 
            // Get the id in the facet of the vertex neighbor of v1 that is not v0           
            v_in_next = v ;
            if( prev_in_facet( f, v ) == from  ){
                grgmesh_debug_assert( is_on_border(f, v) ) ;
                next_in_next = next_in_facet( f, v ) ;
            }
            else {
                grgmesh_debug_assert( is_on_border( f, prev_in_facet(f,v) ) ) ;
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
    void Surface::next_on_border( index_t f, index_t e, index_t& next_f, index_t& next_e ) const {
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
        index_t in0, index_t in1 ) const 
    {
        grgmesh_debug_assert( in0 < vertices_.size() && in1 < vertices_.size() ) ;
        
        // Another possible, probably faster, algorith is to check if the 2 indices 
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for( index_t f = 0; f < nb_cells(); ++f ) {
            bool found = false ;
            index_t prev = surf_vertex_id( f, nb_vertices_in_facet(f)-1 ) ; 
            for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
                index_t p = surf_vertex_id( f, v ) ;
                if( (prev == in0 && p == in1) ||
                    (prev == in1 && p == in0) ) 
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
    index_t Surface::facet_from_model_vertex_ids( index_t i0, index_t i1 ) const {
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
                if( ( prev == i0 && p == i1) ||
                    ( prev == i1 && p == i0 ) ) 
                {
                    edge = prev_in_facet( facet, v ) ;
                    return ;
                }
            }
        }
        else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) return ;
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

                if( p == i0 && next == i1 )
                {
                    edge = v ;
                    return ;
                }
            }
        }
        else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                oriented_edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) return ;
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
    index_t Surface::facet_vertex_id( index_t f, index_t surf_vertex_id_in ) const {
        for( index_t v = 0; v < nb_vertices_in_facet(f); v++ ) {
            if( surf_vertex_id( f, v ) == surf_vertex_id_in ) return v ;
        }
        return NO_ID ;
    }

    /*!
     * @brief Comparator of two vec3
     * 
     * This is a copy, but from where ?
     */
    struct comp_vec3bis {
        bool operator()( const vec3& l, const vec3& r ) const {
            if( l.x != r.x ) return l.x < r.x ;
            if( l.y != r.y ) return l.y < r.y ;
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
        result.resize(0) ;
        for( index_t t = 0; t < nb_cells(); ++t ) {            
            for( index_t v = 0; v < nb_vertices_in_facet(t); v++ ) {
                if( surf_vertex_id( t, v ) == shared_vertex ) {
                    return facets_around_vertex( shared_vertex, result,
                        border_only, t ) ;
                }
            }
        }
        grgmesh_assert_not_reached ;
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

            for( index_t v = 0; v < nb_vertices_in_facet(t); ++v ) {               
                if( surf_vertex_id( t, v ) == P ) 
                {
                    index_t adj_P = adjacent( t, v ) ;
                    index_t prev = prev_in_facet( t, v ) ;
                    index_t adj_prev = adjacent( t, prev ) ;

                    if( adj_P != NO_ADJACENT ){
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
                    }
                    else result.push_back( t ) ;
                    
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
    vec3 Surface::facet_barycenter( index_t f ) const {
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
    double Surface::facet_area( index_t f ) const {
        double result = 0 ;
        for( index_t i = 1; i+1 < nb_vertices_in_facet( f ); i++ ) 
        { 
            result += Utils::triangle_area( 
                vertex( f, 0 ), vertex( f, i ), vertex( f, i+1 ) ) ;
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
    vec3 Surface::facet_normal( index_t f ) const {
        const vec3& p0 = vertex( f, 0 )  ;
        const vec3& p1 = vertex( f, 1 )  ;
        const vec3& p2 = vertex( f, 2 )  ;
        vec3 c0 = cross(p0-p2, p1-p2) ;       
        return normalize( c0 ) ;
    }

    /*!
     * @brief Compute the normal to the surface vertices
     * @details The normal at a point is computed as the mean of the normal
     * to its adjacent facets.
     * 
     * @param[out] normals Coordinates of the normal vectors to the vertices
     */
    void Surface::vertex_normals( std::vector< vec3 >& normals ) const {
        normals.resize( nb_vertices() ) ;
        for( index_t f = 0; f < nb_cells(); f++ ) {
            vec3 normal = facet_normal( f ) ;
            for( index_t p = 0; p < nb_vertices_in_facet( f ); p++ ) {
                index_t id = surf_vertex_id( f, p ) ;
                normals[id] += normal ;
            }
        }
        for( index_t p = 0; p < nb_vertices(); p++ ) {
            normals[p] = normalize( normals[p] ) ;
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


    /*!
     * \todo Check the SurfaceMutator code and comment it 
     * @param id
     * @param p
     *
    void SurfaceMutator::set_vertex( index_t id, const vec3& p ) {
        M_.model_->vertices_[M_.vertices_[id]] = p ;
    }

    /*!
     *
     * @param p
     * @return
     *
    vec3& SurfaceMutator::vertex( index_t p ) const
    {
        return M_.model_->vertices_[ M_.vertices_[p] ] ;
    }


    /*!
     * @brief Cut the Surface along the line 
     * @details First modify to NO_ADJACENT the neighbors the edges that are along the line 
     * and then duplicate the points along this new boundary
     * Corners are not duplicated - maybe they should be in some cases but not in general..
     * 
     * @param[in] L The Line
     */
    void SurfaceMutator::cut_by_line( const Line& L ) {        
        for( index_t i= 0; i+1 < L.nb_vertices(); ++i ) {
            index_t p0 = L.model_vertex_id( i ) ;
            index_t p1 =  (i == L.nb_vertices()-1) ? L.model_vertex_id(0) : L.model_vertex_id( i+1 ) ;

            index_t f = Surface::NO_ID ;
            index_t v = Surface::NO_ID ;
            S_.edge_from_model_vertex_ids(p0, p1, f, v) ;
            grgmesh_debug_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

            index_t f2 = S_.adjacent( f, v ) ;
            index_t v2 = Surface::NO_ID ;
            grgmesh_debug_assert( f2 != Surface::NO_ADJACENT ) ;
            S_.edge_from_model_vertex_ids( p0, p1, f2, v2 ) ;
            grgmesh_debug_assert( v2 != Surface::NO_ID ) ;

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
        S_.oriented_edge_from_model_vertex_ids( L.model_vertex_id( 0 ), L.model_vertex_id( 1 ), f, v ) ;
        grgmesh_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

        index_t id0 = S_.surf_vertex_id( f, v ) ;
        index_t id1 = S_.surf_vertex_id( f, S_.next_in_facet(f,v) ) ;
        
        // Stopping criterion
        index_t last_vertex = L.model_vertex_id( L.nb_vertices()-1 ) ;
        // Hopefully we have all the vertices on the Line.. 
        /// \todo Check that all vertices on the line are recovered
        while( S_.model_vertex_id( id1 ) != last_vertex ) {
            // Get the next vertex on the border 
            // Same algorithm than in determine_line_vertices function
            index_t next_f = Surface::NO_ID ;
            index_t id1_in_next = Surface::NO_ID ;
            index_t next_id1_in_next = Surface::NO_ID ;

            // Get the next facet and next triangle on this boundary
            S_.next_on_border( f, S_.facet_vertex_id(f, id0), S_.facet_vertex_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;
            grgmesh_assert(
                next_f != Surface::NO_ID && id1_in_next != Surface::NO_ID
                    && next_id1_in_next != Surface::NO_ID ) ;
            
            index_t next_id1 = S_.surf_vertex_id( next_f, next_id1_in_next ) ;
            
            // Duplicate the vertex at id1 
            // After having determined the next 1 we can probably get both at the same time 
            // but I am lazy, and we must be careful not to break next_on_border function (Jeanne)
            std::vector< index_t > facets_around_id1 ;
            S_.facets_around_vertex( id1, facets_around_id1, false, f ) ;

            S_.vertices_.push_back( S_.model_vertex_id(id1) ) ;
            grgmesh_debug_assert( S_.nb_vertices() > 0 ) ;
            index_t new_id1 = S_.nb_vertices()-1 ;
            
            for( index_t i = 0; i < facets_around_id1.size(); ++i ){
                index_t cur_f = facets_around_id1[i] ;
                for( index_t cur_v = 0; cur_v < S_.nb_vertices_in_facet( cur_f ) ; cur_v++ )
                {
                    if( S_.surf_vertex_id( cur_f, cur_v ) == id1 ) {
                        S_.facets_[ S_.facet_begin( cur_f ) + cur_v ] = new_id1 ;
                        break ;
                    }
                }
            }
            // Update
            f = next_f ;
            id0 = new_id1 ;
            id1 = next_id1 ;
        }
        /// \todo Check qu'on ne coupe pas complètement la surface, si on a 2 surfaces à la fin c'est la merde
    }




    /*!
     * @brief Compute the size (volume, area, length) of an Element
     * 
     * @param[in] E Element to evaluate
     */
    double BoundaryModelElementMeasure::size( const BoundaryModelElement* E ) {        
        double result = 0. ;

        /// If this element has children sum up their sizes
        for( index_t i = 0; i < E->nb_children(); ++i ){
            result += BoundaryModelElementMeasure::size( &E->child(i) )  ;
        }
        if( result != 0 ) return result ;

        /// Else it is a base element and its size is computed
        
        // If this is a region 
        if( E->element_type() == BoundaryModelElement::BM_REGION ) {
            // Compute the volume if this is a region
            for( index_t i = 0; i < E->nb_boundaries(); i++ ) {
                const Surface& surface = dynamic_cast< const Surface& >( E->boundary( i ) ) ;
                   
                for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                    const vec3& p0 = surface.vertex(t, 0 ) ;
                    for( index_t v = 1; v+1 < surface.nb_vertices_in_facet(t); ++v ){
                        double cur_volume = ( dot( p0,
                            cross( surface.vertex( t, v ), surface.vertex( t, v+1 ) ) ) )
                            / static_cast< double >( 6 ) ;
                        E->side(i) ? result -= cur_volume : result += cur_volume ;
                    }
                }
            }
            return fabs( result ) ;
        }
        else if( E->element_type() == BoundaryModelElement::BM_CORNER ) {
            return 0 ; 
        }
        else if( E->element_type() == BoundaryModelElement::BM_LINE ) {
            const Line* L = dynamic_cast< const Line* >( E ) ;
            grgmesh_assert( L != nil ) ;           
            for( index_t i = 1; i < E->nb_vertices(); ++i ){
                result += GEO::Geom::distance( E->vertex( i ), E->vertex( i-1 ) ) ;
            }
            return result ;
        }
        else if( E->element_type() == BoundaryModelElement::BM_SURFACE ) {
            const Surface* S = dynamic_cast< const Surface* >( E ) ;
            grgmesh_assert( S != nil ) ;

            for( index_t i = 0; i < S->nb_cells(); i++ ) {
                result += S->facet_area( i ) ;
            }
            return result ;
        }
        grgmesh_assert_not_reached ;
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
    vec3 BoundaryModelElementMeasure::barycenter ( 
        const BoundaryModelElement* E, const std::vector< index_t >& cells ) 
    {
        vec3 result(0., 0., 0. ) ;
        double size = 0 ;

        const Line* L = dynamic_cast< const Line* >( E ) ;
        if( L != nil ) {
            for( index_t i = 0; i < cells.size(); ++ i ) {
                result += L->segment_length( cells[i] ) * L->segment_barycenter( cells[i] ) ;
                size   += L->segment_length( cells[i] ) ;
            }
            return size > epsilon ? result/size : result ;
        }
        const Surface* S = dynamic_cast< const Surface* >( E ) ;
        if( S != nil ) {
            for( index_t i = 0; i < cells.size(); ++ i ) {
                result += S->facet_area( cells[i] ) * S->facet_barycenter( cells[i] ) ;
                size   += S->facet_area( cells[i] ) ;
            }
            return size > epsilon ? result/size : result ;
        }        
        grgmesh_assert_not_reached ;
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
            for( index_t i = 1; i < L->nb_vertices(); ++i ){
                // Distance between a vertex and a segment
                const vec3& p0 = L->vertex( i-1 ) ;
                const vec3& p1 = L->vertex( i ) ;

                double distance_pt_2_segment  = FLT_MAX ;
                vec3 c = (p1-p0)/2 ;        
                double half = GEO::distance( p1, c ) ;
                double cp_dot_p0p1 = dot( p-c, p1-p0 ) ;

                if( cp_dot_p0p1 < -half ) distance_pt_2_segment =  GEO::distance( p0, p ) ;
                else if( cp_dot_p0p1 > half ) distance_pt_2_segment = GEO::distance( p1, p ) ;
                else {
                    vec3 projection = c + cp_dot_p0p1*(p1-p0) ;
                    distance_pt_2_segment = GEO::distance( projection, p ) ;
                }    
                result = distance_pt_2_segment < result ? distance_pt_2_segment : result ;
            }        
            return result ;
        }
        const Surface* S = dynamic_cast< const Surface* >( E ) ;
        if( S != nil ) {
            for( index_t i = 0; i < S->nb_cells(); i++ ) {
                for( index_t j = 1; j+1 < S->nb_vertices_in_facet( i ); ++j ) {
                    double cur = GEO::Geom::point_triangle_squared_distance( 
                        p, S->vertex( i, 0 ), S->vertex( i, j ), S->vertex( i, j+1 ) ) ;
                    if( cur < result ) result = cur ;                               
                }
            }
            if( result != FLT_MAX ) result = sqrt( result ) ;
            return result ;           
        }

        const Corner* C = dynamic_cast < const Corner* >( E ) ;
        if( C != nil ) {
            return GEO::distance( C->vertex(), p ) ;
        }

        // If it is not one of the basic types - compute it for the children
        // if any 
        if( E->nb_children() == 0 ) {
            grgmesh_assert_not_reached ;
            return result ;
        }
        else {
            for( index_t i = 0; i < E->nb_children(); ++i ){
                double dist = distance( &E->child(i), p ) ;
                result = ( dist < result ) ? dist : result ;
            }
            return result ;
        }       
    }

}

