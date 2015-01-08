/*! This file is originally part of the Geomodeling plugin of Graphite 
 *  The initial version was developed by Jeanne Pellerin and modified by Arnaud Botella
 *  
 *  This version is modified by Jeanne Pellerin - WIAS
 */

/*! \author Jeanne Pellerin */

#include <grgmesh/boundary_model_element.h>
#include <grgmesh/boundary_model.h>
#include <grgmesh/utils.h>

#include <geogram/basic/geometry_nd.h>
#include <set>
#include <stack>
#include <fstream>

//#define _USE_MATH_DEFINES // Otherwise M_PI not defined interfere with other math stuff ? Jeanne
//#include <math.h>

namespace GRGMesh {


/********************************************************************************************/
/******             BoundaryModelElement implementation     ***********************************/
/********************************************************************************************/

    const BoundaryModelElement& BoundaryModelElement::parent() const
    {
        grgmesh_debug_assert( has_parent() ) ;
        switch( dim() ) {
            case 1:
                return model_->contact( parent_ ) ;
            case 2:
                return model_->one_interface( parent_ ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::boundary( index_t x ) const
    {
        grgmesh_debug_assert( x < nb_boundaries() ) ;
        index_t id = boundaries_.at(x) ;
        switch( dim() ) {
            case 1:
                return model_->corner( id ) ;
            case 2:
                if( has_parent() ) return model_->line( id ) ;
                return model_->contact( id ) ;
            case 3:
                return model_->surface( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::in_boundary( index_t x ) const
    {
        grgmesh_debug_assert( x < nb_in_boundary() ) ;
        index_t id = in_boundary_.at(x) ;
        switch( dim() ) {
            case 0:
                return model_->line( id ) ;
            case 1:
                if( has_parent() ) return model_->surface( id ) ;
                return model_->one_interface( id ) ;
            case 2:
                return model_->region( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::child( index_t x ) const
    {
        grgmesh_debug_assert( !has_parent() && x < nb_children() ) ;
        index_t id = children_.at(x) ;
        switch( dim() ) {
            case 1:
                return model_->line( id ) ;
            case 2:
                return model_->surface( id ) ;
            case 3:
                return model_->region( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

    void BoundaryModelElement::copy_macro_topology(
        const BoundaryModelElement& rhs,
        BoundaryModel& model )
    {
        model_ = &model ;
        name_ = rhs.name_ ;
        id_ = rhs.id_ ;
        dim_ = rhs.dim_ ;
        type_ = rhs.type_ ;
        parent_ = rhs.parent_ ;
        boundaries_ = rhs.boundaries_ ;
        sides_ = rhs.sides_ ;
        in_boundary_ = rhs.in_boundary_ ;
        children_ = rhs.children_ ;
    }
  
     /*! Returns true if this element or one of the element containing it
     *  is on the Volume Of Interest
     *  This info is strored in the type of the element
     */
    bool BoundaryModelElement::is_on_voi() const
    {
        if( type_ == ALL ) 
        {
            for( index_t j = 0; j < nb_in_boundary(); ++j ) {
                GEOL_FEATURE t = in_boundary( j ).type() ;
                if( t == VOI || t == STRATI_VOI || t == FAULT_VOI ){
                    return true ;
                }
            }
        } else if( type_ == VOI || type_ == STRATI_VOI || type_ == FAULT_VOI )
        {
            return true ;
        }
        return false ;
    }
    

/***********************************************************************************************/

    const vec3& Corner::vertex( index_t p ) const
    {
        return model_->vertex( p_ ) ;
    }

    /*void Corner::copy_macro_topology(
        const Corner& rhs,
        BoundaryModel& model )
    {
        BoundaryModelElement::copy_macro_topology( rhs, model ) ;
    }*/
    
    
                

/********************************************************************************************/
/******             Line implementation            ***********************************/
/********************************************************************************************/

    Line::Line( BoundaryModel* model, index_t id ):
        BoundaryModelElement( model, 1, id )
    { 
        boundaries_.resize( 2, nil) ; 
    }

    Line::Line(
        BoundaryModel* model,
        index_t id,
        const std::vector< index_t >& vertices )
        : BoundaryModelElement( model, 1, id ), vertices_( vertices )
    {
    }

    Line::Line(
        BoundaryModel* model,
        index_t id,
        index_t corner0,
        index_t corner1,
        const std::vector< index_t >& vertices
    ):  BoundaryModelElement( model, 1, id ),
        vertices_( vertices )
    {
        boundaries_.push_back( corner0 ) ;
        boundaries_.push_back( corner1 ) ;
    } ;

   /* void Line::copy_macro_topology(
        const Line& rhs,
        BoundaryModel& model )
    {
        BoundaryModelElement::copy_macro_topology( rhs, model ) ;
       // is_inside_border_ = rhs.is_inside_border_ ;
    }*/
 
    /*void Line::add_in_boundary( index_t e )
    {
        for( index_t i = 0; i < nb_in_boundary(); i++ ) {
            if( in_boundary_[i] == e ) {
                assert( !is_inside_border_[i] ) ;
                is_inside_border_[i] = true ;
                return ;
            }
        }
        BoundaryModelElement::add_in_boundary( e ) ;
        is_inside_border_.push_back( false ) ;
    }*/
    
    bool Line::is_inside_border(
        const BoundaryModelElement& surface ) const
    {
        // Find out if this surface is twice in the in_boundary vector
        return std::count( in_boundary_.begin(), in_boundary_.end(), surface.id() ) > 1 ;
    }

    const vec3& Line::vertex( index_t line_vertex_id ) const {
        return model_->vertex( vertices_.at( line_vertex_id ) ) ;
    }

    vec3 Line::segment_barycenter( index_t s ) const {
        return 0.5*( vertex(s) + vertex(s+1) ) ;
    }
    
    double Line::segment_length( index_t s ) const {
        return length( vertex(s+1)-vertex(s) ) ;
    }
    double Line::total_length() const {
        double result = 0 ;
        for( index_t s = 0; s < nb_cells(); s++ ) {
            result += segment_length( s ) ;
        }
        return result ;
    }

   /* bool Line::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }
    signed_index_t Line::find( const vec3& p ) const
    {
        for( index_t i = 0; i < vertices_.size(); ++i ) {
            if( vertex( i ) == p ) return i ;
        }
        return -1 ;
    } */

    
    void LineMutator::set_vertex( index_t id, const vec3& p ) {
        M_.model_->vertices_[M_.vertices_[id]] = p ;
    }

    vec3& LineMutator::vertex( index_t p ) const
    {
        return M_.model_->vertices_[ M_.vertices_[p] ] ;
    }

/********************************************************************************************/
/******             Surface implementation            ***********************************/
/********************************************************************************************/

    const vec3& Surface::vertex( index_t f, index_t v ) const
    {
        grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
        return vertex( surf_vertex_id( f, v) ) ;
    }


    const vec3& Surface::vertex( index_t surf_vertex_id ) const {
        return model_->vertex( vertices_.at(surf_vertex_id) ) ;
    }

    index_t Surface::model_facet_id( index_t f ) const {
        return model_->model_facet( id_, f ) ;
    }

    void Surface::set_first_triangle_as_key()
    {
        // I guess it should'nt be a problem if the first facet is not a triangle
        // that's only a guess (Jeanne)
        key_facet_ = KeyFacet( model_->vertex( vertices_[facets_[0]] ),
            model_->vertex( vertices_[facets_[1]] ),
            model_->vertex( vertices_[facets_[2]] ) ) ;
    }

    /*int Surface::adjcent_in_neighbor( index_t f, index_t e ) const
    {
        signed_index_t adj = adjacent( f, e ) ;
        if( adj == -1 ) return -1 ;

        for( index_t i = 0; i < nb_vertices_in_facet( adj ); i++ ) {
            if( adjacent( adj, i ) == f ) return i ;
        }
        return -1 ;
    }*/


    /**
      * Parcours d'un bord dans une surface - dans un sens ou dans l'autre
      *
      * Arnaud- je dirai que la fonction d'avant �tait bugg�e
      * 
      * Need of two indices in input - so that we are able to go to the next 
      * edge on border in any direction and avoid going back when the next edge on boundary 
      * is in the same facet
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
        //int next_id1 = -1 ;

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

    void Surface::next_on_border( index_t f, index_t e, index_t& next_f, index_t& next_e ) const {
        index_t v = next_in_facet( f, e ) ;
        index_t next_in_next = NO_ID ;
        return next_on_border( f, e, v, next_f, next_e, next_in_next ) ;
    }


    /**
     * Check if the facet has an edge with the given ids in the Surface
     * Returns the id of the vertex in the facet at which start the edge
     * Returns -1 if no edge is found
     *
    signed_index_t Surface::has_edge( signed_index_t f, signed_index_t v0, signed_index_t v1 ) const {

        index_t prev = surf_vertex_id( f, nb_vertices_in_facet(f)-1 ) ; 
        for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
            index_t p = vertex_index( f, v ) ;
            if( (prev == v0 && p == v1) ||
                (prev == v1 && p == v0) ) return prev ;
            prev = p ;
        }
        return -1 ;                               
    }*/

    /** 
     * Check if the surface has an edge starting at v0 and ending at v1
     * Returns the id in the facet at which starts the edge
     * Returns -1 if this edge is not found
     *
    signed_index_t Surface::has_oriented_edge( signed_index_t f, signed_index_t v0, signed_index_t v1 ) const {
        index_t prev = surf_vertex_id( f, nb_vertices_in_facet(f)-1 ) ; 
        for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
            index_t p = surf_vertex_id( f, v ) ;
            if( prev == v0 && p == v1 ) return v ;
            prev = p ;
        }
        return -1 ;
    }*/

    /** 
     * Find the first facet of the surface that has an edge 
     * linking the two vertices (ids in the surface)
     */ 
    index_t Surface::facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const {
        grgmesh_debug_assert( in0 < vertices_.size() && in1 < vertices_.size() ) ;
        
        // Check for all the facets 

        // Sans doute un truc plus rapide - regarder si les deux indices se suivent
        // dans facets_ et ensuite v�rifier si c'est bien la m�me facette

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

    index_t Surface::facet_from_model_vertex_ids( index_t i0, index_t i1 ) const {
        index_t facet = NO_ID ;
        index_t edge = NO_ID ;
        edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
        return facet ;
    }

    /**
     * Get the id of one facet and the corresponding edge 
     * There might be two !! Get only the first
     */
    void Surface::edge_from_model_vertex_ids(
        index_t i0,
        index_t i1,
        index_t& facet,
        index_t& edge ) const
    {
         // Copy from above .. tant pis
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
        // Si on arrive l� on a rien trouv� put facet to -1
        facet = NO_ID ;
            
    }

     /**
     * Get the id of one facet and the corresponding edge 
     * There might be two !! Get only the first
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
        // Si on arrive l� on a rien trouv� put facet to -1
        facet = NO_ID ;
            
    }



   /* bool Surface::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }

    signed_index_t Surface::find( const vec3& p ) const
    {
        for( index_t i = 0; i < vertices_.size(); ++i ) {
            if( model_->vertex( vertices_[i] ) == p ) return i ;
        }
        return -1 ;
    }*/

    /*! Returns the id of a facet that has these two vertices 
     *  if any else returns -1 
     *  
     *  WARNING There might TWO such facets in the surface
     */
    /*int Surface::find_facet( const vec3& p0, const vec3& p1 ) const {
        // There might be several vertices with the same coordinates
        // Test all possible pairs
        
        std::vector< signed_index_t > i0 ;
        std::vector< signed_index_t > i1 ;

        for( index_t v = 0; v < vertices_.size(); ++v ) {
            if( model_->vertex( vertices_[v] ) == p0 ) i0.push_back( v ) ;
            if( model_->vertex( vertices_[v] ) == p1 ) i1.push_back( v ) ;
        }
        signed_index_t t = -1 ;
        for( signed_index_t i = 0; i < i0.size(); ++i ){
            for( signed_index_t j = 0; j < i1.size(); ++j ){
                t = facet_from_surface_vertex_ids( i0[i], i1[j] ) ; 
                if( t != -1 ) return t ;
            }
        }
        return -1 ;
    }*/

    /** 
     * Find the first edge that contains the 2 given vertices (ids in the surface)
     * Return the id of the vertex at which start the edge or -1 if no edge is found
     *
    signed_index_t Surface::find_edge( signed_index_t id0, signed_index_t id1 ) const {
        for( index_t f = 0; f < nb_cells(); ++f ) {
             signed_index_t p = has_edge( f, id0, id1 ) ;
             if( p!= -1 ) return p ;
        }
        return -1 ; 
    }*/
    

    /**
     * Returns the vertex of the facet which id in the surface is the given one
     */
    index_t Surface::facet_vertex_id( index_t t, index_t surf_vertex_id_in ) const {
        for( index_t v = 0; v < nb_vertices_in_facet(t); v++ ) {
            if( surf_vertex_id( t, v ) == surf_vertex_id_in ) return v ;
        }
        return NO_ID ;
    }

/*    signed_index_t Surface::facet_vertex_id( signed_index_t t, const vec3& p ) const {
        for( index_t v = 0; v < nb_vertices_in_facet(t); v++ ) {
            if( vertex( t, v ) == p ) return v ;
        }
        return -1 ;
    }*/

    // this is a copy
    struct comp_vec3bis {
        bool operator()( const vec3& l, const vec3& r ) const {
            if( l.x != r.x ) return l.x < r.x ;
            if( l.y != r.y ) return l.y < r.y ;
            return l.z < r.z ;
        }
    } ;

   

    /*int Surface::edge_id( signed_index_t t, signed_index_t p0, signed_index_t p1 ) const {        
        signed_index_t t_0 = facet_vertex_id( t, p0 ) ;
        signed_index_t t_1 = facet_vertex_id( t, p1 ) ;
       
        if( t_0 > t_1 ) { signed_index_t tmp = t_0 ; t_0 = t_1 ; t_1 = tmp ; }
        
        if     ( t_0 == 0 && t_1 == 1 ) return 2 ;
        else if( t_0 == 0 && t_1 == 2 ) return 1 ;
        else if( t_0 == 1 && t_1 == 2 ) return 0 ;
        else return -1 ;
    }*/


    /**
     * \todo Find a way to make this faster !! It is not that bad actually
     * How ? 
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

    /** Determine the facets sharing the given vertex (id in the surface)
     * 
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

    vec3 Surface::facet_barycenter( index_t f ) const {
        vec3 barycenter( 0., 0., 0. ) ;
        for( index_t i = 0; i < nb_vertices_in_facet( f ); i++ ) {
            barycenter += vertex( f, i ) ;
        }
        return barycenter / nb_vertices_in_facet( f ) ;
    }
     
    double Surface::facet_area( index_t f ) const {
        double result = 0 ;
        for( index_t i = 1; i+1 < nb_vertices_in_facet( f ); i++ ) 
        { 
            result += Utils::triangle_area( 
                vertex( f, 0 ), vertex( f, i ), vertex( f, i+1 ) ) ;
        }
        return result ;
    }

    /**
     * \brief Returns the normal to the triangle made by the first 3 vertices
     * of the facet
     * WARNING : if the facet is not planar calling this has no meaning
     */
    vec3 Surface::facet_normal( index_t f ) const {
        const vec3& p0 = vertex( f, 0 )  ;
        const vec3& p1 = vertex( f, 1 )  ;
        const vec3& p2 = vertex( f, 2 )  ;
        vec3 c0 = cross(p0-p2, p1-p2) ;       
        return normalize( c0 ) ;
    }

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


   
  /*  Box3d Surface::bbox() const
    {
        Box3d result ;
        for( index_t p = 0; p < nb_vertices(); p++ ) {
            result.add_vertex( vertex( p ) ) ;
        }
        return result ;
    }*/

    void SurfaceMutator::set_vertex( index_t id, const vec3& p ) {
        M_.model_->vertices_[M_.vertices_[id]] = p ;
    }

    vec3& SurfaceMutator::vertex( index_t p ) const
    {
        return M_.model_->vertices_[ M_.vertices_[p] ] ;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    double BoundaryModelElementMeasure::size( const BoundaryModelElement* E ) {        
        double result = 0. ;

        // If this element has children sum up their sizes
        for( index_t i = 0; i < E->nb_children(); ++i ){
            result += BoundaryModelElementMeasure::size( &E->child(i) )  ;
        }
        if( result != 0 ) return result ;

        // Else it is a base element
        
        // If this is a region 
        if( E->dim() == 3 ) {
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
        else if( E->dim() == 0 ) {
            return 0 ; 
        }
        else if( E->dim() == 1 ) {
            const Line* L = dynamic_cast< const Line* >( E ) ;
            grgmesh_assert( L != nil ) ;           
            for( index_t i = 1; i < E->nb_vertices(); ++i ){
                result += GEO::Geom::distance( E->vertex( i ), E->vertex( i-1 ) ) ;
            }
            return result ;
        }
        else if( E->dim() == 2 ) {
            const Surface* S = dynamic_cast< const Surface* >( E ) ;
            grgmesh_assert( S != nil ) ;

            for( index_t i = 0; i < S->nb_cells(); i++ ) {
                result += S->facet_area( i ) ;
            }
            return result ;
        }
        grgmesh_assert_not_reached ;
    }        

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

