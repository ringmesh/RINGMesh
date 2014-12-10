/*! This file is originally part of the Geomodeling plugin of Graphite 
 *  The initial version was developed by Jeanne Pellerin and modified by Arnaud Botella
 *  
 *  This version is modified by Jeanne Pellerin - WIAS
 */

/*! \author Jeanne Pellerin */

#include <grgmesh/boundary_model.h>
#include <grgmesh/boundary_model_element.h>
#include <grgmesh/utils.h>
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
        if( !has_parent() ) return nil ;
        switch( dim() ) {
            case 1:
                return model_->contact( parent_ ) ;
            case 2:
                return model_->one_interface( parent_ ) ;
            default:
                assert( false ) ;
                return nil ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::boundary( unsigned int x ) const
    {
        if( x >= nb_boundaries() ) return nil ;
        unsigned int id = boundaries_.at(x) ;
        switch( dim() ) {
            case 1:
                return model_->corner( id ) ;
            case 2:
                if( has_parent() ) return model_->line( id ) ;
                return model_->contact( id ) ;
            case 3:
                return model_->surface( id ) ;
            default:
                assert( false ) ; ;
                return nil ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::in_boundary( unsigned int x ) const
    {
        if( x >= nb_in_boundary() ) return nil ;
        unsigned int id = in_boundary_.at(x) ;
        switch( dim() ) {
            case 0:
                return model_->line( id ) ;
            case 1:
                if( has_parent() ) return model_->surface( id ) ;
                return model_->one_interface( id ) ;
            case 2:
                return model_->region( id ) ;
            default:
                assert( false ) ; ;
                return nil ;
        }
    }

    const BoundaryModelElement& BoundaryModelElement::child( unsigned int x ) const
    {
        if( has_parent() || x >= nb_children() ) return nil ;
        unsigned int id = children_.at(x) ;
        switch( dim() ) {
            case 1:
                return model_->line( id ) ;
            case 2:
                return model_->surface( id ) ;
            case 3:
                return model_->region( id ) ;
            default:
                assert( false ) ; ;
                return nil ;
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


    unsigned int BoundaryModelElement::nb_cells() const {
        unsigned int result = 0 ;
         for( unsigned int i=0; i < nb_children(); ++i ){
            result += child( i ).nb_cells() ;
        }
        return result ;
    }
  
    

/***********************************************************************************************/

    const vec3& Corner::point( unsigned int p ) const
    {
        return model_->point( p_ ) ;
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

    Line::Line( BoundaryModel* model, int id ):
        BoundaryModelElement( model, 1, id )
    { 
        boundaries_.resize( 2, nil) ; 
    }

    Line::Line(
        BoundaryModel* model,
        int id,
        const std::vector< unsigned int >& points )
        : BoundaryModelElement( model, 1, id ), points_( points )
    {
    }

    Line::Line(
        BoundaryModel* model,
        int id,
        unsigned int corner0,
        unsigned int corner1,
        const std::vector< unsigned int >& points
    ):  BoundaryModelElement( model, 1, id ),
        points_( points )
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
 
    /*void Line::add_in_boundary( unsigned int e )
    {
        for( unsigned int i = 0; i < nb_in_boundary(); i++ ) {
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

    const vec3& Line::point( unsigned int line_point_id ) const {
        return model_->point( points_.at( line_point_id ) ) ;
    }

    vec3 Line::segment_barycenter( uint32 s ) const {
        return 0.5*( point(s) + point(s+1) ) ;
    }
    
    double Line::segment_length( uint32 s ) const {
        return length( point(s+1)-point(s) ) ;
    }

   /* bool Line::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }
    int Line::find( const vec3& p ) const
    {
        for( unsigned int i = 0; i < vertices_.size(); ++i ) {
            if( point( i ) == p ) return i ;
        }
        return -1 ;
    } */

    
    void LineMutator::set_point( unsigned int id, const vec3& p ) {
        M_.model_->points_[M_.points_[id]] = p ;
    }

    vec3& LineMutator::point( unsigned int p ) const
    {
        return M_.model_->points_[ M_.points_[p] ] ;
    }

/********************************************************************************************/
/******             Surface implementation            ***********************************/
/********************************************************************************************/

    const vec3& Surface::point( uint32 f, uint32 v ) const
    {
        grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
        return model_->point( points_.at(facets_[facet_begin( f ) + v]) ) ;
    }


    const vec3& Surface::point( uint32 surf_point_id ) const {
        return model_->point( points_.at(surf_point_id) ) ;
    }

    void Surface::set_first_triangle_as_key()
    {
        // I guess it should'nt be a problem if the first facet is not a triangle
        // that's only a guess (Jeanne)
        key_facet_ = KeyFacet( model_->point( points_[facets_[0]] ),
            model_->point( points_[facets_[1]] ),
            model_->point( points_[facets_[2]] ) ) ;
    }

    /*int Surface::adjcent_in_neighbor( unsigned int f, unsigned int e ) const
    {
        int adj = adjacent( f, e ) ;
        if( adj == -1 ) return -1 ;

        for( unsigned int i = 0; i < nb_points_in_facet( adj ); i++ ) {
            if( adjacent( adj, i ) == f ) return i ;
        }
        return -1 ;
    }*/


    /**
      * Parcours d'un bord dans une surface - dans un sens ou dans l'autre
      *
      * Arnaud- je dirai que la fonction d'avant était buggée
      * 
      * Need of two indices in input - so that we are able to go to the next 
      * edge on border in any direction and avoid going back when the next edge on boundary 
      * is in the same facet
      */
    void Surface::next_on_border( 
        uint32 f, uint32 from, uint32 v, 
        int32& next_f, int32& v_in_next, int32& next_in_next ) const 
    {
        grgmesh_assert( v < nb_points_in_facet( f ) ) ;
        grgmesh_assert( is_on_border( f, v ) || is_on_border( f, from ) ) ;

        uint32 V = surf_point_id( f, v ) ;

        // We want the next triangle that is on the boundary and share V
        // If there is no such triangle, the next point on the boundary 
        // is the point of F neighbor of V that is not from 
        

        // Get the facets around the shared point that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)        
        std::vector< uint32 > facets ;  
        int nb_around = facets_around_point( V, facets, true, f ) ;
        grgmesh_assert( nb_around < 3 && nb_around > 0 ) ;

        next_f = facets[0] ;
        //int next_id1 = -1 ;

        if( nb_around == 2 ) {
            if( next_f == f ) next_f = facets[1] ;
            grgmesh_debug_assert( next_f != f ) ;

            // Now get the other point that is on the boundary opposite to p1
            v_in_next = facet_point_id( next_f, V ) ;
            grgmesh_assert( v_in_next != -1 ) ;

            // The edges containing V in next_f are
            // the edge starting at v_in_next and the one ending there
            int prev_v_in_next = prev_in_facet( next_f, v_in_next )  ;

            bool e0_on_boundary = is_on_border( next_f, v_in_next ) ;
            bool e1_on_boundary = is_on_border( next_f, prev_v_in_next ) ;

            // Only one must be on the boundary otherwise there is a corner missing
            grgmesh_assert( e0_on_boundary != e1_on_boundary ) ;

            // From the edge that is on boundary get the next vertex on this boundary
            // If the edge starting at p_in_next is on boundary, new_point is its next
            // If the edge ending at p_in_next is on boundary, new point is its prev
            next_in_next = e0_on_boundary ? next_in_facet( next_f, v_in_next ) : prev_v_in_next ;

        } else if( nb_around == 1 ) {
            // V must be in two border edges of facet f 
            // Get the id in the facet of the point neighbor of v1 that is not v0           
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

    void Surface::next_on_border( uint32 f, uint32 e, int32& next_f, int32& next_e ) const {
        uint32 v = next_in_facet( f, e ) ;
        int32 next_in_next = -1 ;
        return next_on_border( f, e, v, next_f, next_e, next_in_next ) ;
    }


    /**
     * Check if the facet has an edge with the given ids in the Surface
     * Returns the id of the point in the facet at which start the edge
     * Returns -1 if no edge is found
     *
    int Surface::has_edge( int f, int v0, int v1 ) const {

        unsigned int prev = surf_point_id( f, nb_points_in_facet(f)-1 ) ; 
        for( unsigned int v = 0; v < nb_points_in_facet( f ); ++v ) {
            unsigned int p = point_index( f, v ) ;
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
    int Surface::has_oriented_edge( int f, int v0, int v1 ) const {
        unsigned int prev = surf_point_id( f, nb_points_in_facet(f)-1 ) ; 
        for( unsigned int v = 0; v < nb_points_in_facet( f ); ++v ) {
            unsigned int p = surf_point_id( f, v ) ;
            if( prev == v0 && p == v1 ) return v ;
            prev = p ;
        }
        return -1 ;
    }*/

    /** 
     * Find the first facet of the surface that has an edge 
     * linking the two points (ids in the surface)
     */ 
    int Surface::facet_from_surface_point_ids( uint32 in0, uint32 in1 ) const {
        grgmesh_debug_assert( in0 < points_.size() && in1 < points_.size() ) ;
        
        // Check for all the facets 

        // Sans doute un truc plus rapide - regarder si les deux indices se suivent
        // dans facets_ et ensuite vérifier si c'est bien la même facette

        for( unsigned int f = 0; f < nb_cells(); ++f ) {
            int found = -1 ;
            unsigned int prev = surf_point_id( f, nb_points_in_facet(f)-1 ) ; 
            for( unsigned int v = 0; v < nb_points_in_facet( f ); ++v ) {
                unsigned int p = surf_point_id( f, v ) ;
                if( (prev == in0 && p == in1) ||
                    (prev == in1 && p == in0) ) 
                {
                        found = prev ;
                        break ;
                }
                prev = p ;
            }            
            if( found != -1 ) {
                return f ;
            }
        }
        return -1 ; 
    }

    int Surface::facet_from_model_point_ids( uint32 i0, uint32 i1 ) const {
        int facet = -1 ;
        int edge = -1 ;
        edge_from_model_point_ids( i0, i1, facet, edge ) ;
        return facet ;
    }

    /**
     * Get the id of one facet and the corresponding edge 
     * There might be two !! Get only the first
     */
    void Surface::edge_from_model_point_ids( uint32 i0, uint32 i1, int32& facet, int32& edge ) const {
         // Copy from above .. tant pis
        edge = -1 ;
        
        // If a facet is given, look for the edge in this facet only
        if( facet != -1 ) {
            for( unsigned int v = 0; v < nb_points_in_facet( facet ); ++v ) {
                uint32 prev = model_point_id( facet, prev_in_facet( facet, v )  ) ; 
                uint32 p = model_point_id( facet, v ) ;
                if( ( prev == i0 && p == i1) ||
                    ( prev == i1 && p == i0 ) ) 
                {
                    edge = prev_in_facet( facet, v ) ;
                    return ;
                }
            }
        }
        else {
            for( unsigned int f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                edge_from_model_point_ids( i0, i1, facet, edge ) ;
                if( edge != -1 ) return ;
            }
        }
        // Si on arrive là on a rien trouvé put facet to -1
        facet = -1 ; 
            
    }

     /**
     * Get the id of one facet and the corresponding edge 
     * There might be two !! Get only the first
     */
    void Surface::oriented_edge_from_model_point_ids( uint32 i0, uint32 i1, int32& facet, int32& edge ) const {
         // Copy from above .. tant pis
        edge = -1 ;
        
        // If a facet is given, look for the oriented edge in this facet only
        if( facet != -1 ) {
            for( unsigned int v = 0; v < nb_points_in_facet( facet ); ++v ) {
                uint32 p = model_point_id( facet, v ) ;
                uint32 next = model_point_id( facet, next_in_facet( facet, v )  ) ; 

                if( p == i0 && next == i1 )
                {
                    edge = v ;
                    return ;
                }
            }
        }
        else {
            for( unsigned int f = 0; f < nb_cells(); ++f ) {
                facet = f ;
                oriented_edge_from_model_point_ids( i0, i1, facet, edge ) ;
                if( edge != -1 ) return ;
            }
        }
        // Si on arrive là on a rien trouvé put facet to -1
        facet = -1 ; 
            
    }



   /* bool Surface::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }

    int Surface::find( const vec3& p ) const
    {
        for( unsigned int i = 0; i < points_.size(); ++i ) {
            if( model_->point( points_[i] ) == p ) return i ;
        }
        return -1 ;
    }*/

    /*! Returns the id of a facet that has these two points 
     *  if any else returns -1 
     *  
     *  WARNING There might TWO such facets in the surface
     */
    /*int Surface::find_facet( const vec3& p0, const vec3& p1 ) const {
        // There might be several points with the same coordinates
        // Test all possible pairs
        
        std::vector< int > i0 ;
        std::vector< int > i1 ;

        for( unsigned int v = 0; v < points_.size(); ++v ) {
            if( model_->point( points_[v] ) == p0 ) i0.push_back( v ) ;
            if( model_->point( points_[v] ) == p1 ) i1.push_back( v ) ;
        }
        int t = -1 ;
        for( int i = 0; i < i0.size(); ++i ){
            for( int j = 0; j < i1.size(); ++j ){
                t = facet_from_surface_point_ids( i0[i], i1[j] ) ; 
                if( t != -1 ) return t ;
            }
        }
        return -1 ;
    }*/

    /** 
     * Find the first edge that contains the 2 given points (ids in the surface)
     * Return the id of the point at which start the edge or -1 if no edge is found
     *
    int Surface::find_edge( int id0, int id1 ) const {
        for( unsigned int f = 0; f < nb_cells(); ++f ) {
             int p = has_edge( f, id0, id1 ) ;
             if( p!= -1 ) return p ;
        }
        return -1 ; 
    }*/
    

    /**
     * Returns the point of the facet which id in the surface is the given one
     */
    int Surface::facet_point_id( uint32 t, uint32 surf_point_id_in ) const {
        for( unsigned int v = 0; v < nb_points_in_facet(t); v++ ) {
            if( surf_point_id( t, v ) == surf_point_id_in ) return v ;
        }
        return -1 ;
    }

/*    int Surface::facet_point_id( int t, const vec3& p ) const {
        for( unsigned int v = 0; v < nb_points_in_facet(t); v++ ) {
            if( point( t, v ) == p ) return v ;
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

   

    /*int Surface::edge_id( int t, int p0, int p1 ) const {        
        int t_0 = facet_point_id( t, p0 ) ;
        int t_1 = facet_point_id( t, p1 ) ;
       
        if( t_0 > t_1 ) { int tmp = t_0 ; t_0 = t_1 ; t_1 = tmp ; }
        
        if     ( t_0 == 0 && t_1 == 1 ) return 2 ;
        else if( t_0 == 0 && t_1 == 2 ) return 1 ;
        else if( t_0 == 1 && t_1 == 2 ) return 0 ;
        else return -1 ;
    }*/


    /**
     * \todo Find a way to make this faster !! It is not that bad actually
     * How ? 
     */
    int Surface::facets_around_point(
        uint32 shared_point,
        std::vector< uint32 >& result,
        bool border_only ) const
    {
        result.resize(0) ;
        for( unsigned int t = 0; t < nb_cells(); ++t ) {            
            for( unsigned int v = 0; v < nb_points_in_facet(t); v++ ) {
                if( surf_point_id( t, v ) == shared_point ) {
                    return facets_around_point( shared_point, result,
                        border_only, t ) ;
                }
            }
        }
        grgmesh_assert_not_reached ;
        return -1 ;
    }


    /** Determine the facets sharing the given point (id in the surface)
     * 
     */
    int Surface::facets_around_point(
        uint32 P,
        std::vector< uint32 >& result,
        bool border_only,
        uint32 f0 ) const
    {
        result.resize( 0 ) ;

        // Flag the visited facets
        std::vector< int > visited ;
        visited.reserve( 10 ) ;
        
        // Stack of the adjacent facets
        std::stack< int > S ;
        S.push( f0 ) ;
        visited.push_back( f0 ) ;
        
        do {
            int t = S.top() ;
            S.pop() ;         

            for( unsigned int v = 0; v < nb_points_in_facet(t); ++v ) {               
                if( surf_point_id( t, v ) == P ) 
                {
                    int adj_P = adjacent( t, v ) ;                    
                    unsigned int prev = prev_in_facet( t, v ) ;
                    int adj_prev = adjacent( t, prev ) ;

                    if( adj_P != -1 ){
                        // The edge starting at P is not on the boundary
                        if( !Utils::contains( visited, adj_P ) ) {
                            S.push( adj_P ) ;
                            visited.push_back( adj_P ) ;
                        }
                    }                  
                    if( adj_prev != -1 ) {
                        // The edge ending at P is not on the boundary
                        if( !Utils::contains( visited, adj_prev ) ) {
                            S.push( adj_prev ) ;
                            visited.push_back( adj_prev  ) ;
                        }
                    }
                    
                    if( border_only ) {
                        if( adj_P == -1 || adj_prev == -1 ) { 
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

    vec3 Surface::facet_barycenter( int f ) const {
        vec3 barycenter( 0., 0., 0. ) ;
        for( unsigned int i = 0; i < nb_points_in_facet( f ); i++ ) {
            barycenter += point( f, i ) ;
        }
        return barycenter / nb_points_in_facet( f ) ;
    }
     
    double Surface::facet_area( int f ) const {
        double result = 0 ;
        for( unsigned int i = 1; i+1 < nb_points_in_facet( f ); i++ ) 
        { 
            result += Utils::triangle_area( 
                point( f, 0 ), point( f, i ), point( f, i+1 ) ) ;
        }
        return result ;
    }

    /**
     * \brief Returns the normal to the triangle made by the first 3 points
     * of the facet
     * WARNING : if the facet is not planar calling this has no meaning
     */
    vec3 Surface::facet_normal( int f ) const {
        const vec3& p0 = point( f, 0 )  ;
        const vec3& p1 = point( f, 1 )  ;
        const vec3& p2 = point( f, 2 )  ;
        vec3 c0 = cross(p0-p2, p1-p2) ;       
        return normalize( c0 ) ;
    }


    unsigned int Surface::closest_point_in_facet(
        unsigned int f,
        const vec3& v ) const
    {
       unsigned int result = 0 ;
       double dist = DBL_MAX ;
       for( unsigned int p = 0; p < nb_points_in_facet( f ); p++ ) {
           double distance = length2( v - point( f, p ) ) ;
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
        for( unsigned int p = 0; p < nb_points(); p++ ) {
            result.add_point( point( p ) ) ;
        }
        return result ;
    }*/

    void SurfaceMutator::set_point( unsigned int id, const vec3& p ) {
        M_.model_->points_[M_.points_[id]] = p ;
    }

    vec3& SurfaceMutator::point( unsigned int p ) const
    {
        return M_.model_->points_[ M_.points_[p] ] ;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


    vec3 BoundaryModelElementMeasure::barycenter ( 
        const BoundaryModelElement* E, const std::vector< uint32 >& cells ) 
    {
        vec3 result(0., 0., 0. ) ;
        double size = 0 ;

        const Line* L = dynamic_cast< const Line* >( E ) ;
        if( L != nil ) {
            for( unsigned int i = 0; i < cells.size(); ++ i ) {
                result += L->segment_length( cells[i] ) * L->segment_barycenter( cells[i] ) ;
                size   += L->segment_length( cells[i] ) ;
            }
            return size > epsilon ? result/size : result ;
        }
        const Surface* S = dynamic_cast< const Surface* >( E ) ;
        if( S != nil ) {
            for( unsigned int i = 0; i < cells.size(); ++ i ) {
                result += S->facet_area( cells[i] ) * S->facet_barycenter( cells[i] ) ;
                size   += S->facet_area( cells[i] ) ;
            }
            return size > epsilon ? result/size : result ;
        }
        
        grgmesh_assert_not_reached ;
        return result ;
    }



#ifdef TOTOTO


     /*! Return the neighbors of the element = the ones that share one of
     *  its boundaries. Only the boundaries of the given type are considered.
     *  If exclude_voi is true, the negighbors belonging to the VOI are ignored
     */
    std::vector< const BoundaryModelElement* > 
        BoundaryModelElementMeasure::compute_neighbors(
        GEOL_FEATURE through,
        bool exclude_voi
    ) const 
    {
        std::vector< const BoundaryModelElement* > result ;
     
        for( int i = 0; i < boundaries_.size(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            if( through != ALL && b->type() != through ) continue ;

            for( int j = 0; j < b->nb_in_boundary(); ++j ){
                const BoundaryModelElement* bb = b->in_boundary( j ) ;
                if( bb == this ) continue ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                else result.push_back( bb ) ;
            }            
        }
        std::sort( result.begin(), result.end() ) ;
        int size = std::unique(result.begin(), result.end()) - result.begin() ;
        result.resize( size ) ;

        return result ;
    }

    double BoundaryModelElementMeasure::size() const {
        double result = 0. ;
        // If this element has children sum up their sizes
        for( unsigned int i = 0; i < nb_children(); ++i ){
            result += child(i)->size() ;
        }

        if( result == 0 ){
            if( dim_ == 3 ) {
                // Compute the volume if this is a region
                for( unsigned int i = 0; i < nb_boundaries(); i++ ) {
                    const Surface* surface = dynamic_cast< const Surface* >( boundary( i ) ) ;
                   
                    for( unsigned int t = 0; t < surface->nb_cells(); t++ ) {
                        const vec3& p0 = surface->point(t, 0 ) ;
                        for( unsigned int v = 1; v+1 < surface->nb_points_in_facet(t); ++v ){
                            double cur_volume = ( dot( p0,
                                cross( surface->point( t, v ), surface->point( t, v+1 ) ) ) )
                                / static_cast< double >( 6 ) ;
                            sides_[i] ? result -= cur_volume : result += cur_volume ;
                        }
                    }
                }
            }
        }
        return fabs( result ) ;
    }

    

    double BoundaryModelElementMeasure::distance( const vec3& p ) const {
        if( nb_children() == 0 ) return FLT_MAX ;
        else {
            double result = FLT_MAX ;
            for( unsigned int i = 0; i < nb_children(); ++i ){
                double dist = child(i)->distance(p) ;
                result = ( dist < result ) ? dist : result ;
            }
            return result ;
        }
    }
    double BoundaryModelElementMeasure::distance( BoundaryModelElement* e ) const {
        if( nb_children() == 0 ) return FLT_MAX ;
        else {
            double result = FLT_MAX ;
            for( unsigned int i = 0; i < nb_children(); ++i ){
                double dist = child(i)->distance(e) ;
                result = ( dist < result ) ? dist : result ;                
            }
            return result ;
        }
    }

    double BoundaryModelElementMeasure::min_angle( BoundaryModelElement* e ) const {
        if( nb_children() == 0 ) return 999 ;
        if( dim_ == 1 || dim_ == 2 ) return 999 ;

        double result = 180. ;
        for( unsigned int i=0; i < nb_children(); ++i ){
            double cur =  child( i )->min_angle( e )  ;
            result = result < cur ? result : cur ;
        }
        return result ;        
    }

    /*! Mean orientation
     *  For surfaces it is the average Normal
     */
    vec3 BoundaryModelElement::average_orientation() const {
        if( nb_children() == 0 ) return vec3(-99999, -99999, -99999) ;
        if( dim_ == 0 || dim_ == 3 ) return vec3(-99999, -99999, -99999) ;
    
        vec3 result (0,0,0) ;
        double total_size = 0 ;
        for( unsigned int i=0; i < nb_children(); ++i ){
            double s = child( i )->size() ;
            result += s * child( i )->average_orientation() ;
            total_size += s ;
        }
        if(total_size > 10e-30 ) return result/total_size ;
        else return vec3(0.,0.,0.) ;
    }

    /*! Probably not correct 
     */
    double BoundaryModelElement::average_angle_to_y() const {
        vec3 n = average_orientation() ;
        if( n.x >= 0 ) return std::acos( n.y / sqrt( n.x*n.x + n.y*n.y ) )*180./M_PI ;
        else return M_PI - std::acos( n.y / sqrt( n.x*n.x + n.y*n.y ) )*180./M_PI ;
    }

    /*! Probably not correct 
     */
    double BoundaryModelElement::average_dip() const {
        vec3 n = average_orientation() ;
        double phi = std::acos( n.z )*180. / M_PI ;
        if( phi <= 90. ) return 90.-phi ;
        else return phi-90. ;
    }


    /*! What information goes out in the CSV file in which
     *  some complexity will be computed
     */
    void BoundaryModelElement::print_categories( std::ostream& out ) {
        out << "Dimension" << SEP 
            << "Id" << SEP
            << "Name" << SEP
            << "Parent" << SEP
            << "Type" << SEP
            << "Nb boundaries"<< SEP
            << "Nb in boundary of" << SEP 
            << "Number of neighbors" << SEP
            << "Nb children" << SEP
            << "Size" << SEP

            << "Boundary ids" << SEP
            << "In boundary of ids" << SEP
            << "Neighbor ids" << SEP
            << "Children ids " << SEP
            
            << std::endl ;
    }
    
    /* To synchronize with the categories = columns */
    void BoundaryModelElement::print( std::ostream& out ) const {    
    
        out << dim_ << SEP
            << id_ << SEP
            << name_ << SEP
            << parent_ << SEP ;
        
        BoundaryModelMeasure::print_type( out, type_, dim_ ) ;

        out << boundaries_.size() << SEP
            << in_boundary_.size() << SEP ;

        std::vector< const BoundaryModelElement* > n = compute_neighbors( ALL, false ) ;
        out << n.size() << SEP ;

        if( children_.size() == 0 ) out << "-1" << SEP ;
        else out << children_.size() << SEP ;
    
        out << size() << SEP ;
        
        print_vector_contents( boundaries_, out ) ; 
        print_vector_contents( in_boundary_, out ) ; 
        print_vector_contents( n, out ) ; 
        print_vector_contents( children_, out ) ; 

        out << std::endl ;
    }


    struct CompPair {
         bool operator()( const std::pair< double, double >& r, const std::pair< double, double >& l ) const {
            return r.first < l.first ;
        }
    } ;

    /* To get some statistics on a vector of pairs
     * Not very efficient.
     */
    void print_stats( std::ostream& out, std::vector< std::pair< double, double> >& values, 
        double min1 = -DBL_MAX, double min2 = -DBL_MAX, double min3 = -DBL_MAX,
        double max1 =  DBL_MAX, double max2 =  DBL_MAX, double max3 =  DBL_MAX
    ) { 
        if( values.size() == 0 ){
            out << "No values " << SEP ;
            out << "" << SEP
                << "" << SEP
                << "" << SEP  ;            
            if( min1 != -DBL_MAX ) out << "" << SEP ;
            if( min2 != -DBL_MAX ) out << "" << SEP ;
            if( min3 != -DBL_MAX ) out << "" << SEP ;
            if( max1 !=  DBL_MAX ) out << "" << SEP ;
            if( max2 !=  DBL_MAX ) out << "" << SEP ;
            if( max3 !=  DBL_MAX ) out << "" << SEP ;
            return ;
        }
        std::sort( values.begin(), values.end() ) ;

        double W = 0 ;
        double V = 0 ;

        double q1 = 0. ;
        double q2 = 0. ;
        double q3 = 0. ;

        double r1 = 0. ;
        double r2 = 0. ;
        double r3 = 0. ;
        
        for( unsigned int i = 0 ; i < values.size();  ++i ){
            V += values[i].second * values[i].first ;
            W += values[i].second ;

            // Some of the tests are useless if the default values are used
            if( values[i].first < min1 ) q1 += values[i].second ;
            if( values[i].first < min2 ) q2 += values[i].second ;
            if( values[i].first < min3 ) q3 += values[i].second ;

            if( values[i].first > max1 ) r1 += values[i].second ;
            if( values[i].first > max2 ) r2 += values[i].second ;
            if( values[i].first > max3 ) r3 += values[i].second ;
        }
        if( W < 10e-30 ) {
            out << "Nil weights" << SEP << std::endl ;
            out << "" << SEP
                << "" << SEP
                << "" << SEP ;             
            if( min1 != -DBL_MAX ) out << "" << SEP ;
            if( min2 != -DBL_MAX ) out << "" << SEP ;
            if( min3 != -DBL_MAX ) out << "" << SEP ;
            if( max1 !=  DBL_MAX ) out << "" << SEP ;
            if( max2 !=  DBL_MAX ) out << "" << SEP ;
            if( max3 !=  DBL_MAX ) out << "" << SEP ;
            return ;
        }

        double mean = V/W ;

        double SD = 0 ;
        for( unsigned int i = 0 ; i < values.size();  ++i ){
            SD += values[i].second * (values[i].first-mean)*(values[i].first-mean) ;
        }
        SD = sqrt( SD / W ) ;           

        out << values[0].first     << SEP
            << values.back().first << SEP 
            << mean                << SEP
            << SD                  << SEP ;
        if( min1 != -DBL_MAX ) 
            out << q1              << SEP ;
        if( min2 != -DBL_MAX ) 
            out << q2            << SEP ;
        if( min3 != -DBL_MAX ) 
            out << q3            << SEP ;       

        if( max1 !=  DBL_MAX ) 
            out << r1            << SEP ;
        if( max2 !=  DBL_MAX ) 
            out << r2            << SEP ;
        if( max3 !=  DBL_MAX ) 
            out << r3            << SEP ;

        out << W                 << SEP ;

    }
   
    void BoundaryModelElement::print_complexity_categories( std::ostream& out ) {
        out << "ELEMENT"           << SEP 
            << "ID"                << SEP
            << "NAME"              << SEP
            << "BORDER ELEMENTS"   << SEP
            << "INCIDENT ELEMENTS" << SEP
            << "NEIGHBORS"         << SEP 
            << "SIZE"              << SEP
            << "PROJECTED SIZE"    << SEP
            << "ASPECT RATIO"      << SEP
            << "AZIMUT"            << SEP
            << "DIP"               << SEP ;
            // Distance statistics
        out << "MIN DIST"          << SEP
            << "MAX DIST"          << SEP 
            << "W. AV. DIST"       << SEP
            << "W. STD DEV DIST"   << SEP
            << " < 1m"            << SEP 
            << " < 10m"           << SEP 
            << " < 100m"          << SEP 
            << "Total size"      << SEP ;
                 
            // Angle statistics
        out << "MIN ANGLE"         << SEP
            << "MAX ANGLE"         << SEP
            << "W. AV. ANGLE"       << SEP
            << "W. STD DEV ANGLE"   << SEP
            << " < 10 deg"        << SEP 
            << " > 170 deg"       << SEP
            << "Total size"      << SEP
            << std::endl ;
    }
    
    /*! Print the complexity for a region 
     *  \todo consider other possible types 
     */
    void BoundaryModelElement::print_complexity( std::ostream& out ) const {

        out << "Region" << SEP 
            << id_ << SEP
            << name_ << SEP
            << nb_boundary_elements() << SEP
            << "" << SEP
            << compute_neighbors(ALL).size() << SEP 
            << "" << SEP
            << "" << SEP
            << "" << SEP
            << "" << SEP
            << "" << SEP ;

        std::vector< std::pair< double, double > > values ;
  
        compute_distances( values ) ;
        print_stats( out, values, 1., 10., 100. ) ;

        compute_angles( values ) ;
        print_stats( out, values, 10., -DBL_MAX, -DBL_MAX, 170. ) ;
           
        out << std::endl ;
    }

    /*! Returns true if this element or one of the element containing it
     *  is on the Volume Of Interest
     *  This info is strored in the type of the element
     */
    bool BoundaryModelElement::is_on_voi() const
    {
        if( type_ == ALL ) 
        {
            for( int j = 0; j < nb_in_boundary(); ++j ) {
                GEOL_FEATURE t = in_boundary( j )->type() ;
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

    void get_all_boundary_elements(
        const BoundaryModelElement* e,
        std::vector< const BoundaryModelElement* >& result )
    {
        result.clear() ;
        std::set< const BoundaryModelElement* > b1 ;
        for( int i = 0; i < e->nb_boundaries(); ++i ) {
            const BoundaryModelElement* b = e->boundary( i ) ;
            result.push_back( b ) ;
            for( unsigned int j = 0; j < b->nb_boundaries(); b++ ) {
                b1.insert( b->boundary( j ) ) ;
            }
        }

        std::set< const BoundaryModelElement* > b2 ;
        for( std::set< const BoundaryModelElement* >::const_iterator it( b1.begin() );
            it != b1.end(); ++it ) {
            for( unsigned int i = 0; i < ( *it )->nb_boundaries(); i++ ) {
                b2.insert( ( *it )->boundary( i ) ) ;
            }
        }        
        result.insert( result.end(), b1.begin(), b1.end() ) ;
        result.insert( result.end(), b2.begin(), b2.end() ) ;        
    }


    
    /*! Do the 2 elements share a boundary ? 
     */
    bool are_connected(
        const BoundaryModelElement* e1,
        const BoundaryModelElement* e2,
        bool check_all )
    {

        if( !check_all ) {
            // Check only if the elements share a boundary
            for( int i = 0; i < e1->nb_boundaries(); ++i ) {
                const BoundaryModelElement* b1 = e1->boundary( i ) ;
                for( int j = 0; j < e2->nb_boundaries(); ++j ) {
                    if( e2->boundary( i ) == b1 ) return true ;
                }
            }
        }
        else {
            std::vector< const BoundaryModelElement* > b1 ;
            std::vector< const BoundaryModelElement* > b2 ;

            get_all_boundary_elements( e1, b1 ) ;
            get_all_boundary_elements( e2, b2 ) ;
           
            for( int i = 0; i < b1.size(); ++i ){
                if( std::count( b2.begin(), b2.end(), b1[i] ) > 0 ) return true ;
            }
        }
        return false ;
    }

    /*! Approximated distances between boundaries 
     *  Distances are measured from the barycenters of the triangles/segments/ points 
     *  to the triangles/segments/points on other boundaries
     *     
     *  For volumes distances are only measured from the horizons to the other surfaces (because smapling is correct)
     *  For surfaces distances are measured from contacts that are not on the VOI
     * 
     *  Distance between intersecting surface tends toward 0.
     */
    void BoundaryModelElement::compute_distances( std::vector< std::pair< double, double > >& values ) const {
    
        if( dim_ == 0 || dim_ == 1 ) return ;

        // If the element has no or only one boundary get out
        if( boundaries_.size() < 2 ) {
            values.clear() ;
            return ;
        }

        // Allocate room for the values vector
        int values_size = 0 ;
        std::vector< int > nb_simplex ;
        std::vector< const BoundaryModelElement* > from ;
        for( int i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;

            // Skip surfaces that are not horizons( unconf. )  when 
            // the element is a region
            if( dim_ == 3 && b->type() != STRATI ) continue ;
            if( dim_ == 2 && b->is_on_voi() ) continue ;

            from.push_back( b ) ;            
            nb_simplex.push_back( b->nb_cells() ) ;
            values_size += nb_simplex.back() ;
        }
     
        values.resize( values_size ) ;

        int count = 0 ;

        std::vector< const BoundaryModelElement* > to ;
        for( int i = 0; i < from.size(); ++i ){
            const BoundaryModelElement* b = from[i] ;
            
            to.clear() ;
        
            for( int j = 0; j < nb_boundaries(); ++j ){
                const BoundaryModelElement* b1 = boundary( j ) ;

                if( b == b1
                    || (b->is_on_voi() && b1->is_on_voi())
                    || ( dim_ == 3 && b->parent() != nil && b->parent() == b1->parent() )
                  ) continue ;
                else to.push_back( b1 ) ;            
            }

            for( int j = 0; j < nb_simplex[i]; ++j ) {
                double min = DBL_MAX ;
                for( int k = 0; k < to.size(); ++k ) {
                    double cur = b->distance(j, to[k]) ;
                    min = min < cur ? min : cur ;
                }
                values[count+j].first = min ;
                values[count+j].second = b->cell_size( j ) ;
            }

            count += nb_simplex[i] ;
        }
    }

    /*! Fills values with angles measured between the boundaries
     *  A weight (length of a segment or fixed at 1 for corner) is associated to the
     *  angle which is in degree between 0 and 180 
     * 
     *  Angles between two surfaces that have the same parent are ignored
     *  
     *  Returns 999 if the computation fails
     */
    void BoundaryModelElement::compute_angles( std::vector< std::pair< double, double > >& values ) const {

        // It is rather difficult to get the size of the values vector 
        // here. Should not be too big, and greatly inferior than when used for distances just before

        // If the element has no or only one boundary get out
        if( boundaries_.size() < 2 ) {
            values.clear() ;
            return ;
        }

        values.clear() ;
        for( int i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            for( int j = i + 1; j < nb_boundaries(); ++j ) {
                const BoundaryModelElement* b1 = boundary( j ) ;
                bool same_side = false ;
                if( dim_ == 3 && sides_[i] == sides_[j] ) same_side = true ;
                if( are_connected( b, b1, false )
                    && !( dim_ == 3 && b->parent() != nil
                        && b->parent() == b1->parent() ) ) {
                    // Add the values for this contact
                    b->angles( b1, values, same_side ) ;
                }
            }
        }                      
    }

    int BoundaryModelElement::nb_boundary_elements( bool exclude_voi ) const
    {
        unsigned int result = 0 ;
        std::set< const BoundaryModelElement* > b2 ;
        for( int i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            if( exclude_voi && b->is_on_voi() ) continue ;
            if( b->dim() == 0 ){
                // Remove false corners should be removed on closed lines 
                const Corner* c = dynamic_cast< const Corner* >( b ) ;
                if( !c->is_real() ) continue ;
            }
            // Otherwise add an element and check its boundaries
            result++ ;

            for( int j = 0; j < b->nb_boundaries(); ++j ) {
                const BoundaryModelElement* bb = b->boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                if( bb->dim() == 0 ){
                    const Corner* c = dynamic_cast< const Corner* >( bb ) ;
                    if( !c->is_real() ) continue ;
                }

                b2.insert( bb ) ;
            }
        }
        result += b2.size() ;

        std::set< const BoundaryModelElement* > b3 ;
        for( std::set< const BoundaryModelElement* >::iterator it( b2.begin() );
            it != b2.end(); ++it ) {
            for( int j = 0; j < ( *it )->nb_boundaries(); ++j ) {
                const BoundaryModelElement* bb = ( *it )->boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                if( bb->dim() == 0 ) {
                    const Corner* c = dynamic_cast< const Corner* >( bb ) ;
                    if( !c->is_real() ) continue ;
                }
           
                else b3.insert( bb ) ;
            }          
        }
        result += b3.size() ;
        return result ;
    }

    int BoundaryModelElement::nb_incident_elements( bool exclude_voi ) const
    {
        int result = 0 ;
               
        std::set< const BoundaryModelElement* > b2 ;
        for( int i = 0; i < nb_in_boundary(); ++i ) {
            const BoundaryModelElement* b = in_boundary( i ) ;
            if( exclude_voi && b->is_on_voi() ) continue ;
            else {
                result++ ;
                for( int j = 0; j < b->nb_in_boundary(); ++j ) {
                    const BoundaryModelElement* bb = b->in_boundary( j ) ;
                    if( exclude_voi && bb->is_on_voi() ) continue ;
                    else b2.insert( bb ) ;
                }
            }
        }            
        result += b2.size() ;

        std::set< const BoundaryModelElement* > b3 ;
        for( std::set< const BoundaryModelElement* >::iterator it( b2.begin() );
            it != b2.end(); ++it ) {
            for( int j = 0; j < ( *it )->nb_in_boundary(); ++j ) {
                const BoundaryModelElement* bb = ( *it )->in_boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                else b3.insert( bb ) ;
            }           
        }
        result += b3.size() ;
        return result ;
    }

    void Corner::print_complexity( std::ostream& out ) const {
        if( !is_real() ) return ;

        out << "Corner" << SEP 
            <<  id_ << SEP
            << "" << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP

            // All other fields are empty
            << std::endl ;
    }

    double Line::size() const {
        double result = 0. ;
        for( unsigned int i = 1; i < vertices_.size(); ++i ){
            result += length( point( i )-point( i-1 ) ) ;
        }
        // Ca a changé ça
        if( is_closed() ) result += length( model_->point( vertices_.back() )-point( 0 ) );
        return result ;
    }

    double Line::cell_size( int i ) const {
        if( i < vertices_.size()-1 ) return length(point( i+1 ) - point( i)) ;
        else {
            assert( i < vertices_.size() ) ;
            return length( model_->point( vertices_.back() )-point( 0 ) );
        }
    }


    /* Min distance betweeen the given point and the segments of this contact part
     */
    double Line::distance( const vec3& p ) const {
        double result = FLT_MAX ;            
        for( unsigned int i = 1; i < vertices_.size(); ++i ){
            // Distance betweena a point and a segment COPY from smwh else
            const vec3& p0 = point( i-1 ) ;
            const vec3& p1 = point( i ) ;

            double distance_pt_2_segment  = FLT_MAX ;
            vec3 c = (p1-p0)/2 ;        
            double half = ::Geomesh::distance( p1, c ) ;
            double cp_dot_p0p1 = dot( p-c, p1-p0 ) ;

            if( cp_dot_p0p1 < -half ) distance_pt_2_segment =  ::Geomesh::distance( p0, p ) ;
            else if( cp_dot_p0p1 > half ) distance_pt_2_segment = ::Geomesh::distance( p1, p ) ;
            else {
                vec3 projection = c + cp_dot_p0p1*(p1-p0) ;
                distance_pt_2_segment = ::Geomesh::distance( projection, p ) ;
            }    
            result = distance_pt_2_segment < result ? distance_pt_2_segment : result ;
        }
        if( is_closed() ) { // Plus code ligne changé bon 
            // COPY BEUUURKKK !!
            const vec3& p0 = model_->point( vertices_.back() ) ;
            const vec3& p1 = point( 0 ) ;

            double distance_pt_2_segment  = FLT_MAX ;
            vec3 c = (p1-p0)/2 ;        
            double half = ::Geomesh::distance( p1, c ) ;
            double cp_dot_p0p1 = dot( p-c, p1-p0 ) ;

            if( cp_dot_p0p1 < -half ) distance_pt_2_segment = ::Geomesh::distance( p0, p ) ;
            else if( cp_dot_p0p1 > half ) distance_pt_2_segment = ::Geomesh::distance( p1, p ) ;
            else {
                vec3 projection = c + cp_dot_p0p1*(p1-p0) ;
                distance_pt_2_segment = ::Geomesh::distance( projection, p ) ;
            }    
            result = distance_pt_2_segment < result ? distance_pt_2_segment : result ;
        }
        return result ;
    }
    vec3 Line::average_orientation() const {
        vec3 s(0., 0., 0. ) ;
        for( unsigned int i = 1; i < vertices_.size(); ++i ){
            s += point( i-1 ) - point( i );
        }
        if( is_closed() ) { // to update cf Line
            s += point( 0 )-model_->point( vertices_.back() ) ;
        }
        return normalize( s ) ;
    }


    /*! Minimal distance betweeen the vertices of this element to
     *  the triangle of e.
     *  Totally inefficient but not the priority right now
     */
    double Line::distance( BoundaryModelElement* e ) const {
        double result = FLT_MAX ;
        for( unsigned int i = 0; i < vertices_.size(); ++i ) {
            double d = e->distance( point( i ) ) ;
            result = result < d ? result : d ;
        }
        return result ;
    }

    /*! Min distance fron the i-th simplex barycenter to the given elements
     */
    double Line::distance( int s, BoundaryModelElement* to ) const {         
        vec3 centroid ;
        if( s < vertices_.size()-1 ) {
            centroid = ( point( s+1 )+point( s )) * 0.5 ;
        }
        else {
            assert( s == vertices_.size()-1 ) ;
            centroid = (point( 0 )+model_->point( vertices_.back())) * 0.5 ;
        }        
        return to->distance( centroid ) ;        
    }

    /*! Add the measured angle with a weight at 1. for the corners at
     *  which this contact intersect "with"
     */
    void Line::angles( 
        BoundaryModelElement* in, 
        std::vector< std::pair< double, double > >& values,
        bool same_side
    ) const {
        Line* cp = dynamic_cast< Line* >( in ) ;
        if( cp == nil ) return ;


        std::vector< const BoundaryModelElement* > shared ;
        for( int i = 0; i < nb_boundaries(); ++i ){
            for( int j = 0; j < cp->nb_boundaries(); ++j ){
                if( cp->boundaries_[j] == boundaries_[i] ) {
                    shared.push_back( cp->boundary( j ) ) ;
                }
            }
        }

        for( int i = 0; i < shared.size(); ++i ){
            const Corner* c = dynamic_cast< const Corner* >( shared[i] ) ;
            assert( c != nil ) ;

            const vec3& p = c->point() ; 

            vec3 e1 = point( 1 ) - point( 0 ) ;
            if( p != point( 0 ) ){
                assert( p == model_->point( vertices_.back() ) ) ;
                e1 =  point( vertices_.size()-2 )- model_->point( vertices_.back() ) ;
            }
            vec3 e2 = cp->point( 1 ) - cp->point( 0 ) ;
            if( p != cp->point( 0 ) ) {
                assert( p == cp->point( cp->nb_points()-1 ) ) ;
                e2 = cp->point( cp->nb_points()-2 ) - cp->point( cp->nb_points()-1 ) ;
            }
            e1 = normalize( e1 ) ;
            e2 = normalize( e2 ) ;

            double a = std::acos( dot(e1,e2) ) * 180. / M_PI ;
            values.push_back( std::pair< double, double >( a, 1. ) ) ;
        }

    }

    double Line::min_angle( BoundaryModelElement* in ) const {
        double result = 999. ;

        Line* cp = dynamic_cast< Line* >( in ) ;
        if( cp != nil ) {

            std::vector< const BoundaryModelElement* > shared ;
            for( int i = 0; i < nb_boundaries(); ++i ){
                for( int j = 0; j < cp->nb_boundaries(); ++j ){
                    if( cp->boundaries_[j] == boundaries_[i] ) {
                        shared.push_back( cp->boundary( j ) ) ;
                    }
                }
            }
            
            for( int i = 0; i < shared.size(); ++i ){
                const Corner* c = dynamic_cast< const Corner* >( shared[i] ) ;
                assert( c != nil ) ;

                const vec3& p = c->point() ; 

                vec3 e1 = point( 1 ) - point( 0 ) ;
                if( p != point( 0 ) ){
                    assert( p == model_->point( vertices_.back() ) ) ;
                    e1 =  point( vertices_.size()-2 )- model_->point( vertices_.back() ) ;
                }
                vec3 e2 = cp->point( 1 ) - cp->point( 0 ) ;
                if( p != cp->point( 0 ) ) {
                    assert( p == cp->point( cp->nb_points()-1 ) ) ;
                    e2 = cp->point( cp->nb_points()-2 ) - cp->point( cp->nb_points()-1 ) ;
                }
                e1 = normalize( e1 ) ;
                e2 = normalize( e2 ) ;
                
                double a = std::acos( dot(e1,e2) ) * 180. / M_PI ;
                //if( dot(e1,e2) < 0 ) a = std::acos( dot(-e1,e2) ) * 180. / Pi ;

                result = result < a ? result : a ;
            }                
        }         
        return result ;     
    }

    void Line::print_complexity( std::ostream& out ) const {
        out << "Contact"           << SEP 
            << id_                << SEP
            << ""              << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP
            << compute_neighbors(ALL).size() << SEP 
            << size()              << SEP ;

        double d = 0. ;
        if( !is_closed() ){ // to update
            d = length(point(0) - model_->point( vertices_.back())) ;
        }
        else {
            // Compute the max distance between 2 point divided by two
            for( unsigned int i = 0; i < vertices_.size() ; ++i ) {
                for( unsigned int j = 0; j < vertices_.size() ; ++j ) {
                    double dd = length2(point(i)- point(j) ) ;
                    d =  d > dd ? d : dd ;
                }
            }
            d = sqrt(d)/2. ;
        }
        double aspect = size() ;
        if( d > 10e-30 ) aspect /= d ;
        else aspect = FLT_MAX ;

        out << d     << SEP 
            << aspect << SEP               
            << average_angle_to_y() << SEP
            << average_dip()        << SEP


            << std::endl ;
    }

     double Surface::size() const {
        double result = 0. ;
        for( unsigned int i = 0; i < nb_cells(); i++ ) {
            result += triangle_area( point( i, 0 ), point( i, 1 ),
                point( i, 2 ) ) ;
            if( !is_triangle( i ) ) {
                result += triangle_area( point( i, 0 ), point( i, 2 ),
                    point( i, 3 ) ) ;
            }
        }
        return result ;
    }

    double Surface::cell_size( int f ) const
    {
        grgmesh_assert( f < nb_cells() ) ;

        double result = 0 ;
        const vec3& p0 = point( f, 0 ) ;
        for( unsigned int i = 1; i+1 < nb_points_in_facet(f); ++i ){
            result += triangle_area( p0, point(f,i), point(f,i+1) ) ;
        }
        return result ;
    }

    double Surface::distance( const vec3& p ) const {
        double result = FLT_MAX ;
        
        grgmesh_assert_not_reached ;
        // To implement for polygonal facets

        for( unsigned int i = 0; i < nb_cells(); i++ ) {
            double cur_result = point_triangle_squared_distance( p, point( i, 0 ),
                    point( i, 1 ), point( i, 2 ) ) ;
            if( !is_triangle( i ) ) {
                cur_result += point_triangle_squared_distance( p, point( i, 0 ), point( i, 2 ),
                    point( i, 3 ) ) ;
                cur_result /= static_cast< double >( 2.0 ) ;
            }
            result =  cur_result < result ? cur_result : result ;
        }
        if( result != FLT_MAX ) result = sqrt( result ) ;
        return result ;
    }

    vec3 Surface::closest_normal( const vec3& p, vec3& b ) const {
       vec3 result ;
       
       grgmesh_assert_not_reached ;
       // To implement for polygonal facets
       
        double min_dist = FLT_MAX ;
       
       for( unsigned int i = 0; i < nb_cells(); i++ ) {
            double cur_dist = point_triangle_squared_distance( p, point( i, 0 ),
                    point( i, 1 ), point( i, 2 ) ) ;
            if( !is_triangle( i ) ) {
                cur_dist += point_triangle_squared_distance( p, point( i, 0 ), point( i, 2 ),
                    point( i, 3 ) ) ;
                cur_dist /= static_cast< double >( 2.0 ) ;
            }
            if( cur_dist < min_dist) {
                min_dist = cur_dist ;
                result = facet_normal( i ) ;
                b = point(i, 0 ) ;
            }

        }       
        return result ;
    }
   
    /*! Minimal distance betweeen the vertices of this element to
     *  the triangle of e.
     *
     *  Totally inefficient but not the priority right now
     */
    double Surface::distance( BoundaryModelElement* e ) const {
        double result = FLT_MAX ;

        for( unsigned int i = 0; i < points_.size(); ++i ) {
            double d =  e->distance( model_->point( points_[i] ) ) ;
            result = result < d ? result : d ;
        }
        return result ;
    }

    double Surface::distance( int t, BoundaryModelElement* to ) const {
        assert( t < nb_cells() ) ;

        grgmesh_assert_not_reached ;
        // To implement for polygonal facets

        if( is_triangle( t ) ) {
            return to->distance( point( t, 0 ) + point( t, 1 ) + point( t, 2 ) ) / 3. ;
        } else {
            return to->distance(
                point( t, 0 ) + point( t, 1 ) + point( t, 2 ) + point( t, 3 ) ) / 4. ;
        }
    }

    vec3 Surface::facet_normal( int t ) const {
        const vec3& p0 = point( t, 0 )  ;
        const vec3& p1 = point( t, 1 )  ;
        const vec3& p2 = point( t, 2 )  ;
        vec3 c0 = cross(p0-p2, p1-p2) ;
        if( !is_triangle(t) ) {
            grgmesh_assert_not_reached ;
            // To implement for polygonal facet -
            // OK if they are planar, but if they are not.... What is the normal            
        }
        return normalize( c0 ) ;
    }

    vec3 Surface::average_orientation() const {
        double total_a = 0 ;
        vec3 result (0,0,0) ;
        for( unsigned int t = 0; t < nb_cells() ; ++t){
            vec3 n = facet_normal( t ) ;
            double a = cell_size( t ) ;

            result += a*n ;
            total_a += a ;
        }
        if( total_a > 10e-30 ) return result/total_a ;
        else return vec3(0.,0.,0.) ;
    }

    void Surface::angles( 
        BoundaryModelElement* in , 
        std::vector< std::pair< double, double > >& values,
        bool same_side
    ) const {
        if( in == this ) return ;
        Surface* sp = dynamic_cast< Surface* >( in ) ;
        if( sp == nil ) return ;

        std::vector< const BoundaryModelElement* > shared ;
        for( int i = 0; i < nb_boundaries(); ++i ) {
            for( int j = 0; j < sp->nb_boundaries(); ++j ) {
                if( sp->boundaries_[j] == boundaries_[i] ) {
                    shared.push_back( sp->boundary( j ) ) ;
                }
            }
        }

        for( int i = 0; i < shared.size(); ++i ){
            const Line* cp = dynamic_cast< const Line* >( shared[i] ) ;
            assert( cp != nil ) ;

            // Find the triangles sharing each segment 
            for( unsigned int j = 1; j < cp->nb_points(); ++j ){
                const vec3& p0 = cp->point( j-1 ) ;
                const vec3& p1 = cp->point( j ) ;

                int t1 = facet_from_surface_point_ids( p0, p1 ) ;
                int t2 = sp->facet_from_surface_point_ids( p0, p1 ) ;
                if( t1 == -1 || t2 == -1 ) {
                    // Should not happen
                    continue ;
                }
                // Get the angle between these triangles
                double d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
           
                double a = std::acos( d ) * 180 / M_PI ;
                if( same_side ) a = 180. - a ;
                values.push_back( std::pair< double, double >(a, length(p0-p1)) ) ;
            }
            // If the line is closed check the closing segment (Copy)
            if( cp->is_closed() ) { // to update
                const vec3& p0 = cp->point( 0 ) ;
                const vec3& p1 = cp->point( cp->nb_points()-2 ) ;

                int t1 = facet_from_surface_point_ids( p0, p1 ) ;
                int t2 = sp->facet_from_surface_point_ids( p0, p1 ) ;
                if( t1 == -1 || t2 == -1 ) {
                    // Should not happen
                    continue ;
                }

                double d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                 
                double a = std::acos( d ) * 180 / M_PI ;
                if( same_side ) a = 180. - a ;
                values.push_back( std::pair< double, double >( a, length(p0-p1) ) ) ;
            }
        }
    }

    /*! If e is a Surface that share a boundary line with this Surface
     *  returns the min angle between the two
     *  else return 999.
     */
    double Surface::min_angle( BoundaryModelElement* in ) const {
        double result = 999. ;
        if( in == this ) return result ;

        Surface* sp = dynamic_cast< Surface* >( in ) ;
        if( sp != nil ) {
            // Find the contact shared by the two
            std::vector< const BoundaryModelElement* > shared ;
            for( int i = 0; i < nb_boundaries(); ++i ) {
                for( int j = 0; j < sp->nb_boundaries(); ++j ) {
                    if( sp->boundaries_[j] == boundaries_[i] ) {
                        shared.push_back( sp->boundary( j ) ) ;
                    }
                }
            }
            
            for( int i = 0; i < shared.size(); ++i ){
                const Line* cp = dynamic_cast< const Line* >( shared[i] ) ;
                assert( cp != nil ) ;
                
                // Find the triangle sharing each segment 
                 for( unsigned int j = 1; j < cp->nb_points(); ++j ){
                    const vec3& p0 = cp->point( j-1 ) ;
                    const vec3& p1 = cp->point( j ) ;

                    int t1 = facet_from_surface_point_ids( p0, p1 ) ;
                    int t2 = sp->facet_from_surface_point_ids( p0, p1 ) ;
                    if( t1 == -1 || t2 == -1 ) {
                        // Should not happen
                        continue ;
                    }
                    // Get the angle between these triangles
                    double d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                    //if( d < 0 ) d = dot(-triangle_normal(t1), sp->triangle_normal(t2)) ;                   
                    double a = std::acos( d ) * 180 / M_PI ;

                    result = result < a ? result : a ;
                 }
                 // If the line is closed check the closing segment (Copy)
                 if( cp->is_closed() ) { // to update
                    const vec3& p0 = cp->point( 0 ) ;
                    const vec3& p1 = cp->point( cp->nb_points()-2 ) ;

                    int t1 = facet_from_surface_point_ids( p0, p1 ) ;
                    int t2 = sp->facet_from_surface_point_ids( p0, p1 ) ;
                    if( t1 == -1 || t2 == -1 ) {
                        // Should not happen
                        continue ;
                    }

                    double d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                    //if( d < 0 ) d = dot(-triangle_normal(t1), sp->triangle_normal(t2)) ;                   
                    double a = std::acos( d ) * 180 / M_PI ;

                    result = result < a ? result : a ;
                 }
            }                
        }
         
        return result ;
    }

    void Surface::print_complexity( std::ostream& out ) const {
        out << "Surface"           << SEP 
            << id_                << SEP
            << parent()->name()              << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP
            << compute_neighbors(ALL).size() << SEP 
            << size()              << SEP 
            << ""    << SEP
            << ""      << SEP
            
            << average_angle_to_y() << SEP
            << average_dip()        << SEP ;
    

        std::vector< std::pair< double, double > > values ;
  
        compute_distances( values ) ;
        print_stats( out, values, 1., 10., 100. ) ;

        compute_angles( values ) ;
        print_stats( out, values, 10., -DBL_MAX, -DBL_MAX, 170. ) ;
           
        out << std::endl ;

    }

    void Surface::print_mesh( const std::string& filename ) const
    {
        std::ofstream file( filename.c_str(), std::ios::trunc | std::ios::out ) ;
        file << "Surface : " << id() << std::endl ;
        file << "========== Points =========" << std::endl ;
        for( unsigned int p = 0; p < nb_points(); p++ ) {
            file << p << " -> " << point( p ) << std::endl ;
        }
        file << "========== Facets =========" << std::endl ;
        for( unsigned int p = 0; p < nb_cells(); p++ ) {
            file << p << " ->" ;
            for( unsigned int v = 0; v < nb_points_in_facet(p); v++ ) {
                file << " " << surf_point_id( p, v ) ;
            }
            file << std::endl ;
        }
        file << "========== Facet ptr =========" << std::endl ;
        for( unsigned int p = 0; p < nb_cells(); p++ ) {
            file << p << " -> " << facet_begin(p) << " " << facet_end(p) << std::endl ;
        }
        file << "========== Adjacents =========" << std::endl ;
        for( unsigned int p = 0; p < nb_cells(); p++ ) {
            file << p << " ->" ;
            for( unsigned int v = 0; v < nb_points_in_facet(p); v++ ) {
                file << " " << adjacent( p, v ) ;
            }
            file << std::endl ;
        }

    }

    void Surface::point_normal( std::vector< vec3 >& normals ) const
    {
        normals.resize( nb_points() ) ;
        for( unsigned int f = 0; f < nb_cells(); f++ ) {
            vec3 normal = facet_normal( f ) ;
            for( unsigned int p = 0; p < nb_points_in_facet( f ); p++ ) {
                unsigned int id = surf_point_id( f, p ) ;
                normals[id] += normal ;
            }
        }
        for( unsigned int p = 0; p < nb_points(); p++ ) {
            normals[p] = normalize( normals[p] ) ;
        }
    }
#endif
}

