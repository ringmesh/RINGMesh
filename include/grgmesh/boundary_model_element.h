/*[
* Association Scientifique pour la Geologie et ses Applications (ASGA)
* Copyright (c) 1993-2013 ASGA. All Rights Reserved.
*
* This program is a Trade Secret of the ASGA and it is not to be:
* - reproduced, published, or disclosed to other,
* - distributed or displayed,
* - used for purposes or on Sites other than described
*   in the GOCAD Advancement Agreement,
* without the prior written authorization of the ASGA. Licencee
* agrees to attach or embed this Notice on all copies of the program,
* including partial copies or modified versions thereof.
]*/ 
 
/*! \author Jeanne Pellerin */

#ifndef __GRGMESH_BOUNDARY_MODEL_ELEMENT__
#define __GRGMESH_BOUNDARY_MODEL_ELEMENT__

#include <grgmesh/common.h>
#include <grgmesh/vecn.h>
#include <vector> 
#include <map>
#include <string>
#include <algorithm> 
#include <cmath>

namespace GRGMesh {
    class BoundaryModel ;
    class Box3d ;
}


namespace GRGMesh {
/********************************************************************************************/
/******             General stuff               *********************************************/
/********************************************************************************************/   
   
    // I do not know where to put these
    // if in BoundaryModel.h it poses problems at the pre-compilation step with gomgen JP
 
    // For CSV file, but could be changed easily    
    static std::string SEP = "," ;
    static std::string SEP_1 = "   " ;
    
    static float64 EPSILON2 = 10e-6 ;
     

    template< class T > void print_vector_contents( const std::vector< T >& v, std::ostream& out ){
        for( uint32 i =0; i < v.size(); ++i ){
            out << v[i] << SEP_1 ;
        }
        out << SEP ;
    }  

   /* Intermediate structure to build the contacts when loading a
       BoundaryModel from .ml file 
    */
    struct Border {
        Border( int32 part, int32 corner, int32 p0, int32 p1):
            part_id_(part), corner_id_(corner), p0_(p0), p1_(p1) {};
        
        int32 part_id_ ; // Id of the SurfacePart owning this Border
        int32 corner_id_ ; // id of p0 in the BoundaryModel corner vector

        // Ids of the starting corner and second point on the border in the Map
        // to which this Border belong
        int32 p0_ ;
        int32 p1_ ;
    } ;

    static uint32 e2v_for_triangle[3][2] = {
        {2,1},
        {0,2},
        {1,0}
    } ;
    static uint32 e2v_for_quad[4][2] = {
        {1,0},
        {2,1},
        {3,2},
        {0,3}
    } ;


    template< class T > static bool vector_contains( const std::vector< T > cont, const T& value ) {
        for( uint32 i = 0; i < cont.size(); i++ ) {
            if( cont[i] == value ) {
                return true ;
            }
        }
        return false ;
    }

    struct KeyFacet {
        KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 ):
            p0_(p0), p1_(p1), p2_(p2){} ;
        KeyFacet():p0_(), p1_(), p2_() {} ;
        bool is_default() const {
            if( p0_ == vec3() && 
                p1_ == vec3() && 
                p2_ == vec3() ) return true ;
            else {
                grgmesh_assert( p0_ != p1_ && p0_ != p2_ && p1_ != p2_ ) ;
                return false ;
            }
        }
       
    public:
        vec3 p0_ ;
        vec3 p1_ ;
        vec3 p2_ ;
    } ;

       
    enum GEOL_FEATURE {
        ALL,
        STRATI,
        FAULT,
        VOI,
        STRATI_FAULT,
        STRATI_VOI,
        FAULT_VOI
    } ;    
    static GEOL_FEATURE default_type = ALL ;

    /********************************************************************************************/
    /***********       Components of a Model declaration                 ************************/
    /********************************************************************************************/

    static std::vector< vec3 > empty_vector ;
    static std::vector< uint32 > empty_uint_vector ;
    static vec3 dummy = vec3( 0, 0, 0 ) ;
    /*! Base manifold element in a B-Rep model
     */
    class GRGMESH_API BoundaryModelElement {
        friend class BoundaryModelBuilder ;
    public:
        BoundaryModelElement(
            BoundaryModel* model,
            int32 dim = -1,
            int32 id = -1,
            int32 parent = -1,
            const GEOL_FEATURE& type = default_type )
            :
                model_( model ),
                name_( "" ),
                id_( id ),
                dim_( dim ),
                type_( type ),
                parent_( parent )
        {
        }
        virtual ~BoundaryModelElement()
        {
        }

        // Accessors
        const BoundaryModel* model() const { return model_ ; }
        const std::string& name() const { return name_ ; }
        int32 id() const { return id_ ; }
        int32 dim() const { return dim_ ; }
        GEOL_FEATURE type() const { return type_ ; }
        bool is_on_voi() const ;
        bool side( int32 i ) const { return sides_[i] ; }
        bool has_parent() const { return parent_ != -1 ; }
        const BoundaryModelElement* parent() const ;
        int32 parent_id() const { return parent_ ; }
        uint32 nb_boundaries() const { return boundaries_.size() ; }
        uint32 boundary_id( uint32 x ) const { return boundaries_[x] ; }
        const BoundaryModelElement* boundary( uint32 x ) const ;
        uint32 nb_in_boundary() const { return in_boundary_.size() ; }
        uint32 in_boundary_id( uint32 x ) const { return in_boundary_[x] ; }
        const BoundaryModelElement* in_boundary( uint32 x ) const ;
        uint32 nb_children() const { return children_.size() ; }
        uint32 child_id( uint32 x ) const { return children_[x] ; }
        const BoundaryModelElement* child( uint32 x ) const ;
        bool contains( const BoundaryModelElement* in ) const {
            return this == in || std::count( children_.begin(), children_.end(), in->id() ) > 0 ; }
        virtual uint32 nb_simplices() const ;
        virtual const vec3& vertex( uint32 p = 0 ) const {
            grgmesh_assert_not_reached ;
            return dummy ;
        }
        virtual uint32 nb_vertices() const {
            grgmesh_assert_not_reached ;
            return 0 ;
        }
        virtual bool is_triangulated() const {
            grgmesh_debug_assert( dim_ == 3 ) ;
            for( uint32 s = 0; s < nb_boundaries(); s++ ) {
                if( !boundary( s )->is_triangulated() ) return false ;
            }
            return true ;
        }
        
        // Modifiers
        void set_parent( int32 p ){ parent_ = p ; }
        void set_name( const std::string& name ) { name_ = name ; }
        void set_type( GEOL_FEATURE type ) { type_ = type ; } 
        void set_dim( int32 dim ) { dim_ = dim ; }
        void set_id( int32 id ) { id_ = id ; }


        // Functions
        std::vector< const BoundaryModelElement* > compute_neighbors(
            GEOL_FEATURE through, bool exclude_boundaries = true ) const ;
        virtual float64 size () const ;
        virtual float64 simplex_size( int32 /*i*/ ) const { return 0. ; }
        virtual float64 distance( const vec3& p ) const ;
        virtual float64 distance( BoundaryModelElement* e ) const ;
        virtual float64 distance( int32 /*simplex*/, const BoundaryModelElement* /*to*/ ) const { return -1 ; }
        virtual float64 min_angle( BoundaryModelElement* e ) const ;
        virtual void angles(
            const BoundaryModelElement* /*with*/,
            std::vector< std::pair< float64, float64 > >& /*values*/,
            bool /*same_side*/ ) const {} ;
        virtual vec3 average_orientation() const ;      
        float64 average_angle_to_y() const ;
        float64 average_dip() const ;


        static void print_categories( std::ostream& out ) ;
        void print( std::ostream& out ) const ;      

        static void print_complexity_categories( std::ostream& out ) ;
        virtual void print_complexity( std::ostream& out ) const ;

        void compute_distances( std::vector< std::pair< float64, float64 > >& values ) const ;
        void compute_angles( std::vector< std::pair< float64, float64 > >& values ) const ;
        int32 nb_incident_elements( bool exclude_voi = true ) const ;
        int32 nb_boundary_elements( bool exclude_voi = true ) const ;

    protected:
        void copy_macro_topology(
            const BoundaryModelElement& rhs,
            BoundaryModel& model ) ;

        void add_boundary( uint32 b ) { boundaries_.push_back( b ) ; }
        void add_oriented_boundary( uint32 b, bool side ) { boundaries_.push_back(b) ; sides_.push_back(side) ; }
        virtual void add_in_boundary( uint32 e ) { in_boundary_.push_back(e) ; }
        void add_child( uint32 e ){ children_.push_back( e ) ; }
        void change_boundary_side( int32 id ) ;

    protected :
        BoundaryModel* model_ ;
        /// Name of the element if any, default is nothing
        std::string name_ ;
        /// Id of this element in the appropriate vector of the BoundaryModel owning it
        int32 id_ ;
        /// Dimension of the element 0 corner; 1 line; 2 surface; 3 region
        int32 dim_ ;
        /// Geological type for this object, default is ALL 
        GEOL_FEATURE type_ ;

        /// Elements on the boundary whose dim = this->dim - 1
        std::vector< uint32 > boundaries_ ;
        /// Filled for volumes, gives the side + or - of the surface on which the region is
        std::vector< bool > sides_ ; 
        /// Elements on whose boundary this element is dim = this->dim + 1
        std::vector< uint32 > in_boundary_ ;

        /// Id of the geological parent ( geological group of elements to which belong this base element ). If there is is no parent default value is -1.
        int32 parent_ ;
        
        /// Children = the components of this element
        std::vector< uint32 > children_ ;
    } ;


    /*-----------------------------------------------------------------------------------------*/
    /*! Corners of the model
     *  Dimension 0
     *  Store its location p_
     *  Is on the boundaries of contacts, but only the surface parts are known when building it
     */
    class GRGMESH_API Corner : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
    public:
        Corner(
            BoundaryModel* model,
            int32 id = -1,
            uint32 p = 0 )
            : BoundaryModelElement( model, 0, id ), p_( p )
        {
        }
        virtual ~Corner()
        {
        }
        bool is_real() const { return in_boundary_.size() > 1 ; }

        virtual const vec3& vertex( uint32 p = 0 ) const ;
        virtual uint32 nb_vertices() const { return 1 ; }
        virtual float64 size() const { return 0. ; }
        virtual float64 distance( const vec3& p ) const { return (vertex() - p).length() ; }
        virtual float64 distance( BoundaryModelElement* e ) const { return e->distance( vertex() ) ; }
        virtual uint32 nb_simplices() const { return 1 ; }

        virtual void print_complexity( std::ostream& out ) const ;

    private:
        void copy_macro_topology(
            const Corner& rhs,
            BoundaryModel& model ) ;

        void set_corner( uint32 p ) { p_ = p ; }
     private:
        uint32 p_ ;
    };
    /*-----------------------------------------------------------------------------------------*/        
    /*! A part of the contact between 2 interfaces
     * 
     *  boundaries -> 1 or 2 corners
     *  in_boundary -> several surface parts
     */
    class GRGMESH_API ContactPart : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class ContactPartMutator ;
    public:
        ContactPart( BoundaryModel* model, int32 id = -1 ) ;
        ContactPart(
            BoundaryModel* model,
            int32 id,
            const std::vector< uint32 >& points ) ;
        ContactPart(
            BoundaryModel* model,
            int32 id,
            uint32 corner0,
            uint32 corner1,
            const std::vector< uint32 >& points ) ;
        virtual ~ContactPart(){} ;

        bool is_inside_border( int32 in_surface_id ) const { return is_inside_border_[in_surface_id] ; }
        bool is_inside_border( const BoundaryModelElement& e ) const {
            for( uint32 i = 0; i < nb_in_boundary(); i++ ) {
                if( in_boundary_[i] == e.id() ) {
                    return is_inside_border( i ) ;
                }
            }
            return false ;
        }

        virtual const vec3& vertex( uint32 p ) const ;
        virtual uint32 nb_vertices() const { return vertices_.size() ; }
        
        float64 resolution( uint32 p ) const { return resolution_[p] ; }
        void set_resolutions( const std::vector< float64 >& resolutions ) {
            resolution_ = resolutions ;
        }

        bool contains( const vec3& p ) const ;
        int32 find( const vec3& p ) const ;
        int32 point_id( const vec3& p ) const ;
        bool is_closed () const { 
            return (boundaries_[0]!= nil ) && (boundaries_[0] == boundaries_[1]) ; 
        }  

        virtual float64 size() const ;
        virtual float64 simplex_size( int32 i ) const ;
        virtual float64 distance( const vec3& p ) const ;
        virtual float64 distance( BoundaryModelElement* e ) const ;
        virtual float64 distance( int32 simplex, BoundaryModelElement* to ) const ;
        virtual float64 min_angle( BoundaryModelElement* e ) const ;
        virtual vec3 average_orientation() const ;      
        virtual uint32 nb_simplices() const {
            if( is_closed() ) return vertices_.size() ; 
            else return vertices_.size()-1 ; } 
        virtual void angles(
            BoundaryModelElement* with,
            std::vector< std::pair< float64, float64 > >& values,
            bool same_side = false ) const ;
        
        virtual void print_complexity( std::ostream& out ) const ;

    private:
        void copy_macro_topology(
            const ContactPart& rhs,
            BoundaryModel& model ) ;
        void set_vertices( const std::vector< uint32 >& vertices ) {
            vertices_ = vertices ;
        }
        virtual void add_in_boundary( uint32 e ) ;
        void set_is_inside_border( bool x ) { is_inside_border_.push_back(x) ; }
    private:
        /// In case of a closed line, the last vertex (equal to the first) is not stored.
        std::vector< uint32 > vertices_ ;
        std::vector< bool > is_inside_border_ ;
        std::vector< float64 > resolution_ ;
    } ;

    class GRGMESH_API ContactPartMutator {
    public:
        ContactPartMutator( ContactPart& M )
            : M_( M )
        {
        }
        ContactPartMutator( const ContactPart& M )
            : M_( const_cast< ContactPart& >( M ) )
        {
        }
        vec3& point( uint32 p ) const ;
        void set_point( uint32 id, const float64* p ) {
            set_point( id, vec3( p ) ) ;
        }
        void set_point( uint32 id, const vec3& p ) ;
        std::vector< uint32 >& points() const { return M_.vertices_ ; }
        std::vector< float64 >& resolution() const { return M_.resolution_ ; }

        void clear()
        {
            M_.vertices_.clear() ;
            M_.resolution_.clear() ;
        }

    private:
        ContactPart& M_ ;
    };
    /*-----------------------------------------------------------------------------------------*/
    
    struct FacetEdge {
        FacetEdge() : facet_ ( -1 ), edge_( -1 )
        {
        }
        FacetEdge( int32 t ) : facet_ ( t ), edge_( -1 )
        {
        }
        FacetEdge( int32 t, int32 e ) : facet_ ( t ), edge_( e )
        {
        }
        bool operator==( const FacetEdge& rhs ) const {
            return ( facet_ == rhs.facet_ ) && ( edge_ == rhs.edge_ ) ;
        }
        bool operator!=( const FacetEdge& rhs ) const {
            return ( facet_ != rhs.facet_ ) || ( edge_ != rhs.edge_ ) ;
        }
        int32 facet_ ;
        int32 edge_ ;
    } ;

    /*! A SurfacePart dimension 
     * 
     * boundaries_ -> contact parts
     * in_boundary_ -> 1 or 2 regions
     *
     */
    class GRGMESH_API SurfacePart : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class SurfacePartMutator ;
    public:
        SurfacePart(
            BoundaryModel* model,
            int32 id = -1,
            int32 parent = -1,
            const GEOL_FEATURE& type = default_type )
            :
                BoundaryModelElement( model, 2, id, parent, type ),
                is_triangulated_( true )
        {
        }
        virtual ~SurfacePart(){ } ;
       
        void set_key_facet( const KeyFacet& key ) { key_facet_ = key ; }
        
        virtual float64 size() const ;
        virtual float64 simplex_size( int32 i ) const ;
        virtual float64 distance( const vec3& p ) const ;
        virtual float64 distance( BoundaryModelElement* e ) const ;
        virtual float64 distance( int32 simplex, BoundaryModelElement* to ) const ;
        virtual float64 min_angle( BoundaryModelElement* e ) const ;
        virtual void angles( 
            BoundaryModelElement* with,
            std::vector< std::pair< float64, float64 > >& values,
            bool same_side = false ) const ;
        virtual vec3 average_orientation() const ;
        virtual uint32 nb_simplices() const { return facets_.empty() ? 0 : facet_ptr_.size() - 1 ; }

        const KeyFacet& key_facet() const
        {
            return key_facet_ ;
        }
        void set_first_triangle_as_key() ;
        virtual void print_complexity( std::ostream& out ) const ;
        void print_mesh( const std::string& filename ) const ;
        
        virtual bool is_triangulated() const { return is_triangulated_ ; }
        virtual uint32 nb_vertices() const { return points_.size() ; }
        uint32 size_of_facets() const { return facets_.size() ; }
        uint32 facet_begin( int32 f ) const { return facet_ptr_[f] ; }
        uint32 facet_end( int32 f ) const { return facet_ptr_[f+1] ; }
        uint32 nb_points_in_facet( int32 f ) const {
            return facet_end( f ) - facet_begin( f ) ; }
        const vec3& point( int32 f, int32 v ) const ;
        uint32 point_index( int32 f, int32 v ) const { return facets_[facet_begin(f)+v] ; }
        virtual const vec3& vertex( uint32 v ) const ;
        int32 adjacent( int32 f, int32 e ) const { return adjacent_[facet_begin(f)+e] ; }
        int32 adjcent_in_neighbor( uint32 f, uint32 e ) const ;
        bool is_on_border( int32 t, int32 e ) const { return adjacent( t, e ) == -1 ; }
        bool is_on_border( int32 t ) const {
            return is_on_border( t, 0 ) || is_on_border( t, 1 ) || is_on_border( t, 2 ) ;
        }
        FacetEdge next_on_border( const FacetEdge& te ) const ;
        uint32 edge_vertex( int32 f, int32 e, int32 v ) const {
            if( is_triangle( f ) ) {
                return e2v_for_triangle[e][v] ;
            } else {
                return e2v_for_quad[e][v] ;
            }
        }
        bool is_triangle( int32 f ) const { return nb_points_in_facet( f ) == 3 ; }
        bool is_resolution_set() const { return resolution_.size() == nb_vertices() ; }
        bool is_U_set() const { return U_.size() == nb_simplices() ; }
        bool is_V_set() const { return V_.size() == nb_simplices() ; }
        bool is_W_set() const { return W_.size() == nb_simplices() ; }
        bool is_UVW_set() const { return is_U_set() && is_V_set() && is_W_set() ; }

        int32 find_triangle( const vec3& p0, const vec3& p1, const vec3& p2 ) const ;
        int32 find_triangle( int32 id0, int32 id1 ) const ;
        int32 find_triangle( const vec3& p0, const vec3& p1 ) const ;
        int32 point_id( int32 t, int32 p0 ) const ;
        int32 point_id( int32 t, const vec3& p ) const ;
        int32 edge_id( int32 t, int32 p0, int32 p1 ) const ;
        vec3 facet_normal( int32 t ) const ;
        void point_normal( std::vector< vec3 >& normals ) const ;
        Box3d bbox() const ;
        bool contains( const vec3& p ) const ;
        int32 find( const vec3& p ) const ;

        int32 triangles_around_point(
            int32 shared_point,
            std::vector<int>& triangles, 
            bool border_only ) const ;
        vec3 barycenter( int32 f ) const ;
        int32 triangles_around_point_with_hint(
            int32 shared_point,
            std::vector< int32 >& result,
            bool border_only,
            int32 triangle_hint ) const ;
        uint32 closest_point_in_facet( uint32 f, const float64* point ) const {
            return closest_point_in_facet( f, vec3( point ) ) ;
        }
        uint32 closest_point_in_facet( uint32 f, const vec3& point ) const ;

        float64 resolution( int32 p ) const { return resolution_[p] ; }
        float64 resolution( int32 f, int32 v ) const { return resolution_[facets_[facet_begin(f)+v]] ; }
        float64 facet_resolution( uint32 f ) const ;
        const vec3& U( int32 f ) const { return U_[f] ; }
        const vec3& V( int32 f ) const { return V_[f] ; }
        const vec3& W( int32 f ) const { return W_[f] ; }
    private:
        int32& adjacent( int32 f, int32 e ) { return adjacent_[facet_begin(f)+e] ; }

        void compute_is_triangulated() {
            for( uint32 f = 0; f < nb_simplices(); f++ ) {
                if( !is_triangle( f ) ) {
                    is_triangulated_ = false ;
                    break ;
                }
            }
        }
        void compute_adjacent_facets() ;
        bool key_facet_orientation() ;
        int32 has_edge( int32 trgl, int32 id0, int32 id1 ) const ;

        bool same_point( int32 i, int32 j ) const {
            if( i == j ) return true ;
            else if( points_[i] == points_[j] ) return true ;
            else return false ;
        }


        bool set_points_and_facets(
            const std::vector< uint32 >& points,
            const std::vector< uint32>& facets,
            const std::vector< uint32>& facet_ptr,
            bool compute_adjacent = true
        ){
            points_ = points ;
            facets_ = facets ;
            facet_ptr_ = facet_ptr ;
            compute_is_triangulated() ;
            if( compute_adjacent ) compute_adjacent_facets() ;
            return key_facet_orientation() ;
        }
        void set_adjacent_facets( const std::vector<int32>& in ){
            grgmesh_assert( in.size() == facets_.size() ) ;
            adjacent_.resize(0) ;
            adjacent_ = in ;
        }

    private:
        KeyFacet key_facet_ ;

        std::vector< uint32 > points_ ;
        std::vector< uint32 > facets_ ;
        std::vector< uint32 > facet_ptr_ ;

        /// adjacent triangles one per edge, -1 if on the border
        std::vector< int32 > adjacent_ ;

        bool is_triangulated_ ;

        //Attributes
        std::vector< float64 > resolution_ ;
        std::vector< vec3 > U_, V_, W_ ;

    };  

    class GRGMESH_API SurfacePartMutator {
    public:
        SurfacePartMutator( SurfacePart& M )
            : M_( M )
        {
        }
        SurfacePartMutator( const SurfacePart& M )
            : M_( const_cast< SurfacePart& >( M ) )
        {
        }
        void set_point( uint32 id, const float64* p ) {
            set_point( id, vec3( p ) ) ;
        }
        void set_point( uint32 id, const vec3& p ) ;
        vec3& point( uint32 p ) const ;
        std::vector< uint32 >& points() const { return M_.points_ ; }
        std::vector< uint32 >& facets() const { return M_.facets_ ; }
        std::vector< uint32 >& facet_ptr() const { return M_.facet_ptr_ ; }
        std::vector< int32 >& adjacents() const { return M_.adjacent_ ; }

        std::vector< float64 >& resolution() const { return M_.resolution_ ; }
        std::vector< vec3 >& U() const { return M_.U_ ; }
        std::vector< vec3 >& V() const { return M_.V_ ; }
        std::vector< vec3 >& W() const { return M_.W_ ; }

        void clear()
        {
            M_.points_.clear() ;
            M_.facets_.clear() ;
            M_.facet_ptr_.clear() ;
            M_.adjacent_.clear() ;
            M_.resolution_.clear() ;
            M_.U_.clear() ;
            M_.V_.clear() ;
            M_.W_.clear() ;
        }


    private:
        SurfacePart& M_ ;
    };
}
#endif

