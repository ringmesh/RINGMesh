/*! File part of phD work of Jeanne Pellerin ? copyright ASGA ? 
 *  
 * 
 *  Developed with Arnaud Botella
 *
 *  Will GRGMesh go open source ?
 */ 

#ifndef __GRGMESH_BOUNDARY_MODEL_ELEMENT__
#define __GRGMESH_BOUNDARY_MODEL_ELEMENT__

#include <grgmesh/common.h>
#include <grgmesh/grgmesh_assert.h>

#include <vector> 
#include <map>
#include <string>
#include <algorithm> 
#include <cmath>
#include <cassert>

namespace GRGMesh {
    class BoundaryModel ;
}


namespace GRGMesh {


    /********************************************************************************************/
    /******             General stuff               *********************************************/
    /********************************************************************************************/   
    
   /**
     * \brief A triangle that set the orientation of one Surface's facets
     * 
     */
    struct KeyFacet {
        KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 ):
            p0_(p0), p1_(p1), p2_(p2){} ;

        KeyFacet():p0_(), p1_(), p2_() {} ;
    
        bool is_default() const {
            if( p0_.x == 1 && p0_.y == -1 && p0_.z == -1 &&
                p1_.x == 1 && p1_.y == -1 && p1_.z == -1  && 
                p2_.x == 1 && p2_.y == -1 && p2_.z == -1  ) return true ;
            else {
                // la compilation ne passe pas �a m'�nerve
                /// \todo put this assert back
                //assert( !(p0_==p1_) && !(p0_==p2_) && !(p1_==p2_) ) ;
                return false ;
            }
        }
    public:
        vec3 p0_ ;
        vec3 p1_ ;
        vec3 p2_ ;
    } ;

    /**
     * Types for BoundaryModelElement 
     * \todo Read all types, this is not sufficient
     */
    enum GEOL_FEATURE {
        ALL,
        STRATI,
        FAULT,
        VOI,
        STRATI_FAULT,
        STRATI_VOI,
        FAULT_VOI
    } ;
    /// Default type is all
    static GEOL_FEATURE default_type = ALL ;

    /********************************************************************************************/
    /***********       Components of a Model declaration                 ************************/
    /********************************************************************************************/

    /**
     * \brief Generic class describing the elements of a BoundaryModel
    */
    class GRGMESH_API BoundaryModelElement {
        friend class BoundaryModelBuilder ;

    public:
        const static BoundaryModelElement dummy_element = BoundaryModelElement( nil, 0 ) ;
        const static index_t NO_ID = index_t( -1 ) ;

        BoundaryModelElement(
            BoundaryModel* model, index_t dim, index_t id = NO_ID, index_t parent = NO_ID,
            GEOL_FEATURE type = default_type )
            : model_( model ),  name_( "" ), id_( id ), dim_( dim ), type_( type ),
            parent_( parent )
        {
        }
        
        virtual ~BoundaryModelElement() { }

        const BoundaryModel& model() const { return *model_ ; }
        
        const std::string& name() const { return name_ ; }
        index_t id() const { return id_ ; }
        index_t dim() const { return dim_ ; }
        GEOL_FEATURE type() const { return type_ ; }

        bool is_on_voi() const ;
        bool side( index_t i ) const { return sides_[i] ; }

        bool has_parent() const { return parent_ != NO_ID ; }
        const BoundaryModelElement& parent() const ;
        index_t parent_id() const { return parent_ ; }
        
        index_t nb_boundaries() const { return boundaries_.size() ; }
        index_t boundary_id( index_t x ) const { return boundaries_[x] ; }
        const BoundaryModelElement& boundary( index_t x ) const ;
        
        index_t nb_in_boundary() const { return in_boundary_.size() ; }
        index_t in_boundary_id( index_t x ) const { return in_boundary_[x] ; }
        const BoundaryModelElement& in_boundary( index_t x ) const ;
        
        index_t nb_children() const { return children_.size() ; }
        index_t child_id( index_t x ) const { return children_[x] ; }
        const BoundaryModelElement& child( index_t x ) const ;
        

        virtual index_t nb_cells() const {
            grgmesh_assert_not_reached ;  return 0 ;
        }
        virtual index_t nb_points() const {
            grgmesh_assert_not_reached ;  return 0 ;
        }
        virtual index_t model_point_id( index_t p = 0 ) const {
            grgmesh_assert_not_reached ; return 0 ;
        }
  
        virtual const vec3& point( index_t p = 0 ) const {
            grgmesh_assert_not_reached ; return dummy_vec3 ;
        }
  
        
        
        /*virtual bool is_triangulated() const {
            grgmesh_assert( dim_ == 3 ) ;
            for( index_t s = 0; s < nb_boundaries(); s++ ) {
                if( !boundary( s ).is_triangulated() ) return false ;
            }
            return true ;
        }*/
        
    protected:
        void copy_macro_topology(
            const BoundaryModelElement& rhs, BoundaryModel& model ) ;
        
        void set_parent( index_t p ){ parent_ = p ; }
        void set_name( const std::string& name ) { name_ = name ; }
        void set_type( GEOL_FEATURE type ) { type_ = type ; } 
        void set_dim( index_t dim ) { dim_ = dim ; }
        void set_id( index_t id ) { id_ = id ; }
        
        void add_boundary( index_t b ) { boundaries_.push_back( b ) ; }
        void add_boundary( index_t b, bool side ) { boundaries_.push_back(b) ; sides_.push_back(side) ; }
        virtual void add_in_boundary( index_t e ) { in_boundary_.push_back(e) ; }
        void add_child( index_t e ){ children_.push_back( e ) ; }
       

    protected :
        /// \todo Can we have something else than a POINTER ?? default constructor is needed...
        BoundaryModel* model_ ;

        /// Name of the element, empty string if none nothing
        std::string name_ ;

        /// Id of this element in the appropriate vector of the BoundaryModel owning it
        index_t id_ ;
        
        /// Dimension of the element 0 corner; 1 line; 2 surface; 3 region
        index_t dim_ ;
        
        /// Geological type for this object, default is ALL 
        GEOL_FEATURE type_ ;

        /// Elements on the boundary of this element - their dimension is dim_-1
        std::vector< index_t > boundaries_ ;

        /// Flag on which side of the boundary is this element 
        /// Filled for volumetric regions only + (true) or - (false)
        std::vector< bool > sides_ ; 
        
        /// Elements in which boundary this element is - their dimension is dim_+1
        std::vector< index_t > in_boundary_ ;

        /// Index of the parent (group of elements to which belong this).
        /// Default value is NO_ID.
        index_t parent_ ;

        /// The group elements making up this one, empty for basic elements
        std::vector< index_t > children_ ;
    } ;

    /*-----------------------------------------------------------------------------------------*/
    /*! Corners of the model
    *  Dimension 0
    */
    class GRGMESH_API Corner : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
    public:
        Corner(
            BoundaryModel* model,
            signed_index_t id = -1,
            index_t p = 0 )
            : BoundaryModelElement( model, 0, id ), p_( p )
        {
        }
        virtual ~Corner()
        {
        }
        bool is_real() const { return in_boundary_.size() > 1 ; }

        virtual const vec3& point( index_t p = 0 ) const ;
        
        virtual index_t nb_cells() const { return 1 ; }
        virtual index_t nb_points() const { return 1 ; }
        virtual index_t model_point_id( index_t p = 0 ) const { return p_ ; } 
        
    private:
        //void copy_macro_topology( const Corner& rhs, BoundaryModel& model ) ;

        void set_point( index_t p ) { p_ = p ; }
    private:
        index_t p_ ;
    };
    /*-----------------------------------------------------------------------------------------*/        
    /*! A part of the contact between 2 interfaces
    * 
    *  boundaries -> 1 or 2 corners
    *  in_boundary -> several surface parts
    */
    class GRGMESH_API Line: public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class LineMutator ;
    public:
        Line( BoundaryModel* model, index_t id = NO_ID ) ;
        Line(
            BoundaryModel* model,
            index_t id,
            const std::vector< index_t >& points ) ;
        Line(
            BoundaryModel* model,
            index_t id,
            index_t corner0,
            index_t corner1,
            const std::vector< index_t >& points ) ;
        virtual ~Line(){} ;

   
        virtual index_t nb_cells() const {            
            return points_.size()-1 ; 
        }
        // If the line is closed the last point is equal to the first one
        virtual index_t nb_points() const {             
            return points_.size() ; 
        }
        virtual index_t model_point_id( index_t p ) const {
            return points_.at(p) ;
        }

        //bool contains( const vec3& p ) const ;
        //int find( const vec3& p ) const ;
        //int facet_point_id( const vec3& p ) const ;
        bool is_closed () const { 
            return (boundaries_[0]!= nil ) && (boundaries_[0] == boundaries_[1]) ; 
        }  

        // Returns true if this line is twice in the boundary of the surface e
        bool is_inside_border( const BoundaryModelElement& e ) const ;
            
        virtual const vec3& point( index_t line_point_id ) const ;
        inline vec3 segment_barycenter( index_t s ) const ;

        inline double segment_length( index_t s ) const ;
            

    private:
        /*void copy_macro_topology(
            const Line& rhs,
            BoundaryModel& model ) ;*/
        void set_vertices( const std::vector< index_t >& vertices ) {
            points_ = vertices ;
        }
        //virtual void add_in_boundary( index_t e ) ;
        //void set_is_inside_border( bool x ) { is_inside_border_.push_back(x) ; }
    private:
        /// In case of a closed line, the last vertex (equal to the first) is not stored.
        /// ATTENTION �A A CHANG� SINON C'EST TROP CHIANT, et en fait c'�tait pas bon
        std::vector< index_t > points_ ;
                 
        //std::vector< bool > is_inside_border_ ;
    } ;

    class GRGMESH_API LineMutator {
    public:
        LineMutator( Line& M )
            : M_( M )
        {
        }
        LineMutator( const Line& M )
            : M_( const_cast< Line& >( M ) )
        {
        }
        void set_point( index_t id, const vec3& p ) ;
        vec3& point( index_t p ) const ;
        std::vector< index_t >& points() const { return M_.points_ ; }

        void clear()
        {
            M_.points_.clear() ;
        }

    private:
        Line& M_ ;
    };
    /*-----------------------------------------------------------------------------------------*/
       
    /**
    * \brief A polygonal manifold surface part of a BoundaryModel
    * 
    * This is a BoundaryModelElement of dimension 2.
    * Its boundaries_ is made of Line - It is in_boundary_ of 1 or 2 Region
    *
    */
    class GRGMESH_API Surface : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class SurfaceMutator ;
    public:
        const static index_t NO_ADJACENT = index_t( -1 ) ;
        Surface(
            BoundaryModel* model,
            index_t id = NO_ID,
            index_t parent = NO_ID,
            const GEOL_FEATURE& type = default_type )
            :
            BoundaryModelElement( model, 2, id, parent, type ),
            is_triangulated_( false )
        {
        }
        virtual ~Surface(){} ;
       
        virtual index_t nb_cells() const { return facets_.empty() ? 0 : facet_ptr_.size() - 1 ; }
        virtual index_t nb_points() const { return points_.size() ; }               
        virtual index_t model_point_id( index_t p ) const {
            return points_[p] ;
        }
        
        //bool key_facet_orientation() ;
        const KeyFacet& key_facet() const { return key_facet_ ; }       

        bool is_triangulated() const { return is_triangulated_ ; }
        
        index_t facet_begin( index_t f ) const { return facet_ptr_.at(f) ; }
        index_t facet_end( index_t f ) const { return facet_ptr_.at(f+1) ; }
        index_t nb_points_in_facet( index_t f ) const { return facet_end( f ) - facet_begin( f ) ; }
        index_t next_in_facet( index_t f, index_t v ) const { 
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            if( v != nb_points_in_facet(f)-1 ) return v+1 ;
            else return 0 ;
        }
        index_t prev_in_facet( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            if( v > 0 ) return v-1 ;
            else return nb_points_in_facet(f)-1 ;
        }            
        
        const vec3& point( index_t f, index_t v ) const ;

        /** 
         * Access to the ids of the points in the model
         * The same point can appear several times - boundaries inside 
         * Necessary for export
         * Should be inline - but linking fails because implementation is in cpp
         * and cannot be there because including boundary_model.h results in circular includes
         */        
        virtual const vec3& point( index_t surf_point_id ) const ;
        
        /** Returns the id of point \param v in facet \param f in this surface */
        index_t surf_point_id( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return facets_[facet_begin(f)+v] ; 
        }
        /** Returns the id of point  \param v, in facet \param f in the parent BoundaryModel */ 
        index_t model_point_id( index_t f, index_t v ) const { 
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return points_[ surf_point_id( f, v ) ] ;
        }     

        index_t surf_point_id( index_t model_point_id ) const {
            for( index_t i = 0; i < points_.size() ; ++i ){
                if ( points_[i] == model_point_id ) return i ;
            }
            return NO_ID ;
        }
        
        /** Returns the id of the adjacent facet of \param f in this surface along the edge starting at \param v */
        index_t adjacent( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return adjacent_[facet_begin(f)+v] ; 
        }
        //int adjacent_in_neighbor( index_t f, index_t e ) const ;
        
        bool is_on_border( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return adjacent( f, v ) == -1 ; 
        }
        bool is_on_border( index_t f ) const {
            for( index_t adj = 0; adj < nb_points_in_facet(f); adj++ ) {
                if( is_on_border( f, adj ) ) return true ;
            }
            return false ;
        }
        
        /**
         * For a bidirectional border traversal 
         * From the input facet f, get the facet that share point v and 
         * get the ids of point v and of the following point in this next facet
         * next facet may be the same, from is required to avoid going back 
         */
        void next_on_border(
            index_t f, index_t from, index_t v,
            index_t& next_f, index_t& v_in_next, index_t& to  ) const ;
        
        /**
         * One directional border traversal
         */
        void next_on_border(
            index_t f,
            index_t e,
            index_t& next_f,
            index_t& next_e ) const ;

        //FacetEdge next_on_border( const FacetEdge& te ) const ;
       
        bool is_triangle( index_t f ) const { return nb_points_in_facet( f ) == 3 ; }

//        signed_index_t find_facet( signed_index_t id0, signed_index_t id1 ) const ;
        index_t facet_from_model_point_ids( index_t i0, index_t i1 ) const ;
        void edge_from_model_point_ids(
            index_t i0,
            index_t i1,
            index_t& f,
            index_t& e ) const ;
        void oriented_edge_from_model_point_ids(
            index_t i0,
            index_t i1,
            index_t& facet,
            index_t& edge ) const ;

        
        //int find_facet( const vec3& p0, const vec3& p1 ) const ; // TO DO A VIRER !!
        //int find_edge( signed_index_t id0, signed_index_t id1 ) const ;

        /** Returns the id of the point with surface id p0 in the given facet */
        index_t facet_point_id( index_t t, index_t surf_point_id ) const ;
        /** Returns the id of the points with the given coordinates in the given facet */
     //   signed_index_t facet_point_id( signed_index_t t, const vec3& p ) const ;
        //int edge_id( signed_index_t t, signed_index_t p0, signed_index_t p1 ) const ;

        //Box3d bbox() const ;
        //bool contains( const vec3& p ) const ;
        //int find( const vec3& p ) const ;

        index_t facets_around_point(
            index_t surf_point_id, 
            std::vector< index_t >& result, 
            bool border_only ) const ;
        
        vec3 facet_barycenter( index_t f ) const ;
        double facet_area( index_t f ) const ;
        vec3 facet_normal( index_t f ) const ;
        void point_normals( std::vector< vec3 >& normals ) const ;
        
        index_t facets_around_point(
            index_t surf_point_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet ) const ;

        index_t closest_point_in_facet( index_t f, const vec3& point ) const ;

    private:
        void set_key_facet( const KeyFacet& key ) { key_facet_ = key ; }
        void set_first_triangle_as_key() ;
        
        // On en a besoin � la construction
        // On ne peut pas utiliser model_point ids - car on ne diff�rencie pas
        // alors les arr�tes sur des bords internes = 1 bord - mais points identiques sur
        // triangles en face dans la m�me surfaace
        index_t facet_from_surface_point_ids( index_t i0, index_t i1 ) const ;


        void set_adjacent( index_t f, index_t e, index_t adjacent ) {
            adjacent_[facet_begin(f)+e] = adjacent ;
        }

        void compute_is_triangulated() {
            for( index_t f = 0; f < nb_cells(); f++ ) {
                if( !is_triangle( f ) ) {
                    is_triangulated_ = false ;
                    return ;
                }
            }
            is_triangulated_ = true ;
        }       
        //int has_edge( signed_index_t f, signed_index_t v0, signed_index_t v1 ) const ;
        //int has_oriented_edge( signed_index_t f, signed_index_t v0, signed_index_t v1 ) const ;

        /*bool same_point( signed_index_t i, signed_index_t j ) const {
            if( i == j ) return true ;
            else if( points_[i] == points_[j] ) return true ;
            else return false ;
        }*/


        void set_geometry(
            const std::vector< index_t >& points,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr )
        {
            // Are these copies parallelized ?
            points_ = points ;
            facets_ = facets ;
            facet_ptr_ = facet_ptr ;

            compute_is_triangulated() ;
        }

        void set_adjacent( const std::vector< index_t >& adjacent ){
            grgmesh_assert( adjacent.size() == facets_.size() ) ;
            adjacent_ = adjacent ;
        }

       
    private:
        KeyFacet key_facet_ ;

        std::vector< index_t > points_ ;
        std::vector< index_t > facets_ ;
        std::vector< index_t > facet_ptr_ ;

    
        // The adjacent facet is given for each vertex of each facet for the edge
        // starting at this vertex.
        // NO_ADJACENT if the edge is along a Line (Surface boundary)
        std::vector< index_t > adjacent_ ;

        bool is_triangulated_ ;
        
    };  


    class GRGMESH_API SurfaceMutator {
    public:
        SurfaceMutator( Surface& M )
            : M_( M )
        {
        }
        SurfaceMutator( const Surface& M )
            : M_( const_cast< Surface& >( M ) )
        {
        }
        void set_point( index_t id, const vec3& p ) ;
        vec3& point( index_t p ) const ;
        std::vector< index_t >& points() const { return M_.points_ ; }
        std::vector< index_t >& facets() const { return M_.facets_ ; }
        std::vector< index_t >& facet_ptr() const { return M_.facet_ptr_ ; }
        std::vector< index_t >& adjacents() const { return M_.adjacent_ ; }

        void clear()
        {
            M_.points_.clear() ;
            M_.facets_.clear() ;
            M_.facet_ptr_.clear() ;
            M_.adjacent_.clear() ;
        }

    private:
        Surface& M_ ;
    };



    /**
     * \brief Class to answer geometrical request on BoundaryModelElement
     */
    class GRGMESH_API BoundaryModelElementMeasure {
    public:
        static double size( const BoundaryModelElement* E ) ;
        static double cell_size( const BoundaryModelElement* E, index_t cell ) ;
        static double distance( const BoundaryModelElement* from,  const vec3& p ) ;
        static double distance( const BoundaryModelElement* from, const BoundaryModelElement* to ) ;
        static vec3 barycenter ( const BoundaryModelElement* E ) ;
        static vec3 barycenter ( const BoundaryModelElement* E, const std::vector< index_t >& cells ) ;
    } ;

} // namespace

#endif

